#!/usr/bin/env python3
"""Spatial structure analysis for cloud mask products.

Uses numpy vectorized operations for speed (no scipy dependency).
Analyzes a center crop to keep computation tractable.

Usage:
    python scripts/spatial_analysis.py --date 20220803 --orbit 0740
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def label_connected_components(binary: np.ndarray) -> tuple[np.ndarray, int]:
    """Union-Find based connected component labeling (4-connectivity)."""
    h, w = binary.shape
    parent = np.arange(h * w, dtype=np.int32)
    rank = np.zeros(h * w, dtype=np.int32)
    labels = np.zeros(h * w, dtype=np.int32)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx == ry:
            return
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1

    flat = binary.ravel()
    # First pass: union neighbors
    for i in range(h):
        for j in range(w):
            idx = i * w + j
            if flat[idx] == 0:
                continue
            # Left neighbor
            if j > 0 and flat[idx - 1] == 1:
                union(idx, idx - 1)
            # Up neighbor
            if i > 0 and flat[idx - w] == 1:
                union(idx, idx - w)

    # Second pass: assign labels
    label_map = {}
    current_label = 0
    for idx in range(h * w):
        if flat[idx] == 0:
            continue
        root = find(idx)
        if root not in label_map:
            current_label += 1
            label_map[root] = current_label
        labels[idx] = label_map[root]

    return labels.reshape(h, w), current_label


def numpy_boundary_length(binary: np.ndarray) -> int:
    """Compute boundary length using numpy shifts."""
    # A pixel is on the boundary if it's 1 and any neighbor is 0 or edge
    h, w = binary.shape
    # Pad with zeros
    padded = np.pad(binary, 1, mode='constant', constant_values=0)
    # Check 4 neighbors
    up = padded[:-2, 1:-1]
    down = padded[2:, 1:-1]
    left = padded[1:-1, :-2]
    right = padded[1:-1, 2:]
    # Boundary: binary=1 AND (any neighbor=0)
    is_boundary = (binary == 1) & ((up == 0) | (down == 0) | (left == 0) | (right == 0))
    return int(np.sum(is_boundary))


def numpy_neighbor_agreement(cm: np.ndarray) -> float:
    """Compute fraction of interior pixels where >=3 of 4 neighbors agree."""
    h, w = cm.shape
    padded = np.pad(cm, 1, mode='edge')
    up = padded[:-2, 1:-1]
    down = padded[2:, 1:-1]
    left = padded[1:-1, :-2]
    right = padded[1:-1, 2:]
    agree = ((up == cm).astype(int) + (down == cm).astype(int) +
             (left == cm).astype(int) + (right == cm).astype(int))
    # Interior pixels only
    interior = np.ones_like(cm, dtype=bool)
    interior[0, :] = False
    interior[-1, :] = False
    interior[:, 0] = False
    interior[:, -1] = False
    return float(np.mean(agree[interior] >= 3))


def numpy_salt_pepper(cm: np.ndarray) -> dict:
    """Detect salt-and-pepper noise: isolated pixels different from all 4 neighbors."""
    h, w = cm.shape
    padded = np.pad(cm, 1, mode='edge')
    up = padded[:-2, 1:-1]
    down = padded[2:, 1:-1]
    left = padded[1:-1, :-2]
    right = padded[1:-1, 2:]
    # Isolated: all 4 neighbors differ from center
    all_diff = (up != cm) & (down != cm) & (left != cm) & (right != cm)
    # Also check: at least 3 neighbors differ
    n_diff = ((up != cm).astype(int) + (down != cm).astype(int) +
              (left != cm).astype(int) + (right != cm).astype(int))
    isolated = int(np.sum(all_diff))
    near_isolated = int(np.sum(n_diff >= 3))
    return {
        'fully_isolated': isolated,
        'near_isolated': near_isolated,
        'fully_isolated_pct': 100.0 * isolated / cm.size,
        'near_isolated_pct': 100.0 * near_isolated / cm.size,
    }


def compute_region_stats(labels: np.ndarray, n_regions: int) -> dict:
    """Compute region size statistics from labeled array."""
    if n_regions == 0:
        return {'count': 0, 'mean': 0, 'median': 0, 'max': 0, 'p10': 0, 'p90': 0,
                'n_iso': 0, 'iso_pct': 0}

    # Count pixels per label
    sizes = np.bincount(labels.ravel())[1:]  # skip label 0
    sizes = sizes[sizes > 0]

    return {
        'count': n_regions,
        'mean': float(np.mean(sizes)),
        'median': float(np.median(sizes)),
        'max': int(np.max(sizes)),
        'p10': int(np.percentile(sizes, 10)) if len(sizes) > 0 else 0,
        'p90': int(np.percentile(sizes, 90)) if len(sizes) > 0 else 0,
        'n_iso': int(np.sum(sizes == 1)),
        'iso_pct': 100.0 * np.sum(sizes == 1) / max(n_regions, 1),
    }


def compute_spatial_metrics(cm: np.ndarray) -> dict:
    """Compute comprehensive spatial structure metrics."""
    total = cm.size
    cloudy = (cm == 0).astype(np.int8)
    clear = (cm == 3).astype(np.int8)
    prob_cloudy = (cm == 1).astype(np.int8)
    prob_clear = (cm == 2).astype(np.int8)

    # Cloud distribution
    cloudy_pct = 100.0 * np.sum(cloudy) / total
    clear_pct = 100.0 * np.sum(clear) / total
    prob_cloudy_pct = 100.0 * np.sum(prob_cloudy) / total
    prob_clear_pct = 100.0 * np.sum(prob_clear) / total

    # Connected regions (cloudy)
    logger.info("  Labeling cloudy connected regions...")
    t0 = time.time()
    cloud_labels, n_cloud = label_connected_components(cloudy)
    logger.info(f"  Cloudy regions: {n_cloud} ({time.time()-t0:.1f}s)")

    # Connected regions (clear)
    logger.info("  Labeling clear connected regions...")
    t0 = time.time()
    clear_labels, n_clear = label_connected_components(clear)
    logger.info(f"  Clear regions: {n_clear} ({time.time()-t0:.1f}s)")

    # Region stats
    cloud_stats = compute_region_stats(cloud_labels, n_cloud)
    clear_stats = compute_region_stats(clear_labels, n_clear)

    # Boundary length
    logger.info("  Computing boundary length...")
    t0 = time.time()
    boundary = numpy_boundary_length(cloudy)
    logger.info(f"  Boundary: {boundary} ({time.time()-t0:.1f}s)")

    # Neighbor agreement
    logger.info("  Computing neighbor agreement...")
    t0 = time.time()
    agreement = numpy_neighbor_agreement(cm)
    logger.info(f"  Agreement: {agreement:.4f} ({time.time()-t0:.1f}s)")

    # Salt-and-pepper noise
    logger.info("  Detecting salt-and-pepper noise...")
    t0 = time.time()
    sp = numpy_salt_pepper(cm)
    logger.info(f"  Done ({time.time()-t0:.1f}s)")

    return {
        'total_pixels': total,
        'cloudy_pct': cloudy_pct,
        'clear_pct': clear_pct,
        'prob_cloudy_pct': prob_cloudy_pct,
        'prob_clear_pct': prob_clear_pct,
        'cloud_fraction': cloudy_pct + prob_cloudy_pct,
        'cloud_regions': cloud_stats,
        'clear_regions': clear_stats,
        'boundary_length': boundary,
        'neighbor_agreement': agreement,
        'salt_pepper': sp,
    }


def print_metrics(metrics: dict):
    """Print metrics in readable format."""
    sp = metrics['salt_pepper']
    cr = metrics['cloud_regions']
    cl = metrics['clear_regions']

    print("\n" + "=" * 60)
    print("  Spatial Structure Analysis")
    print("=" * 60)

    print(f"\n  --- Cloud Distribution ---")
    print(f"  Total pixels:       {metrics['total_pixels']:>12,}")
    print(f"  Cloudy:             {metrics['cloudy_pct']:>11.2f}%")
    print(f"  Prob Cloudy:        {metrics['prob_cloudy_pct']:>11.2f}%")
    print(f"  Prob Clear:         {metrics['prob_clear_pct']:>11.2f}%")
    print(f"  Clear:              {metrics['clear_pct']:>11.2f}%")
    print(f"  Cloud Fraction:     {metrics['cloud_fraction']:>11.2f}%")

    print(f"\n  --- Connected Regions (Cloudy) ---")
    print(f"  Number of regions:  {cr['count']:>12,}")
    print(f"  Isolated (size=1):  {cr['n_iso']:>12,} ({cr['iso_pct']:.1f}%)")
    print(f"  Mean size:          {cr['mean']:>12.0f} px")
    print(f"  Median size:        {cr['median']:>12.0f} px")
    print(f"  Max size:           {cr['max']:>12,} px")
    print(f"  P10 / P90:          {cr['p10']:>6} / {cr['p90']:>6} px")

    print(f"\n  --- Connected Regions (Clear) ---")
    print(f"  Number of regions:  {cl['count']:>12,}")
    print(f"  Isolated (size=1):  {cl['n_iso']:>12,} ({cl['iso_pct']:.1f}%)")
    print(f"  Mean size:          {cl['mean']:>12.0f} px")
    print(f"  Median size:        {cl['median']:>12.0f} px")
    print(f"  Max size:           {cl['max']:>12,} px")

    print(f"\n  --- Spatial Quality ---")
    print(f"  Boundary length:    {metrics['boundary_length']:>12,}")
    print(f"  Neighbor agreement: {metrics['neighbor_agreement']:>11.4f}")
    print(f"  Fully isolated px:  {sp['fully_isolated']:>12,} ({sp['fully_isolated_pct']:.3f}%)")
    print(f"  Near-isolated (3/4):{sp['near_isolated']:>12,} ({sp['near_isolated_pct']:.3f}%)")
    print("=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Spatial analysis of cloud mask")
    parser.add_argument("--date", default="20220803", help="Date (YYYYMMDD)")
    parser.add_argument("--orbit", default="0740", help="Orbit time tag (HHMM)")
    parser.add_argument("--crop", type=int, default=800,
                        help="Crop size for analysis (0=full swath, slow)")
    args = parser.parse_args()

    # Run Fortran backend
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

    mersi_dir = Path('/data/Data_yuq/mersi') / args.date
    l1b_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{args.date}_{args.orbit}_1000M_MS.HDF"
    geo_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{args.date}_{args.orbit}_GEO1K_MS.HDF"
    nwp_path = Path('/data/nwp') / args.date / 'ORG' / f'gfs0p25_41L_{args.date}_06_00'

    logger.info("Reading data...")
    pxldat = read_l1b_data(str(l1b_path))
    geo = read_geo_data(str(geo_path))
    nwp = read_nwp_binary(str(nwp_path))
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    logger.info(f"Swath: {n_elem} x {n_line}")

    # Import and run native backend
    import os
    os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))
    from fy3_cloudmask.algorithm.native_backend import process_swath_native

    logger.info("Running Fortran backend...")
    t0 = time.time()
    result = process_swath_native(
        ref_vis=np.ascontiguousarray(pxldat[:, :, :19].astype(np.float32)),
        tbb_ir=np.ascontiguousarray(pxldat[:, :, 19:].astype(np.float32)),
        lat=np.ascontiguousarray(geo['lat'].astype(np.float32)),
        lon=np.ascontiguousarray(geo['lon'].astype(np.float32)),
        satzen=np.ascontiguousarray(geo['vza'].astype(np.float32)),
        solzen=np.ascontiguousarray(geo['sza'].astype(np.float32)),
        relaz=np.ascontiguousarray(np.zeros_like(geo['sza']).astype(np.float32)),
        glint=np.ascontiguousarray(geo['glint_angle'].astype(np.float32)),
        sfctmp=np.ascontiguousarray(nwp_interp['tsfc'].astype(np.float32)),
        pmsl=np.ascontiguousarray(nwp_interp['pmsl'].astype(np.float32)),
        uwind=np.ascontiguousarray(nwp_interp['u_wind'].astype(np.float32)),
        vwind=np.ascontiguousarray(nwp_interp['v_wind'].astype(np.float32)),
        tpw=np.ascontiguousarray(nwp_interp['tpw'].astype(np.float32)),
        elev=np.ascontiguousarray(geo['elevation'].astype(np.float32)),
        eco=np.ascontiguousarray(geo['eco_type'].astype(np.int8)),
        lsf=np.ascontiguousarray(geo['lsf'].astype(np.int8)),
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    t_f90 = time.time() - t0
    logger.info(f"Fortran: {t_f90:.1f}s")

    cm = result['cloud_mask']
    conf = result['confidence']

    # Crop for spatial analysis (Union-Find is O(n) but slow in Python)
    if args.crop > 0 and cm.shape[0] > args.crop and cm.shape[1] > args.crop:
        ci, cj = cm.shape[0] // 2, cm.shape[1] // 2
        half = args.crop // 2
        cm_crop = cm[ci-half:ci+half, cj-half:cj+half]
        logger.info(f"Cropped to {args.crop}x{args.crop} center for spatial analysis")
    else:
        cm_crop = cm

    logger.info("Computing spatial metrics...")
    metrics = compute_spatial_metrics(cm_crop)

    # Add confidence stats (full swath)
    valid = (cm >= 0) & (cm <= 3)
    conf_valid = conf[valid]
    metrics['conf_mean'] = float(np.mean(conf_valid))
    metrics['conf_std'] = float(np.std(conf_valid))
    metrics['conf_median'] = float(np.median(conf_valid))

    print_metrics(metrics)

    # Save
    out_dir = Path(__file__).resolve().parent.parent / 'output' / 'validation'
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f'spatial_{args.date}_{args.orbit}.txt'
    with open(out_path, 'w') as f:
        for key, val in metrics.items():
            if isinstance(val, dict):
                for k2, v2 in val.items():
                    f.write(f"{key}.{k2}: {v2}\n")
            else:
                f.write(f"{key}: {val}\n")
    logger.info(f"Saved to {out_path}")


if __name__ == '__main__':
    main()
