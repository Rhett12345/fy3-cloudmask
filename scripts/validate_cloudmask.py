#!/usr/bin/env python3
"""Automated validation of cloud mask products.

Computes statistical metrics for cloud mask quality assessment.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path
from collections import Counter

import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def flood_fill_label(mask: np.ndarray) -> tuple[np.ndarray, int]:
    """Simple connected-component labeling using flood fill.

    Args:
        mask: 2D boolean array.

    Returns:
        Tuple of (labeled_array, n_labels).
    """
    h, w = mask.shape
    labeled = np.zeros((h, w), dtype=np.int32)
    label = 0

    for i in range(h):
        for j in range(w):
            if mask[i, j] and labeled[i, j] == 0:
                label += 1
                # BFS
                queue = [(i, j)]
                labeled[i, j] = label
                while queue:
                    ci, cj = queue.pop(0)
                    for di, dj in [(-1,0),(1,0),(0,-1),(0,1)]:
                        ni, nj = ci+di, cj+dj
                        if 0 <= ni < h and 0 <= nj < w and mask[ni, nj] and labeled[ni, nj] == 0:
                            labeled[ni, nj] = label
                            queue.append((ni, nj))

    return labeled, label


def compute_cloud_mask_metrics(cm_array: np.ndarray, conf_array: np.ndarray) -> dict:
    """Compute cloud mask quality metrics.

    Args:
        cm_array: (n_elem, n_line) cloud mask values (0-3).
        conf_array: (n_elem, n_line) confidence values.

    Returns:
        Dictionary of metrics.
    """
    valid = (cm_array >= 0) & (cm_array <= 3)
    n_valid = np.sum(valid)
    if n_valid == 0:
        return {}

    cm = cm_array[valid]
    conf = conf_array[valid]

    # Basic distribution
    counts = Counter(cm)
    total = len(cm)

    cloudy_pct = counts.get(0, 0) / total * 100
    prob_cloudy_pct = counts.get(1, 0) / total * 100
    prob_clear_pct = counts.get(2, 0) / total * 100
    clear_pct = counts.get(3, 0) / total * 100

    # Cloud fraction (cloudy + prob_cloudy)
    cloud_fraction = (counts.get(0, 0) + counts.get(1, 0)) / total * 100

    # Connected region analysis (for cloud pixels)
    cloud_mask = (cm_array == 0) | (cm_array == 1)

    # Boundary length estimation
    grad_x = np.abs(np.diff(cloud_mask.astype(int), axis=0))
    grad_y = np.abs(np.diff(cloud_mask.astype(int), axis=1))
    boundary_length = int(np.sum(grad_x) + np.sum(grad_y))

    # Confidence statistics (nan-aware)
    conf_mean = float(np.nanmean(conf))
    conf_std = float(np.nanstd(conf))
    conf_median = float(np.nanmedian(conf))

    # Confidence distribution
    conf_hist, _ = np.histogram(conf, bins=[0, 0.5, 0.66, 0.95, 0.99, 1.01])

    return {
        'total_pixels': total,
        'cloudy_pct': cloudy_pct,
        'prob_cloudy_pct': prob_cloudy_pct,
        'prob_clear_pct': prob_clear_pct,
        'clear_pct': clear_pct,
        'cloud_fraction': cloud_fraction,
        'boundary_length': boundary_length,
        'conf_mean': conf_mean,
        'conf_std': conf_std,
        'conf_median': conf_median,
        'conf_hist_0_05': int(conf_hist[0]),
        'conf_hist_05_066': int(conf_hist[1]),
        'conf_hist_066_095': int(conf_hist[2]),
        'conf_hist_095_099': int(conf_hist[3]),
        'conf_hist_099_1': int(conf_hist[4]),
    }


def compare_metrics(f90_metrics: dict, py_metrics: dict) -> dict:
    """Compare metrics between Fortran and Python.

    Args:
        f90_metrics: Fortran metrics.
        py_metrics: Python metrics.

    Returns:
        Dictionary of differences.
    """
    diff = {}
    for key in f90_metrics:
        if key in py_metrics and isinstance(f90_metrics[key], (int, float)):
            f90_val = f90_metrics[key]
            py_val = py_metrics[key]
            diff[f'{key}_f90'] = f90_val
            diff[f'{key}_py'] = py_val
            if py_val != 0:
                diff[f'{key}_ratio'] = f90_val / py_val
            diff[f'{key}_diff'] = f90_val - py_val
    return diff


def print_metrics(label: str, metrics: dict):
    """Print metrics in a formatted way."""
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  Total pixels:       {metrics.get('total_pixels', 0):>12,}")
    print(f"  Cloudy:             {metrics.get('cloudy_pct', 0):>11.2f}%")
    print(f"  Prob Cloudy:        {metrics.get('prob_cloudy_pct', 0):>11.2f}%")
    print(f"  Prob Clear:         {metrics.get('prob_clear_pct', 0):>11.2f}%")
    print(f"  Confident Clear:    {metrics.get('clear_pct', 0):>11.2f}%")
    print(f"  Cloud Fraction:     {metrics.get('cloud_fraction', 0):>11.2f}%")
    print(f"  Boundary Length:    {metrics.get('boundary_length', 0):>12,}")
    print(f"  Conf Mean:          {metrics.get('conf_mean', 0):>11.4f}")
    print(f"  Conf Std:           {metrics.get('conf_std', 0):>11.4f}")
    print(f"  Conf Median:        {metrics.get('conf_median', 0):>11.4f}")
    print(f"  Conf [0, 0.5):      {metrics.get('conf_hist_0_05', 0):>12,}")
    print(f"  Conf [0.5, 0.66):   {metrics.get('conf_hist_05_066', 0):>12,}")
    print(f"  Conf [0.66, 0.95):  {metrics.get('conf_hist_066_095', 0):>12,}")
    print(f"  Conf [0.95, 0.99):  {metrics.get('conf_hist_095_099', 0):>12,}")
    print(f"  Conf [0.99, 1.0]:   {metrics.get('conf_hist_099_1', 0):>12,}")


def main():
    parser = argparse.ArgumentParser(description="Validate cloud mask products")
    parser.add_argument("--date", default="20220803", help="Date (YYYYMMDD)")
    parser.add_argument("--orbit", default="0740", help="Orbit time tag (HHMM)")
    args = parser.parse_args()

    date_str = args.date
    time_tag = args.orbit

    # Find files
    mersi_dir = Path('/data/Data_yuq/mersi') / date_str
    l1b_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF"
    geo_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_GEO1K_MS.HDF"
    nwp_path = Path('/data/nwp') / date_str / 'ORG' / f'gfs0p25_41L_{date_str}_06_00'

    if not l1b_path.exists() or not geo_path.exists():
        logger.error("Input files not found")
        sys.exit(1)

    # Read data
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

    logger.info("Reading data...")
    pxldat = read_l1b_data(str(l1b_path))
    geo = read_geo_data(str(geo_path))
    nwp = read_nwp_binary(str(nwp_path))
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    logger.info(f"Swath: {n_elem} x {n_line}")

    # Run Fortran backend
    logger.info("Running Fortran backend...")
    t0 = time.time()
    result_f90 = process_swath_native(
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
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    t_f90 = time.time() - t0
    logger.info(f"Fortran: {t_f90:.1f}s")

    cm_f90 = result_f90['cloud_mask']
    conf_f90 = result_f90['confidence']

    # Compute metrics
    logger.info("Computing Fortran metrics...")
    f90_metrics = compute_cloud_mask_metrics(cm_f90, conf_f90)
    print_metrics("Fortran Cloud Mask Metrics", f90_metrics)

    # Save metrics to file
    output_dir = Path(__file__).resolve().parent.parent / 'output' / 'validation'
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics_file = output_dir / f'metrics_{date_str}_{time_tag}.txt'
    with open(metrics_file, 'w') as f:
        f.write(f"Date: {date_str}\n")
        f.write(f"Orbit: {time_tag}\n")
        f.write(f"Swath: {n_elem} x {n_line}\n")
        f.write(f"Fortran time: {t_f90:.1f}s\n\n")
        for key, value in f90_metrics.items():
            f.write(f"{key}: {value}\n")

    logger.info(f"Metrics saved to {metrics_file}")

    # Per-region analysis
    logger.info("Performing per-region analysis...")
    lat = geo['lat']

    # Polar region (|lat| > 60)
    polar_mask = np.abs(lat) > 60
    if np.any(polar_mask):
        polar_cm = cm_f90[polar_mask]
        polar_conf = conf_f90[polar_mask]
        polar_counts = Counter(polar_cm)
        polar_total = len(polar_cm)
        print(f"\n  Polar region (|lat| > 60):")
        print(f"    Pixels: {polar_total:,}")
        print(f"    Cloudy: {polar_counts.get(0, 0) / polar_total * 100:.2f}%")
        print(f"    Prob Cloudy: {polar_counts.get(1, 0) / polar_total * 100:.2f}%")
        print(f"    Prob Clear: {polar_counts.get(2, 0) / polar_total * 100:.2f}%")
        print(f"    Confident Clear: {polar_counts.get(3, 0) / polar_total * 100:.2f}%")

    # Non-polar region
    non_polar_mask = np.abs(lat) <= 60
    if np.any(non_polar_mask):
        non_polar_cm = cm_f90[non_polar_mask]
        non_polar_counts = Counter(non_polar_cm)
        non_polar_total = len(non_polar_cm)
        print(f"\n  Non-polar region (|lat| <= 60):")
        print(f"    Pixels: {non_polar_total:,}")
        print(f"    Cloudy: {non_polar_counts.get(0, 0) / non_polar_total * 100:.2f}%")
        print(f"    Prob Cloudy: {non_polar_counts.get(1, 0) / non_polar_total * 100:.2f}%")
        print(f"    Prob Clear: {non_polar_counts.get(2, 0) / non_polar_total * 100:.2f}%")
        print(f"    Confident Clear: {non_polar_counts.get(3, 0) / non_polar_total * 100:.2f}%")


if __name__ == "__main__":
    main()
