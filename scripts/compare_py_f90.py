#!/usr/bin/env python3
"""Compare Python vs Fortran cloud mask backends on a subset of pixels.

Usage:
    python scripts/compare_py_f90.py
    python scripts/compare_py_f90.py --date 20220803 --orbit 0740 --sample 50000
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_swath
from fy3_cloudmask.config import load_config

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

MERSI_ROOT = Path('/data/Data_yuq/mersi')
NWP_ROOT = Path('/data/nwp')

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp


def run_fortran_backend(pxldat, geo, nwp_interp):
    """Run Fortran native backend on full swath."""
    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
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
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    t_f90 = time.time() - t0
    return result, t_f90


def run_python_subset(pxldat, geo, nwp_interp, indices, thresholds):
    """Run Python backend on a subset of pixel indices."""
    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]

    # Build full arrays for swath call (need full arrays for 3x3 extraction)
    lsf_swath = np.zeros((n_elem, n_line), dtype=np.int32)  # will be overridden per-pixel
    sst_swath = np.zeros((n_elem, n_line), dtype=np.float64)

    cm_results = np.full(len(indices), -1, dtype=np.int32)
    conf_results = np.zeros(len(indices), dtype=np.float64)

    # Run pixel-by-pixel on subset
    from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_pixel, _extract_3x3

    t0 = time.time()
    for idx, (i, j) in enumerate(indices):
        eco_val = int(geo['eco_type'][i, j])
        lat_val = geo['lat'][i, j]
        lon_val = geo['lon'][i, j]
        sza_val = geo['sza'][i, j]
        vza_val = geo['vza'][i, j]
        elev_val = geo['elevation'][i, j]

        # Simple land/sea classification from eco_type
        # Matches Fortran: lsf=0 for water, lsf=2 for coast, lsf=1 for land
        if eco_val == 0:
            lsf_val = 0  # water
        elif eco_val == 14:
            lsf_val = 2  # coast
        else:
            lsf_val = 1  # land

        # Extract 3x3 neighborhoods
        indat_3x3_11um = _extract_3x3(pxldat[:, :, 23], i, j, n_elem, n_line)
        indat_3x3_vis = _extract_3x3(pxldat[:, :, 2], i, j, n_elem, n_line)

        result = run_cloud_mask_pixel(
            pxldat=pxldat[i, j, :],
            lat=lat_val,
            lon=lon_val,
            elevation=elev_val,
            lsf=lsf_val,
            sza=sza_val,
            vza=vza_val,
            glint_angle=geo['glint_angle'][i, j],
            eco_type=eco_val,
            snow_mask_val=0,
            sst=0.0,
            nwp_sfctmp=nwp_interp['tsfc'][i, j],
            nwp_pmsl=nwp_interp['pmsl'][i, j],
            nwp_u_wind=nwp_interp['u_wind'][i, j],
            nwp_v_wind=nwp_interp['v_wind'][i, j],
            nwp_precip_water=nwp_interp['tpw'][i, j],
            sensor_id=21,  # FY-3D
            bt_clr=np.zeros(7, dtype=np.float64),
            thresholds=thresholds,
            indat_3x3_11um=indat_3x3_11um,
            indat_3x3_vis=indat_3x3_vis,
        )
        cm_results[idx] = result.cloud_mask
        conf_results[idx] = result.confidence

    t_py = time.time() - t0
    return cm_results, conf_results, t_py


def main():
    parser = argparse.ArgumentParser(description="Compare Python vs Fortran cloud mask")
    parser.add_argument("--date", default="20220803", help="Date (YYYYMMDD)")
    parser.add_argument("--orbit", default="0740", help="Orbit time tag (HHMM)")
    parser.add_argument("--sample", type=int, default=50000, help="Number of pixels to sample")
    args = parser.parse_args()

    date_str = args.date
    time_tag = args.orbit

    # Find files
    mersi_dir = MERSI_ROOT / date_str
    l1b_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF"
    geo_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_GEO1K_MS.HDF"

    if not l1b_path.exists() or not geo_path.exists():
        logger.error("Input files not found")
        sys.exit(1)

    # Find NWP
    obs_hour = int(time_tag[:2])
    nwp_hours = [0, 3, 6, 9, 12, 15, 18, 21]
    best_hh = max(h for h in nwp_hours if h <= obs_hour) if obs_hour >= min(nwp_hours) else max(nwp_hours)
    nwp_path = NWP_ROOT / date_str / 'ORG' / f'gfs0p25_41L_{date_str}_{best_hh:02d}_00'
    if not nwp_path.exists():
        logger.error(f"NWP not found: {nwp_path}")
        sys.exit(1)

    logger.info(f"L1b: {l1b_path}")
    logger.info(f"GEO: {geo_path}")
    logger.info(f"NWP: {nwp_path}")

    # Read data
    logger.info("Reading L1b...")
    pxldat = read_l1b_data(str(l1b_path))
    logger.info("Reading GEO...")
    geo = read_geo_data(str(geo_path))
    logger.info("Reading NWP...")
    nwp = read_nwp_binary(str(nwp_path))
    logger.info("Interpolating NWP...")
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    logger.info(f"Swath size: {n_elem} x {n_line} = {n_elem*n_line:,} pixels")

    if not is_native_available():
        logger.error("Native backend not available!")
        sys.exit(1)

    # Run Fortran on full swath
    logger.info("Running Fortran backend (full swath)...")
    result_f90, t_f90 = run_fortran_backend(pxldat, geo, nwp_interp)
    logger.info(f"  Fortran: {t_f90:.1f}s")

    cm_f90 = result_f90['cloud_mask']
    conf_f90 = result_f90['confidence']

    # Sample pixels for Python comparison
    np.random.seed(42)
    n_sample = min(args.sample, n_elem * n_line)
    flat_indices = np.random.choice(n_elem * n_line, size=n_sample, replace=False)
    indices_2d = [(idx // n_line, idx % n_line) for idx in flat_indices]

    # Load thresholds
    import yaml
    thresholds_path = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
    with open(thresholds_path) as f:
        thresholds = yaml.safe_load(f)

    logger.info(f"Running Python backend ({n_sample:,} pixels)...")
    cm_py, conf_py, t_py = run_python_subset(pxldat, geo, nwp_interp, indices_2d, thresholds)
    logger.info(f"  Python: {t_py:.1f}s ({n_sample/t_py:.0f} pix/s)")

    # Extract Fortran results at sampled locations
    cm_f90_sampled = np.array([cm_f90[i, j] for i, j in indices_2d])
    conf_f90_sampled = np.array([conf_f90[i, j] for i, j in indices_2d])

    # Compare
    valid = (cm_py >= 0) & (cm_py <= 3) & (cm_f90_sampled >= 0) & (cm_f90_sampled <= 3)
    n_valid = int(np.sum(valid))

    if n_valid == 0:
        logger.error("No valid pixels for comparison!")
        sys.exit(1)

    agree = int(np.sum(cm_py[valid] == cm_f90_sampled[valid]))
    agree_rate = agree / n_valid

    conf_corr = float(np.corrcoef(conf_py[valid], conf_f90_sampled[valid])[0, 1]) if n_valid > 1 else 0.0
    conf_diff = np.abs(conf_py[valid] - conf_f90_sampled[valid])

    cat_labels = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}
    py_dist = {c: int(np.sum(cm_py[valid] == c)) for c in range(4)}
    f90_dist = {c: int(np.sum(cm_f90_sampled[valid] == c)) for c in range(4)}

    # Confusion matrix
    confusion = {}
    for py_c in range(4):
        for f90_c in range(4):
            n = int(np.sum((cm_py[valid] == py_c) & (cm_f90_sampled[valid] == f90_c)))
            if n > 0:
                confusion[f"{cat_labels[py_c]}->{cat_labels[f90_c]}"] = n

    # Full swath Fortran distribution
    f90_valid = (cm_f90 >= 0) & (cm_f90 <= 3)
    f90_full_dist = {c: int(np.sum(cm_f90[f90_valid] == c)) for c in range(4)}
    f90_full_total = int(np.sum(f90_valid))

    # Print results
    print(f"\n{'='*70}")
    print(f"Python vs Fortran Comparison: {date_str} orbit {time_tag}")
    print(f"{'='*70}")
    print(f"  Swath:            {n_elem} x {n_line} = {n_elem*n_line:,} pixels")
    print(f"  Sampled:          {n_sample:,} pixels")
    print(f"  Valid compared:   {n_valid:,}")
    print(f"  Agreement:        {agree:,} / {n_valid:,} = {agree_rate:.4f} ({agree_rate*100:.1f}%)")
    print(f"  Conf correlation: {conf_corr:.4f}")
    print(f"  Conf diff (mean): {np.mean(conf_diff):.4f}")
    print(f"  Conf diff (std):  {np.std(conf_diff):.4f}")

    print(f"\n  Python distribution (sampled):")
    for c in range(4):
        v = py_dist[c]
        pct = 100.0 * v / n_valid if n_valid > 0 else 0
        print(f"    {cat_labels[c]:<20} {v:>10,} ({pct:>5.1f}%)")

    print(f"\n  Fortran distribution (sampled):")
    for c in range(4):
        v = f90_dist[c]
        pct = 100.0 * v / n_valid if n_valid > 0 else 0
        print(f"    {cat_labels[c]:<20} {v:>10,} ({pct:>5.1f}%)")

    print(f"\n  Fortran distribution (full swath):")
    for c in range(4):
        v = f90_full_dist[c]
        pct = 100.0 * v / f90_full_total if f90_full_total > 0 else 0
        print(f"    {cat_labels[c]:<20} {v:>10,} ({pct:>5.1f}%)")

    print(f"\n  Confusion matrix (Python -> Fortran):")
    for k, v in sorted(confusion.items()):
        print(f"    {k:<40} {v:>10,}")

    print(f"\n  Timing:")
    print(f"    Python (sampled):  {t_py:.1f}s ({n_sample/t_py:.0f} pix/s)")
    print(f"    Fortran (full):    {t_f90:.1f}s ({n_elem*n_line/t_f90/1e6:.1f} Mpix/s)")
    est_py_full = t_py * (n_elem * n_line) / n_sample
    print(f"    Est. Python full:  {est_py_full:.0f}s ({est_py_full/t_f90:.0f}x slower)")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
