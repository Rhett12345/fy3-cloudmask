#!/usr/bin/env python3
"""Compare Python vs Fortran on a strip of lines with full spatial context.

Processes lines 800-1000 (200 lines x 2048 pixels = ~410K pixels) with both
backends and compares pixel-by-pixel.

Usage:
    python scripts/compare_strip.py
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_swath

import yaml
thresholds_path = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
with open(thresholds_path) as f:
    thresholds = yaml.safe_load(f)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp


def main():
    l1b = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
    geo_p = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF'
    nwp_p = '/data/nwp/20220803/ORG/gfs0p25_41L_20220803_06_00'

    print("Reading data...")
    pxldat = read_l1b_data(l1b)
    geo = read_geo_data(geo_p)
    nwp = read_nwp_binary(nwp_p)
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    print(f"Full swath: {n_elem} x {n_line}")

    # Select a strip of lines
    j_start, j_end = 800, 1000
    n_strip = j_end - j_start
    print(f"Strip: lines {j_start}-{j_end} ({n_strip} lines, {n_elem*n_strip:,} pixels)")

    # --- Fortran on full swath ---
    print("\nRunning Fortran (full swath)...")
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
        lsf=np.ascontiguousarray(geo['lsf'].astype(np.int8)),
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    t_f90 = time.time() - t0
    print(f"  Fortran: {t_f90:.1f}s")

    cm_f90_full = result_f90['cloud_mask']
    conf_f90_full = result_f90['confidence']

    # Extract strip from Fortran results
    cm_f90 = cm_f90_full[:, j_start:j_end]
    conf_f90 = conf_f90_full[:, j_start:j_end]

    # --- Python on strip ---
    # Prepare strip data
    pxldat_strip = pxldat[:, j_start:j_end, :]
    lat_strip = geo['lat'][:, j_start:j_end]
    lon_strip = geo['lon'][:, j_start:j_end]
    elev_strip = geo['elevation'][:, j_start:j_end]
    eco_strip = geo['eco_type'][:, j_start:j_end]
    sza_strip = geo['sza'][:, j_start:j_end]
    vza_strip = geo['vza'][:, j_start:j_end]
    glint_strip = geo['glint_angle'][:, j_start:j_end]
    sfctmp_strip = nwp_interp['tsfc'][:, j_start:j_end]
    pmsl_strip = nwp_interp['pmsl'][:, j_start:j_end]
    uwind_strip = nwp_interp['u_wind'][:, j_start:j_end]
    vwind_strip = nwp_interp['v_wind'][:, j_start:j_end]
    tpw_strip = nwp_interp['tpw'][:, j_start:j_end]

    # Derive lsf from eco_type (same as Fortran GEO reader)
    # eco_type: 0=water, others=land
    lsf_strip = np.where(eco_strip == 0, 0, 1).astype(np.int32)

    print(f"\nRunning Python (strip: {n_elem}x{n_strip}={n_elem*n_strip:,} pixels)...")
    t0 = time.time()
    cm_bit, qa_bit, cm_py, conf_py = run_cloud_mask_swath(
        pxldat_swath=pxldat_strip.astype(np.float64),
        lat_swath=lat_strip,
        lon_swath=lon_strip,
        elevation_swath=elev_strip,
        lsf_swath=lsf_strip,
        sza_swath=sza_strip,
        vza_swath=vza_strip,
        glint_angle_swath=glint_strip,
        eco_type_swath=eco_strip.astype(np.int32),
        snow_mask_swath=np.zeros((n_elem, n_strip), dtype=np.int32),
        sst_swath=np.zeros((n_elem, n_strip), dtype=np.float64),
        nwp_sfctmp_swath=sfctmp_strip,
        nwp_pmsl_swath=pmsl_strip,
        nwp_u_wind_swath=uwind_strip,
        nwp_v_wind_swath=vwind_strip,
        nwp_precip_water_swath=tpw_strip,
        bt_clr_swath=np.zeros((n_elem, n_strip, 7), dtype=np.float64),
        sensor_id=21,
        thresholds=thresholds,
    )
    t_py = time.time() - t0
    print(f"  Python: {t_py:.1f}s ({n_elem*n_strip/t_py:.0f} pix/s)")

    # --- Compare ---
    valid = (cm_py >= 0) & (cm_py <= 3) & (cm_f90 >= 0) & (cm_f90 <= 3)
    n_valid = int(np.sum(valid))
    n_total = cm_py.size

    agree = int(np.sum(cm_py[valid] == cm_f90[valid]))
    agree_rate = agree / n_valid if n_valid > 0 else 0

    conf_corr = float(np.corrcoef(conf_py[valid], conf_f90[valid])[0, 1]) if n_valid > 1 else 0.0
    conf_diff = np.abs(conf_py[valid] - conf_f90[valid])

    cat_labels = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}
    py_dist = {c: int(np.sum(cm_py[valid] == c)) for c in range(4)}
    f90_dist = {c: int(np.sum(cm_f90[valid] == c)) for c in range(4)}

    confusion = {}
    for py_c in range(4):
        for f90_c in range(4):
            n = int(np.sum((cm_py[valid] == py_c) & (cm_f90[valid] == f90_c)))
            if n > 0:
                confusion[f"{cat_labels[py_c]}->{cat_labels[f90_c]}"] = n

    print(f"\n{'='*70}")
    print(f"Python vs Fortran Comparison (strip lines {j_start}-{j_end})")
    print(f"{'='*70}")
    print(f"  Total pixels:     {n_total:,}")
    print(f"  Valid pixels:     {n_valid:,}")
    print(f"  Agreement:        {agree:,} / {n_valid:,} = {agree_rate:.4f} ({agree_rate*100:.1f}%)")
    print(f"  Conf correlation: {conf_corr:.4f}")
    print(f"  Conf diff (mean): {np.mean(conf_diff):.4f}")

    print(f"\n  Python distribution:")
    for c in range(4):
        v = py_dist[c]
        pct = 100.0 * v / n_valid if n_valid > 0 else 0
        print(f"    {cat_labels[c]:<20} {v:>10,} ({pct:>5.1f}%)")

    print(f"\n  Fortran distribution:")
    for c in range(4):
        v = f90_dist[c]
        pct = 100.0 * v / n_valid if n_valid > 0 else 0
        print(f"    {cat_labels[c]:<20} {v:>10,} ({pct:>5.1f}%)")

    print(f"\n  Confusion matrix (Python -> Fortran):")
    for k, v in sorted(confusion.items()):
        print(f"    {k:<40} {v:>10,}")

    print(f"\n  Timing:")
    print(f"    Python (strip):  {t_py:.1f}s")
    print(f"    Fortran (full):  {t_f90:.1f}s")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
