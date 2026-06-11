#!/usr/bin/env python3
"""Analyze FY-3D data distributions to understand threshold needs.

Reads L1b + GEO data and computes key test value distributions by surface type.
Output: histograms of BTD 11-12, BTD 11-4, 0.64um refl, GEMI, 1.38um refl.
"""

import os
import sys
import math
import numpy as np
from pathlib import Path

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_fortran_only import read_l1b_data, read_geo_data


def main():
    l1b_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
    geo_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF'

    print("Reading L1b data...")
    pxldat = read_l1b_data(l1b_path)  # (nElem, nLine, 25)

    print("Reading GEO data...")
    geo = read_geo_data(geo_path)

    # Read LandSeaMask from GEO file
    import h5py
    with h5py.File(geo_path, 'r') as f:
        lsf_raw = f['Geolocation/LandSeaMask'][:].astype(np.int32)
    geo['lsf'] = lsf_raw.T

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    print(f"Swath: {n_elem} x {n_line}")

    # Band indices (MERSI-II)
    # VIS: 0.64um=idx2, 0.86um=idx3, 1.38um=idx18 (0-based)
    # IR: 3.8um=idx19, 8.5um=idx22, 11um=idx23, 12um=idx24
    ref064 = pxldat[:, :, 2]   # 0.64um
    ref086 = pxldat[:, :, 3]   # 0.86um
    ref138 = pxldat[:, :, 18]  # 1.38um
    ir38 = pxldat[:, :, 19]    # 3.8um
    ir11 = pxldat[:, :, 23]    # 11um
    ir12 = pxldat[:, :, 24]    # 12um

    lat = geo['lat']
    lon = geo['lon']
    sza = geo['sza']
    lsf = geo['lsf']

    # Classify: day/night, land/water
    is_day = sza < 85.0
    is_land = (lsf == 1)
    is_water = (lsf == 0)
    is_coast = (lsf == 2)

    # Valid data mask
    vis_valid = (ref064 > -99.0) & (ref064 <= 2.3)
    ir_valid = (ir11 > 0.0) & (ir11 < 1000.0) & (ir12 > 0.0) & (ir12 < 1000.0)
    ir38_valid = (ir38 > 0.0) & (ir38 < 1000.0)

    # Compute test values
    btd_11_12 = np.where(ir_valid, ir11 - ir12, np.nan)
    btd_11_4 = np.where(ir_valid & ir38_valid, ir11 - ir38, np.nan)

    # GEMI ratio
    s1 = ref064 * 100.0
    s2 = ref086 * 100.0
    etan = 2.0 * (s2 - s1) + 1.5 * s2 + 0.5 * s1
    etad = s2 + s1 + 0.5
    with np.errstate(divide='ignore', invalid='ignore'):
        eta = np.where(etad > 0, etan / etad, 0.0)
        vrat = np.where(s1 < 100.0,
                        eta * (1.0 - 0.25 * eta) - ((s1 - 0.125) / (100.0 - s1)),
                        0.0)
    vrat = np.where(vis_valid & (ref086 > -99.0), vrat, np.nan)

    # Surface type masks
    land_day = is_land & is_day & vis_valid & ir_valid
    water_day = is_water & is_day & vis_valid & ir_valid
    land_nite = is_land & ~is_day & ir_valid

    print("\n" + "=" * 70)
    print("Data Distribution Analysis")
    print("=" * 70)

    def print_stats(name, data, thresholds=None):
        """Print statistics for a test value array."""
        valid = data[~np.isnan(data)]
        if len(valid) == 0:
            print(f"  {name}: no valid data")
            return
        print(f"  {name}:")
        print(f"    N={len(valid):,}  mean={np.mean(valid):.3f}  std={np.std(valid):.3f}")
        print(f"    min={np.min(valid):.3f}  p1={np.percentile(valid,1):.3f}  "
              f"p5={np.percentile(valid,5):.3f}  p10={np.percentile(valid,10):.3f}")
        print(f"    p25={np.percentile(valid,25):.3f}  p50={np.percentile(valid,50):.3f}  "
              f"p75={np.percentile(valid,75):.3f}")
        print(f"    p90={np.percentile(valid,90):.3f}  p95={np.percentile(valid,95):.3f}  "
              f"p99={np.percentile(valid,99):.3f}  max={np.max(valid):.3f}")
        if thresholds:
            mid = thresholds[1]
            lo_pct = 100.0 * np.sum(valid < thresholds[2]) / len(valid)
            mid_pct = 100.0 * np.sum((valid >= thresholds[2]) & (valid < thresholds[1])) / len(valid)
            hi_pct = 100.0 * np.sum((valid >= thresholds[1]) & (valid < thresholds[0])) / len(valid)
            above_pct = 100.0 * np.sum(valid >= thresholds[0]) / len(valid)
            print(f"    Threshold analysis (lo={thresholds[2]}, mid={thresholds[1]}, hi={thresholds[0]}):")
            print(f"      < hicut (clear, conf=1.0): {lo_pct:.1f}%")
            print(f"      hicut-mid (transition):    {mid_pct:.1f}%")
            print(f"      mid-locut (cloudy):        {hi_pct:.1f}%")
            print(f"      > locut (cloud, conf=0.0): {above_pct:.1f}%")

    # ================================================================
    # BTD 11-12 analysis
    # ================================================================
    print("\n--- BTD 11-12um (thin cirrus test) ---")
    print_stats("Land Day", btd_11_12[land_day])
    print_stats("Water Day", btd_11_12[water_day])

    # ================================================================
    # BTD 11-4 analysis
    # ================================================================
    print("\n--- BTD 11-4um (fog/low cloud test) ---")
    print_stats("Land Day", btd_11_4[land_day],
                thresholds=[-14.0, -12.0, -10.0])  # current Fortran thresholds
    print_stats("Water Day", btd_11_4[water_day],
                thresholds=[-10.0, -8.0, -6.0])

    # With proposed new thresholds
    print("\n  With PROPOSED thresholds [-34, -14, -4]:")
    print_stats("Land Day (proposed)", btd_11_4[land_day],
                thresholds=[-34.0, -14.0, -4.0])

    # ================================================================
    # 0.64um reflectance analysis
    # ================================================================
    print("\n--- 0.64um Reflectance ---")
    print_stats("Land Day", ref064[land_day],
                thresholds=[0.24, 0.20, 0.16])
    print_stats("Water Day", ref064[water_day])

    # ================================================================
    # GEMI ratio analysis
    # ================================================================
    print("\n--- GEMI Ratio (vrat) ---")
    print_stats("Land Day", vrat[land_day],
                thresholds=[1.80, 1.85, 1.90])

    # ================================================================
    # 1.38um reflectance analysis
    # ================================================================
    print("\n--- 1.38um Reflectance (thin cirrus NIR) ---")
    print_stats("Land Day", ref138[land_day],
                thresholds=[0.04, 0.035, 0.03])

    # ================================================================
    # Simulate confidence for land_day
    # ================================================================
    print("\n" + "=" * 70)
    print("Simulated Confidence Distribution (Land Day)")
    print("=" * 70)

    from fy3_cloudmask.algorithm.confidence import conf_test

    mask = land_day & ~np.isnan(btd_11_12)
    n = np.sum(mask)
    print(f"  Valid pixels: {n:,}")

    # BTD 11-12 confidence (using tview threshold ~ 2.0 as approximation)
    dfthrsh_approx = 2.0
    locut = dfthrsh_approx + 0.3 * dfthrsh_approx
    hicut = dfthrsh_approx - 1.25
    c_btd1112 = np.array([conf_test(v, locut, dfthrsh_approx, hicut, 1.0)
                          for v in btd_11_12[mask]])

    # 0.64um confidence
    ref064_thr = [0.24, 0.20, 0.16, 1.0]
    c_ref064 = np.array([conf_test(v, ref064_thr[0], ref064_thr[1], ref064_thr[2], ref064_thr[3])
                         for v in ref064[mask]])

    # GEMI confidence
    vrat_thr = [1.80, 1.85, 1.90, 1.0]
    vrat_masked = np.where(np.isnan(vrat[mask]), 0.0, vrat[mask])
    c_vrat = np.array([conf_test(v, vrat_thr[0], vrat_thr[1], vrat_thr[2], vrat_thr[3])
                       for v in vrat_masked])

    # 1.38um confidence
    ref138_thr = [0.04, 0.035, 0.03, 1.0]
    c_ref138 = np.array([conf_test(v, ref138_thr[0], ref138_thr[1], ref138_thr[2], ref138_thr[3])
                         for v in ref138[mask]])

    # Group minimums (matching Fortran LandDay logic)
    # Group 1 (PFMFT/NFMFT): commented out, cmin1=1.0 always
    cmin1 = np.ones(n)
    # Group 2 (BTD 11-12 only, BTD 11-4 commented out)
    cmin2 = c_btd1112
    # Group 3 (0.64um + GEMI)
    cmin3 = np.minimum(c_ref064, c_vrat)
    # Group 4 (1.38um)
    cmin4 = c_ref138

    # Geometric mean
    # Count active groups
    ng1 = np.ones(n)  # always "active" with cmin1=1.0
    ng2 = np.ones(n)  # BTD 11-12 always runs
    ng3 = np.ones(n)  # VIS always runs for land_day
    ng4 = np.ones(n)  # NIR always runs
    groups = 4.0

    pre_conf = cmin1 * cmin2 * cmin3 * cmin4
    conf_sim = np.power(pre_conf, 1.0 / groups)

    print(f"\n  Simulated confidence stats:")
    print(f"    mean={np.mean(conf_sim):.4f}  std={np.std(conf_sim):.4f}")
    print(f"    p1={np.percentile(conf_sim,1):.4f}  p5={np.percentile(conf_sim,5):.4f}  "
          f"p10={np.percentile(conf_sim,10):.4f}")
    print(f"    p25={np.percentile(conf_sim,25):.4f}  p50={np.percentile(conf_sim,50):.4f}  "
          f"p75={np.percentile(conf_sim,75):.4f}")
    print(f"    p90={np.percentile(conf_sim,90):.4f}  p95={np.percentile(conf_sim,95):.4f}  "
          f"p99={np.percentile(conf_sim,99):.4f}")

    print(f"\n  Cloud mask distribution:")
    print(f"    Cloudy (conf<=0.66):      {100*np.sum(conf_sim<=0.66)/n:.2f}%")
    print(f"    Prob Cloudy (0.66-0.95):  {100*np.sum((conf_sim>0.66)&(conf_sim<=0.95))/n:.2f}%")
    print(f"    Prob Clear (0.95-0.99):   {100*np.sum((conf_sim>0.95)&(conf_sim<=0.99))/n:.2f}%")
    print(f"    Confident Clear (>0.99):  {100*np.sum(conf_sim>0.99)/n:.2f}%")

    # Show per-group contribution
    print(f"\n  Per-group mean confidence:")
    print(f"    Group 1 (PFMFT/NFMFT, disabled): mean={np.mean(cmin1):.4f}")
    print(f"    Group 2 (BTD 11-12):             mean={np.mean(cmin2):.4f}")
    print(f"    Group 3 (VIS 0.64+GEMI):         mean={np.mean(cmin3):.4f}")
    print(f"    Group 4 (NIR 1.38):              mean={np.mean(cmin4):.4f}")

    # Show how many pixels have each group at 0.0
    print(f"\n  Pixels with group conf=0.0:")
    print(f"    Group 1: {100*np.sum(cmin1==0)/n:.2f}%")
    print(f"    Group 2: {100*np.sum(cmin2==0)/n:.2f}%")
    print(f"    Group 3: {100*np.sum(cmin3==0)/n:.2f}%")
    print(f"    Group 4: {100*np.sum(cmin4==0)/n:.2f}%")
    print(f"    Any group=0: {100*np.sum((cmin1==0)|(cmin2==0)|(cmin3==0)|(cmin4==0))/n:.2f}%")

    # Confidence histogram
    print(f"\n  Confidence histogram:")
    bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.66, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]
    for i in range(len(bins)-1):
        count = np.sum((conf_sim >= bins[i]) & (conf_sim < bins[i+1]))
        if count > 0:
            print(f"    [{bins[i]:.2f}, {bins[i+1]:.2f}): {count:>10,} ({100*count/n:5.2f}%)")


if __name__ == '__main__':
    main()
