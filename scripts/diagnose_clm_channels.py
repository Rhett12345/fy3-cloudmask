#!/usr/bin/env python3
"""Diagnose CLM results against channel reflectance/brightness temperature.

Checks physical consistency:
- Cloudy pixels should have higher VIS reflectance (visible channels)
- Cloudy pixels should have lower 11um BT (cold tops)
- Clear sky BT should be close to surface temperature
- BTD (11-12um) patterns should match CLM

Usage:
    python scripts/diagnose_clm_channels.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import h5py
import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

# MERSI-II channel indices (0-based in pxldat array)
# VIS channels: 0-18 (ref)
# IR channels: 19-24 (tbb)
# Key channels:
#   CH1 (0.47um) = idx 0
#   CH2 (0.55um) = idx 1
#   CH3 (0.65um) = idx 2
#   CH5 (0.86um) = idx 4
#   CH6 (1.38um) = idx 5
#   CH20 (3.7um) = idx 19
#   CH22 (8.6um) = idx 21
#   CH23 (11um)  = idx 22 (actually 22 in 0-based, channel 23)
#   CH24 (12um)  = idx 23
#   CH25 (10.8um) = idx 24

# From constants.py: IR_OFFSET = 19
# IR bands in pxldat[:,:,19:] correspond to bands 20-25
# So: pxldat[:,:,19] = band20 (3.7um), pxldat[:,:,20] = band21 (4.0um)
#     pxldat[:,:,21] = band22 (8.6um), pxldat[:,:,22] = band23 (11um)
#     pxldat[:,:,23] = band24 (12um), pxldat[:,:,24] = band25 (10.8um)

BT_37_IDX = 19   # 3.7um
BT_86_IDX = 21   # 8.6um
BT_11_IDX = 22   # 11um
BT_12_IDX = 23   # 12um
BT_108_IDX = 24   # 10.8um
REF_065_IDX = 2   # 0.65um
REF_086_IDX = 4   # 0.86um
REF_138_IDX = 5   # 1.38um
REF_213_IDX = 16  # 2.13um

CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}


def load_clm(h5_path):
    """Load cloud mask from HDF5 file."""
    with h5py.File(h5_path, 'r') as f:
        # Try different key patterns
        for key in ['CloudMask', 'cloud_mask', 'CLM', 'clm', 'cloud_mask_result',
                    'Cloud_Mask_1km/Cloud_Mask_Value', 'cm']:
            if key in f:
                return np.array(f[key])
        # Print available keys
        print(f"Available keys in {h5_path}:")
        f.visit(lambda k: print(f"  {k}: {f[k].shape if hasattr(f[k], 'shape') else 'group'}"))
        return None


def main():
    l1b_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
    geo_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF'
    nwp_path = '/data/nwp/20220803/ORG/gfs0p25_41L_20220803_06_00'

    old_clm_path = '/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_20220803_0740_CLM_CLA.h5'
    new_clm_path = '/tmp/fresh_f90_0740.h5'

    print("Reading L1B data...")
    pxldat = read_l1b_data(l1b_path)
    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    print(f"  Swath: {n_elem} x {n_line}")

    print("Reading GEO data...")
    geo = read_geo_data(geo_path)

    print("Reading NWP data...")
    nwp = read_nwp_binary(nwp_path)
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    # Extract key channels
    bt_11um = pxldat[:, :, BT_11_IDX]  # 11um BT
    bt_12um = pxldat[:, :, BT_12_IDX]  # 12um BT
    bt_37um = pxldat[:, :, BT_37_IDX]  # 3.7um BT
    bt_86um = pxldat[:, :, BT_86_IDX]  # 8.6um BT
    ref_065 = pxldat[:, :, REF_065_IDX]  # 0.65um reflectance
    ref_086 = pxldat[:, :, REF_086_IDX]  # 0.86um reflectance
    ref_138 = pxldat[:, :, REF_138_IDX]  # 1.38um reflectance

    btd_11_12 = bt_11um - bt_12um  # 11-12um BTD
    btd_11_37 = bt_11um - bt_37um  # 11-3.7um BTD

    # Surface temperature from NWP
    tsfc = nwp_interp['tsfc']

    # Load CLM results
    print("\nLoading cloud mask results...")
    clm_old = load_clm(old_clm_path)
    clm_new = load_clm(new_clm_path)

    if clm_old is None or clm_new is None:
        print("ERROR: Could not load cloud mask data!")
        return

    # Ensure same shape
    print(f"  Old CLM shape: {clm_old.shape}")
    print(f"  New CLM shape: {clm_new.shape}")

    # Transpose if needed (HDF5 might be stored differently)
    if clm_old.shape == (n_line, n_elem):
        clm_old = clm_old.T
    if clm_new.shape == (n_line, n_elem):
        clm_new = clm_new.T

    # Flatten for statistics
    total = n_elem * n_line

    # Mask invalid data
    valid_bt = (bt_11um > 100) & (bt_11um < 350)
    valid_ref = (ref_065 >= 0) & (ref_065 <= 2.0)

    print("\n" + "=" * 70)
    print("Channel Statistics by Cloud Mask Category")
    print("=" * 70)

    for label, clm_data in [("OLD (5/26)", clm_old), ("NEW (after fix)", clm_new)]:
        print(f"\n--- {label} ---")

        for cm_val in range(4):
            mask = (clm_data == cm_val) & valid_bt
            n = np.sum(mask)
            if n == 0:
                continue

            pct = 100.0 * n / total
            print(f"\n  {CM_LABELS[cm_val]}: {n:,} pixels ({pct:.1f}%)")
            print(f"    11um BT:   mean={np.mean(bt_11um[mask]):.1f}K, std={np.std(bt_11um[mask]):.1f}K, "
                  f"range=[{np.min(bt_11um[mask]):.1f}, {np.max(bt_11um[mask]):.1f}]")
            print(f"    12um BT:   mean={np.mean(bt_12um[mask]):.1f}K, std={np.std(bt_12um[mask]):.1f}K")
            print(f"    BTD 11-12: mean={np.mean(btd_11_12[mask]):.2f}K, std={np.std(btd_11_12[mask]):.2f}K")
            print(f"    BTD 11-37: mean={np.mean(btd_11_37[mask]):.2f}K, std={np.std(btd_11_37[mask]):.2f}K")
            print(f"    BTD 11-86: mean={np.mean(bt_11um[mask]-bt_86um[mask]):.2f}K")

            mask_ref = mask & valid_ref
            if np.sum(mask_ref) > 0:
                print(f"    0.65um:    mean={np.mean(ref_065[mask_ref]):.4f}, std={np.std(ref_065[mask_ref]):.4f}")
                print(f"    0.86um:    mean={np.mean(ref_086[mask_ref]):.4f}, std={np.std(ref_086[mask_ref]):.4f}")
            if np.sum(mask & (ref_138 >= 0)) > 0:
                print(f"    1.38um:    mean={np.mean(ref_138[mask & (ref_138>=0)]):.4f}")

            # BT - Tsfc difference
            mask_sfc = mask & (tsfc > 100)
            if np.sum(mask_sfc) > 0:
                print(f"    BT11-Tsfc: mean={np.mean(bt_11um[mask_sfc]-tsfc[mask_sfc]):.1f}K")

    # Physical consistency checks
    print("\n" + "=" * 70)
    print("Physical Consistency Checks")
    print("=" * 70)

    for label, clm_data in [("OLD", clm_old), ("NEW", clm_new)]:
        print(f"\n--- {label} ---")

        # Check 1: Cloudy should have lower 11um BT than clear
        cloudy_bt = bt_11um[(clm_data == 0) & valid_bt]
        clear_bt = bt_11um[(clm_data == 3) & valid_bt]
        if len(cloudy_bt) > 0 and len(clear_bt) > 0:
            print(f"  11um BT: cloudy mean={np.mean(cloudy_bt):.1f}K, clear mean={np.mean(clear_bt):.1f}K, "
                  f"diff={np.mean(clear_bt)-np.mean(cloudy_bt):.1f}K {'OK' if np.mean(clear_bt) > np.mean(cloudy_bt) else 'WARNING: inverted!'}")

        # Check 2: Cloudy should have higher VIS reflectance
        cloudy_ref = ref_065[(clm_data == 0) & valid_ref]
        clear_ref = ref_065[(clm_data == 3) & valid_ref]
        if len(cloudy_ref) > 0 and len(clear_ref) > 0:
            print(f"  0.65um:  cloudy mean={np.mean(cloudy_ref):.4f}, clear mean={np.mean(clear_ref):.4f}, "
                  f"{'OK' if np.mean(cloudy_ref) > np.mean(clear_ref) else 'WARNING: inverted!'}")

        # Check 3: BTD 11-12 should be smaller for thin cirrus
        # (higher BTD = more likely clear in some tests)

        # Check 4: Night clear should have BT close to Tsfc
        sza = geo['sza']
        night_mask = sza > 90
        if np.sum(night_mask & (clm_data == 3) & valid_bt) > 0:
            night_clear_diff = bt_11um[night_mask & (clm_data == 3) & valid_bt] - tsfc[night_mask & (clm_data == 3) & valid_bt]
            print(f"  Night clear BT-Tsfc: mean={np.mean(night_clear_diff):.1f}K (should be ~0)")

    # Check agreement vs channel difference
    print("\n" + "=" * 70)
    print("Where OLD and NEW disagree - channel analysis")
    print("=" * 70)

    disagree = (clm_old != clm_new) & valid_bt
    n_disagree = np.sum(disagree)
    print(f"\nDisagreeing pixels: {n_disagree:,} ({100*n_disagree/total:.1f}%)")

    # For disagreeing pixels, which channels differ most?
    if n_disagree > 0:
        print(f"\n  Disagreeing pixels channel stats:")
        print(f"    11um BT:   mean={np.mean(bt_11um[disagree]):.1f}K")
        print(f"    BTD 11-12: mean={np.mean(btd_11_12[disagree]):.2f}K")
        print(f"    0.65um:    mean={np.mean(ref_065[disagree & valid_ref]):.4f}")

        # Transition analysis: what transitions happen?
        print(f"\n  Transition analysis:")
        for old_val in range(4):
            for new_val in range(4):
                if old_val == new_val:
                    continue
                trans_mask = (clm_old == old_val) & (clm_new == new_val) & valid_bt
                n_trans = np.sum(trans_mask)
                if n_trans > 100:
                    print(f"    {CM_LABELS[old_val]:20s} -> {CM_LABELS[new_val]:20s}: {n_trans:>10,} pixels")
                    print(f"      11um BT: {np.mean(bt_11um[trans_mask]):.1f}K, BTD 11-12: {np.mean(btd_11_12[trans_mask]):.2f}K")

    # NFMFT / PFMFT consistency check
    print("\n" + "=" * 70)
    print("NFMFT/PFMFT Test Value Distribution")
    print("=" * 70)

    # PFMFT: max(BT11, BT12) - BT37, should be high for clouds
    # NFMFT: BT11 - BT12, should be negative for some clouds
    pfmft = np.maximum(bt_11um, bt_12um) - bt_37um
    nfmft = bt_11um - bt_12um  # same as btd_11_12

    for label, clm_data in [("OLD", clm_old), ("NEW", clm_new)]:
        print(f"\n--- {label} ---")
        for cm_val in range(4):
            mask = (clm_data == cm_val) & valid_bt
            if np.sum(mask) == 0:
                continue
            print(f"  {CM_LABELS[cm_val]:20s}: PFMFT mean={np.mean(pfmft[mask]):.2f}K, NFMFT mean={np.mean(nfmft[mask]):.2f}K")


if __name__ == "__main__":
    main()
