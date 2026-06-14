#!/usr/bin/env python3
"""Check confidence distribution and group test contributions."""
import os, sys
from pathlib import Path
import h5py
import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data

l1b_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
old_path = '/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_20220803_0740_CLM_CLA.h5'
new_path = '/tmp/fresh_f90_0740.h5'

print("Reading L1B data...")
pxldat = read_l1b_data(l1b_path)
n_elem, n_line = pxldat.shape[0], pxldat.shape[1]

# Load confidence
with h5py.File(old_path, 'r') as f:
    conf_old = np.array(f['Cloud_Mask_1km/Confidence'])
    print(f"Old confidence shape: {conf_old.shape}")

with h5py.File(new_path, 'r') as f:
    conf_new = np.array(f['conf'])
    print(f"New confidence shape: {conf_new.shape}")

# Transpose if needed
if conf_old.shape == (n_line, n_elem):
    conf_old = conf_old.T
if conf_new.shape == (n_line, n_elem):
    conf_new = conf_new.T

# Load cloud mask
with h5py.File(old_path, 'r') as f:
    cm_old = np.array(f['Cloud_Mask_1km/Cloud_Mask_Value'])

with h5py.File(new_path, 'r') as f:
    cm_new = np.array(f['cm'])

if cm_old.shape == (n_line, n_elem):
    cm_old = cm_old.T
if cm_new.shape == (n_line, n_elem):
    cm_new = cm_new.T

valid_bt = (pxldat[:,:,22] > 100) & (pxldat[:,:,22] < 350)

print("\n" + "="*70)
print("Confidence Distribution by Cloud Mask Category")
print("="*70)

CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}

for label, cm_data, conf_data in [("OLD", cm_old, conf_old), ("NEW", cm_new, conf_new)]:
    print(f"\n--- {label} ---")
    for cm_val in range(4):
        mask = (cm_data == cm_val) & valid_bt
        n = np.sum(mask)
        if n == 0:
            continue
        print(f"  {CM_LABELS[cm_val]:20s}: n={n:>10,}, conf mean={np.mean(conf_data[mask]):.4f}, "
              f"std={np.std(conf_data[mask]):.4f}, min={np.min(conf_data[mask]):.4f}, max={np.max(conf_data[mask]):.4f}")

print("\n" + "="*70)
print("Overall Confidence Statistics")
print("="*70)

for label, conf_data in [("OLD", conf_old), ("NEW", conf_new)]:
    valid = valid_bt
    print(f"\n--- {label} ---")
    print(f"  Mean: {np.mean(conf_data[valid]):.4f}")
    print(f"  Std:  {np.std(conf_data[valid]):.4f}")
    print(f"  Min:  {np.min(conf_data[valid]):.4f}")
    print(f"  Max:  {np.max(conf_data[valid]):.4f}")

    # Histogram
    bins = [0, 0.01, 0.1, 0.3, 0.5, 0.66, 0.95, 0.99, 1.01]
    hist, _ = np.histogram(conf_data[valid], bins=bins)
    print(f"\n  Confidence histogram:")
    for i in range(len(bins)-1):
        print(f"    [{bins[i]:.2f}, {bins[i+1]:.2f}): {hist[i]:>10,} ({100*hist[i]/np.sum(hist):.1f}%)")

# Check where confidence is very low
print("\n" + "="*70)
print("Low Confidence Pixels (conf < 0.01)")
print("="*70)

for label, cm_data, conf_data in [("OLD", cm_old, conf_old), ("NEW", cm_new, conf_new)]:
    low_conf = (conf_data < 0.01) & valid_bt
    n_low = np.sum(low_conf)
    print(f"\n--- {label}: {n_low:,} pixels with conf < 0.01 ---")
    if n_low > 0:
        for cm_val in range(4):
            n = np.sum((cm_data == cm_val) & low_conf)
            if n > 0:
                print(f"  {CM_LABELS[cm_val]:20s}: {n:>10,}")
