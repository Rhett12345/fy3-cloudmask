#!/usr/bin/env python
"""Analyze cloud mask output: statistics, spatial distribution, and RGB correlation."""

import sys
import numpy as np
import h5py
from pathlib import Path

H5_PATH = '/data/Data_yuq/fy3_cloud/20220803/FY3D_MERSI_20220803_0740_CLM_CLA.h5'

with h5py.File(H5_PATH, 'r') as f:
    print("=== HDF5 datasets ===")
    for key in f.keys():
        ds = f[key]
        print(f"  {key}: shape={ds.shape}, dtype={ds.dtype}")

    cm = f['cloud_mask'][:]
    conf = f['confidence'][:]
    qa = f['qa_bits'][:]
    testbits = f['testbits'][:]
    cla = f['cloud_amount'][:]

n_total = cm.size
print(f"\n{'='*60}")
print(f"Cloud Mask Statistics (total {n_total:,} pixels)")
print(f"{'='*60}")

labels = {0: 'Cloudy', 1: 'Prob Cloudy', 2: 'Prob Clear', 3: 'Confident Clear'}
for val, name in labels.items():
    count = np.sum(cm == val)
    pct = 100 * count / n_total
    mask = cm == val
    if mask.any():
        mean_conf = conf[mask].mean()
        print(f"  {name:20s}: {count:>10,} ({pct:6.2f}%)  mean_conf={mean_conf:.4f}")
    else:
        print(f"  {name:20s}: {count:>10,} ({pct:6.2f}%)")

print(f"\n  Mean confidence: {conf.mean():.4f}")
print(f"  Std confidence:  {conf.std():.4f}")

# Confidence distribution by CM class
print(f"\n{'='*60}")
print("Confidence distribution by CM class")
print(f"{'='*60}")
for val, name in labels.items():
    mask = cm == val
    if mask.any():
        c = conf[mask]
        print(f"  {name:20s}: min={c.min():.4f}  p25={np.percentile(c,25):.4f}  "
              f"median={np.median(c):.4f}  p75={np.percentile(c,75):.4f}  max={c.max():.4f}")

# Confidence histogram
print(f"\n{'='*60}")
print("Confidence histogram (10 bins)")
print(f"{'='*60}")
bins = np.linspace(0, 1, 11)
hist, edges = np.histogram(conf, bins=bins)
for i in range(len(hist)):
    print(f"  [{edges[i]:.2f}, {edges[i+1]:.2f}): {hist[i]:>10,} ({100*hist[i]/n_total:5.2f}%)")

# Spatial distribution: check if cloud/clear are geographically coherent
print(f"\n{'='*60}")
print("Spatial coherence check (cross-track profile)")
print(f"{'='*60}")
# cm shape is (n_elem, n_line) = (2048, 2000)
n_elem, n_line = cm.shape
for line_idx in [0, n_line//4, n_line//2, 3*n_line//4, n_line-1]:
    row = cm[:, line_idx]
    cloudy_pct = 100 * np.sum(row == 0) / n_elem
    clear_pct = 100 * np.sum(row == 3) / n_elem
    print(f"  Line {line_idx:4d}: Cloudy={cloudy_pct:5.1f}%  Clear={clear_pct:5.1f}%")

# Along-track profile
print(f"\n  Along-track profile (column center):")
col_idx = n_elem // 2
col = cm[col_idx, :]
cloudy_pct = 100 * np.sum(col == 0) / n_line
clear_pct = 100 * np.sum(col == 3) / n_line
print(f"  Col {col_idx}: Cloudy={cloudy_pct:5.1f}%  Clear={clear_pct:5.1f}%")

# Check testbits utilization
print(f"\n{'='*60}")
print("Testbits utilization (which tests fired)")
print(f"{'='*60}")
# testbits shape: (n_elem, n_line, 6)
if testbits.ndim == 3:
    tb = testbits.reshape(-1, 6)
else:
    tb = testbits.reshape(-1, 6)
# Count how many pixels have each bit set
for byte_idx in range(6):
    for bit in range(8):
        bit_pos = byte_idx * 8 + bit
        count = np.sum((tb[:, byte_idx] >> bit) & 1)
        if count > 0:
            pct = 100 * count / n_total
            print(f"  Bit {bit_pos:2d} (byte {byte_idx}, bit {bit}): {count:>10,} ({pct:5.2f}%)")

# Check for suspicious patterns
print(f"\n{'='*60}")
print("Sanity checks")
print(f"{'='*60}")

# 1. All pixels same class?
for val, name in labels.items():
    count = np.sum(cm == val)
    if count == n_total:
        print(f"  WARNING: ALL pixels are {name}! Likely algorithm failure.")
    elif count == 0:
        print(f"  WARNING: ZERO pixels are {name}! Check thresholds.")

# 2. Confidence vs CM consistency
for val, name in labels.items():
    mask = cm == val
    if mask.any():
        mean_c = conf[mask].mean()
        if val == 0 and mean_c > 0.5:
            print(f"  WARNING: Cloudy pixels have high mean confidence ({mean_c:.4f})")
        if val == 3 and mean_c < 0.5:
            print(f"  WARNING: Clear pixels have low mean confidence ({mean_c:.4f})")

# 3. Cloud amount statistics
print(f"\n  Cloud amount (5km): shape={cla.shape}")
valid_cla = cla[cla >= 0]
if valid_cla.size > 0:
    print(f"  Range: [{valid_cla.min():.2f}, {valid_cla.max():.2f}]")
    print(f"  Mean: {valid_cla.mean():.2f}")
    print(f"  Std:  {valid_cla.std():.2f}")

print("\nDone.")
