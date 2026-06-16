#!/usr/bin/env python3
"""Debug confidence differences between Python and Fortran."""
import h5py
import numpy as np

with h5py.File('/tmp/python_ref_0740.h5', 'r') as f:
    cm_py = np.array(f['cm'])
    conf_py = np.array(f['conf'])

with h5py.File('/tmp/fresh_f90_0740.h5', 'r') as f:
    cm_f90 = np.array(f['cm'])
    conf_f90 = np.array(f['conf'])

CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}

print("Confidence distribution by cloud mask category:")
print(f"{'Category':20s} | {'Py conf mean':>12s} | {'F90 conf mean':>13s} | {'Count':>10s}")
print("-" * 60)

for c in range(4):
    py_mask = cm_py == c
    f90_mask = cm_f90 == c
    py_conf = conf_py[py_mask] if np.sum(py_mask) > 0 else [0]
    f90_conf = conf_f90[f90_mask] if np.sum(f90_mask) > 0 else [0]
    print(f"Python {CM_LABELS[c]:15s} | {np.mean(py_conf):>12.4f} | {np.mean(f90_conf):>13.4f} | {np.sum(py_mask):>10,}")

# Check key transitions
for py_c, f90_c, label in [(1, 0, 'prob_cloudy->cloudy'), (0, 3, 'cloudy->confident_clear'),
                            (1, 3, 'prob_cloudy->confident_clear')]:
    mask = (cm_py == py_c) & (cm_f90 == f90_c)
    if np.sum(mask) > 0:
        print(f"\n{label} ({np.sum(mask):,} pixels):")
        print(f"  Python conf: mean={np.mean(conf_py[mask]):.4f}, range=[{np.min(conf_py[mask]):.4f}, {np.max(conf_py[mask]):.4f}]")
        print(f"  Fortran conf: mean={np.mean(conf_f90[mask]):.4f}, range=[{np.min(conf_f90[mask]):.4f}, {np.max(conf_f90[mask]):.4f}]")

# Confidence histogram
for name, conf in [("Python", conf_py), ("Fortran", conf_f90)]:
    print(f"\n{name} confidence histogram:")
    hist, bins = np.histogram(conf, bins=[0, 0.01, 0.5, 0.66, 0.95, 0.99, 1.01])
    for i in range(len(hist)):
        print(f"  [{bins[i]:.2f}, {bins[i+1]:.2f}): {hist[i]:>10,} ({100*hist[i]/conf.size:.1f}%)")

# Where does Fortran give conf=0 (cloudy) but Python gives high conf?
mask = (conf_f90 < 0.01) & (conf_py > 0.66)
print(f"\nF90 conf<0.01 but Py conf>0.66: {np.sum(mask):,} pixels")
if np.sum(mask) > 0:
    print(f"  Python conf mean: {np.mean(conf_py[mask]):.4f}")
    print(f"  Python CM distribution: cloudy={np.sum(cm_py[mask]==0)}, prob_cloudy={np.sum(cm_py[mask]==1)}, prob_clear={np.sum(cm_py[mask]==2)}, confident_clear={np.sum(cm_py[mask]==3)}")

# Where does Fortran give conf=1 (confident_clear) but Python gives low conf?
mask = (conf_f90 > 0.99) & (conf_py < 0.66)
print(f"\nF90 conf>0.99 but Py conf<0.66: {np.sum(mask):,} pixels")
if np.sum(mask) > 0:
    print(f"  Python conf mean: {np.mean(conf_py[mask]):.4f}")
    print(f"  Python CM distribution: cloudy={np.sum(cm_py[mask]==0)}, prob_cloudy={np.sum(cm_py[mask]==1)}, prob_clear={np.sum(cm_py[mask]==2)}, confident_clear={np.sum(cm_py[mask]==3)}")
