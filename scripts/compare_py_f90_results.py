#!/usr/bin/env python3
"""Compare Python and Fortran cloud mask results."""
import numpy as np
import h5py

# Load results
with h5py.File('/tmp/python_ref_0740.h5', 'r') as f:
    cm_py = np.array(f['cm'])
    conf_py = np.array(f['conf'])

with h5py.File('/tmp/fresh_f90_0740.h5', 'r') as f:
    cm_f90 = np.array(f['cm'])
    conf_f90 = np.array(f['conf'])

CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}

print("=" * 70)
print("Python vs Fortran Comparison")
print("=" * 70)

# Overall agreement
agree = np.sum(cm_py == cm_f90)
total = cm_py.size
print(f"\nOverall agreement: {agree:,} / {total:,} = {100*agree/total:.1f}%")

# Confusion matrix
print("\nConfusion Matrix (Python -> Fortran):")
print(f"{'':20s} | {'F90 cloudy':>12s} | {'F90 prob_cloudy':>15s} | {'F90 prob_clear':>14s} | {'F90 confident':>13s}")
print("-" * 80)
for py_val in range(4):
    row = []
    for f90_val in range(4):
        n = np.sum((cm_py == py_val) & (cm_f90 == f90_val))
        row.append(n)
    print(f"{'Python ' + CM_LABELS[py_val]:20s} | {row[0]:>12,} | {row[1]:>15,} | {row[2]:>14,} | {row[3]:>13,}")

# Distribution comparison
print("\nDistribution Comparison:")
print(f"{'Category':20s} | {'Python':>10s} | {'Fortran':>10s} | {'Diff':>10s}")
print("-" * 55)
for c in range(4):
    py_n = np.sum(cm_py == c)
    f90_n = np.sum(cm_f90 == c)
    print(f"{CM_LABELS[c]:20s} | {py_n:>10,} | {f90_n:>10,} | {f90_n-py_n:>+10,}")

# Confidence comparison
print("\nConfidence Statistics:")
print(f"{'':20s} | {'Python':>10s} | {'Fortran':>10s}")
print("-" * 45)
print(f"{'Mean':20s} | {np.mean(conf_py):>10.4f} | {np.mean(conf_f90):>10.4f}")
print(f"{'Std':20s} | {np.std(conf_py):>10.4f} | {np.std(conf_f90):>10.4f}")
print(f"{'Min':20s} | {np.min(conf_py):>10.4f} | {np.min(conf_f90):>10.4f}")
print(f"{'Max':20s} | {np.max(conf_py):>10.4f} | {np.max(conf_f90):>10.4f}")

# Where do they disagree most?
disagree = cm_py != cm_f90
print(f"\nDisagreeing pixels: {np.sum(disagree):,} ({100*np.sum(disagree)/total:.1f}%)")

# Check specific transitions
print("\nTop transitions (Python -> Fortran):")
transitions = {}
for py_val in range(4):
    for f90_val in range(4):
        if py_val == f90_val:
            continue
        n = np.sum((cm_py == py_val) & (cm_f90 == f90_val))
        if n > 0:
            transitions[f"{CM_LABELS[py_val]}->{CM_LABELS[f90_val]}"] = n

for k, v in sorted(transitions.items(), key=lambda x: -x[1])[:10]:
    print(f"  {k:30s}: {v:>10,}")
