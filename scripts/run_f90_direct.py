#!/usr/bin/env python3
"""Run Fortran backend directly on 0740 orbit and save results."""
import os, sys, time
from pathlib import Path
import h5py
import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

l1b = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
geo_p = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF'
nwp_p = '/data/nwp/20220803/ORG/gfs0p25_41L_20220803_06_00'

print("Reading data...")
pxldat = read_l1b_data(l1b)
geo = read_geo_data(geo_p)
nwp = read_nwp_binary(nwp_p)
nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
print(f"Swath: {n_elem} x {n_line} = {n_elem*n_line:,} pixels")

if not is_native_available():
    print("ERROR: Native backend not available!")
    sys.exit(1)

print("\nRunning Fortran backend...")
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
print(f"  Fortran: {t_f90:.1f}s ({n_elem*n_line/t_f90/1e6:.1f} Mpix/s)")

cm = result['cloud_mask']
conf = result['confidence']

# Replace NaN values with 0.0 (cloudy)
nan_mask = np.isnan(conf)
if np.sum(nan_mask) > 0:
    print(f"\nWarning: {np.sum(nan_mask)} NaN confidence values replaced with 0.0")
    conf[nan_mask] = 0.0
    cm[nan_mask] = 0  # Set to cloudy

# Statistics
CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}
print("\nCloud Mask Distribution:")
for c in range(4):
    n = np.sum(cm == c)
    print(f"  {CM_LABELS[c]:20s}: {n:>10,} ({100*n/cm.size:.1f}%)")

print(f"\nConfidence Statistics:")
print(f"  Mean: {np.mean(conf):.4f}")
print(f"  Std:  {np.std(conf):.4f}")

# Save
out_path = '/tmp/fresh_f90_0740.h5'
print(f"\nSaving to {out_path}...")
with h5py.File(out_path, 'w') as f:
    f.create_dataset('cm', data=cm)
    f.create_dataset('conf', data=conf)
    f.create_dataset('lat', data=geo['lat'])
    f.create_dataset('lon', data=geo['lon'])
print("Done!")
