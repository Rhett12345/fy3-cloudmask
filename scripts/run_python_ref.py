#!/usr/bin/env python3
"""Run Python reference implementation and compare with Fortran."""
import os, sys, time
from pathlib import Path
import h5py
import numpy as np
import yaml

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_swath

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

# Load thresholds
thresholds_path = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
with open(thresholds_path) as f:
    thresholds = yaml.safe_load(f)

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

# Derive lsf from eco_type
eco = geo['eco_type']
lsf = np.where(eco == 0, 0, 1).astype(np.int32)

print("\nRunning Python reference implementation...")
t0 = time.time()
cm_bit, qa_bit, cm_py, conf_py = run_cloud_mask_swath(
    pxldat_swath=pxldat.astype(np.float64),
    lat_swath=geo['lat'],
    lon_swath=geo['lon'],
    elevation_swath=geo['elevation'],
    lsf_swath=lsf,
    sza_swath=geo['sza'],
    vza_swath=geo['vza'],
    glint_angle_swath=geo['glint_angle'],
    eco_type_swath=eco.astype(np.int32),
    snow_mask_swath=np.zeros((n_elem, n_line), dtype=np.int32),
    sst_swath=np.zeros((n_elem, n_line), dtype=np.float64),
    nwp_sfctmp_swath=nwp_interp['tsfc'],
    nwp_pmsl_swath=nwp_interp['pmsl'],
    nwp_u_wind_swath=nwp_interp['u_wind'],
    nwp_v_wind_swath=nwp_interp['v_wind'],
    nwp_precip_water_swath=nwp_interp['tpw'],
    bt_clr_swath=np.zeros((n_elem, n_line, 7), dtype=np.float64),
    sensor_id=21,
    thresholds=thresholds,
)
t_py = time.time() - t0
print(f"  Python: {t_py:.1f}s ({n_elem*n_line/t_py:.0f} pix/s)")

# Statistics
CM_LABELS = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}
print("\nPython Cloud Mask Distribution:")
for c in range(4):
    n = np.sum(cm_py == c)
    print(f"  {CM_LABELS[c]:20s}: {n:>10,} ({100*n/cm_py.size:.1f}%)")

print(f"\nPython Confidence Statistics:")
print(f"  Mean: {np.mean(conf_py):.4f}")
print(f"  Std:  {np.std(conf_py):.4f}")

# Save
out_path = '/tmp/python_ref_0740.h5'
print(f"\nSaving to {out_path}...")
with h5py.File(out_path, 'w') as f:
    f.create_dataset('cm', data=cm_py)
    f.create_dataset('conf', data=conf_py)
print("Done!")
