#!/usr/bin/env python3
"""Process 2020-03-08 all 3 orbits with Python backend, operational + recalibration.

Outputs to 20200308_py/ (separate from Fortran 20200308/).
"""
import os, sys, time
from pathlib import Path

import h5py
import numpy as np
import yaml
import cfgrib
from scipy.interpolate import RegularGridInterpolator

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_swath
from fy3_cloudmask.io.recalibration import RecalibrationManager

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data

# ---- Paths ----
NWP_GRIB = '/data/nwp/20200308/ORG/fnl_20200308_12_00.grib2'
THRESHOLDS_PATH = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
RECAL_DIR = '/home/liusy2020/yuq/cloudmask/fy3d_recali'
OUT_DIR = '/data/Data_yuq/fy3_cloud/20200308_py'

ORBITS = [
    ('1345', '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_1000M_MS.HDF',
            '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_GEO1K_MS.HDF'),
]

os.makedirs(OUT_DIR, exist_ok=True)

# ---- Load thresholds ----
with open(THRESHOLDS_PATH) as f:
    thresholds = yaml.safe_load(f)

# ---- Load recalibration ----
recal_mgr = RecalibrationManager(RECAL_DIR)
cal0, cal1 = recal_mgr.load_coefficients('20200308')
print(f"Recalibration coeffs: cal0={cal0}, cal1={cal1}")

# ---- Load NWP once (same GRIB for all 3 orbits) ----
print("\nLoading NWP from GRIB2...")
t0 = time.time()
ds = cfgrib.open_dataset(NWP_GRIB, backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
grib_lat = ds.latitude.values
grib_lon = ds.longitude.values
tsfc_k = ds['t'].values

# Fill missing fields with zeros
for key in ['u', 'v', 'prmsl', 'pwat']:
    if key not in ds:
        ds[key] = (['latitude', 'longitude'], np.zeros_like(tsfc_k))

u_wind = ds['u'].values
v_wind = ds['v'].values
pmsl_hpa = ds['prmsl'].values / 100.0
tpw_kgm2 = ds['pwat'].values
print(f"  GRIB grid: lat {grib_lat.shape}, lon {grib_lon.shape}")
print(f"  NWP loaded in {time.time()-t0:.1f}s")


def interp_nwp(pixel_lat, pixel_lon):
    """Bilinearly interpolate NWP fields to pixel grid."""
    lon_360 = np.where(pixel_lon < 0, pixel_lon + 360, pixel_lon)

    if grib_lat[0] > grib_lat[-1]:
        lat_asc = grib_lat[::-1]
        tsfc_asc = tsfc_k[::-1, :]
        u_asc = u_wind[::-1, :]
        v_asc = v_wind[::-1, :]
        pmsl_asc = pmsl_hpa[::-1, :]
        tpw_asc = tpw_kgm2[::-1, :]
    else:
        lat_asc = grib_lat
        tsfc_asc = tsfc_k
        u_asc = u_wind
        v_asc = v_wind
        pmsl_asc = pmsl_hpa
        tpw_asc = tpw_kgm2

    def _interp(field_asc):
        interp = RegularGridInterpolator(
            (lat_asc, grib_lon), field_asc, bounds_error=False, fill_value=None)
        points = np.stack([pixel_lat.ravel(), lon_360.ravel()], axis=-1)
        return interp(points).reshape(pixel_lat.shape).astype(np.float64)

    return {
        'sfctmp': _interp(tsfc_asc),
        'u_wind': _interp(u_asc),
        'v_wind': _interp(v_asc),
        'pmsl': _interp(pmsl_asc),
        'tpw': _interp(tpw_asc),
    }


def process_orbit(name, l1b_path, geo_path, recal_cal0, recal_cal1, calib_label):
    """Process one orbit with given calibration and return cm, conf, lat, lon."""
    print(f"\n  [{calib_label}] Reading L1b/GEO...")
    t0 = time.time()
    pxldat = read_l1b_data(l1b_path, recal_cal0=recal_cal0, recal_cal1=recal_cal1)
    geo = read_geo_data(geo_path)
    print(f"    Read done in {time.time()-t0:.1f}s")

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    print(f"    Swath: {n_elem} x {n_line} = {n_elem*n_line:,} pixels")

    print(f"    Interpolating NWP...")
    t0 = time.time()
    nwp = interp_nwp(geo['lat'], geo['lon'])
    print(f"    NWP interp done in {time.time()-t0:.1f}s")

    lsf = geo['lsf'].astype(np.int32)

    print(f"    Running Python backend...")
    t0 = time.time()
    cm_bit, qa_bit, cm, conf = run_cloud_mask_swath(
        pxldat_swath=pxldat,
        lat_swath=geo['lat'],
        lon_swath=geo['lon'],
        elevation_swath=geo['elevation'],
        lsf_swath=lsf,
        sza_swath=geo['sza'],
        vza_swath=geo['vza'],
        glint_angle_swath=geo['glint_angle'],
        eco_type_swath=geo['eco_type'].astype(np.int32),
        snow_mask_swath=np.zeros((n_elem, n_line), dtype=np.int32),
        sst_swath=np.zeros((n_elem, n_line), dtype=np.float64),
        nwp_sfctmp_swath=nwp['sfctmp'],
        nwp_pmsl_swath=nwp['pmsl'],
        nwp_u_wind_swath=nwp['u_wind'],
        nwp_v_wind_swath=nwp['v_wind'],
        nwp_precip_water_swath=nwp['tpw'],
        bt_clr_swath=np.zeros((n_elem, n_line, 7), dtype=np.float64),
        sensor_id=21,
        thresholds=thresholds,
    )
    elapsed = time.time() - t0
    print(f"    Python backend: {elapsed:.1f}s ({n_elem*n_line/elapsed/1e6:.2f} Mpix/s)")

    # NaN cleanup
    nan_mask = np.isnan(conf)
    if np.sum(nan_mask) > 0:
        print(f"    Warning: {np.sum(nan_mask)} NaN confidence values -> 0")
        conf[nan_mask] = 0.0
        cm[nan_mask] = 0

    return cm, conf, geo['lat'], geo['lon']


# ---- Process all orbits ----
CM_LABELS = {0: 'cloudy', 1: 'prob_cld', 2: 'prob_clr', 3: 'clear'}

for name, l1b_path, geo_path in ORBITS:
    print(f"\n{'='*60}")
    print(f"Orbit {name}")
    print(f"{'='*60}")

    for calib_type, cal0_arr, cal1_arr, suffix in [
        ('operational', None, None, ''),
        ('recal', cal0, cal1, '_recal'),
    ]:
        cm, conf, lat, lon = process_orbit(
            name, l1b_path, geo_path, cal0_arr, cal1_arr, calib_type)

        # Stats
        print(f"    Cloud distribution:")
        for c in range(4):
            n = np.sum(cm == c)
            print(f"      {CM_LABELS[c]:>12s}: {n:>10,} ({100*n/cm.size:.1f}%)")
        print(f"    Confidence: mean={np.mean(conf):.4f}, std={np.std(conf):.4f}")

        # Save
        out_file = os.path.join(OUT_DIR, f'FY3D_MERSI_20200308_{name}_CLM_CLA{suffix}.h5')
        print(f"    Saving to {out_file}")
        with h5py.File(out_file, 'w') as f:
            f.create_dataset('cm', data=cm)
            f.create_dataset('conf', data=conf)
            f.create_dataset('lat', data=lat)
            f.create_dataset('lon', data=lon)
            f.attrs['backend'] = 'python'
            f.attrs['orbit'] = f'2020-03-08 {name[:2]}:{name[2:]}:00 UTC'
            f.attrs['calibration'] = calib_type

print(f"\n{'='*60}")
print(f"All orbits x 2 calibrations = {len(ORBITS)*2} files saved to 20200308_py/")
print(f"{'='*60}")
