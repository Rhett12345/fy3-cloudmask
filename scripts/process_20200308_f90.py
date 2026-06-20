#!/usr/bin/env python3
"""Process 2020-03-08 3 orbits with Fortran native backend (operational + recalibration).

Outputs CLM HDF5 to /data/Data_yuq/fy3_cloud/20200308_f90/
"""

import os, sys, time
from pathlib import Path

import h5py
import numpy as np
import cfgrib
from scipy.interpolate import RegularGridInterpolator

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
from fy3_cloudmask.io.recalibration import RecalibrationManager

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data

RECAL_DIR = '/home/liusy2020/yuq/cloudmask/fy3d_recali'
NWP_GRIB = '/data/nwp/20200308/ORG/fnl_20200308_12_00.grib2'
OUT_DIR = '/data/Data_yuq/fy3_cloud/20200308_f90'

ORBITS = [
    ('1345', '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_1000M_MS.HDF',
            '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_GEO1K_MS.HDF'),
    ('1435', '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1435_1000M_MS.HDF',
            '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1435_GEO1K_MS.HDF'),
    ('1525', '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1525_1000M_MS.HDF',
            '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1525_GEO1K_MS.HDF'),
]

os.makedirs(OUT_DIR, exist_ok=True)

if not is_native_available():
    print("ERROR: Native backend not available!")
    sys.exit(1)

# ---- Load recalibration ----
recal_mgr = RecalibrationManager(RECAL_DIR)
cal0, cal1 = recal_mgr.load_coefficients('20200308')
print(f"Recalibration coeffs: cal0={cal0}, cal1={cal1}")

# ---- Load NWP from GRIB2 ----
print(f"\nLoading NWP from {NWP_GRIB}...")
t0 = time.time()
ds = cfgrib.open_dataset(NWP_GRIB, backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
grib_lat = ds.latitude.values
grib_lon = ds.longitude.values
tsfc_k = ds['t'].values
for key in ['u', 'v', 'prmsl', 'pwat']:
    if key not in ds:
        ds[key] = (['latitude', 'longitude'], np.zeros_like(tsfc_k))
u_wind = ds['u'].values
v_wind = ds['v'].values
pmsl_hpa = ds['prmsl'].values / 100.0
tpw_kgm2 = ds['pwat'].values
print(f"  GRIB grid: lat {grib_lat.shape}, lon {grib_lon.shape}, loaded in {time.time()-t0:.1f}s")


def interp_nwp(pixel_lat, pixel_lon):
    """Bilinear interpolation of NWP fields to pixel grid."""
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


def run_orbit(name, l1b_path, geo_path, cal0_arr, cal1_arr, cal_type):
    out_file = os.path.join(OUT_DIR, f'FY3D_MERSI_20200308_{name}_CLM_CLA.h5')
    if 'recal' in cal_type:
        out_file = os.path.join(OUT_DIR, f'FY3D_MERSI_20200308_{name}_CLM_CLA_recal.h5')
    if os.path.exists(out_file):
        print(f"  [SKIP] {name} {cal_type} (exists)")
        return out_file

    print(f"  Reading L1b/GEO...")
    t0 = time.time()
    pxldat = read_l1b_data(l1b_path, recal_cal0=cal0_arr, recal_cal1=cal1_arr)
    geo = read_geo_data(geo_path)
    print(f"    Read done in {time.time()-t0:.1f}s")

    print(f"  Interpolating NWP...")
    t0 = time.time()
    nwp_int = interp_nwp(geo['lat'], geo['lon'])
    print(f"    Done in {time.time()-t0:.1f}s")

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]

    # Compute btclr (clear-sky BT estimate) from NWP surface temperature
    # btclr(5)=BT11_clr, btclr(6)=BT12_clr; other channels left as 0
    btclr_arr = np.zeros((n_elem, n_line, 7), dtype=np.float32)
    btclr_arr[:, :, 4] = np.maximum(nwp_int['sfctmp'], 275.0)  # BT11 clear > 270K
    btclr_arr[:, :, 5] = btclr_arr[:, :, 4] - 1.5  # typical 11-12um BTD ~1.5K
    btclr_arr = np.ascontiguousarray(btclr_arr)

    print(f"  Running Fortran backend ({n_elem}x{n_line})...")
    t0 = time.time()
    result = process_swath_native(
        ref_vis=np.ascontiguousarray(pxldat[:, :, :19].astype(np.float32)),
        tbb_ir=np.ascontiguousarray(pxldat[:, :, 19:].astype(np.float32)),
        lat=np.ascontiguousarray(geo['lat'].astype(np.float32)),
        lon=np.ascontiguousarray(geo['lon'].astype(np.float32)),
        satzen=np.ascontiguousarray(geo['vza'].astype(np.float32)),
        solzen=np.ascontiguousarray(geo['sza'].astype(np.float32)),
        relaz=np.ascontiguousarray(geo['relaz'].astype(np.float32)),
        glint=np.ascontiguousarray(geo['glint_angle'].astype(np.float32)),
        sfctmp=np.ascontiguousarray(nwp_int['sfctmp'].astype(np.float32)),
        pmsl=np.ascontiguousarray(nwp_int['pmsl'].astype(np.float32)),
        uwind=np.ascontiguousarray(nwp_int['u_wind'].astype(np.float32)),
        vwind=np.ascontiguousarray(nwp_int['v_wind'].astype(np.float32)),
        tpw=np.ascontiguousarray(nwp_int['tpw'].astype(np.float32)),
        elev=np.ascontiguousarray(geo['elevation'].astype(np.float32)),
        eco=np.ascontiguousarray(geo['eco_type'].astype(np.int8)),
        lsf=np.ascontiguousarray(geo['lsf'].astype(np.int8)),
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=btclr_arr,
        n_elem=n_elem, n_line=n_line,
    )
    elapsed = time.time() - t0
    print(f"    Fortran: {elapsed:.1f}s ({n_elem*n_line/elapsed/1e6:.1f} Mpix/s)")

    cm = result['cloud_mask']
    conf = result['confidence']

    nan_mask = np.isnan(conf)
    if np.sum(nan_mask) > 0:
        print(f"    Warning: {np.sum(nan_mask)} NaN conf -> 0")
        conf[nan_mask] = 0.0
        cm[nan_mask] = 0

    for c in range(4):
        n = int(np.sum(cm == c))
        print(f"    class {c}: {n:>10,} ({100*n/cm.size:.1f}%)")
    print(f"    Confidence: mean={np.mean(conf):.4f}")

    print(f"  Saving to {out_file}")
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('cm', data=cm)
        f.create_dataset('conf', data=conf)
        f.create_dataset('lat', data=geo['lat'])
        f.create_dataset('lon', data=geo['lon'])
    return out_file


# ---- Process all orbits ----
for name, l1b_path, geo_path in ORBITS:
    print(f"\n{'='*50}")
    print(f"Orbit {name}")
    print(f"{'='*50}")
    for cal_type, cal0_arr, cal1_arr in [
        ('operational', None, None),
        ('recal', cal0, cal1),
    ]:
        run_orbit(name, l1b_path, geo_path, cal0_arr, cal1_arr, cal_type)

print(f"\n{'='*50}")
print("All orbits processed!")
