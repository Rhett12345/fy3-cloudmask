#!/usr/bin/env python
"""Backend speed comparison + cloud mask visualization.

Runs cloud mask on a full FY-3D orbit with Python/Numba backend, times it,
and generates RGB + cloud mask side-by-side visualizations.

If the native C++/Fortran backend is functional, also compares speed.

Usage:
    python scripts/backend_compare_and_viz.py
"""

import os
import sys
import time
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # 让scripts目录可导入
from find_matched_files import find_best_triplet

import h5py
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm import run_cloud_mask_swath
from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native, get_backend_info
from fy3_cloudmask.output import compute_cloud_amount, write_combined_product

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# ============================================================
# Config
# ============================================================
DATE = '20220803'  # ← 只需改这一个日期

_triplet = find_best_triplet(DATE)
if _triplet is None:
    raise FileNotFoundError(f'找不到 {DATE} 的完整 L1B+GEO+NWP 文件组合')
L1B_FILE = str(_triplet['l1b'])
GEO_FILE = str(_triplet['geo'])
NWP_FILE = str(_triplet['nwp'])
# L1B_FILE = '/data/Data_yuq/mersi/20230606//FY3D_MERSI_GBAL_L1_20230606_1440_1000M_MS.HDF'
# GEO_FILE = '/data/Data_yuq/mersi/20230606//FY3D_MERSI_GBAL_L1_20230606_1440_GEO1K_MS.HDF'
# NWP_FILE = '/data/nwp/20230606/ORG/gfs0p25_41L_20230606_09_00'
THRESHOLDS_FILE = '/home/liusy2020/yuq/cloudmask/fy3_cloudmask/config/thresholds/mersi_ii3d_v8.yaml'
OUTPUT_DIR = '/data/Data_yuq/fy3_cloud'
OMP_NUM_THREADS = 8

os.environ['OMP_NUM_THREADS'] = str(OMP_NUM_THREADS)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Planck constants
H_PLANCK = 6.62606876e-34
C_LIGHT = 2.99792458e+08
K_BOLTZMANN = 1.3806503e-23
C1_PLANCK = 2.0 * H_PLANCK * C_LIGHT * C_LIGHT
C2_PLANCK = H_PLANCK * C_LIGHT / K_BOLTZMANN
IR_WAVENUMBERS = np.array([2643.4359, 2471.654, 1382.621, 1168.182, 933.364, 836.941])
TCS_FY3D = np.array([0.9992917440, 0.9994814177, 0.9989956900, 0.9997135336, 0.9980397975, 0.9983777125])
TCI_FY3D = np.array([0.50718071650, 0.3493280160, 0.40925130837, 0.1014073981, 0.57633464244, 0.4317181810])


import re

def extract_datetime_from_filename(filepath):
    basename = os.path.basename(filepath)
    # 匹配类似 _20230606_1440_ 的模式
    match = re.search(r'_(\d{8})_(\d{4})_', basename)
    if match:
        date_str, time_str = match.groups()
        return f"{date_str}_{time_str}"    # 例如 "20230606_1440"
    # 如果文件名匹配不上，回退到文件属性中的观测时间
    return get_datetime_from_hdf5(filepath)
# ============================================================
# Data Readers
# ============================================================
def read_l1b_data(l1b_path):
    with h5py.File(l1b_path, 'r') as f:
        vis_250 = f['Data/EV_250_Aggr.1KM_RefSB'][:].astype(np.float64)
        vis_1km = f['Data/EV_1KM_RefSB'][:].astype(np.float64)
        ir_250 = f['Data/EV_250_Aggr.1KM_Emissive'][:].astype(np.float64)
        ir_1km = f['Data/EV_1KM_Emissive'][:].astype(np.float64)
        vis_cal = f['Calibration/VIS_Cal_Coeff'][:]
        esd = f.attrs['EarthSun Distance Ratio']
        slope_250 = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Slope']
        intercept_250 = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Intercept']
        slope_1km = f['Data/EV_1KM_Emissive'].attrs['Slope']
        intercept_1km = f['Data/EV_1KM_Emissive'].attrs['Intercept']

    n_pixel, n_line = 2048, 2000
    pxldat = np.zeros((n_pixel, n_line, 25), dtype=np.float64)
    esd2 = esd * esd

    for b in range(4):
        c0, c1, c2 = vis_cal[b]
        dn = vis_250[b]
        refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        refl = np.clip(refl, -0.1, 2.0)
        pxldat[:, :, b] = refl.T

    for b in range(15):
        c0, c1, c2 = vis_cal[b + 4]
        dn = vis_1km[b]
        refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        refl = np.clip(refl, -0.1, 2.0)
        pxldat[:, :, b + 4] = refl.T

    # for b in range(4):
    #     radiance_mw = (ir_1km[b] + intercept_1km[b]) * slope_1km[b]
    #     radiance_mw = np.maximum(radiance_mw, 0.01)
    #     wvn = IR_WAVENUMBERS[b]
    #     vs = 100.0 * wvn
    #     bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs**3 / (1e-5 * radiance_mw) + 1.0)
    #     bt = (bt_raw - TCI_FY3D[b]) / TCS_FY3D[b]
    #     bt = np.clip(bt, 150.0, 350.0)
    #     pxldat[:, :, b + 19] = bt.T
    #
    # for b in range(2):
    #     radiance_mw = (ir_250[b] + intercept_250[b]) * slope_250[b]
    #     radiance_mw = np.maximum(radiance_mw, 0.01)
    #     wvn = IR_WAVENUMBERS[b + 4]
    #     vs = 100.0 * wvn
    #     bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs**3 / (1e-5 * radiance_mw) + 1.0)
    #     bt = (bt_raw - TCI_FY3D[b + 4]) / TCS_FY3D[b + 4]
    #     bt = np.clip(bt, 150.0, 350.0)
    #     pxldat[:, :, b + 23] = bt.T
    # EV_250_Aggr.1KM_Emissive → pxldat[23]=11μm, pxldat[24]=12μm
    # 先读250m IR（band24=11μm, band25=12μm），用IR_WAVENUMBERS[4]和[5]
    for b in range(2):
        radiance_mw = (ir_250[b] + intercept_250[b]) * slope_250[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b + 4]  # [4]=933→11μm, [5]=836→12μm
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs ** 3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b + 4]) / TCS_FY3D[b + 4]
        bt = np.clip(bt, 150.0, 350.0)
        bt[bt <= 150.05] = -999.0
        pxldat[:, :, b + 23] = bt.T  # → [23]=11μm, [24]=12μm

    # EV_1KM_Emissive → pxldat[19]=3.8μm, [20]=4.05μm, [21]=7.3μm, [22]=8.5μm
    # 用IR_WAVENUMBERS[0:4]
    for b in range(4):
        radiance_mw = (ir_1km[b] + intercept_1km[b]) * slope_1km[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b]  # [0]=2643→3.8μm, [1]=2471→4.05μm, [2]=1382→7.3μm, [3]=1168→8.5μm
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs ** 3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b]) / TCS_FY3D[b]
        bt = np.clip(bt, 150.0, 350.0)
        bt[bt <= 150.05] = -999.0
        pxldat[:, :, b + 19] = bt.T  # → [19]=3.8μm, [20]=4.05μm, [21]=7.3μm, [22]=8.5μm

    return pxldat


def read_geo_data(geo_path):
    with h5py.File(geo_path, 'r') as f:
        lat = f['Geolocation/Latitude'][:].astype(np.float64)
        lon = f['Geolocation/Longitude'][:].astype(np.float64)
        dem = f['Geolocation/DEM'][:].astype(np.float64)
        lsf = f['Geolocation/LandSeaMask'][:].astype(np.int32)
        eco = f['Geolocation/LandCover'][:].astype(np.int32)
        sza_raw = f['Geolocation/SolarZenith'][:].astype(np.float64)
        vza_raw = f['Geolocation/SensorZenith'][:].astype(np.float64)
        saa_raw = f['Geolocation/SolarAzimuth'][:].astype(np.float64)
        vaa_raw = f['Geolocation/SensorAzimuth'][:].astype(np.float64)

    sza = sza_raw / 100.0
    vza = vza_raw / 100.0
    sza_rad = np.radians(sza)
    vza_rad = np.radians(vza)
    daa_rad = np.radians(saa_raw / 100.0 - vaa_raw / 100.0)
    cos_glint = np.cos(sza_rad) * np.cos(vza_rad) + np.sin(sza_rad) * np.sin(vza_rad) * np.cos(daa_rad)
    cos_glint = np.clip(cos_glint, -1.0, 1.0)
    glint_angle = np.degrees(np.arccos(cos_glint))

    return {
        'lat': lat.T, 'lon': lon.T, 'elevation': dem.T,
        'lsf': lsf.T, 'eco_type': eco.T, 'sza': sza.T,
        'vza': vza.T, 'glint_angle': glint_angle.T,
    }


# def read_nwp_binary(nwp_path, nvar=283):
#     nlon, nlat = 1440, 721
#     data = np.fromfile(nwp_path, dtype=np.float32)
#     expected = nlon * nlat * nvar
#     if data.size < expected:
#         raise ValueError(f"NWP file too small: {data.size} < {expected}")
#     arr = data[:expected].reshape(nlon, nlat, nvar).transpose(2, 1, 0)
#     arr_shifted = np.empty_like(arr)
#     arr_shifted[:, :, :720] = arr[:, :, 720:]
#     arr_shifted[:, :, 720:] = arr[:, :, :720]
#     return {
#         'lon': np.linspace(-180, 179.75, nlon),
#         'lat': np.linspace(-90, 90, nlat),
#         'tsfc': arr_shifted[2], 'pmsl': arr_shifted[1],
#         'u_wind': arr_shifted[7], 'v_wind': arr_shifted[8], 'tpw': arr_shifted[9],
#     }

def read_nwp_binary(nwp_path, nvar=283):
    nlon, nlat = 1440, 721
    data = np.fromfile(nwp_path, dtype=np.float32)
    expected = nlon * nlat * nvar
    if data.size < expected:
        raise ValueError(f"NWP file too small: {data.size} < {expected}")
    # 正确顺序是 (nvar, nlat, nlon)
    arr = data[:expected].reshape(nvar, nlat, nlon)
    # 经度从0~360 → -180~180，把后720列移到前面
    arr_shifted = np.empty_like(arr)
    arr_shifted[:, :, :720] = arr[:, :, 720:]
    arr_shifted[:, :, 720:] = arr[:, :, :720]
    return {
        'lon': np.linspace(-180, 179.75, nlon),
        'lat': np.linspace(-90, 90, nlat),
        'tsfc': arr_shifted[2],   # 地表温度 mean≈282K ✓
        'pmsl': arr_shifted[1],   # 海平面气压 mean≈101047Pa ✓
        'u_wind': arr_shifted[7], # U风 ✓
        'v_wind': arr_shifted[8], # V风 ✓
        'tpw': arr_shifted[9],    # 可降水量 ✓
    }


def interpolate_nwp(nwp, lat_px, lon_px):
    nwp_lat = nwp['lat']
    nwp_lon = nwp['lon']
    lat_idx = np.searchsorted(nwp_lat, lat_px.ravel()).clip(0, len(nwp_lat) - 1)
    lon_idx = np.searchsorted(nwp_lon, lon_px.ravel()).clip(0, len(nwp_lon) - 1)
    shape = lat_px.shape
    fields = {}
    for key in ['tsfc', 'pmsl', 'u_wind', 'v_wind', 'tpw']:
        fields[key] = nwp[key][lat_idx, lon_idx].reshape(shape)
    return fields


# ============================================================
# Visualization
# ============================================================
def make_rgb(pxldat):
    """True-color RGB: R=0.64um, G=0.55um, B=0.47um."""
    r = pxldat[:, :, 2]
    g = pxldat[:, :, 1]
    b = pxldat[:, :, 0]

    def stretch(band, plow=2, phigh=98):
        vmin, vmax = np.percentile(band[band > 0], [plow, phigh])
        band_clip = np.clip(band, vmin, vmax)
        return (band_clip - vmin) / (vmax - vmin + 1e-10)

    rgb = np.stack([stretch(r), stretch(g), stretch(b)], axis=-1)
    return np.clip(rgb, 0, 1)


def save_visualizations(rgb, cm_tmp, py_elapsed, native_elapsed, speedup, n_total, output_subdir, dt_str):
    """Generate RGB + Cloud Mask visualizations as PNG files."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # ---- Full-res standalone cloud mask ----
    step = 2
    cm_sub = cm_tmp[::step, ::step]
    cloudy_pct = 100 * np.sum(cm_tmp == 0) / n_total
    probc_pct = 100 * np.sum(cm_tmp == 1) / n_total
    probcl_pct = 100 * np.sum(cm_tmp == 2) / n_total
    clear_pct = 100 * np.sum(cm_tmp == 3) / n_total

    fig, ax = plt.subplots(figsize=(16, 16))
    cmap = matplotlib.colors.ListedColormap(['#e62020', '#ff9420', '#60a0ff', '#208020'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ax.imshow(np.flipud(cm_sub.transpose(1, 0)), origin='lower',
              cmap=cmap, norm=norm, interpolation='nearest')
    # ax.set_title(f'FY-3D MERSI-II Cloud Mask\n2023-06-06 05:00 UTC | Python backend: {py_elapsed:.0f}s',
    #              fontsize=14, fontweight='bold')
    # ax.set_xlabel('Along-track pixel (subsampled 1:2)')
    # ax.set_ylabel('Cross-track pixel (subsampled 1:2)')
    legend_patches = [
        mpatches.Patch(color='#e62020', label=f'Cloudy           {cloudy_pct:.1f}%'),
        mpatches.Patch(color='#ff9420', label=f'Prob Cloudy      {probc_pct:.1f}%'),
        mpatches.Patch(color='#60a0ff', label=f'Prob Clear       {probcl_pct:.1f}%'),
        mpatches.Patch(color='#208020', label=f'Confident Clear  {clear_pct:.1f}%'),
    ]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=10, framealpha=0.9)
    # cloudmask_path = os.path.join(OUTPUT_DIR, 'FY3D_20230606_0500_CloudMask.png')
    cloudmask_path = os.path.join(output_subdir, f'FY3D_{dt_str}_CloudMask.png')
    fig.savefig(cloudmask_path, dpi=120, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    logger.info(f"  Cloud mask: {cloudmask_path} ({os.path.getsize(cloudmask_path)/1024:.0f} KB)")

    # ---- Side-by-side: RGB + Cloud Mask ----
    step2 = 4
    rgb_sub = rgb[::step2, ::step2, :]
    cm_sub2 = cm_tmp[::step2, ::step2]

    fig2, axes = plt.subplots(1, 2, figsize=(22, 11))
    speedup_str = f' | Native: {native_elapsed:.0f}s ({speedup:.1f}x speedup)' if native_elapsed else ''
    # fig2.suptitle(f'FY-3D MERSI-II Cloud Mask 2023-06-06 05:00 UTC'
    #               f' | Python: {py_elapsed:.0f}s{speedup_str}',
    #               fontsize=14, fontweight='bold')

    axes[0].imshow(np.flipud(rgb_sub.transpose(1, 0, 2)), origin='lower')
    # axes[0].set_title('True Color Composite (R:0.64μm G:0.55μm B:0.47μm)')
    # axes[0].set_xlabel('Along-track pixel')
    # axes[0].set_ylabel('Cross-track pixel')

    axes[1].imshow(np.flipud(cm_sub2.transpose(1, 0)), origin='lower',
                   cmap=cmap, norm=norm, interpolation='nearest')
    # axes[1].set_title('Cloud Mask  (Red=Cloudy Orange=ProbCloudy Blue=ProbClear Green=Clear)')
    # axes[1].set_xlabel('Along-track pixel')
    axes[1].legend(handles=legend_patches, loc='upper right', fontsize=8, framealpha=0.9)

    plt.tight_layout()
    # sidebyside_path = os.path.join(OUTPUT_DIR, 'FY3D_20230606_0500_RGB_CloudMask.png')
    sidebyside_path = os.path.join(output_subdir, f'FY3D_{dt_str}_RGB_CloudMask.png')
    fig2.savefig(sidebyside_path, dpi=120, bbox_inches='tight', facecolor='white')
    plt.close(fig2)
    logger.info(f"  Side-by-side: {sidebyside_path} ({os.path.getsize(sidebyside_path)/1024:.0f} KB)")


# ================================================================
# Main
# ================================================================
def main():
    print("=" * 70)
    print("FY-3D Cloud Mask — Backend Comparison + Visualization")
    print("=" * 70)
    print(f"  OMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS', 'not set')}")
    print(f"  Backend info: {get_backend_info()}")
    print(f"  Output dir: {OUTPUT_DIR}")
    print()

    # ============================================================
    # 动态提取日期时间并创建输出目录
    # ============================================================
    L1B_FILE = str(_triplet['l1b'])
    GEO_FILE = str(_triplet['geo'])
    NWP_FILE = str(_triplet['nwp'])

    dt_str = extract_datetime_from_filename(L1B_FILE)  # 如 "20230606_1440"
    output_subdir = os.path.join(OUTPUT_DIR, dt_str[:8])  # 按日期分子目录
    os.makedirs(output_subdir, exist_ok=True)
    print(f"Output datetime: {dt_str} -> {output_subdir}")
    # ============================================================

    # --- Load data ---
    logger.info("Loading L1b data...")
    t0 = time.time()
    pxldat = read_l1b_data(L1B_FILE)
    logger.info(f"  L1b loaded ({time.time()-t0:.1f}s), shape={pxldat.shape}")
    logger.info("Loading GEO data...")
    t0 = time.time()
    geo = read_geo_data(GEO_FILE)
    logger.info(f"  GEO loaded ({time.time()-t0:.1f}s)")

    logger.info("Loading NWP data...")
    t0 = time.time()
    nwp = read_nwp_binary(NWP_FILE)
    logger.info(f"  NWP loaded ({time.time()-t0:.1f}s)")

    lat, lon = geo['lat'], geo['lon']
    n_elem, n_line = lat.shape

    logger.info("Interpolating NWP to swath...")
    t0 = time.time()
    nwp_interp = interpolate_nwp(nwp, lat, lon)
    logger.info(f"  NWP interpolated ({time.time()-t0:.1f}s)")

    # Prepare common inputs
    nwp_sfctmp = nwp_interp['tsfc']
    nwp_pmsl = nwp_interp['pmsl'] / 100.0
    nwp_u = nwp_interp['u_wind']
    nwp_v = nwp_interp['v_wind']
    nwp_pw = nwp_interp['tpw']

    bt_clr = np.zeros((n_elem, n_line, 7), dtype=np.float64)
    bt_clr[:, :, 4] = nwp_sfctmp
    bt_clr[:, :, 5] = nwp_sfctmp - 1.0
    bt_clr[:, :, 0] = nwp_sfctmp - 10.0
    bt_clr[:, :, 1] = nwp_sfctmp - 20.0
    bt_clr[:, :, 2] = nwp_sfctmp - 10.0
    bt_clr[:, :, 3] = nwp_sfctmp - 5.0

    snow_mask = np.zeros((n_elem, n_line), dtype=np.int32)
    sst = np.zeros((n_elem, n_line), dtype=np.float64)

    with open(THRESHOLDS_FILE) as f:
        thresholds = yaml.safe_load(f)

    n_total = n_elem * n_line
    print(f"\n  Swath: {n_elem} x {n_line} = {n_total:,} pixels\n")

    # ================================================================
    # 1. Python/Numba backend (always works)
    # ================================================================
    print("=" * 70)
    print("BACKEND 1: Python + Numba JIT (pure Python)")
    print("=" * 70)
    t0 = time.time()
    py_bitarray, py_qa, py_cm, py_conf = run_cloud_mask_swath(
        pxldat_swath=pxldat,
        lat_swath=lat, lon_swath=lon,
        elevation_swath=geo['elevation'],
        lsf_swath=geo['lsf'],
        sza_swath=geo['sza'],
        vza_swath=geo['vza'],
        glint_angle_swath=geo['glint_angle'],
        eco_type_swath=geo['eco_type'],
        snow_mask_swath=snow_mask, sst_swath=sst,
        nwp_sfctmp_swath=nwp_sfctmp, nwp_pmsl_swath=nwp_pmsl,
        nwp_u_wind_swath=nwp_u, nwp_v_wind_swath=nwp_v,
        nwp_precip_water_swath=nwp_pw,
        bt_clr_swath=bt_clr,
        sensor_id=21, thresholds=thresholds,
    )
    py_elapsed = time.time() - t0
    py_rate = n_total / py_elapsed
    cloudy_pct = 100 * np.sum(py_cm == 0) / n_total
    pcloudy_pct = 100 * np.sum(py_cm == 1) / n_total
    pclear_pct = 100 * np.sum(py_cm == 2) / n_total
    clear_pct = 100 * np.sum(py_cm == 3) / n_total

    print(f"  Time:       {py_elapsed:.1f}s")
    print(f"  Rate:       {py_rate:,.0f} pixels/s ({n_total/py_elapsed/1e6:.2f} MPixels/s)")
    print(f"  Cloudy:      {np.sum(py_cm==0):>8,} ({cloudy_pct:5.1f}%)")
    print(f"  Prob Cloudy: {np.sum(py_cm==1):>8,} ({pcloudy_pct:5.1f}%)")
    print(f"  Prob Clear:  {np.sum(py_cm==2):>8,} ({pclear_pct:5.1f}%)")
    print(f"  Conf Clear:  {np.sum(py_cm==3):>8,} ({clear_pct:5.1f}%)")
    print(f"  Mean conf:   {py_conf.mean():.4f}")

    # ================================================================
    # 2. Native C++/Fortran backend
    # ================================================================
    native_elapsed = None
    speedup = None
    native_cm = None
    native_ok = False

    print()
    print("=" * 70)
    print(f"BACKEND 2: C++/Fortran + OpenMP (OMP_NUM_THREADS={OMP_NUM_THREADS})")
    print("=" * 70)

    if not is_native_available():
        print("  SKIP: Native backend .so not available")
    else:
        try:
            ref_vis = pxldat[:, :, :19].astype(np.float32)
            tbb_ir = pxldat[:, :, 19:25].astype(np.float32)
            lat_f32 = lat.astype(np.float32)
            lon_f32 = lon.astype(np.float32)
            satzen = geo['vza'].astype(np.float32)
            solzen = geo['sza'].astype(np.float32)
            relaz = np.zeros_like(lat_f32)
            glint = geo['glint_angle'].astype(np.float32)
            sfctmp = np.ascontiguousarray(nwp_sfctmp.astype(np.float32))
            pmsl = np.ascontiguousarray(nwp_pmsl.astype(np.float32))
            uwind = np.ascontiguousarray(nwp_u.astype(np.float32))
            vwind = np.ascontiguousarray(nwp_v.astype(np.float32))
            tpw = np.ascontiguousarray(nwp_pw.astype(np.float32))
            elev = geo['elevation'].astype(np.float32)
            eco = geo['eco_type'].astype(np.int8)
            snmask = snow_mask.astype(np.int8)
            btclr_f32 = bt_clr.astype(np.float32)

            print(f"  Calling process_swath_native({n_elem}, {n_line}) ...")
            t0 = time.time()
            native_result = process_swath_native(
                ref_vis, tbb_ir, lat_f32, lon_f32,
                satzen, solzen, relaz, glint,
                sfctmp, pmsl, uwind, vwind, tpw,
                elev, eco, snmask, btclr_f32,
                n_elem, n_line,
            )
            native_elapsed = time.time() - t0
            native_rate = n_total / native_elapsed
            native_cm = native_result['cloud_mask']
            native_conf = native_result['confidence']
            native_ok = True

            nc_cloudy = 100 * np.sum(native_cm == 0) / n_total
            nc_clear = 100 * np.sum(native_cm == 3) / n_total
            print(f"  Time:       {native_elapsed:.1f}s")
            print(f"  Rate:       {native_rate:,.0f} pixels/s ({n_total/native_elapsed/1e6:.2f} MPixels/s)")
            print(f"  Cloudy:     {np.sum(native_cm==0):>8,} ({nc_cloudy:.1f}%)")
            print(f"  Clear:      {np.sum(native_cm==3):>8,} ({nc_clear:.1f}%)")
        except Exception as e:
            print(f"  ERROR: {e}")
            print("  Native backend crashed — likely libgfortran ABI mismatch.")
            print("  The pre-compiled .so needs gfortran >= 10, system has older libgfortran.")
            print("  To fix: rebuild with 'cd ext && ./build.sh --install' using conda gfortran.")

    # ================================================================
    # 3. Comparison
    # ================================================================
    print()
    print("=" * 70)
    print("COMPARISON")
    print("=" * 70)
    if native_ok:
        speedup = py_elapsed / native_elapsed
        agreement = np.mean(py_cm == native_cm)
        conf_mae = np.mean(np.abs(py_conf - native_conf))
        print(f"  Python:     {py_elapsed:.1f}s ({py_rate:,.0f} px/s)")
        print(f"  Native:     {native_elapsed:.1f}s ({native_rate:,.0f} px/s)")
        print(f"  Speedup:    {speedup:.1f}x")
        print(f"  Agreement:  {agreement*100:.2f}% pixel-level")
        print(f"  Conf MAE:   {conf_mae:.6f}")
    else:
        print(f"  Python only: {py_elapsed:.1f}s ({py_rate:,.0f} px/s)")
        print("  Native backend not functional — no comparison available.")

    # ================================================================
    # 4. Visualization
    # ================================================================
    print()
    print("=" * 70)
    print("VISUALIZATION")
    print("=" * 70)

    logger.info("Generating RGB composite...")
    rgb = make_rgb(pxldat)
    logger.info(f"  RGB shape: {rgb.shape}, range: [{rgb.min():.2f}, {rgb.max():.2f}]")

    # 对齐云掩码与RGB的行列方向
    logger.info(f"  py_cm shape before fix: {py_cm.shape}")
    if py_cm.shape != rgb.shape[:2]:
        logger.info(f"  Shape mismatch! Transposing py_cm: {py_cm.shape} -> {py_cm.T.shape}")
        py_cm = py_cm.T
        py_conf = py_conf.T
        py_qa = py_qa.T
        py_bitarray = py_bitarray.T
    logger.info(f"  py_cm shape after fix: {py_cm.shape}")

    logger.info("Computing 5km cloud amount...")
    qa_proc = (py_cm >= 0).astype(np.int32)
    cloud_amount, cloud_amount_qa = compute_cloud_amount(py_cm, qa_proc)

    # 统计各CM类别的平均RGB亮度
    print("\n  各CM类别平均RGB_R亮度：")
    r_band = pxldat[:, :, 2]
    for c, name in [(0, 'Cloudy'), (1, 'ProbCloudy'), (2, 'ProbClear'), (3, 'Clear')]:
        mask = py_cm == c
        if mask.any():
            print(f"  {name:12s}: 平均RGB_R={r_band[mask].mean():.3f}  像素数={mask.sum():>8,}")

    print("=========================\n")

    save_visualizations(rgb, py_cm, py_elapsed, native_elapsed, speedup, n_total, output_subdir, dt_str)

    # ================================================================
    # 5. Write HDF5
    # ================================================================
    logger.info("Writing HDF5 output...")
    box_size = 5
    n_elem_5km = n_elem // box_size
    n_line_5km = n_line // box_size
    lon_5km = lon[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]
    lat_5km = lat[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]

    # h5_path = os.path.join(OUTPUT_DIR, 'FY3D_MERSI_20230606_0500_CLM_CLA.h5')
    h5_path = os.path.join(output_subdir, f'FY3D_MERSI_{dt_str}_CLM_CLA.h5')
    write_combined_product(
        h5_path,
        py_bitarray, py_qa, py_cm, py_conf,
        cloud_amount, cloud_amount_qa,
        lon, lat, lon_5km, lat_5km,
        attributes={
            'input_l1b': L1B_FILE,
            'input_geo': GEO_FILE,
            'input_nwp': NWP_FILE,
            'backend': 'Python+Numba',
            'python_time_s': f'{py_elapsed:.1f}',
            'native_time_s': f'{native_elapsed:.1f}' if native_elapsed else 'N/A',
            'speedup': f'{speedup:.1f}x' if speedup else 'N/A',
        },
    )
    logger.info(f"  HDF5: {h5_path}")

    # ================================================================
    # Summary
    # ================================================================
    print()
    print("=" * 70)
    print("OUTPUT FILES")
    print("=" * 70)
    for f in sorted(Path(output_subdir).glob("*.png")):
        size_kb = f.stat().st_size / 1024
        print(f"  {f}  ({size_kb:.0f} KB)")
    for f in sorted(Path(output_subdir).glob("*.h5")):
        size_mb = f.stat().st_size / 1024 / 1024
        print(f"  {f}  ({size_mb:.1f} MB)")

    print()
    print("Done.")


if __name__ == '__main__':
    main()