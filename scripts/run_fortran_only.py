#!/usr/bin/env python3
"""Run Fortran-only cloud mask on multiple orbits for stability/performance testing.

Processes 5-7 dates (~15-21 orbits) using ONLY the native C++/Fortran backend.
Reports per-phase timing, stability metrics, and CPU performance.

Usage:
    python scripts/run_fortran_only.py
    python scripts/run_fortran_only.py --omp-threads 4
    python scripts/run_fortran_only.py --output /data/Data_yuq/fy3_cloud/validation_fortran
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# ---- Planck constants for BT conversion ----
H_PLANCK = 6.62606876e-34
C_LIGHT = 2.99792458e+08
K_BOLTZMANN = 1.3806503e-23
C1_PLANCK = 2.0 * H_PLANCK * C_LIGHT * C_LIGHT
C2_PLANCK = H_PLANCK * C_LIGHT / K_BOLTZMANN
IR_WAVENUMBERS = np.array([2643.4359, 2471.654, 1382.621, 1168.182, 909.458, 836.941])
TCS_FY3D = np.array([0.9992917440, 0.9994814177, 0.9989956900, 0.9997135336, 0.9980397975, 0.9983777125])
TCI_FY3D = np.array([0.50718071650, 0.3493280160, 0.40925130837, 0.1014073981, 0.57633464244, 0.4317181810])

# ---- Confirmed test dates with complete L1B+GEO+NWP triplets ----
CONFIRMED_DATES = [
    ("20220803", "summer", "Mid-latitude summer, 3 morning orbits"),
    ("20220816", "summer", "Mid-latitude summer, 3 orbits"),
    ("20230102", "winter", "Winter cold land, 3 early morning orbits"),
    ("20230403", "spring", "Spring transitional, 3 orbits"),
    ("20240102", "winter", "Early winter, 3 morning orbits"),
    ("20240402", "spring", "Spring, 3 orbits"),
    ("20250302", "spring", "Spring snow melt, 3 midday orbits"),
]


def read_l1b_data(
    l1b_path: str,
    recal_cal0: np.ndarray | None = None,
    recal_cal1: np.ndarray | None = None,
) -> np.ndarray:
    """Read L1b HDF5 and convert to reflectance/BT. Returns (2048, 2000, 25).

    Parameters
    ----------
    l1b_path : str
        Path to L1b HDF5 file.
    recal_cal0 : ndarray, shape (7,), optional
        Recalibration intercept for bands 0-6. If None, uses onboard calibration.
    recal_cal1 : ndarray, shape (7,), optional
        Recalibration slope for bands 0-6. If None, uses onboard calibration.
    """
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

    use_recal = recal_cal0 is not None and recal_cal1 is not None

    # Bands 0-3: 250m channels (from vis_250)
    for b in range(4):
        dn = vis_250[b]
        if use_recal and b < 7:
            refl = recal_cal0[b] + recal_cal1[b] * dn
        else:
            c0, c1, c2 = vis_cal[b]
            refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        pxldat[:, :, b] = np.clip(refl, -0.1, 2.0).T

    # Bands 4-18: 1km channels (from vis_1km)
    for b in range(15):
        dn = vis_1km[b]
        if use_recal and (b + 4) < 7:
            refl = recal_cal0[b + 4] + recal_cal1[b + 4] * dn
        else:
            c0, c1, c2 = vis_cal[b + 4]
            refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        pxldat[:, :, b + 4] = np.clip(refl, -0.1, 2.0).T

    for b in range(2):
        radiance_mw = (ir_250[b] + intercept_250[b]) * slope_250[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b + 4]
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs ** 3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b + 4]) / TCS_FY3D[b + 4]
        bt = np.clip(bt, 150.0, 350.0)
        bt[bt <= 150.05] = -999.0
        pxldat[:, :, b + 23] = bt.T

    for b in range(4):
        radiance_mw = (ir_1km[b] + intercept_1km[b]) * slope_1km[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b]
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs ** 3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b]) / TCS_FY3D[b]
        bt = np.clip(bt, 150.0, 350.0)
        bt[bt <= 150.05] = -999.0
        pxldat[:, :, b + 19] = bt.T

    return pxldat


def read_geo_data(geo_path: str) -> dict:
    """Read GEO HDF5."""
    with h5py.File(geo_path, 'r') as f:
        lat = f['Geolocation/Latitude'][:].astype(np.float64)
        lon = f['Geolocation/Longitude'][:].astype(np.float64)
        dem = f['Geolocation/DEM'][:].astype(np.float64)
        eco = f['Geolocation/LandCover'][:].astype(np.int32)
        lsf_raw = f['Geolocation/LandSeaMask'][:].astype(np.int8)
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
    glint_angle = np.degrees(np.arccos(np.clip(cos_glint, -1.0, 1.0)))

    return {
        'lat': lat.T, 'lon': lon.T, 'elevation': dem.T,
        'eco_type': eco.T, 'sza': sza.T, 'vza': vza.T,
        'glint_angle': glint_angle.T, 'lsf': lsf_raw.T,
    }


def read_nwp_binary(nwp_path: str, nvar: int = 283) -> dict:
    """Read NWP binary file. Supports both gfs0p25 (0.25-deg) and fnl (1-deg) formats."""
    data = np.fromfile(nwp_path, dtype=np.float32)

    # Detect format by file size
    # gfs0p25: 1440x721x283 = 293,821,920 values (1.1GB)
    # fnl: 360x181x179 = 11,633,640 values (44.5MB)
    if data.size >= 293821920:
        # gfs0p25 format (0.25-degree, 283 variables)
        nlon, nlat = 1440, 721
        nvar_actual = 283
        arr = data[:nlon * nlat * nvar_actual].reshape(nvar_actual, nlat, nlon)
        arr_shifted = np.empty_like(arr)
        arr_shifted[:, :, :720] = arr[:, :, 720:]
        arr_shifted[:, :, 720:] = arr[:, :, :720]
        lon = np.linspace(-180, 179.75, nlon)
        lat = np.linspace(-90, 90, nlat)
        # Variable indices for gfs0p25: tsfc=2, pmsl=1, u_wind=7, v_wind=8, tpw=9
        return {
            'lon': lon, 'lat': lat,
            'tsfc': arr_shifted[2], 'pmsl': arr_shifted[1],
            'u_wind': arr_shifted[7], 'v_wind': arr_shifted[8], 'tpw': arr_shifted[9],
        }
    else:
        # fnl format (1-degree, 179 variables)
        nlon, nlat = 360, 181
        nvar_actual = data.size // (nlon * nlat)
        if nvar_actual < 10:
            raise ValueError(f"NWP file too small for fnl format: {data.size} values")
        arr = data[:nvar_actual * nlat * nlon].reshape(nvar_actual, nlat, nlon)
        # fnl format: longitude starts at 0, shift to -180..180
        arr_shifted = np.empty_like(arr)
        arr_shifted[:, :, :180] = arr[:, :, 180:]
        arr_shifted[:, :, 180:] = arr[:, :, :180]
        lon = np.linspace(-180, 179, nlon)
        lat = np.linspace(-90, 90, nlat)
        # Variable indices for fnl: pmsl=1, tsfc=2, u_wind=4, v_wind=5, tpw=9
        return {
            'lon': lon, 'lat': lat,
            'tsfc': arr_shifted[2], 'pmsl': arr_shifted[1],
            'u_wind': arr_shifted[4], 'v_wind': arr_shifted[5], 'tpw': arr_shifted[9],
        }


def interpolate_nwp(nwp: dict, lat_px: np.ndarray, lon_px: np.ndarray) -> dict:
    """Nearest-neighbor interpolation of NWP to pixel locations."""
    nwp_lat = nwp['lat']
    nwp_lon = nwp['lon']
    lat_idx = np.searchsorted(nwp_lat, lat_px.ravel()).clip(0, len(nwp_lat) - 1)
    lon_idx = np.searchsorted(nwp_lon, lon_px.ravel()).clip(0, len(nwp_lon) - 1)
    shape = lat_px.shape
    fields = {}
    for key in ['tsfc', 'pmsl', 'u_wind', 'v_wind', 'tpw']:
        fields[key] = nwp[key][lat_idx, lon_idx].reshape(shape)
    return fields


def run_single_orbit(
    l1b_path: str,
    geo_path: str,
    nwp_path: str,
    output_dir: str,
    recal_cal0: np.ndarray | None = None,
    recal_cal1: np.ndarray | None = None,
    suffix: str = "",
) -> dict:
    """Process one orbit with Fortran-only backend. Returns timing/stats dict.

    Parameters
    ----------
    l1b_path, geo_path, nwp_path : str
        Input file paths.
    output_dir : str
        Output directory.
    recal_cal0, recal_cal1 : ndarray, shape (7,), optional
        Recalibration coefficients for bands 0-6. If None, uses onboard calibration.
    """
    orbit_tag = Path(l1b_path).stem.split('_1000M')[0].split('_')[-2:]
    orbit_tag = '_'.join(orbit_tag)
    date_dir = orbit_tag[:8]
    orbit_dir = os.path.join(output_dir, date_dir)
    os.makedirs(orbit_dir, exist_ok=True)

    logger.info(f"Processing orbit: {orbit_tag}")

    # Phase 1: Data reading
    t_io_start = time.time()

    logger.info("  Reading L1b data...")
    t0 = time.time()
    pxldat = read_l1b_data(l1b_path, recal_cal0=recal_cal0, recal_cal1=recal_cal1)
    t_l1b = time.time() - t0

    logger.info("  Reading GEO data...")
    t0 = time.time()
    geo = read_geo_data(geo_path)
    t_geo = time.time() - t0

    logger.info("  Reading NWP data...")
    t0 = time.time()
    nwp = read_nwp_binary(nwp_path)
    t_nwp_read = time.time() - t0

    logger.info("  Interpolating NWP...")
    t0 = time.time()
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])
    t_nwp_interp = time.time() - t0

    t_io_total = time.time() - t_io_start

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    n_total = n_elem * n_line

    # Phase 2: Prepare arrays for Fortran
    t_prep_start = time.time()
    ref_vis = np.ascontiguousarray(pxldat[:, :, :19].astype(np.float32))
    tbb_ir = np.ascontiguousarray(pxldat[:, :, 19:].astype(np.float32))
    sfctmp = np.ascontiguousarray(nwp_interp['tsfc'].astype(np.float32))
    pmsl = np.ascontiguousarray(nwp_interp['pmsl'].astype(np.float32))
    uwind = np.ascontiguousarray(nwp_interp['u_wind'].astype(np.float32))
    vwind = np.ascontiguousarray(nwp_interp['v_wind'].astype(np.float32))
    tpw = np.ascontiguousarray(nwp_interp['tpw'].astype(np.float32))
    t_prep = time.time() - t_prep_start

    # Phase 3: Fortran compute
    logger.info("  Running Fortran native backend...")
    t0 = time.time()
    result = process_swath_native(
        ref_vis=ref_vis, tbb_ir=tbb_ir,
        lat=np.ascontiguousarray(geo['lat'].astype(np.float32)),
        lon=np.ascontiguousarray(geo['lon'].astype(np.float32)),
        satzen=np.ascontiguousarray(geo['vza'].astype(np.float32)),
        solzen=np.ascontiguousarray(geo['sza'].astype(np.float32)),
        relaz=np.ascontiguousarray(np.zeros_like(geo['sza']).astype(np.float32)),
        glint=np.ascontiguousarray(geo['glint_angle'].astype(np.float32)),
        sfctmp=sfctmp, pmsl=pmsl, uwind=uwind, vwind=vwind, tpw=tpw,
        elev=np.ascontiguousarray(geo['elevation'].astype(np.float32)),
        eco=np.ascontiguousarray(geo['eco_type'].astype(np.int8)),
        lsf=np.ascontiguousarray(geo['lsf'].astype(np.int8)),
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    t_fortran = time.time() - t0

    # Phase 4: Save output
    t_save_start = time.time()
    h5_path = os.path.join(orbit_dir, f'FY3D_MERSI_{orbit_tag}_CLM_CLA{suffix}.h5')
    with h5py.File(h5_path, 'w') as f:
        grp = f.create_group('Cloud_Mask_1km')
        grp.create_dataset('Latitude',  data=geo['lat'].astype(np.float32))
        grp.create_dataset('Longitude', data=geo['lon'].astype(np.float32))
        grp.create_dataset('Cloud_Mask_Value', data=result['cloud_mask'].astype(np.int8))
        grp.create_dataset('Confidence', data=result['confidence'].astype(np.float32))
        tb = result.get('cm_bitarray', result.get('testbits', np.zeros((n_elem, n_line, 6), dtype=np.uint8)))
        qa = result.get('qa_bitarray', result.get('qa_bits', np.zeros((n_elem, n_line, 10), dtype=np.uint8)))
        grp.create_dataset('TestBits', data=tb.astype(np.uint8))
        grp.create_dataset('QA_Bits',  data=qa.astype(np.uint8))
    t_save = time.time() - t_save_start

    t_total = t_io_total + t_prep + t_fortran + t_save

    # Stability analysis
    cm = result['cloud_mask']
    conf = result['confidence']
    valid = (cm >= 0) & (cm <= 3)
    n_valid = int(np.sum(valid))

    cm_dist = {}
    for cat in range(4):
        cm_dist[int(cat)] = int(np.sum(cm[valid] == cat)) if n_valid > 0 else 0

    conf_valid = conf[valid] if n_valid > 0 else np.array([0.0])
    conf_mean = float(np.mean(conf_valid))
    conf_std = float(np.std(conf_valid))

    # Check for anomalous patterns
    fill_pixels = int(np.sum(cm < 0)) + int(np.sum(cm > 3))
    stability = "OK" if n_valid > 0 and fill_pixels < n_total * 0.5 else "WARN"

    stats = {
        'orbit_tag': orbit_tag,
        'l1b': l1b_path,
        'geo': geo_path,
        'nwp': nwp_path,
        'n_elem': n_elem,
        'n_line': n_line,
        'n_total': n_total,
        'n_valid': n_valid,
        'fill_pixels': fill_pixels,
        'stability': stability,
        'timing': {
            'l1b_read': t_l1b,
            'geo_read': t_geo,
            'nwp_read': t_nwp_read,
            'nwp_interp': t_nwp_interp,
            'io_total': t_io_total,
            'array_prep': t_prep,
            'fortran_compute': t_fortran,
            'save': t_save,
            'total': t_total,
        },
        'cloud_mask_distribution': cm_dist,
        'confidence': {'mean': conf_mean, 'std': conf_std},
        'throughput_mpix_per_sec': (n_total / 1e6) / t_fortran if t_fortran > 0 else 0,
    }

    logger.info(f"  Done in {t_total:.1f}s (Fortran: {t_fortran:.1f}s, "
                f"IO: {t_io_total:.1f}s, throughput: {stats['throughput_mpix_per_sec']:.1f} Mpix/s)")
    logger.info(f"  Valid: {n_valid}/{n_total}, stability: {stability}")
    logger.info(f"  CM distribution: cloudy={cm_dist.get(0,0)}, prob_cloudy={cm_dist.get(1,0)}, "
                f"prob_clear={cm_dist.get(2,0)}, clear={cm_dist.get(3,0)}")

    return stats


def main():
    parser = argparse.ArgumentParser(description="Fortran-only batch cloud mask")
    parser.add_argument("--output", default="/data/Data_yuq/fy3_cloud/validation_fortran",
                        help="Output directory")
    parser.add_argument("--omp-threads", type=int, default=8, help="OpenMP threads")
    parser.add_argument("--max-orbits", type=int, default=3, help="Max orbits per date")
    parser.add_argument("--recal-dir", default=None,
                        help="Path to recalibration data directory (e.g., ../fy3d_recali)")
    parser.add_argument("--date", default=None,
                        help="Process specific date (YYYYMMDD) instead of hardcoded list")
    args = parser.parse_args()

    os.environ['OMP_NUM_THREADS'] = str(args.omp_threads)

    if not is_native_available():
        logger.error("Native backend NOT available! Build with: cd ext && ./build.sh --install")
        sys.exit(1)

    logger.info(f"Native backend available. OMP_NUM_THREADS={args.omp_threads}")

    # Setup recalibration if requested
    recal_mgr = None
    if args.recal_dir:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))
        from fy3_cloudmask.io.recalibration import RecalibrationManager
        recal_mgr = RecalibrationManager(args.recal_dir)
        logger.info(f"Recalibration enabled: {args.recal_dir}")

    # Import find_matched_files to get triplets
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from find_matched_files import find_matched_triplets

    os.makedirs(args.output, exist_ok=True)
    all_stats = []

    # Use provided date or hardcoded list
    if args.date:
        dates_to_process = [(args.date, "custom", f"Custom date {args.date}")]
    else:
        dates_to_process = CONFIRMED_DATES

    for date_str, category, description in dates_to_process:
        logger.info(f"\n{'='*60}")
        logger.info(f"Date: {date_str} ({category}) - {description}")
        logger.info(f"{'='*60}")

        # Load recalibration coefficients for this date
        recal_cal0, recal_cal1 = None, None
        if recal_mgr is not None:
            try:
                recal_cal0, recal_cal1 = recal_mgr.load_coefficients(date_str)
                logger.info(f"  Recalibration loaded for {date_str}")
            except FileNotFoundError as e:
                logger.warning(f"  Recalibration not found for {date_str}: {e}")

        triplets = find_matched_triplets(date_str)
        if not triplets:
            logger.warning(f"  No triplets found for {date_str}, skipping")
            continue

        selected = triplets[:args.max_orbits]
        for t in selected:
            l1b = str(t['l1b'])
            orbit_tag = '_'.join(Path(l1b).stem.split('_1000M')[0].split('_')[-2:])
            date_orbit_dir = os.path.join(args.output, date_str)

            # Onboard calibration
            onboard_h5 = os.path.join(date_orbit_dir, f'FY3D_MERSI_{orbit_tag}_CLM_CLA.h5')
            if os.path.exists(onboard_h5):
                logger.info(f"  [SKIP] {orbit_tag} onboard (exists)")
            else:
                try:
                    stats = run_single_orbit(
                        l1b_path=l1b,
                        geo_path=str(t['geo']),
                        nwp_path=str(t['nwp']),
                        output_dir=args.output,
                        suffix="",
                    )
                    stats['date'] = date_str
                    stats['category'] = category
                    stats['description'] = description
                    stats['recalibration'] = False
                    all_stats.append(stats)
                except Exception as e:
                    logger.error(f"  Failed (onboard): {e}")
                    import traceback
                    traceback.print_exc()

            # Recalibration
            recal_h5 = os.path.join(date_orbit_dir, f'FY3D_MERSI_{orbit_tag}_CLM_CLA_recal.h5')
            if recal_cal0 is None:
                pass
            elif os.path.exists(recal_h5):
                logger.info(f"  [SKIP] {orbit_tag} recal (exists)")
            else:
                try:
                    stats = run_single_orbit(
                        l1b_path=l1b,
                        geo_path=str(t['geo']),
                        nwp_path=str(t['nwp']),
                        output_dir=args.output,
                        recal_cal0=recal_cal0,
                        recal_cal1=recal_cal1,
                        suffix="_recal",
                    )
                    stats['date'] = date_str
                    stats['category'] = category
                    stats['description'] = description
                    stats['recalibration'] = True
                    all_stats.append(stats)
                except Exception as e:
                    logger.error(f"  Failed (recal): {e}")
                    import traceback
                    traceback.print_exc()

    # Print summary table
    print(f"\n{'='*100}")
    print(f"{'Orbit':<25} {'Category':<10} {'Total(s)':>9} {'F90(s)':>8} {'IO(s)':>7} "
          f"{'Mpix/s':>8} {'Valid%':>7} {'Stability':>10}")
    print(f"{'='*100}")

    for s in all_stats:
        t = s['timing']
        valid_pct = 100.0 * s['n_valid'] / s['n_total'] if s['n_total'] > 0 else 0
        print(f"{s['orbit_tag']:<25} {s['category']:<10} "
              f"{t['total']:>8.1f} {t['fortran_compute']:>7.1f} {t['io_total']:>6.1f} "
              f"{s['throughput_mpix_per_sec']:>7.1f} {valid_pct:>6.1f}% {s['stability']:>10}")

    # Performance summary
    if all_stats:
        fortran_times = [s['timing']['fortran_compute'] for s in all_stats]
        total_times = [s['timing']['total'] for s in all_stats]
        throughputs = [s['throughput_mpix_per_sec'] for s in all_stats]
        print(f"\n{'='*100}")
        print("PERFORMANCE SUMMARY")
        print(f"{'='*100}")
        print(f"  Orbits processed:  {len(all_stats)}")
        print(f"  Fortran compute:   mean={np.mean(fortran_times):.1f}s, "
              f"min={np.min(fortran_times):.1f}s, max={np.max(fortran_times):.1f}s")
        print(f"  Total (incl IO):   mean={np.mean(total_times):.1f}s, "
              f"min={np.min(total_times):.1f}s, max={np.max(total_times):.1f}s")
        print(f"  Throughput:        mean={np.mean(throughputs):.1f} Mpix/s, "
              f"max={np.max(throughputs):.1f} Mpix/s")
        print(f"  IO overhead:       mean={np.mean([s['timing']['io_total'] for s in all_stats]):.1f}s "
              f"({100*np.mean([s['timing']['io_total'] for s in all_stats])/np.mean(total_times):.0f}% of total)")

        # Stability check
        unstable = [s for s in all_stats if s['stability'] != 'OK']
        if unstable:
            print(f"\n  WARNING: {len(unstable)} orbits have stability issues!")
            for s in unstable:
                print(f"    {s['orbit_tag']}: {s['stability']}")
        else:
            print(f"\n  All {len(all_stats)} orbits processed stably.")

        # Per-surface-type breakdown from cloud mask distribution
        print(f"\n{'='*100}")
        print("CLOUD MASK DISTRIBUTION (aggregated)")
        print(f"{'='*100}")
        total_cm = {0: 0, 1: 0, 2: 0, 3: 0}
        for s in all_stats:
            for k, v in s['cloud_mask_distribution'].items():
                total_cm[int(k)] += v
        total_valid = sum(total_cm.values())
        if total_valid > 0:
            labels = {0: 'cloudy', 1: 'prob_cloudy', 2: 'prob_clear', 3: 'confident_clear'}
            for k in range(4):
                pct = 100.0 * total_cm[k] / total_valid
                print(f"  {labels[k]:<20} {total_cm[k]:>12,} ({pct:>5.1f}%)")


if __name__ == "__main__":
    main()
