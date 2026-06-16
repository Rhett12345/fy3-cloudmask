#!/usr/bin/env python3
"""Debug isolated polar pixels: trace test group confidences.

Usage:
    python scripts/debug_polar_isolated.py
"""

import os
import sys
import time
from pathlib import Path

import numpy as np

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import process_swath_native
from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_pixel, _extract_3x3
from fy3_cloudmask.algorithm.tests.polar_day import polar_day_land
from fy3_cloudmask.algorithm.surface_classifier import classify_pixel_surface
from fy3_cloudmask.config import load_config

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp


def main():
    date_str = "20220803"
    time_tag = "0740"

    mersi_dir = Path('/data/Data_yuq/mersi') / date_str
    l1b_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF"
    geo_path = mersi_dir / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_GEO1K_MS.HDF"
    nwp_path = Path('/data/nwp') / date_str / 'ORG' / f'gfs0p25_41L_{date_str}_06_00'

    logger.info("Reading data...")
    pxldat = read_l1b_data(str(l1b_path))
    geo = read_geo_data(str(geo_path))
    nwp = read_nwp_binary(str(nwp_path))
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    logger.info(f"Swath: {n_elem} x {n_line}")

    # Run Fortran backend
    logger.info("Running Fortran backend...")
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
    logger.info(f"Fortran: {time.time()-t0:.1f}s")

    cm = result['cloud_mask']
    conf = result['confidence']

    # Find isolated cloudy pixels (cm=0, all 4 neighbors != 0)
    cloudy = (cm == 0)
    padded = np.pad(cm, 1, mode='edge')
    up = padded[:-2, 1:-1]
    down = padded[2:, 1:-1]
    left = padded[1:-1, :-2]
    right = padded[1:-1, 2:]
    all_diff = (up != cm) & (down != cm) & (left != cm) & (right != cm)
    isolated = cloudy & all_diff

    # Filter to polar region (lat > 60)
    polar_mask = geo['lat'] > 60.0
    isolated_polar = isolated & polar_mask

    # Filter to interior (not edge)
    interior = np.ones_like(cm, dtype=bool)
    interior[0, :] = False
    interior[-1, :] = False
    interior[:, 0] = False
    interior[:, -1] = False
    isolated_polar_interior = isolated_polar & interior

    n_isolated = int(np.sum(isolated_polar_interior))
    logger.info(f"Isolated polar interior cloudy pixels: {n_isolated}")

    # Get coordinates of isolated pixels
    iso_coords = np.argwhere(isolated_polar_interior)
    if len(iso_coords) == 0:
        logger.info("No isolated polar pixels found!")
        return

    # Sample some pixels
    np.random.seed(42)
    n_sample = min(20, len(iso_coords))
    sample_idx = np.random.choice(len(iso_coords), size=n_sample, replace=False)
    sample_coords = iso_coords[sample_idx]

    # Load Python thresholds
    import yaml
    thresholds_path = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
    with open(thresholds_path) as f:
        thresholds = yaml.safe_load(f)

    print(f"\n{'='*80}")
    print(f"Debug: Isolated Polar Cloudy Pixels (Fortran)")
    print(f"{'='*80}")
    print(f"Total isolated polar interior: {n_isolated}")
    print(f"Sampling {n_sample} pixels for detailed analysis\n")

    for idx, (i, j) in enumerate(sample_coords):
        lat_val = geo['lat'][i, j]
        lon_val = geo['lon'][i, j]
        sza_val = geo['sza'][i, j]
        vza_val = geo['vza'][i, j]
        eco_val = int(geo['eco_type'][i, j])
        lsf_val = int(geo['lsf'][i, j])
        elev_val = geo['elevation'][i, j]
        sfctmp_val = nwp_interp['tsfc'][i, j]

        cm_val = cm[i, j]
        conf_val = conf[i, j]

        # Get neighbor cm values
        neighbor_cm = [cm[i-1, j], cm[i+1, j], cm[i, j-1], cm[i, j+1]]

        # Classify pixel
        flags, _ = classify_pixel_surface(
            pxldat[i, j, :], lat_val, lon_val, elev_val, lsf_val,
            sza_val, vza_val, geo['glint_angle'][i, j], eco_val, 0,
            0.0,  # sst
            nwp_interp['tsfc'][i, j], nwp_interp['pmsl'][i, j],
            nwp_interp['u_wind'][i, j], nwp_interp['v_wind'][i, j],
            nwp_interp['tpw'][i, j],
            21,  # sensor_id
            thresholds
        )

        # Run Python pixel-level to get test details
        indat_3x3_11um = _extract_3x3(pxldat[:, :, 23], i, j, n_elem, n_line)
        indat_3x3_vis = _extract_3x3(pxldat[:, :, 2], i, j, n_elem, n_line)

        py_result = run_cloud_mask_pixel(
            pxldat=pxldat[i, j, :],
            lat=lat_val, lon=lon_val, elevation=elev_val,
            lsf=lsf_val, sza=sza_val, vza=vza_val,
            glint_angle=geo['glint_angle'][i, j],
            eco_type=eco_val, snow_mask_val=0, sst=0.0,
            nwp_sfctmp=sfctmp_val,
            nwp_pmsl=nwp_interp['pmsl'][i, j],
            nwp_u_wind=nwp_interp['u_wind'][i, j],
            nwp_v_wind=nwp_interp['v_wind'][i, j],
            nwp_precip_water=nwp_interp['tpw'][i, j],
            sensor_id=21,
            bt_clr=np.zeros(7, dtype=np.float64),
            thresholds=thresholds,
            indat_3x3_11um=indat_3x3_11um,
            indat_3x3_vis=indat_3x3_vis,
        )

        # Key band values
        masv66 = pxldat[i, j, 2]   # 0.64um
        masv88 = pxldat[i, j, 3]   # 0.86um
        masv138 = pxldat[i, j, 18]  # 1.38um (after swap, band 19 -> index 18)
        masir11 = pxldat[i, j, 23]  # 11um
        masir12 = pxldat[i, j, 24]  # 12um
        masir4 = pxldat[i, j, 19]   # 3.8um

        print(f"--- Pixel {idx+1}: ({i}, {j}) ---")
        print(f"  Location:  lat={lat_val:.2f}, lon={lon_val:.2f}")
        print(f"  Geometry:  sza={sza_val:.1f}, vza={vza_val:.1f}")
        print(f"  Surface:   eco={eco_val}, lsf={lsf_val}, elev={elev_val:.0f}m")
        print(f"  Temp:      sfctmp={sfctmp_val:.1f}K")
        print(f"  Flags:     polar={flags.polar}, land={flags.land}, snow={flags.snow}, "
              f"ice={flags.ice}, day={flags.day}")
        print(f"  Fortran:   cm={cm_val}, conf={conf_val:.4f}")
        print(f"  Python:    cm={py_result.cloud_mask}, conf={py_result.confidence:.4f}")
        print(f"  Neighbors: cm={neighbor_cm} (agreement={sum(n == cm_val for n in neighbor_cm)}/4)")
        print(f"  Bands:")
        print(f"    0.64um:  {masv66:.4f}  (threshold 50%: 0.18)")
        print(f"    0.86um:  {masv88:.4f}")
        print(f"    1.38um:  {masv138:.4f}  (threshold 50%: 0.18)")
        print(f"    3.8um:   {masir4:.2f}K")
        print(f"    11um:    {masir11:.2f}K")
        print(f"    12um:    {masir12:.2f}K")
        print(f"    11-12:   {masir11-masir12:.4f}K")
        print()

    # Summary statistics for all isolated polar pixels
    iso_lat = geo['lat'][isolated_polar_interior]
    iso_lon = geo['lon'][isolated_polar_interior]
    iso_conf = conf[isolated_polar_interior]
    iso_eco = geo['eco_type'][isolated_polar_interior]
    iso_lsf = geo['lsf'][isolated_polar_interior]
    iso_sza = geo['sza'][isolated_polar_interior]
    iso_v66 = pxldat[:, :, 2][isolated_polar_interior]
    iso_v138 = pxldat[:, :, 18][isolated_polar_interior]

    print(f"\n{'='*80}")
    print(f"Summary Statistics for All {n_isolated} Isolated Polar Pixels")
    print(f"{'='*80}")
    print(f"  Latitude:    {np.min(iso_lat):.1f} - {np.max(iso_lat):.1f} (mean {np.mean(iso_lat):.1f})")
    print(f"  Longitude:   {np.min(iso_lon):.1f} - {np.max(iso_lon):.1f}")
    print(f"  SZA:         {np.min(iso_sza):.1f} - {np.max(iso_sza):.1f} (mean {np.mean(iso_sza):.1f})")
    print(f"  Confidence:  {np.min(iso_conf):.4f} - {np.max(iso_conf):.4f} (mean {np.mean(iso_conf):.4f})")
    print(f"  0.64um ref:  {np.min(iso_v66):.4f} - {np.max(iso_v66):.4f} (mean {np.mean(iso_v66):.4f})")
    print(f"  1.38um ref:  {np.min(iso_v138):.4f} - {np.max(iso_v138):.4f} (mean {np.mean(iso_v138):.4f})")

    # Eco type distribution
    eco_unique, eco_counts = np.unique(iso_eco, return_counts=True)
    print(f"  Eco types:")
    for e, c in zip(eco_unique, eco_counts):
        print(f"    {e}: {c} ({100*c/n_isolated:.1f}%)")

    # LSF distribution
    lsf_unique, lsf_counts = np.unique(iso_lsf, return_counts=True)
    print(f"  LSF values:")
    for l, c in zip(lsf_unique, lsf_counts):
        print(f"    {l}: {c} ({100*c/n_isolated:.1f}%)")

    # Confidence distribution
    conf_bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.66]
    conf_hist, _ = np.histogram(iso_conf, bins=conf_bins)
    print(f"  Confidence distribution:")
    for k in range(len(conf_hist)):
        print(f"    [{conf_bins[k]:.2f}, {conf_bins[k+1]:.2f}): {conf_hist[k]} ({100*conf_hist[k]/n_isolated:.1f}%)")

    # Check how many have high visible reflectance
    high_ref = iso_v66 > 0.18
    print(f"  Pixels with 0.64um > 0.18: {int(np.sum(high_ref))} ({100*np.sum(high_ref)/n_isolated:.1f}%)")
    high_ref2 = iso_v66 > 0.22
    print(f"  Pixels with 0.64um > 0.22: {int(np.sum(high_ref2))} ({100*np.sum(high_ref2)/n_isolated:.1f}%)")

    print(f"\n{'='*80}")


if __name__ == '__main__':
    main()
