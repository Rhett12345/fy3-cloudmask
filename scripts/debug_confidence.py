#!/usr/bin/env python3
"""Debug confidence calculation for FY-3D cloud mask.

Processes a small subset of pixels and prints intermediate confidence values
to understand why confidence is becoming binary (0.0 or 1.0).
"""

import os
import sys
import numpy as np
from pathlib import Path

# Setup paths
os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_pixel, classify_pixel_surface
from fy3_cloudmask.algorithm.confidence import conf_test, encode_confidence
from fy3_cloudmask.constants import BAD_DATA, IR_11, IR_12, IR_38, BAND_064, BAND_138

# Import data reading functions
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp

import yaml
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def load_thresholds():
    """Load thresholds from YAML file."""
    thresholds_path = Path(__file__).resolve().parent.parent / 'config' / 'thresholds' / 'mersi_ii3d_v8.yaml'
    if thresholds_path.exists():
        with open(thresholds_path) as f:
            return yaml.safe_load(f)
    else:
        logger.warning(f"Threshold file not found: {thresholds_path}")
        return {}


def debug_pixel_confidence(pxldat, lat, lon, elevation, lsf, sza, vza, glint_angle,
                          eco_type, snow_mask_val, sst, nwp_sfctmp, nwp_pmsl,
                          nwp_u_wind, nwp_v_wind, nwp_precip_water, bt_clr, thresholds):
    """Debug confidence calculation for a single pixel."""
    # Classify surface type
    flags, pxldat = classify_pixel_surface(
        pxldat, lat, lon, elevation, lsf, sza, vza, glint_angle,
        eco_type, snow_mask_val, sst, nwp_sfctmp, nwp_pmsl,
        nwp_u_wind, nwp_v_wind, nwp_precip_water, 21, thresholds,
    )

    print(f"  Surface flags: day={flags.day}, polar={flags.polar}, land={flags.land}, "
          f"snow={flags.snow}, ice={flags.ice}, water={flags.water}, "
          f"desert={flags.desert}, coast={flags.coast}")

    if flags.bad_value or flags.bad_geo:
        print("  Bad value/geo - returning confidence=0.0")
        return 0.0, 0, 0

    # Get thresholds for polar day land
    thr = thresholds.get('land_day_polar', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    print(f"  Threshold keys available: land_day_polar={bool(thr)}, pfmft={bool(pfmft_thr)}, nfmft={bool(nfmft_thr)}")

    # Extract band values
    masv66 = pxldat[BAND_064]
    masv188 = pxldat[BAND_138]
    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    print(f"  Band values: VIS0.64={masv66:.4f}, VIS1.38={masv188:.4f}, "
          f"IR3.8={masir4:.2f}, IR11={masir11:.2f}, IR12={masir12:.2f}")

    # Initialize confidence values
    cmin1 = 1.0  # Group 1: IR forward model
    cmin2 = 1.0  # Group 2: BTD tests
    cmin3 = 1.0  # Group 3: Visible tests
    cmin4 = 1.0  # Group 4: NIR test

    ngtests = [0, 0, 0, 0]

    # GROUP 1: PFMFT + NFMFT
    print("\n  GROUP 1: IR Forward Model Tests")

    # PFMFT test
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            masir11 < pfmft_thr.get('bt_11_max', 310.0) and
            (bt_clr[4] - bt_clr[5]) > pfmft_thr.get('btd_min', 0.0)):
        if masir11 > 270.0 and bt_clr[4] > 270.0:
            tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5]) * (masir11 - 260.0) / (bt_clr[4] - 260.0)
        else:
            tv11_12 = masir11 - masir12

        if is_cold_sfc == 1:
            pfmft_cold = pfmft_thr.get('cold', [2.0, 1.5, 1.0, 1.0])
            c1 = conf_test(tv11_12, pfmft_cold[0], pfmft_cold[1], pfmft_cold[2], pfmft_cold[3])
        else:
            pfmft_land = pfmft_thr.get('land', [4.0, 3.5, 3.0, 1.0])
            c1 = conf_test(tv11_12, pfmft_land[0], pfmft_land[1], pfmft_land[2], pfmft_land[3])

        print(f"    PFMFT: tv11_12={tv11_12:.4f}, confidence={c1:.4f}")
        cmin1 = min(cmin1, c1)
        ngtests[0] += 1
    else:
        print(f"    PFMFT: skipped (masir11={masir11:.2f}, masir12={masir12:.2f}, bt_clr_diff={bt_clr[4]-bt_clr[5]:.2f})")

    # NFMFT test
    nfmft_max = nfmft_thr.get('max_threshold', 1.50)
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            (masir11 - masir12) <= nfmft_max):
        tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5])
        nfmft_land = nfmft_thr.get('land', [-23.0, -22.5, -22.0, 1.0])
        c2 = conf_test(tv11_12, nfmft_land[0], nfmft_land[1], nfmft_land[2], nfmft_land[3])

        print(f"    NFMFT: tv11_12={tv11_12:.4f}, confidence={c2:.4f}")
        cmin1 = min(cmin1, c2)
        ngtests[0] += 1
    else:
        print(f"    NFMFT: skipped (masir11-masir12={masir11-masir12:.4f} > {nfmft_max})")

    # GROUP 2: BTD tests
    print("\n  GROUP 2: BTD Tests")

    # 11-12um BTD test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and vza > 0.0:
        masdf1 = masir11 - masir12
        import math
        cosvza = math.cos(vza * math.pi / 180.0)
        schi = 1.0 / cosvza if abs(cosvza) > 1e-6 else 99.0

        from fy3_cloudmask.algorithm.spatial import tview
        diftemp = tview(schi, masir11)
        btd_11_12 = thr.get('btd_11_12', [3.0])
        dfthrsh = btd_11_12[0] if (diftemp < 0.1 or abs(schi - 99.0) < 0.0001) else diftemp

        locut = dfthrsh
        midpt = dfthrsh - 0.3 * dfthrsh
        hicut = midpt - 1.25

        c5 = conf_test(masdf1, locut, midpt, hicut, 1.0)
        print(f"    BTD 11-12: masdf1={masdf1:.4f}, dfthrsh={dfthrsh:.4f}, "
              f"locut={locut:.4f}, midpt={midpt:.4f}, hicut={hicut:.4f}, confidence={c5:.4f}")
        cmin2 = min(cmin2, c5)
        ngtests[1] += 1
    else:
        print(f"    BTD 11-12: skipped (masir11={masir11:.2f}, masir12={masir12:.2f}, vza={vza:.2f})")

    # 11-4um BTD test
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        mas11_4 = masir11 - masir4
        btd_11_4 = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])
        c3 = conf_test(mas11_4, btd_11_4[0], btd_11_4[1], btd_11_4[2], btd_11_4[3])
        print(f"    BTD 11-4: mas11_4={mas11_4:.4f}, thresholds={btd_11_4}, confidence={c3:.4f}")
        cmin2 = min(cmin2, c3)
        ngtests[1] += 1
    else:
        print(f"    BTD 11-4: skipped (masir11={masir11:.2f}, masir4={masir4:.2f})")

    # GROUP 3: Visible tests
    print("\n  GROUP 3: Visible Tests")
    visusd = flags.visusd
    if visusd:
        # 0.64um reflectance test
        VIS_VALID_MIN = -99.0
        VIS_VALID_MAX = 2.3
        if masv66 > VIS_VALID_MIN and masv66 <= VIS_VALID_MAX:
            ref064 = thr.get('ref064', [0.24, 0.20, 0.16, 1.0])
            c4 = conf_test(masv66, ref064[0], ref064[1], ref064[2], ref064[3])
            print(f"    VIS 0.64: masv66={masv66:.4f}, thresholds={ref064}, confidence={c4:.4f}")
            cmin3 = min(cmin3, c4)
            ngtests[2] += 1
        else:
            print(f"    VIS 0.64: skipped (masv66={masv66:.4f} out of range)")
    else:
        print(f"    VIS tests: skipped (visusd={visusd})")

    # GROUP 4: NIR test (1.38um)
    print("\n  GROUP 4: NIR Test (1.38um)")
    hi_elev = flags.hi_elev
    if not hi_elev and masv188 > VIS_VALID_MIN and masv188 <= VIS_VALID_MAX:
        ref138 = thr.get('ref138', [0.04, 0.035, 0.03, 1.0])
        c7 = conf_test(masv188, ref138[0], ref138[1], ref138[2], ref138[3])
        print(f"    NIR 1.38: masv188={masv188:.4f}, thresholds={ref138}, confidence={c7:.4f}")
        cmin4 = min(cmin4, c7)
        ngtests[3] += 1
    else:
        print(f"    NIR 1.38: skipped (masv188={masv188:.4f}, hi_elev={hi_elev})")

    # Final confidence calculation
    print("\n  Final Confidence Calculation:")
    print(f"    Group minimums: cmin1={cmin1:.4f}, cmin2={cmin2:.4f}, cmin3={cmin3:.4f}, cmin4={cmin4:.4f}")
    print(f"    Tests per group: {ngtests}")

    active_groups = sum(1 for g in ngtests if g > 0)
    if active_groups > 0:
        product = 1.0
        for i in range(4):
            if ngtests[i] > 0:
                product *= [cmin1, cmin2, cmin3, cmin4][i]
        confdnc = product ** (1.0 / active_groups)
    else:
        confdnc = 1.0

    print(f"    Active groups: {active_groups}, Final confidence: {confdnc:.4f}")

    # Encode confidence
    bit1, bit2 = encode_confidence(confdnc)
    print(f"    Encoded bits: ({bit1}, {bit2})")

    if confdnc > 0.99:
        cm_class = "Confident Clear"
    elif confdnc > 0.95:
        cm_class = "Prob Clear"
    elif confdnc > 0.66:
        cm_class = "Prob Cloudy"
    else:
        cm_class = "Cloudy"

    print(f"    Cloud mask class: {cm_class}")

    return confdnc, active_groups, sum(ngtests)


def main():
    """Main debug function."""
    # Data paths
    l1b_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_1000M_MS.HDF'
    geo_path = '/data/Data_yuq/mersi/20220803/FY3D_MERSI_GBAL_L1_20220803_0740_GEO1K_MS.HDF'
    nwp_path = '/data/nwp/20220803/ORG/gfs0p25_41L_20220803_06_00'

    if not os.path.exists(l1b_path):
        logger.error(f"L1b file not found: {l1b_path}")
        return

    logger.info("Loading thresholds...")
    thresholds = load_thresholds()

    logger.info("Reading L1b data...")
    pxldat_full = read_l1b_data(l1b_path)

    logger.info("Reading GEO data...")
    geo = read_geo_data(geo_path)

    # Read LandSeaMask from GEO file
    import h5py
    with h5py.File(geo_path, 'r') as f:
        lsf_raw = f['Geolocation/LandSeaMask'][:].astype(np.int32)
    geo['lsf'] = lsf_raw.T  # Transpose to match other arrays

    logger.info("Reading NWP data...")
    nwp = read_nwp_binary(nwp_path)

    logger.info("Interpolating NWP...")
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat_full.shape[0], pxldat_full.shape[1]
    logger.info(f"Swath dimensions: {n_elem} x {n_line}")

    # Debug a sample of pixels
    sample_pixels = [
        (1024, 1000),  # Center of swath
        (1024, 500),   # Quarter from start
        (1024, 1500),  # Quarter from end
        (512, 1000),   # Left side
        (1536, 1000),  # Right side
    ]

    logger.info(f"\n{'='*80}")
    logger.info("Debugging confidence calculation for sample pixels")
    logger.info(f"{'='*80}")

    for i, (col, row) in enumerate(sample_pixels):
        logger.info(f"\nPixel {i+1}: ({col}, {row})")
        logger.info(f"{'='*60}")

        pxldat = pxldat_full[col, row, :]
        lat = geo['lat'][col, row]
        lon = geo['lon'][col, row]
        elevation = geo['elevation'][col, row]
        lsf = int(geo['lsf'][col, row])
        sza = geo['sza'][col, row]
        vza = geo['vza'][col, row]
        glint_angle = geo['glint_angle'][col, row]
        eco_type = int(geo['eco_type'][col, row])
        snow_mask_val = 0
        sst = 0.0
        nwp_sfctmp = nwp_interp['tsfc'][col, row]
        nwp_pmsl = nwp_interp['pmsl'][col, row]
        nwp_u_wind = nwp_interp['u_wind'][col, row]
        nwp_v_wind = nwp_interp['v_wind'][col, row]
        nwp_precip_water = nwp_interp['tpw'][col, row]
        bt_clr = np.zeros(7, dtype=np.float32)  # Placeholder

        logger.info(f"  Location: lat={lat:.2f}, lon={lon:.2f}, elev={elevation:.0f}m")
        logger.info(f"  Geometry: sza={sza:.2f}, vza={vza:.2f}")
        logger.info(f"  NWP: sfctmp={nwp_sfctmp:.2f}K, pmsl={nwp_pmsl:.2f}hPa")
        logger.info(f"  Land/Sea flag: {lsf}")

        confdnc, active_groups, total_tests = debug_pixel_confidence(
            pxldat, lat, lon, elevation, lsf, sza, vza, glint_angle,
            eco_type, snow_mask_val, sst, nwp_sfctmp, nwp_pmsl,
            nwp_u_wind, nwp_v_wind, nwp_precip_water, bt_clr, thresholds
        )

    # Summary statistics
    logger.info(f"\n{'='*80}")
    logger.info("Summary Statistics")
    logger.info(f"{'='*80}")

    # Process a larger sample to get statistics
    n_sample = 1000
    logger.info(f"Processing {n_sample} random pixels for statistics...")

    confidences = []
    for _ in range(n_sample):
        col = np.random.randint(0, n_elem)
        row = np.random.randint(0, n_line)

        pxldat = pxldat_full[col, row, :]
        lat = geo['lat'][col, row]
        lon = geo['lon'][col, row]
        elevation = geo['elevation'][col, row]
        lsf = int(geo['lsf'][col, row])
        sza = geo['sza'][col, row]
        vza = geo['vza'][col, row]
        glint_angle = geo['glint_angle'][col, row]
        eco_type = int(geo['eco_type'][col, row])
        snow_mask_val = 0
        sst = 0.0
        nwp_sfctmp = nwp_interp['tsfc'][col, row]
        nwp_pmsl = nwp_interp['pmsl'][col, row]
        nwp_u_wind = nwp_interp['u_wind'][col, row]
        nwp_v_wind = nwp_interp['v_wind'][col, row]
        nwp_precip_water = nwp_interp['tpw'][col, row]
        bt_clr = np.zeros(7, dtype=np.float32)

        try:
            confdnc, _, _ = debug_pixel_confidence(
                pxldat, lat, lon, elevation, lsf, sza, vza, glint_angle,
                eco_type, snow_mask_val, sst, nwp_sfctmp, nwp_pmsl,
                nwp_u_wind, nwp_v_wind, nwp_precip_water, bt_clr, thresholds
            )
            confidences.append(confdnc)
        except Exception as e:
            logger.warning(f"Error processing pixel ({col}, {row}): {e}")

    if confidences:
        confidences = np.array(confidences)
        logger.info(f"\nConfidence Statistics:")
        logger.info(f"  Total pixels processed: {len(confidences)}")
        logger.info(f"  Mean confidence: {np.mean(confidences):.4f}")
        logger.info(f"  Std confidence: {np.std(confidences):.4f}")
        logger.info(f"  Min confidence: {np.min(confidences):.4f}")
        logger.info(f"  Max confidence: {np.max(confidences):.4f}")

        # Distribution by class
        logger.info(f"\nCloud Mask Class Distribution:")
        logger.info(f"  Cloudy (conf <= 0.66): {np.sum(confidences <= 0.66)} ({100*np.sum(confidences <= 0.66)/len(confidences):.2f}%)")
        logger.info(f"  Prob Cloudy (0.66 < conf <= 0.95): {np.sum((confidences > 0.66) & (confidences <= 0.95))} ({100*np.sum((confidences > 0.66) & (confidences <= 0.95))/len(confidences):.2f}%)")
        logger.info(f"  Prob Clear (0.95 < conf <= 0.99): {np.sum((confidences > 0.95) & (confidences <= 0.99))} ({100*np.sum((confidences > 0.95) & (confidences <= 0.99))/len(confidences):.2f}%)")
        logger.info(f"  Confident Clear (conf > 0.99): {np.sum(confidences > 0.99)} ({100*np.sum(confidences > 0.99)/len(confidences):.2f}%)")

        # Histogram
        logger.info(f"\nConfidence Histogram:")
        bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.66, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]
        for i in range(len(bins)-1):
            count = np.sum((confidences >= bins[i]) & (confidences < bins[i+1]))
            if count > 0:
                logger.info(f"  [{bins[i]:.2f}, {bins[i+1]:.2f}): {count} ({100*count/len(confidences):.2f}%)")


if __name__ == '__main__':
    main()
