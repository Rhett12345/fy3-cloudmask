#!/usr/bin/env python3
"""Test recalibration integration with FY-3D data.

Processes FY-3D orbits with both onboard and recalibration coefficients,
then compares cloud mask results.

Usage:
    python scripts/test_recalibration.py
    python scripts/test_recalibration.py --date 20220803
    python scripts/test_recalibration.py --recal-dir ../fy3d_recali
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

# Setup paths
os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
from fy3_cloudmask.io.recalibration import RecalibrationManager
from fy3_cloudmask.output.writer import write_combined_product
from fy3_cloudmask.output.cloud_amount import compute_cloud_amount_with_coords

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# Data paths
MERSI_ROOT = Path('/data/Data_yuq/mersi')
NWP_ROOT = Path('/data/nwp')
RECAL_BASE = Path('/home/liusy2020/yuq/cloudmask/fy3d_recali')

# NWP config
NWP_PATTERN = 'gfs0p25_41L_{date}_{hh}_00'
NWP_HOURS = [0, 3, 6, 9, 12, 15, 18, 21]

# Import functions from run_fortran_only
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_fortran_only import read_l1b_data, read_geo_data, read_nwp_binary, interpolate_nwp


def find_nwp(date_str: str, obs_hour: int) -> Path | None:
    """Find the nearest NWP file for the given observation hour."""
    best_hh = None
    for hh in sorted(NWP_HOURS, reverse=True):
        if hh <= obs_hour:
            best_hh = hh
            break
    if best_hh is None:
        from datetime import timedelta
        prev_date = (datetime.strptime(date_str, '%Y%m%d') - timedelta(days=1)).strftime('%Y%m%d')
        best_hh = max(NWP_HOURS)
        date_str = prev_date

    fname = NWP_PATTERN.format(date=date_str, hh=f'{best_hh:02d}')
    nwp_path = NWP_ROOT / date_str / 'ORG' / fname
    return nwp_path if nwp_path.exists() else None


def find_test_orbits(date_str: str) -> list[dict]:
    """Find FY-3D test orbits for a given date."""
    mersi_dir = MERSI_ROOT / date_str
    if not mersi_dir.exists():
        logger.warning(f"MERSI directory not found: {mersi_dir}")
        return []

    orbits = []
    for l1b in sorted(mersi_dir.glob('FY3D_MERSI_GBAL_L1_*_1000M_MS.HDF')):
        stem = l1b.stem
        parts = stem.split('_')
        time_tag = parts[-3]  # e.g., '0245'

        # Find matching GEO file
        geo_pattern = stem.replace('1000M', 'GEO1K')
        geo = mersi_dir / (geo_pattern + '.HDF')
        if not geo.exists():
            logger.warning(f"GEO not found for {l1b.name}")
            continue

        # Parse observation hour
        obs_time = datetime.strptime(f"{date_str}_{time_tag}", '%Y%m%d_%H%M')
        nwp = find_nwp(date_str, obs_time.hour)
        if nwp is None:
            logger.warning(f"NWP not found for {l1b.name}")
            continue

        orbits.append({
            'l1b': l1b,
            'geo': geo,
            'nwp': nwp,
            'time_tag': time_tag,
        })

    return orbits


def process_and_save(
    l1b_path: str,
    geo_path: str,
    nwp_path: str,
    output_path: str,
    recal_cal0: np.ndarray | None = None,
    recal_cal1: np.ndarray | None = None,
    attributes: dict | None = None,
) -> dict:
    """Process one orbit and save HDF5 results. Returns stats dict."""
    logger.info(f"  Reading L1b data...")
    t0 = time.time()
    pxldat = read_l1b_data(l1b_path, recal_cal0=recal_cal0, recal_cal1=recal_cal1)
    t_l1b = time.time() - t0

    logger.info(f"  Reading GEO data...")
    geo = read_geo_data(geo_path)

    logger.info(f"  Reading NWP data...")
    nwp = read_nwp_binary(nwp_path)

    logger.info(f"  Interpolating NWP...")
    nwp_interp = interpolate_nwp(nwp, geo['lat'], geo['lon'])

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]

    logger.info(f"  Running Fortran native backend...")
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
    t_fortran = time.time() - t0

    # Extract results
    cm = result['cloud_mask']
    conf = result['confidence']
    cm_bitarray = result.get('cm_bitarray', np.zeros((n_elem, n_line, 6), dtype=np.uint8)).astype(np.uint8)
    qa_bitarray = result.get('qa_bitarray', np.zeros((n_elem, n_line, 10), dtype=np.uint8)).astype(np.uint8)

    # Compute 5km cloud amount
    qa_processed = np.ones((n_elem, n_line), dtype=np.int32)
    cloud_amount, cloud_amount_qa, lon_5km, lat_5km = compute_cloud_amount_with_coords(
        cm.astype(np.int32), qa_processed, geo['lon'], geo['lat'],
    )

    # Build attributes
    if attributes is None:
        attributes = {}
    attributes.update({
        'input_l1b': l1b_path,
        'input_geo': geo_path,
        'input_nwp': nwp_path,
    })

    # Save HDF5
    write_combined_product(
        filepath=output_path,
        cm_bitarray=cm_bitarray,
        cm_qa_bitarray=qa_bitarray,
        cm_tmp=cm.astype(np.int32),
        confidence=conf.astype(np.float64),
        cloud_amount=cloud_amount,
        cloud_amount_qa=cloud_amount_qa,
        lon=geo['lon'],
        lat=geo['lat'],
        lon_5km=lon_5km,
        lat_5km=lat_5km,
        attributes=attributes,
    )

    # Compute stats
    valid = (cm >= 0) & (cm <= 3)
    n_valid = int(np.sum(valid))
    n_total = n_elem * n_line

    cm_dist = {}
    for cat in range(4):
        cm_dist[int(cat)] = int(np.sum(cm[valid] == cat)) if n_valid > 0 else 0

    conf_valid = conf[valid] if n_valid > 0 else np.array([0.0])

    stats = {
        'n_elem': n_elem, 'n_line': n_line,
        'n_valid': n_valid,
        'cloud_mask_distribution': cm_dist,
        'confidence_mean': float(np.mean(conf_valid)),
        'confidence_std': float(np.std(conf_valid)),
        'timing_l1b': t_l1b,
        'timing_fortran': t_fortran,
    }

    return stats


def compare_results(onboard_path: str, recal_path: str) -> dict:
    """Compare onboard vs recalibration results from HDF5 files."""
    import h5py
    with h5py.File(onboard_path, 'r') as f:
        cm_on = f['Cloud_Mask_1km/Cloud_Mask_Value'][:]
        conf_on = f['Cloud_Mask_1km/Confidence'][:]
    with h5py.File(recal_path, 'r') as f:
        cm_re = f['Cloud_Mask_1km/Cloud_Mask_Value'][:]
        conf_re = f['Cloud_Mask_1km/Confidence'][:]

    valid = (cm_on >= 0) & (cm_on <= 3) & (cm_re >= 0) & (cm_re <= 3)
    n_valid = int(np.sum(valid))
    n_total = cm_on.size

    # Agreement rate
    agree = np.sum(cm_on[valid] == cm_re[valid]) if n_valid > 0 else 0
    agree_rate = agree / n_valid if n_valid > 0 else 0.0

    # Confidence difference
    conf_diff = np.abs(conf_on[valid] - conf_re[valid]) if n_valid > 0 else np.array([0.0])

    # Per-category changes
    changes = {}
    for old_cat in range(4):
        for new_cat in range(4):
            if old_cat != new_cat:
                n = int(np.sum((cm_on[valid] == old_cat) & (cm_re[valid] == new_cat)))
                if n > 0:
                    changes[f"{old_cat}->{new_cat}"] = n

    return {
        'n_total': n_total,
        'n_valid': n_valid,
        'agreement_rate': float(agree_rate),
        'confidence_diff_mean': float(np.mean(conf_diff)),
        'confidence_diff_max': float(np.max(conf_diff)),
        'category_changes': changes,
    }


def main():
    parser = argparse.ArgumentParser(description="Test recalibration with FY-3D data")
    parser.add_argument("--date", default="20220803",
                        help="Date to process (YYYYMMDD)")
    parser.add_argument("--recal-dir", default=str(RECAL_BASE),
                        help="Path to recalibration data directory")
    parser.add_argument("--output", default="/tmp/recal_test",
                        help="Output directory")
    parser.add_argument("--max-orbits", type=int, default=2,
                        help="Max orbits to process")
    args = parser.parse_args()

    if not is_native_available():
        logger.error("Native backend NOT available! Build with: cd ext && ./build.sh --install")
        sys.exit(1)

    # Limit OMP threads to avoid race conditions
    if 'OMP_NUM_THREADS' not in os.environ:
        os.environ['OMP_NUM_THREADS'] = '1'

    # Output goes to {output}/{date}/
    date_output = os.path.join(args.output, args.date)
    os.makedirs(date_output, exist_ok=True)

    # Find test orbits
    orbits = find_test_orbits(args.date)
    if not orbits:
        logger.error(f"No test orbits found for {args.date}")
        sys.exit(1)

    orbits = orbits[:args.max_orbits]
    logger.info(f"Found {len(orbits)} test orbits for {args.date}")

    # Load recalibration coefficients
    recal_mgr = RecalibrationManager(args.recal_dir)
    recal_cal0, recal_cal1 = recal_mgr.load_coefficients(args.date)
    logger.info(f"Recalibration loaded for {args.date}")
    logger.info(f"  cal0: {recal_cal0}")
    logger.info(f"  cal1: {recal_cal1}")

    results = []

    for orbit in orbits:
        time_tag = orbit['time_tag']
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing orbit: {time_tag}")
        logger.info(f"  L1b: {orbit['l1b'].name}")

        # Process with onboard calibration
        onboard_path = os.path.join(date_output, f"FY3D_MERSI_{args.date}_{time_tag}_CLM_CLA_onboard.h5")
        logger.info(f"  [1/2] Onboard calibration...")
        stats_on = process_and_save(
            str(orbit['l1b']), str(orbit['geo']), str(orbit['nwp']),
            onboard_path,
            attributes={'backend': 'Fortran+OpenMP (onboard calibration)'},
        )

        # Process with recalibration
        recal_path = os.path.join(date_output, f"FY3D_MERSI_{args.date}_{time_tag}_CLM_CLA_recal.h5")
        logger.info(f"  [2/2] Recalibration...")
        stats_re = process_and_save(
            str(orbit['l1b']), str(orbit['geo']), str(orbit['nwp']),
            recal_path,
            recal_cal0=recal_cal0,
            recal_cal1=recal_cal1,
            attributes={'backend': 'Fortran+OpenMP (recalibration)'},
        )

        # Compare
        comparison = compare_results(onboard_path, recal_path)

        result = {
            'orbit': time_tag,
            'onboard': stats_on,
            'recalibration': stats_re,
            'comparison': comparison,
        }
        results.append(result)

        logger.info(f"  Comparison:")
        logger.info(f"    Agreement rate: {comparison['agreement_rate']:.4f}")
        logger.info(f"    Confidence diff (mean): {comparison['confidence_diff_mean']:.4f}")
        logger.info(f"    Category changes: {comparison['category_changes']}")

    # Save summary
    summary_path = os.path.join(date_output, 'comparison_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary table
    print(f"\n{'='*80}")
    print(f"Recalibration Test Summary - {args.date}")
    print(f"{'='*80}")
    print(f"{'Orbit':<10} {'Agreement':>10} {'ConfDiff':>10} {'OnClear%':>10} {'ReClear%':>10} {'Changes':>15}")
    print(f"{'-'*80}")

    for r in results:
        on = r['onboard']
        re = r['recalibration']
        comp = r['comparison']
        on_clear = on['cloud_mask_distribution'].get(3, 0) / on['n_valid'] * 100 if on['n_valid'] > 0 else 0
        re_clear = re['cloud_mask_distribution'].get(3, 0) / re['n_valid'] * 100 if re['n_valid'] > 0 else 0
        changes = comp['category_changes']
        change_str = ', '.join(f"{k}:{v}" for k, v in changes.items()) if changes else "none"
        print(f"{r['orbit']:<10} {comp['agreement_rate']:>10.4f} {comp['confidence_diff_mean']:>10.4f} "
              f"{on_clear:>10.2f} {re_clear:>10.2f} {change_str:>15}")

    print(f"{'='*80}")
    print(f"Full results saved to: {date_output}")


if __name__ == '__main__':
    main()
