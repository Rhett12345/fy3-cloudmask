#!/usr/bin/env python3
"""Compare FY-3D CLM with MYD35 cloud mask via lat/lon collocation.

Reads MYD35 HDF granules and FY-3D CLM output for matching time windows,
collocates by nearest-neighbor lat/lon, and reports confusion matrix.
"""

import sys
from pathlib import Path
import h5py
import numpy as np
from pyhdf.SD import SD, SDC

# --- Paths ---
FY3_CLM_DIR = Path('/data/Data_yuq/fy3_cloud/20200308_f90')
MYD35_DIR = Path('/data/Data_yuq/aqua_modis/MYD35_L2/20200308')

# FY-3D orbit times (UTC) and corresponding MYD35 time windows
ORBITS = {
    '1345': {'fy3_time': '1345', 'myd35_start': '1325', 'myd35_end': '1440'},
    '1435': {'fy3_time': '1435', 'myd35_start': '1415', 'myd35_end': '1530'},
    '1525': {'fy3_time': '1525', 'myd35_start': '1505', 'myd35_end': '1620'},
}


def read_myd35_cloud_mask(hdf_path):
    """Read Cloud_Mask and lat/lon from MYD35 HDF4."""
    sd = SD(str(hdf_path), SDC.READ)
    cm = sd.select('Cloud_Mask')[:]  # (nlines, npixels) uint16
    lat = sd.select('Latitude')[:]
    lon = sd.select('Longitude')[:]
    sd.end()

    # MYD35 Cloud_Mask bit0 = cloud mask flag (0=determined)
    # bit1-2 = 00=cloudy, 01=uncertain_clear, 10=prob_clear, 11=conf_clear
    # Extract bits 1-2: (cm >> 1) & 3
    cm_4class = (cm >> 1) & 3  # 0=cloudy, 1=uncertain, 2=prob_clear, 3=conf_clear
    return cm_4class.astype(np.int8), lat, lon



def read_fy3_clm(clm_path):
    """Read FY-3D CLM output."""
    with h5py.File(clm_path, 'r') as f:
        cm = f['cm'][:]
        conf = f['conf'][:]
        lat = f['lat'][:]
        lon = f['lon'][:]
    return cm, conf, lat, lon


def collocate_nn(lat1, lon1, lat2, lon2, max_dist_deg=0.1):
    """Find nearest-neighbor matches between two grids.

    Returns indices in grid2 for each pixel in grid1.
    """
    from scipy.spatial import cKDTree
    # Convert to 3D Cartesian
    d2r = np.pi / 180.0
    x1 = np.cos(lat1 * d2r) * np.cos(lon1 * d2r)
    y1 = np.cos(lat1 * d2r) * np.sin(lon1 * d2r)
    z1 = np.sin(lat1 * d2r)

    x2 = np.cos(lat2 * d2r) * np.cos(lon2 * d2r)
    y2 = np.cos(lat2 * d2r) * np.sin(lon2 * d2r)
    z2 = np.sin(lat2 * d2r)

    tree = cKDTree(np.column_stack([x2.ravel(), y2.ravel(), z2.ravel()]))
    pts = np.column_stack([x1.ravel(), y1.ravel(), z1.ravel()])
    dist, idx = tree.query(pts, k=1)

    # chord distance → angular distance
    ang_dist = 2 * np.arcsin(np.clip(dist / 2, 0, 1))
    ang_dist_deg = ang_dist * 180.0 / np.pi

    valid = ang_dist_deg < max_dist_deg
    return idx, ang_dist_deg, valid


def confusion_stats(cm_fy3, cm_myd35, valid_mask):
    """Compute confusion matrix and agreement stats."""
    fy3 = cm_fy3.ravel()[valid_mask]
    myd = cm_myd35.ravel()[valid_mask]

    # Collapse to binary: 0,1=cloudy, 2,3=clear
    fy3_bin = (fy3 >= 2).astype(int)
    myd_bin = (myd >= 2).astype(int)

    n = len(fy3)
    if n == 0:
        return {}

    # 4-class confusion
    conf4 = np.zeros((4, 4), dtype=int)
    for i in range(4):
        for j in range(4):
            conf4[i, j] = np.sum((fy3 == i) & (myd == j))

    # Binary confusion: rows=FY3, cols=MYD35
    # FY3 cloudy (0) vs clear (1), MYD35 cloudy (0) vs clear (1)
    tp = np.sum((fy3_bin == 1) & (myd_bin == 1))  # both clear
    tn = np.sum((fy3_bin == 0) & (myd_bin == 0))  # both cloudy
    fp = np.sum((fy3_bin == 1) & (myd_bin == 0))  # FY3 clear, MYD cloudy
    fn = np.sum((fy3_bin == 0) & (myd_bin == 1))  # FY3 cloudy, MYD clear

    total = tp + tn + fp + fn
    accuracy = (tp + tn) / total * 100 if total > 0 else 0
    pod_cloudy = tn / (tn + fp) * 100 if (tn + fp) > 0 else 0  # MYD cloudy → FY3 cloudy
    pod_clear = tp / (tp + fn) * 100 if (tp + fn) > 0 else 0    # MYD clear → FY3 clear

    # FY3 cloud fraction vs MYD35 cloud fraction
    fy3_cf = np.sum(fy3_bin == 0) / n * 100
    myd_cf = np.sum(myd_bin == 0) / n * 100

    return {
        'n_matched': n,
        'conf4': conf4,
        'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn,
        'accuracy': accuracy,
        'pod_cloudy': pod_cloudy,
        'pod_clear': pod_clear,
        'fy3_cloud_fraction': fy3_cf,
        'myd35_cloud_fraction': myd_cf,
        'fy3_4class_pct': [np.sum(fy3 == i) / n * 100 for i in range(4)],
        'myd35_4class_pct': [np.sum(myd == i) / n * 100 for i in range(4)],
    }


def main():
    for orbit_tag, cfg in ORBITS.items():
        # Find FY-3D files
        fy3_name = f'FY3D_MERSI_20200308_{orbit_tag}_CLM_CLA'
        clm_path = FY3_CLM_DIR / f'{fy3_name}.h5'

        if not clm_path.exists():
            print(f"[{orbit_tag}] CLM file not found: {clm_path}")
            continue

        print(f"\n{'='*60}")
        print(f"  Orbit {orbit_tag}")
        print(f"{'='*60}")

        # Read FY-3D
        fy3_cm, fy3_conf, fy3_lat, fy3_lon = read_fy3_clm(clm_path)
        print(f"  FY-3D: {fy3_cm.shape}, cloud fraction: {np.sum(fy3_cm <= 1) / fy3_cm.size * 100:.1f}%")

        # Find matching MYD35 granules
        myd35_files = sorted(MYD35_DIR.glob('MYD35_L2.A2020068.*.hdf'))
        matching = []
        for fp in myd35_files:
            # Extract HHMM from filename: MYD35_L2.A2020068.HHMM...
            parts = fp.name.split('.')
            if len(parts) >= 3:
                hhmm = parts[2]
                if cfg['myd35_start'] <= hhmm <= cfg['myd35_end']:
                    matching.append(fp)

        if not matching:
            print(f"  No matching MYD35 granules found in [{cfg['myd35_start']}, {cfg['myd35_end']}]")
            continue

        print(f"  Matching MYD35 granules: {len(matching)}")

        # Read all MYD35 granules, flatten into single arrays
        myd_cm_list = []
        myd_lat_list = []
        myd_lon_list = []
        for fp in matching:
            try:
                cm, lat, lon = read_myd35_cloud_mask(fp)
                myd_cm_list.append(cm.ravel())
                myd_lat_list.append(lat.ravel())
                myd_lon_list.append(lon.ravel())
            except Exception as e:
                print(f"  Warning: failed to read {fp.name}: {e}")

        if not myd_cm_list:
            continue

        myd_cm_all = np.concatenate(myd_cm_list)
        myd_lat_all = np.concatenate(myd_lat_list)
        myd_lon_all = np.concatenate(myd_lon_list)
        print(f"  MYD35 combined: {myd_cm_all.shape}")

        # Collocate: subsample FY-3D for speed (every 4th pixel)
        fy3_lat_sub = fy3_lat[::4, ::4]
        fy3_lon_sub = fy3_lon[::4, ::4]
        fy3_cm_sub = fy3_cm[::4, ::4]
        fy3_conf_sub = fy3_conf[::4, ::4]
        print(f"  FY-3D subsampled: {fy3_cm_sub.shape}")

        idx, dist, valid = collocate_nn(fy3_lat_sub, fy3_lon_sub, myd_lat_all, myd_lon_all, max_dist_deg=0.05)
        n_valid = np.sum(valid)
        print(f"  Collocated pixels (within 0.05 deg): {n_valid:,} ({n_valid / valid.size * 100:.1f}%)")

        if n_valid < 100:
            print("  Too few matches, skipping.")
            continue

        # Get matched MYD35 values
        myd_matched = myd_cm_all.ravel()[idx[valid]]
        fy3_matched = fy3_cm_sub.ravel()[valid]

        # Compute stats
        stats = confusion_stats(fy3_matched, myd_matched, np.ones(len(fy3_matched), dtype=bool))

        print(f"\n  --- Binary Comparison (0,1=cloudy / 2,3=clear) ---")
        print(f"  Matched pixels:     {stats['n_matched']:,}")
        print(f"  FY-3D cloud frac:   {stats['fy3_cloud_fraction']:.1f}%")
        print(f"  MYD35 cloud frac:   {stats['myd35_cloud_fraction']:.1f}%")
        print(f"  Overall Accuracy:   {stats['accuracy']:.2f}%")
        print(f"  POD cloudy:         {stats['pod_cloudy']:.2f}%")
        print(f"  POD clear:          {stats['pod_clear']:.2f}%")
        print(f"  TP(both clear)={stats['tp']:,}  TN(both cloudy)={stats['tn']:,}")
        print(f"  FP(FY3 clear,MYD cloudy)={stats['fp']:,}  FN(FY3 cloudy,MYD clear)={stats['fn']:,}")

        print(f"\n  --- 4-Class Distribution ---")
        labels = ['Cloudy', 'Prob Cloudy', 'Prob Clear', 'Conf Clear']
        print(f"  {'Class':<16s} {'FY-3D':>8s} {'MYD35':>8s}")
        for i in range(4):
            print(f"  {labels[i]:<16s} {stats['fy3_4class_pct'][i]:>7.1f}% {stats['myd35_4class_pct'][i]:>7.1f}%")

        print(f"\n  --- 4-Class Confusion Matrix (rows=FY3D, cols=MYD35) ---")
        print(f"  {'':>12s} {'M:Cloudy':>8s} {'M:PCloudy':>8s} {'M:PClear':>8s} {'M:Clear':>8s}")
        fy3_labels = ['F:Cloudy', 'F:PCloudy', 'F:PClear', 'F:Clear']
        for i in range(4):
            row = stats['conf4'][i]
            print(f"  {fy3_labels[i]:>12s} {row[0]:>8d} {row[1]:>8d} {row[2]:>8d} {row[3]:>8d}")


if __name__ == '__main__':
    main()
