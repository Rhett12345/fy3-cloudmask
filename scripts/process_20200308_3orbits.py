#!/usr/bin/env python3
"""Process all 3 FY-3D MERSI orbits for 2020-03-08 and compare with MYD35."""
import os, sys, time
from pathlib import Path
import h5py
import numpy as np
import cfgrib
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
from pyhdf.SD import SD, SDC

os.environ['FY3_CODE_ROOT'] = str(Path(__file__).resolve().parent.parent / 'coeff') + '/'
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'src'))

from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
from run_fortran_only import read_l1b_data, read_geo_data

# ---- Orbits ----
ORBITS = [
    {
        'name': '1345',
        'l1b': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_1000M_MS.HDF',
        'geo': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1345_GEO1K_MS.HDF',
        'nwp': '/data/nwp/20200308/ORG/fnl_20200308_12_00.grib2',
        'myd35': '/data/Data_yuq/aqua_modis/MYD35_L2/20200308/MYD35_L2.A2020068.1345.061.2020069151204.hdf',
    },
    {
        'name': '1435',
        'l1b': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1435_1000M_MS.HDF',
        'geo': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1435_GEO1K_MS.HDF',
        'nwp': '/data/nwp/20200308/ORG/fnl_20200308_12_00.grib2',
        'myd35': '/data/Data_yuq/aqua_modis/MYD35_L2/20200308/MYD35_L2.A2020068.1435.061.2020069151007.hdf',
    },
    {
        'name': '1525',
        'l1b': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1525_1000M_MS.HDF',
        'geo': '/data/Data_yuq/mersi/20200308/FY3D_MERSI_GBAL_L1_20200308_1525_GEO1K_MS.HDF',
        'nwp': '/data/nwp/20200308/ORG/fnl_20200308_12_00.grib2',
        'myd35': '/data/Data_yuq/aqua_modis/MYD35_L2/20200308/MYD35_L2.A2020068.1525.061.2020069151715.hdf',
    },
]

OUT_DIR = '/data/Data_yuq/fy3_cloud/20200308'
os.makedirs(OUT_DIR, exist_ok=True)

if not is_native_available():
    print("ERROR: Native backend not available!")
    sys.exit(1)


def load_nwp(grib_path, pixel_lat, pixel_lon):
    """Load NWP from GRIB2 and interpolate to pixel grid."""
    ds = cfgrib.open_dataset(grib_path, backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}})
    grib_lat = ds.latitude.values
    grib_lon = ds.longitude.values
    tsfc_k = ds['t'].values

    lon_360 = np.where(pixel_lon < 0, pixel_lon + 360, pixel_lon)

    # reverse lat if descending
    if grib_lat[0] > grib_lat[-1]:
        lat_asc = grib_lat[::-1]
        field_asc = tsfc_k[::-1, :]
    else:
        lat_asc = grib_lat
        field_asc = tsfc_k

    interp = RegularGridInterpolator((lat_asc, grib_lon), field_asc, bounds_error=False, fill_value=None)
    points = np.stack([pixel_lat.ravel(), lon_360.ravel()], axis=-1)
    sfctmp = interp(points).reshape(pixel_lat.shape).astype(np.float32)

    # Fill other fields with zeros
    uwind = np.zeros_like(sfctmp)
    vwind = np.zeros_like(sfctmp)
    pmsl = np.full_like(sfctmp, 1013.25)
    tpw = np.zeros_like(sfctmp)

    return sfctmp, uwind, vwind, pmsl, tpw


def process_orbit(orbit):
    """Process one orbit and return cloud mask + geolocation."""
    name = orbit['name']
    print(f"\n{'='*60}")
    print(f"Processing orbit {name}")
    print(f"{'='*60}")

    t0 = time.time()
    pxldat = read_l1b_data(orbit['l1b'])
    geo = read_geo_data(orbit['geo'])
    print(f"  Read L1b/GEO: {time.time()-t0:.1f}s")

    n_elem, n_line = pxldat.shape[0], pxldat.shape[1]
    print(f"  Swath: {n_elem} x {n_line} = {n_elem*n_line:,} pixels")

    # NWP
    t0 = time.time()
    pixel_lat = geo['lat']
    sfctmp, uwind, vwind, pmsl, tpw = load_nwp(orbit['nwp'], pixel_lat, geo['lon'])
    print(f"  NWP interp: {time.time()-t0:.1f}s (sfctmp: {sfctmp.min():.0f}-{sfctmp.max():.0f}K)")

    # Fortran
    t0 = time.time()
    result = process_swath_native(
        ref_vis=np.ascontiguousarray(pxldat[:, :, :19].astype(np.float32)),
        tbb_ir=np.ascontiguousarray(pxldat[:, :, 19:].astype(np.float32)),
        lat=np.ascontiguousarray(pixel_lat.astype(np.float32)),
        lon=np.ascontiguousarray(geo['lon'].astype(np.float32)),
        satzen=np.ascontiguousarray(geo['vza'].astype(np.float32)),
        solzen=np.ascontiguousarray(geo['sza'].astype(np.float32)),
        relaz=np.ascontiguousarray(np.zeros_like(geo['sza'], dtype=np.float32)),
        glint=np.ascontiguousarray(geo['glint_angle'].astype(np.float32)),
        sfctmp=np.ascontiguousarray(sfctmp),
        pmsl=np.ascontiguousarray(pmsl),
        uwind=np.ascontiguousarray(uwind),
        vwind=np.ascontiguousarray(vwind),
        tpw=np.ascontiguousarray(tpw),
        elev=np.ascontiguousarray(geo['elevation'].astype(np.float32)),
        eco=np.ascontiguousarray(geo['eco_type'].astype(np.int8)),
        lsf=np.ascontiguousarray(geo['lsf'].astype(np.int8)),
        snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
        btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
        n_elem=n_elem, n_line=n_line,
    )
    print(f"  Fortran: {time.time()-t0:.1f}s")

    cm = result['cloud_mask']
    conf = result['confidence']
    cm[np.isnan(conf)] = 0
    conf[np.isnan(conf)] = 0.0

    return cm, conf, pixel_lat, geo['lon']


def compare_myd35(cm, lat, lon, myd35_path):
    """Compare FY3D cloud mask with MYD35."""
    myd35 = SD(myd35_path, SDC.READ)
    myd35_cm_raw = myd35.select('Cloud_Mask').get().astype(np.int32)
    byte0 = myd35_cm_raw[0, :, :]
    myd35_cloud_1km = (byte0 >> 1) & 0x03
    myd35_lat_5km = myd35.select('Latitude').get()  # (406, 270)
    myd35_lon_5km = myd35.select('Longitude').get()

    # Aggregate MYD35 to 5km mode
    myd35_cm_5km = np.zeros((406, 270), dtype=np.int32)
    for i in range(406):
        i1, i2 = i*5, min((i+1)*5, 2030)
        for j in range(270):
            j1, j2 = j*5, min((j+1)*5, 1354)
            block = myd35_cloud_1km[i1:i2, j1:j2]
            vals, counts = np.unique(block, return_counts=True)
            myd35_cm_5km[i, j] = vals[np.argmax(counts)]

    # KDTree matching
    myd35_lat_f = myd35_lat_5km.ravel()
    myd35_lon_f = myd35_lon_5km.ravel()
    myd35_cm_f = myd35_cm_5km.ravel()

    fy3_lat_f = lat.ravel()
    fy3_lon_f = lon.ravel()
    fy3_cm_f = cm.ravel()

    lat_min = max(myd35_lat_f.min(), fy3_lat_f.min())
    lat_max = min(myd35_lat_f.max(), fy3_lat_f.max())

    valid = (np.abs(myd35_lat_f) > 1e-6) & (myd35_lat_f >= lat_min) & (myd35_lat_f <= lat_max)
    fy3_ok = (fy3_lat_f >= lat_min) & (fy3_lat_f <= lat_max)

    myd35_points = np.stack([myd35_lat_f[valid], myd35_lon_f[valid]], axis=-1)
    fy3_points = np.stack([fy3_lat_f[fy3_ok], fy3_lon_f[fy3_ok]], axis=-1)

    tree = cKDTree(fy3_points)
    dist, idx = tree.query(myd35_points, k=1)
    good = dist < 0.025

    fy3_match = fy3_cm_f[fy3_ok][idx[good]]
    myd35_match = myd35_cm_f[valid][good]

    try:
        myd35.end()
    except Exception:
        pass

    return fy3_match, myd35_match


def analyze_noise(cm):
    """Count isolated pixels."""
    center = cm[1:-1, 1:-1]
    nbrs = np.stack([
        cm[:-2, :-2], cm[:-2, 1:-1], cm[:-2, 2:],
        cm[1:-1, :-2],                 cm[1:-1, 2:],
        cm[2:, :-2],   cm[2:, 1:-1],  cm[2:, 2:],
    ], axis=-1)
    same = np.sum(nbrs == center[:, :, None], axis=-1)
    return {k: int(np.sum(same == k)) for k in range(9)}, int(np.sum(same == 0))


# ---- Process all orbits ----
all_results = {}
for orb in ORBITS:
    name = orb['name']
    cm, conf, lat, lon = process_orbit(orb)

    # Noise analysis
    noise_dist, n_isolated = analyze_noise(cm)
    total_px = cm.size

    # MYD35 comparison
    fy3_match, myd35_match = compare_myd35(cm, lat, lon, orb['myd35'])

    # Stats
    labels = ['cloudy', 'prob_cld', 'prob_clr', 'clear']
    cm_dist = {i: int(np.sum(cm == i)) for i in range(4)}

    # Binary agreement
    fy3_clear = fy3_match >= 2
    myd35_clear = myd35_match >= 2
    agreement = np.mean(fy3_clear == myd35_clear) * 100

    n = len(fy3_match)
    fy3_cld_pct = np.mean(~fy3_clear) * 100
    myd35_cld_pct = np.mean(~myd35_clear) * 100

    # Save
    out_file = os.path.join(OUT_DIR, f'FY3D_MERSI_20200308_{name}_CLM_CLA.h5')
    with h5py.File(out_file, 'w') as f:
        f.create_dataset('cm', data=cm)
        f.create_dataset('conf', data=conf)
        f.create_dataset('lat', data=lat)
        f.create_dataset('lon', data=lon)
        f.attrs['version'] = 'v3.2.2'
        f.attrs['orbit'] = f'2020-03-08 {name[:2]}:{name[2:]}:00 UTC'

    all_results[name] = {
        'cm_dist': cm_dist,
        'total': total_px,
        'noise_dist': noise_dist,
        'n_isolated': n_isolated,
        'myd35_n': n,
        'myd35_agreement': agreement,
        'fy3_cld_pct': fy3_cld_pct,
        'myd35_cld_pct': myd35_cld_pct,
        'fy3_match': fy3_match,
        'myd35_match': myd35_match,
    }

# ---- Summary ----
print(f"\n{'='*60}")
print(f"SUMMARY: All 3 Orbits")
print(f"{'='*60}")

for name in ['1345', '1435', '1525']:
    r = all_results[name]
    print(f"\n--- Orbit {name} ---")
    print(f"  Total pixels: {r['total']:,}")
    print(f"  Isolated (0/8): {r['n_isolated']:,} ({100*r['n_isolated']/r['total']:.04f}%)")
    print(f"  Cloud dist: cloudy={r['cm_dist'][0]:,}, prob_cld={r['cm_dist'][1]:,}, prob_clr={r['cm_dist'][2]:,}, clear={r['cm_dist'][3]:,}")
    print(f"  MYD35 matched: {r['myd35_n']:,}")
    print(f"  MYD35 agreement (binary): {r['myd35_agreement']:.1f}%")
    print(f"  FY3D cloud: {r['fy3_cld_pct']:.1f}%, MYD35 cloud: {r['myd35_cld_pct']:.1f}%")

# Combined MYD35 confusion matrix
print(f"\n--- Combined 3-Orbit MYD35 Confusion Matrix ---")
labels = ['cloudy', 'prob_cld', 'prob_clr', 'clear']
all_fy3 = np.concatenate([all_results[n]['fy3_match'] for n in ['1345', '1435', '1525']])
all_myd35 = np.concatenate([all_results[n]['myd35_match'] for n in ['1345', '1435', '1525']])
print(f"  Total matched: {len(all_fy3):,}")
all_agree = np.mean((all_fy3 >= 2) == (all_myd35 >= 2)) * 100
print(f"  Overall binary agreement: {all_agree:.1f}%")

print(f"\n{'':>18}", end="")
for l in labels:
    print(f"{'MYD35 '+l:>16}", end="")
print()
for i, li in enumerate(labels):
    row = np.sum(all_fy3 == i)
    print(f"{'FY3D '+li:>18}", end="")
    for j in range(4):
        n_ij = np.sum((all_fy3 == i) & (all_myd35 == j))
        pct = 100 * n_ij / row if row > 0 else 0
        print(f"{n_ij:>8,} {pct:>5.1f}%", end="")
    print()

print("\nDone!")
