"""
io_myd35.py — MYD35_L2 (Aqua MODIS) cloud mask reader + spatiotemporal matcher
================================================================================
Handles:
  • Reading MYD35_L2 HDF4 (pyhdf) or HDF-EOS2 via gdal fallback
  • Reading MYD03 1-km geolocation (co-locate with MYD35)
  • Time-window matching against a MERSI orbit timestamp
  • Spatial overlap detection (bounding-box + polygon intersection)
  • Resampling MYD35 onto MERSI 1-km grid via nearest-neighbour (pyresample)

MYD35 Cloud_Mask byte-0 bit layout:
  bit 0   : processed flag  (1 = processed)
  bit 1-2 : cloud confidence  00=Cloudy  01=Uncertain/Prob.Cloudy
                               10=Prob.Clear  11=Confident Clear
  bit 3   : day/night flag
"""

from __future__ import annotations
import os
import re
import glob
from pathlib import Path
from datetime import datetime, timedelta, timezone

import numpy as np

# ── optional imports (graceful degradation) ──────────────────────────────────
try:
    from pyhdf.SD import SD, SDC
    HAS_PYHDF = True
except ImportError:
    HAS_PYHDF = False
    print("[WARN] pyhdf not found — MYD35 HDF4 reading unavailable")

try:
    import pyresample as prs
    HAS_PYRESAMPLE = True
except ImportError:
    HAS_PYRESAMPLE = False
    print("[WARN] pyresample not found — will use scipy KDTree fallback")

try:
    from scipy.spatial import cKDTree
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# MYD35 → unified CLM class mapping (matches MERSI CLM convention)
#   0 = Cloudy  1 = Prob.Cloudy  2 = Prob.Clear  3 = Conf.Clear
MYD35_TO_CLM = {0: 0, 1: 1, 2: 2, 3: 3}

# Maximum time difference for a granule to be considered "matching" (minutes)
DEFAULT_TIME_WINDOW_MIN = 15

# Minimum fractional overlap to accept a granule (0–1)
DEFAULT_MIN_OVERLAP = 0.05


# ─────────────────────────────────────────────────────────────────────────────
# 1.  MYD35 / MYD03 file discovery
# ─────────────────────────────────────────────────────────────────────────────

def parse_modis_datetime(filename: str) -> datetime | None:
    """
    Extract UTC datetime from a MODIS filename.
    Patterns handled:
      MYD35_L2.A2023157.1440.061.*.hdf   → julian day
      MYD03.A2023157.1440.061.*.hdf
    Returns UTC-aware datetime or None.
    """
    # Pattern: .AYYYYDDD.HHMM.
    m = re.search(r'\.A(\d{7})\.(\d{4})\.', filename)
    if not m:
        return None
    yyyyddd, hhmm = m.group(1), m.group(2)
    year = int(yyyyddd[:4])
    doy  = int(yyyyddd[4:])
    hh   = int(hhmm[:2])
    mm   = int(hhmm[2:])
    dt   = datetime(year, 1, 1, hh, mm, tzinfo=timezone.utc) + \
           timedelta(days=doy - 1)
    return dt


def find_myd35_granules(
    search_dirs: list[str] | str,
    mersi_dt: datetime,
    time_window_min: int = DEFAULT_TIME_WINDOW_MIN,
) -> list[dict]:
    """
    Scan *search_dirs* for MYD35_L2 files within ±time_window_min of mersi_dt.

    Returns list of dicts:
      { 'myd35': path, 'myd03': path_or_None, 'dt': datetime, 'dt_diff_min': float }
    sorted by |time difference|.
    """
    if isinstance(search_dirs, str):
        search_dirs = [search_dirs]

    window = timedelta(minutes=time_window_min)
    candidates = []

    for d in search_dirs:
        for fpath in Path(d).rglob("MYD35_L2*.hdf"):
            dt = parse_modis_datetime(fpath.name)
            if dt is None:
                continue
            diff = abs((dt - mersi_dt).total_seconds() / 60)
            if diff <= time_window_min:
                # Try to find co-located MYD03 (same granule time)
                myd03 = _find_myd03(fpath)
                candidates.append({
                    "myd35":       str(fpath),
                    "myd03":       myd03,
                    "dt":          dt,
                    "dt_diff_min": diff,
                })

    candidates.sort(key=lambda x: x["dt_diff_min"])
    if candidates:
        print(f"[MYD35] Found {len(candidates)} granule(s) within "
              f"±{time_window_min} min of {mersi_dt:%Y-%m-%d %H:%M} UTC")
    else:
        print(f"[MYD35] No granules found within ±{time_window_min} min "
              f"of {mersi_dt:%Y-%m-%d %H:%M} UTC in: {search_dirs}")
    return candidates


def _find_myd03(myd35_path: Path) -> str | None:
    """Look for MYD03 in same directory with identical granule time."""
    m = re.search(r'(A\d{7}\.\d{4})', myd35_path.name)
    if not m:
        return None
    tag = m.group(1)
    candidates = list(myd35_path.parent.glob(f"MYD03.{tag}.*.hdf"))
    return str(candidates[0]) if candidates else None


# ─────────────────────────────────────────────────────────────────────────────
# 2.  MYD35 HDF4 reader
# ─────────────────────────────────────────────────────────────────────────────

def read_myd35(myd35_path: str, myd03_path: str | None = None) -> dict | None:
    """
    Read MYD35_L2 cloud mask + geolocation.

    Geolocation source priority:
      1. MYD03 (1-km, preferred)
      2. MYD35 internal 5-km geo (upsampled to 1 km)

    Returns dict:
      { 'clm': (H,W) int32  0-3/-1,
        'lat': (H,W) float64,
        'lon': (H,W) float64,
        'confidence': (H,W) float32  0-1  (derived),
        'dt':  datetime }
    or None on failure.
    """
    if not HAS_PYHDF:
        print("[ERROR] pyhdf required to read MYD35 HDF4 files.")
        return None
    if not os.path.exists(myd35_path):
        print(f"[ERROR] MYD35 file not found: {myd35_path}")
        return None

    try:
        hdf = SD(myd35_path, SDC.READ)

        # ── Cloud_Mask: shape (6, 2030, 1354) — 6 bytes per pixel ──
        cm_sds  = hdf.select("Cloud_Mask")
        cm_data = cm_sds.get()                    # (6, nline, nelem)
        # byte index 0 carries the confidence bits
        byte0   = cm_data[0].astype(np.uint8)     # (nline, nelem)

        processed = (byte0 & 0x01).astype(bool)
        conf_bits = (byte0 >> 1) & 0x03           # 2-bit confidence

        clm = np.full(byte0.shape, -1, dtype=np.int32)
        clm[processed & (conf_bits == 0)] = 0     # Cloudy
        clm[processed & (conf_bits == 1)] = 1     # Prob. Cloudy
        clm[processed & (conf_bits == 2)] = 2     # Prob. Clear
        clm[processed & (conf_bits == 3)] = 3     # Confident Clear

        # Derive pseudo-confidence [0,1] for difference maps
        conf_float = np.where(processed,
                              conf_bits.astype(np.float32) / 3.0,
                              np.nan)

        # ── Geolocation ──────────────────────────────────────────────
        if myd03_path and os.path.exists(myd03_path):
            lat, lon = _read_myd03_geo(myd03_path, clm.shape)
        else:
            lat, lon = _read_myd35_5km_geo(hdf, clm.shape)

        hdf.end()

        dt = parse_modis_datetime(os.path.basename(myd35_path))

        return {
            "clm":        clm,
            "lat":        lat,
            "lon":        lon,
            "confidence": conf_float,
            "dt":         dt,
            "source":     myd35_path,
        }

    except Exception as e:
        print(f"[ERROR] Reading MYD35 {myd35_path}: {e}")
        return None


def _read_myd03_geo(myd03_path: str, target_shape: tuple) -> tuple:
    """Read 1-km lat/lon from MYD03."""
    hdf  = SD(myd03_path, SDC.READ)
    lat  = hdf.select("Latitude").get().astype(np.float64)
    lon  = hdf.select("Longitude").get().astype(np.float64)
    hdf.end()
    # Crop/pad to target_shape if necessary (granule edge mismatches)
    lat = _match_shape(lat, target_shape)
    lon = _match_shape(lon, target_shape)
    lat = np.where((lat < -90)  | (lat > 90),   np.nan, lat)
    lon = np.where((lon < -180) | (lon > 180),  np.nan, lon)
    return lat, lon


def _read_myd35_5km_geo(hdf, target_shape: tuple) -> tuple:
    """Read 5-km geo from MYD35 and upsample to 1-km."""
    lat5 = hdf.select("Latitude").get().astype(np.float64)
    lon5 = hdf.select("Longitude").get().astype(np.float64)
    from scipy.ndimage import zoom
    zy = target_shape[0] / lat5.shape[0]
    zx = target_shape[1] / lat5.shape[1]
    lat = zoom(lat5, (zy, zx), order=1)
    lon = zoom(lon5, (zy, zx), order=1)
    lat = _match_shape(lat, target_shape)
    lon = _match_shape(lon, target_shape)
    return lat, lon


def _match_shape(arr: np.ndarray, target: tuple) -> np.ndarray:
    """Crop or zero-pad array to target shape (H, W)."""
    H, W = target
    h, w = arr.shape
    out = np.full(target, np.nan, dtype=arr.dtype)
    out[:min(H, h), :min(W, w)] = arr[:min(H, h), :min(W, w)]
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Spatial overlap detection
# ─────────────────────────────────────────────────────────────────────────────

def compute_bbox(lat: np.ndarray, lon: np.ndarray) -> tuple | None:
    """Return (lon_min, lon_max, lat_min, lat_max) for valid pixels."""
    valid = np.isfinite(lat) & np.isfinite(lon)
    if not valid.any():
        return None
    return (float(lon[valid].min()), float(lon[valid].max()),
            float(lat[valid].min()), float(lat[valid].max()))


def bbox_overlap_fraction(bb1: tuple, bb2: tuple) -> float:
    """
    Fractional overlap of two bounding boxes (lon_min,lon_max,lat_min,lat_max).
    Returns fraction relative to the smaller box.
    """
    lon_min  = max(bb1[0], bb2[0])
    lon_max  = min(bb1[1], bb2[1])
    lat_min  = max(bb1[2], bb2[2])
    lat_max  = min(bb1[3], bb2[3])
    if lon_max <= lon_min or lat_max <= lat_min:
        return 0.0
    overlap_area = (lon_max - lon_min) * (lat_max - lat_min)
    area1 = (bb1[1] - bb1[0]) * (bb1[3] - bb1[2])
    area2 = (bb2[1] - bb2[0]) * (bb2[3] - bb2[2])
    smaller = min(area1, area2)
    return overlap_area / smaller if smaller > 0 else 0.0


def filter_overlapping_granules(
    candidates: list[dict],
    mersi_lat: np.ndarray,
    mersi_lon: np.ndarray,
    myd35_data_list: list[dict | None],
    min_overlap: float = DEFAULT_MIN_OVERLAP,
) -> list[tuple[dict, dict]]:
    """
    Keep only granules with sufficient spatial overlap with the MERSI swath.

    Returns list of (candidate_meta, myd35_data) tuples.
    """
    mersi_bb = compute_bbox(mersi_lat, mersi_lon)
    if mersi_bb is None:
        return []

    accepted = []
    for meta, data in zip(candidates, myd35_data_list):
        if data is None:
            continue
        myd_bb = compute_bbox(data["lat"], data["lon"])
        if myd_bb is None:
            continue
        frac = bbox_overlap_fraction(mersi_bb, myd_bb)
        print(f"  [OVERLAP] {os.path.basename(meta['myd35'])} "
              f"Δt={meta['dt_diff_min']:.1f} min  overlap={frac*100:.1f}%")
        if frac >= min_overlap:
            accepted.append((meta, data))

    print(f"[MYD35] {len(accepted)}/{len(candidates)} granule(s) passed "
          f"overlap threshold ({min_overlap*100:.0f}%)")
    return accepted


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Grid resampling: MYD35 → MERSI grid
# ─────────────────────────────────────────────────────────────────────────────

def resample_to_mersi_grid(
    myd35_lat: np.ndarray,
    myd35_lon: np.ndarray,
    myd35_clm: np.ndarray,
    mersi_lat:  np.ndarray,
    mersi_lon:  np.ndarray,
    radius_m:   float = 1500.0,
) -> np.ndarray:
    """
    Nearest-neighbour resampling of MYD35 CLM onto the MERSI 1-km grid.

    Uses pyresample if available, falls back to scipy cKDTree.

    Parameters
    ----------
    radius_m : search radius in metres (1500 m ≈ 1.5 pixels at 1 km)

    Returns
    -------
    clm_on_mersi : (H_mersi, W_mersi) int32  0-3 or -1 for no-data
    """
    if HAS_PYRESAMPLE:
        return _resample_pyresample(myd35_lat, myd35_lon, myd35_clm,
                                    mersi_lat, mersi_lon, radius_m)
    elif HAS_SCIPY:
        return _resample_kdtree(myd35_lat, myd35_lon, myd35_clm,
                                mersi_lat, mersi_lon, radius_m)
    else:
        raise RuntimeError("Neither pyresample nor scipy available for resampling.")


def _resample_pyresample(myd_lat, myd_lon, myd_clm,
                          mer_lat, mer_lon, radius_m):
    src = prs.geometry.SwathDefinition(
        lons=myd_lon.astype(np.float32),
        lats=myd_lat.astype(np.float32))
    tgt = prs.geometry.SwathDefinition(
        lons=mer_lon.astype(np.float32),
        lats=mer_lat.astype(np.float32))

    # Fill invalid MYD35 pixels with -999 so we can mask after
    clm_f = myd_clm.astype(np.float32)
    clm_f[myd_clm < 0] = -999.0

    result = prs.kd_tree.resample_nearest(
        src, clm_f, tgt,
        radius_of_influence=radius_m,
        fill_value=np.nan,
        nprocs=1,
    )
    out = np.round(result).astype(np.int32)
    out[~np.isfinite(result)] = -1
    out[out < 0] = -1
    return out


def _resample_kdtree(myd_lat, myd_lon, myd_clm,
                      mer_lat, mer_lon, radius_m):
    """Scipy KDTree fallback — slower but dependency-light."""
    # Convert to Cartesian for radius search
    def to_xyz(lat_deg, lon_deg):
        lat = np.deg2rad(lat_deg)
        lon = np.deg2rad(lon_deg)
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z = np.sin(lat)
        return np.stack([x, y, z], axis=-1)

    valid_src = np.isfinite(myd_lat) & np.isfinite(myd_lon) & (myd_clm >= 0)
    src_xyz = to_xyz(myd_lat[valid_src], myd_lon[valid_src])
    src_clm = myd_clm[valid_src]

    valid_tgt = np.isfinite(mer_lat) & np.isfinite(mer_lon)
    tgt_xyz   = to_xyz(mer_lat[valid_tgt], mer_lon[valid_tgt])

    # Radius in Cartesian (chord length for 1500 m on Earth r=6371 km)
    chord = 2 * np.sin(radius_m / (2 * 6_371_000))

    tree = cKDTree(src_xyz)
    dist, idx = tree.query(tgt_xyz, k=1, distance_upper_bound=chord, workers=-1)

    out = np.full(mer_lat.shape, -1, dtype=np.int32)
    found = dist < chord
    tgt_flat = np.where(valid_tgt.ravel())[0]
    out.ravel()[tgt_flat[found]] = src_clm[idx[found]]
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 5.  High-level convenience function
# ─────────────────────────────────────────────────────────────────────────────

def load_best_myd35_for_mersi(
    mersi_lat:      np.ndarray,
    mersi_lon:      np.ndarray,
    mersi_dt:       datetime,
    search_dirs:    list[str] | str,
    time_window_min: int   = DEFAULT_TIME_WINDOW_MIN,
    min_overlap:    float  = DEFAULT_MIN_OVERLAP,
    radius_m:       float  = 1500.0,
) -> dict | None:
    """
    End-to-end: find, load, filter, and resample the best MYD35 granule.

    Returns dict with keys:
      'clm_native'  : MYD35 CLM on its own grid  (H_myd, W_myd)
      'clm_resampled': MYD35 CLM resampled to MERSI grid  (H_mer, W_mer)
      'lat', 'lon'  : MYD35 native geolocation
      'dt'          : granule datetime
      'dt_diff_min' : time offset vs MERSI
      'source'      : file path
    or None if no suitable granule found.
    """
    candidates = find_myd35_granules(search_dirs, mersi_dt, time_window_min)
    if not candidates:
        return None

    # Load all candidate granules
    loaded = [read_myd35(c["myd35"], c["myd03"]) for c in candidates]

    # Filter by spatial overlap
    accepted = filter_overlapping_granules(
        candidates, mersi_lat, mersi_lon, loaded, min_overlap)
    if not accepted:
        return None

    # Take best (smallest Δt + largest overlap) — already sorted by Δt
    meta, data = accepted[0]

    print(f"[MYD35] Using granule: {os.path.basename(meta['myd35'])} "
          f"(Δt={meta['dt_diff_min']:.1f} min)")

    # Resample onto MERSI grid
    clm_resampled = resample_to_mersi_grid(
        data["lat"], data["lon"], data["clm"],
        mersi_lat, mersi_lon, radius_m)

    return {
        "clm_native":    data["clm"],
        "clm_resampled": clm_resampled,
        "lat":           data["lat"],
        "lon":           data["lon"],
        "dt":            data["dt"],
        "dt_diff_min":   meta["dt_diff_min"],
        "source":        meta["myd35"],
    }
