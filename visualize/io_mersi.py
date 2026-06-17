"""
io_mersi.py — FY-3D MERSI-II L1B + CLM data readers
=====================================================
Clean, modular version of the I/O functions from visualize_clm_nature.py.

Public API
----------
load_mersi_l1b(l1b_path)
    → dict with keys: rgb (H,W,3), lat, lon

load_clm_hdf5(path)
    → dict with keys: clm (H,W), lat, lon

parse_mersi_datetime(filename)
    → datetime (UTC-aware)

find_l1b_for_clm(clm_path, mersi_root)
    → path string or None
"""

from __future__ import annotations
import os
import re
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import h5py


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Filename utilities
# ─────────────────────────────────────────────────────────────────────────────

def parse_mersi_datetime(filename: str) -> datetime | None:
    """
    Extract UTC datetime from FY-3D MERSI filename.
    Pattern: ..._YYYYMMDD_HHMM_...
    Returns UTC-aware datetime or None.
    """
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(filename))
    if not m:
        return None
    d, t = m.group(1), m.group(2)
    try:
        return datetime(int(d[:4]), int(d[4:6]), int(d[6:8]),
                        int(t[:2]),  int(t[2:]),
                        tzinfo=timezone.utc)
    except ValueError:
        return None


def find_l1b_for_clm(
    clm_path: str,
    mersi_root: str | Path = "/data/Data_yuq/mersi",
) -> str | None:
    """
    Locate the MERSI L1B 1-km HDF file corresponding to a CLM file.
    Searches <mersi_root>/<YYYYMMDD>/FY3D_MERSI_GBAL_L1_<date>_<time>_1000M_MS.HDF
    """
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    p = (Path(mersi_root) / date_str /
         f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF")
    return str(p) if p.exists() else None


# ─────────────────────────────────────────────────────────────────────────────
# 2.  MERSI L1B reader  →  RGB + geolocation
# ─────────────────────────────────────────────────────────────────────────────

def load_mersi_l1b(l1b_path: str) -> dict | None:
    """
    Read MERSI-II 1-km L1B HDF and produce a display-ready RGB array.

    Uses bands B1/B2/B3 (EV_250_Aggr.1KM_RefSB channels 0,1,2)
    and the VIS calibration coefficients.

    Returns
    -------
    dict with:
      'rgb'  : (H, W, 3) float32, values clipped to [0, 1]
      'lat'  : (H, W) float64
      'lon'  : (H, W) float64
      'path' : input file path
    or None on failure.
    """
    if not l1b_path or not os.path.exists(l1b_path):
        print(f"[WARN] L1B not found: {l1b_path}")
        return None

    try:
        with h5py.File(l1b_path, "r") as f:
            vis_250 = f["Data/EV_250_Aggr.1KM_RefSB"][:].astype(np.float64)
            vis_cal = f["Calibration/VIS_Cal_Coeff"][:]
            esd     = float(np.squeeze(f.attrs.get("EarthSun Distance Ratio", 1.0)))

            # Geolocation — try several common paths
            lat, lon = _read_l1b_geo(f)

        esd2   = esd ** 2
        n_bands, n_line, n_elem = vis_250.shape
        rgb = np.zeros((n_line, n_elem, 3), dtype=np.float32)

        for i in range(3):                        # R, G, B → bands 0, 1, 2
            c0, c1, c2 = vis_cal[i]
            dn   = vis_250[i]
            refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
            rgb[:, :, i] = refl.astype(np.float32)

        rgb = _stretch_rgb(rgb)

        # Transpose to match CLM lat/lon axis order (elem, line) ← (line, elem)
        rgb = np.transpose(rgb, (1, 0, 2))

        return {"rgb": rgb, "lat": lat, "lon": lon, "path": l1b_path}

    except Exception as e:
        print(f"[ERROR] Loading L1B {l1b_path}: {e}")
        return None


def _read_l1b_geo(f: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    """Try multiple common geolocation paths inside L1B HDF."""
    # Ordered by likelihood
    paths = [
        ("Geolocation/Latitude",    "Geolocation/Longitude"),
        ("Data/Latitude",           "Data/Longitude"),
        ("Navigation/Latitude",     "Navigation/Longitude"),
    ]
    for lp, lonp in paths:
        if lp in f and lonp in f:
            lat = f[lp][:].astype(np.float64)
            lon = f[lonp][:].astype(np.float64)
            lat = np.where((lat < -90)  | (lat > 90),   np.nan, lat)
            lon = np.where((lon < -180) | (lon > 180),  np.nan, lon)
            return lat, lon
    raise KeyError("Geolocation datasets not found in L1B HDF")


def _stretch_rgb(rgb: np.ndarray, lo_pct: float = 2, hi_pct: float = 98) -> np.ndarray:
    """Per-channel percentile stretch → [0, 1]."""
    out = np.zeros_like(rgb)
    for ch in range(3):
        band  = rgb[:, :, ch]
        valid = band > 0
        if not valid.any():
            continue
        p2, p98 = np.percentile(band[valid], [lo_pct, hi_pct])
        if p98 <= p2:
            continue
        out[:, :, ch] = np.clip((band - p2) / (p98 - p2 + 1e-10), 0, 1)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 3.  MERSI CLM reader  (HDF5)
# ─────────────────────────────────────────────────────────────────────────────

def load_clm_hdf5(path: str) -> dict | None:
    """
    Read MERSI CLM from HDF5.  Handles two storage conventions:

      Flat:   /cm, /conf, /lat, /lon  (root-level datasets)
      Group:  Cloud_Mask_1km/  (grouped datasets)

    Returns
    -------
    dict with keys: clm (H,W int32), lat, lon, path
    or None on failure.
    """
    if not os.path.exists(path):
        print(f"[ERROR] CLM file not found: {path}")
        return None
    try:
        with h5py.File(path, "r") as f:
            # Auto-detect structure: flat (/cm, /lat, /lon) vs grouped (Cloud_Mask_1km/)
            if "cm" in f:
                # Flat structure
                lat = f["lat"][:].astype(np.float64)
                lon = f["lon"][:].astype(np.float64)
                clm = f["cm"][:].astype(np.int32)
                # If cm is all zeros, try deriving from conf
                if np.all(clm == 0) and "conf" in f:
                    derived = _clm_from_confidence(f["conf"][:].astype(np.float64))
                    if np.any(derived >= 0):
                        clm = derived
            elif "Cloud_Mask_1km" in f:
                # Original grouped structure
                grp = f["Cloud_Mask_1km"]
                lat = grp["Latitude"][:].astype(np.float64)
                lon = grp["Longitude"][:].astype(np.float64)
                clm = grp["Cloud_Mask_Value"][:].astype(np.int32)
                if np.all(clm == 0) and "Cloud_Mask" in grp:
                    decoded = _decode_bitmask(grp["Cloud_Mask"][:])
                    if np.any(decoded >= 0):
                        clm = decoded
                if np.all(clm == 0) and "Confidence" in grp:
                    derived = _clm_from_confidence(grp["Confidence"][:])
                    if np.any(derived >= 0):
                        clm = derived
            else:
                print(f"[ERROR] Unknown HDF5 structure in {path}")
                return None

        lat = np.where((lat < -90)  | (lat > 90),   np.nan, lat)
        lon = np.where((lon < -180) | (lon > 180),  np.nan, lon)
        clm = np.where((clm < 0)   | (clm > 3),     -1,    clm)

        return {"clm": clm, "lat": lat, "lon": lon, "path": path}

    except Exception as e:
        print(f"[ERROR] Loading CLM {path}: {e}")
        return None

def _decode_bitmask(cm_raw: np.ndarray) -> np.ndarray:
    """Decode 3-byte (or more) CLM bitmask array → class 0-3."""
    byte0     = cm_raw[:, :, 0].astype(np.uint8)
    processed = (byte0 & 0x01).astype(bool)
    b1        = (byte0 >> 1) & 1
    b2        = (byte0 >> 2) & 1
    result    = np.full(byte0.shape, -1, dtype=np.int32)
    result[processed & (b2 == 0) & (b1 == 0)] = 0
    result[processed & (b2 == 0) & (b1 == 1)] = 1
    result[processed & (b2 == 1) & (b1 == 0)] = 2
    result[processed & (b2 == 1) & (b1 == 1)] = 3
    return result


def _clm_from_confidence(conf: np.ndarray) -> np.ndarray:
    """Bin float confidence [0,1] into 4 CLM classes."""
    clm   = np.full(conf.shape, -1, dtype=np.int32)
    valid = np.isfinite(conf)
    clm[valid & (conf < 0.33)]               = 0
    clm[valid & (conf >= 0.33) & (conf < 0.66)] = 1
    clm[valid & (conf >= 0.66) & (conf < 0.90)] = 2
    clm[valid & (conf >= 0.90)]                  = 3
    return clm


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Distribution printer  (moved here from visualize_clm_nature.py)
# ─────────────────────────────────────────────────────────────────────────────

CLM_LABEL = {0: "Cloudy", 1: "Prob. Cloudy", 2: "Prob. Clear", 3: "Conf. Clear"}


def print_clm_distribution(label: str, clm: np.ndarray) -> None:
    """Print pixel category distribution to stdout."""
    valid_mask  = (clm >= 0) & (clm <= 3)
    total_valid = int(valid_mask.sum())
    total_all   = clm.size
    sep = "═" * 62
    print(sep)
    print(f"  CLM Distribution — {label}")
    print(f"  Total pixels : {total_all:>10,}  |  Valid : {total_valid:>10,}")
    print("─" * 62)
    print(f"  {'Category':<20} {'Class':>5}  {'Count':>10}  {'Ratio':>7}")
    print("─" * 62)
    for v in range(4):
        cnt = int((clm == v).sum())
        pct = 100.0 * cnt / total_valid if total_valid > 0 else 0.0
        bar = "█" * int(pct / 2)
        print(f"  {CLM_LABEL[v]:<20} {v:>5}  {cnt:>10,}  {pct:>6.2f}%  {bar}")
    invalid = int((clm < 0).sum())
    inv_pct = 100.0 * invalid / total_all if total_all > 0 else 0.0
    print("─" * 62)
    print(f"  {'Invalid/unprocessed':<20} {'–':>5}  {invalid:>10,}  {inv_pct:>6.2f}%")
    print(sep + "\n")
