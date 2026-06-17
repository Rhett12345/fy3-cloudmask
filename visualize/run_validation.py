"""
run_validation.py — FY-3D MERSI-II × MYD35 cloud mask validation pipeline
===========================================================================
Orchestrates the full workflow:
  1. Load MERSI CLM (recal + onboard) and L1B RGB
  2. Parse MERSI orbit datetime
  3. Find, load, filter, and resample matching MYD35 granule
  4. Build Figure 1  (MERSI recal vs onboard)
  5. Build Figure 2  (MERSI vs MYD35 truth)
  6. Print & optionally save validation statistics

CLI usage
---------
# Single orbit from HDF5 files:
python run_validation.py \\
    --recal    /data/.../YYYYMMDD_HHMM_CLM_CLA_recal.h5 \\
    --onboard  /data/.../YYYYMMDD_HHMM_CLM_CLA.h5 \\
    --myd35_dir /data/myd35/ \\
    --output_dir ./output/

# Batch mode (auto-discovers CLM pairs):
python run_validation.py \\
    --data_dir  /data/Data_yuq/fy3_cloud/recal_test/ \\
    --myd35_dir /data/myd35/ \\
    --output_dir ./output/ \\
    --step 4

Python API usage
----------------
from run_validation import run_single_orbit

stats = run_single_orbit(
    recal_path   = "path/recal_CLM.h5",
    onboard_path = "path/onboard_CLM.h5",
    myd35_dirs   = ["/data/myd35/"],
    output_dir   = "./output",
)
"""

from __future__ import annotations
import argparse
import os
import re
from pathlib import Path

import numpy as np

from io_mersi import (
    load_clm_hdf5, load_mersi_l1b, find_l1b_for_clm,
    parse_mersi_datetime, print_clm_distribution,
)
from io_myd35 import load_best_myd35_for_mersi
from figure_1 import make_figure1_from_arrays
from figure_2 import make_figure2


# ─────────────────────────────────────────────────────────────────────────────
# Single-orbit runner
# ─────────────────────────────────────────────────────────────────────────────

def run_single_orbit(
    recal_path:      str,
    onboard_path:    str,
    myd35_dirs:      list[str],
    output_dir:      str       = "./output",
    mersi_root:      str       = "/data/Data_yuq/mersi",
    step:            int       = 4,
    time_window_min: int       = 15,
    min_overlap:     float     = 0.05,
    skip_fig1:       bool      = False,
    skip_fig2:       bool      = False,
) -> dict:
    """
    Full pipeline for one orbit.

    Returns
    -------
    dict with keys: 'fig1_output', 'fig2_output', 'stats' (or None if skipped)
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 1. Extract orbit tag for output filenames ────────────────────
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(onboard_path))
    orbit_tag = f"{m.group(1)}_{m.group(2)}" if m else "orbit"
    date_str  = ""
    if m:
        d, t = m.group(1), m.group(2)
        date_str = f"{d[:4]}-{d[4:6]}-{d[6:8]}  {t[:2]}:{t[2:]} UTC"

    print(f"\n{'═'*65}")
    print(f"  Orbit: {orbit_tag}   ({date_str})")
    print(f"{'═'*65}")

    # ── 2. Load MERSI CLM ────────────────────────────────────────────
    recal_data   = load_clm_hdf5(recal_path)
    onboard_data = load_clm_hdf5(onboard_path)
    if recal_data is None or onboard_data is None:
        print("[ERROR] CLM loading failed — skipping orbit.")
        return {}

    lat         = recal_data["lat"]
    lon         = recal_data["lon"]
    recal_clm   = recal_data["clm"]
    onboard_clm = onboard_data["clm"]

    # ── 3. Load MERSI L1B RGB ────────────────────────────────────────
    l1b_path = find_l1b_for_clm(recal_path, mersi_root)
    l1b_data = load_mersi_l1b(l1b_path) if l1b_path else None
    rgb      = l1b_data["rgb"] if l1b_data else None
    if rgb is None:
        print("[WARN] No L1B RGB — panels will show grey footprint.")

    # ── 4. Figure 1: MERSI recal vs onboard ─────────────────────────
    fig1_out = str(out_dir / f"fig1_{orbit_tag}.png")
    if not skip_fig1:
        make_figure1_from_arrays(
            recal_clm   = recal_clm,
            onboard_clm = onboard_clm,
            lat         = lat,
            lon         = lon,
            rgb         = rgb,
            output      = fig1_out,
            date_str    = date_str,
            step        = step,
        )
    else:
        print(f"[SKIP] Figure 1")

    # ── 5. MYD35 matching ────────────────────────────────────────────
    stats = None
    fig2_out = None

    if skip_fig2:
        print("[SKIP] Figure 2")
    else:
        mersi_dt = parse_mersi_datetime(onboard_path)
        if mersi_dt is None:
            print("[WARN] Could not parse MERSI datetime — skipping Figure 2.")
        elif not myd35_dirs:
            print("[WARN] No MYD35 search directories provided — skipping Figure 2.")
        else:
            myd35_data = load_best_myd35_for_mersi(
                mersi_lat       = lat,
                mersi_lon       = lon,
                mersi_dt        = mersi_dt,
                search_dirs     = myd35_dirs,
                time_window_min = time_window_min,
                min_overlap     = min_overlap,
            )

            if myd35_data is None:
                print("[WARN] No matching MYD35 granule found — Figure 2 skipped.")
            else:
                fig2_out = str(out_dir / f"fig2_{orbit_tag}.png")
                stats = make_figure2(
                    mersi_lat      = lat,
                    mersi_lon      = lon,
                    recal_clm      = recal_clm,
                    onboard_clm    = onboard_clm,
                    myd35_data     = myd35_data,
                    mersi_rgb      = rgb,
                    output         = fig2_out,
                    mersi_date_str = date_str,
                    step           = step,
                )

    return {
        "fig1_output": fig1_out,
        "fig2_output": fig2_out,
        "stats":       stats,
        "orbit_tag":   orbit_tag,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Batch runner
# ─────────────────────────────────────────────────────────────────────────────

def find_orbit_pairs(data_dir: str) -> list[tuple[str, str, str, str]]:
    """
    Scan data_dir for onboard/recal CLM pairs.
    Returns list of (onboard_path, recal_path, date_str, time_str).
    """
    pairs = []
    for onboard in sorted(Path(data_dir).rglob("*_CLM_CLA.h5")):
        if "_recal" in onboard.name:
            continue
        recal = onboard.with_name(
            onboard.name.replace("_CLM_CLA.h5", "_CLM_CLA_recal.h5"))
        if not recal.exists():
            print(f"[WARN] Missing recal for {onboard.name}")
            continue
        m = re.search(r'(\d{8})_(\d{4})', onboard.name)
        if not m:
            continue
        pairs.append((str(onboard), str(recal), m.group(1), m.group(2)))
    return pairs


def run_batch(
    data_dir:        str,
    myd35_dirs:      list[str],
    output_dir:      str   = "./output",
    mersi_root:      str   = "/data/Data_yuq/mersi",
    step:            int   = 4,
    time_window_min: int   = 15,
    min_overlap:     float = 0.05,
    overwrite:       bool  = False,
) -> list[dict]:
    """Run the full pipeline for all orbits in data_dir."""
    pairs = find_orbit_pairs(data_dir)
    if not pairs:
        print(f"[ERROR] No CLM pairs found in {data_dir}")
        return []

    print(f"[BATCH] Found {len(pairs)} orbit pair(s)")
    all_results = []

    for onboard, recal, date, time in pairs:
        orbit_tag = f"{date}_{time}"
        out_dir   = Path(output_dir)
        fig1_out  = out_dir / f"fig1_{orbit_tag}.png"
        fig2_out  = out_dir / f"fig2_{orbit_tag}.png"

        if not overwrite and fig1_out.exists() and fig2_out.exists():
            print(f"[SKIP] {orbit_tag} (outputs exist)")
            continue

        try:
            result = run_single_orbit(
                recal_path      = recal,
                onboard_path    = onboard,
                myd35_dirs      = myd35_dirs,
                output_dir      = str(out_dir),
                mersi_root      = mersi_root,
                step            = step,
                time_window_min = time_window_min,
                min_overlap     = min_overlap,
                skip_fig1       = (not overwrite and fig1_out.exists()),
                skip_fig2       = (not overwrite and fig2_out.exists()),
            )
            all_results.append(result)
        except Exception as e:
            print(f"[ERROR] {orbit_tag}: {e}")
            import traceback; traceback.print_exc()

    print(f"\n[DONE] Processed {len(all_results)} orbit(s) → {output_dir}/")
    return all_results


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="FY-3D MERSI-II × MYD35 cloud mask validation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single orbit
  python run_validation.py \\
      --recal    ./data/20230606_1440_CLM_CLA_recal.h5 \\
      --onboard  ./data/20230606_1440_CLM_CLA.h5 \\
      --myd35_dir /data/myd35/2023/ \\
      --output_dir ./output/

  # Batch
  python run_validation.py \\
      --data_dir   /data/Data_yuq/fy3_cloud/recal_test/ \\
      --myd35_dir  /data/myd35/ \\
      --output_dir ./output/ \\
      --step 4 --time_window 15
        """,
    )
    # Input modes
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--data_dir",  help="Batch: directory with CLM pairs")
    g.add_argument("--onboard",   help="Single: onboard CLM HDF5 path")

    parser.add_argument("--recal",       help="Recalibration CLM HDF5 path (single mode)")
    parser.add_argument("--myd35_dir",   nargs="+", default=[],
                        metavar="DIR",   help="Directory/ies to search for MYD35 files")
    parser.add_argument("--output_dir",  default="./output")
    parser.add_argument("--mersi_root",  default="/data/Data_yuq/mersi")
    parser.add_argument("--step",        type=int,   default=4,
                        help="Subsampling step (4 = every 4th pixel)")
    parser.add_argument("--time_window", type=int,   default=15,
                        help="MYD35 temporal search window in minutes (default 15)")
    parser.add_argument("--min_overlap", type=float, default=0.05,
                        help="Minimum spatial overlap fraction (default 0.05)")
    parser.add_argument("--overwrite",   action="store_true",
                        help="Overwrite existing output files")
    parser.add_argument("--skip_fig1",   action="store_true")
    parser.add_argument("--skip_fig2",   action="store_true")

    args = parser.parse_args()

    if args.data_dir:
        run_batch(
            data_dir        = args.data_dir,
            myd35_dirs      = args.myd35_dir,
            output_dir      = args.output_dir,
            mersi_root      = args.mersi_root,
            step            = args.step,
            time_window_min = args.time_window,
            min_overlap     = args.min_overlap,
            overwrite       = args.overwrite,
        )
    else:
        if not args.recal:
            parser.error("--recal is required in single-orbit mode")
        run_single_orbit(
            recal_path      = args.recal,
            onboard_path    = args.onboard,
            myd35_dirs      = args.myd35_dir,
            output_dir      = args.output_dir,
            mersi_root      = args.mersi_root,
            step            = args.step,
            time_window_min = args.time_window,
            min_overlap     = args.min_overlap,
            skip_fig1       = args.skip_fig1,
            skip_fig2       = args.skip_fig2,
        )
