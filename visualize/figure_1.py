"""
figure_1.py — MERSI-II CLM comparison figure  (4-panel, Nature style)
======================================================================
Panel layout:
  (a) RGB true colour (MERSI L1B)
  (b) Recalibration CLM
  (c) Onboard CLM
  (d) Class agreement map  (Recal − Onboard)

v2 bug-fix notes
----------------
- Per-class % stats used to be drawn as a floating box INSIDE the map
  axes. On long thin polar swaths the box overflowed the panel and
  visually collided with the colorbar text (see earlier screenshots).
  They are now rendered as a one-line caption BELOW each panel via
  plot_utils.panel_caption(), which always has its own dedicated space.
- Figure canvas is sized per-panel (3.2in × 2.8in) rather than forced
  into a fixed 174mm total width, giving every map enough room.

Usage
-----
# From HDF5 files:
from figure_1 import make_figure1_from_files
make_figure1_from_files(
    recal_path   = "path/to/recal_CLM.h5",
    onboard_path = "path/to/onboard_CLM.h5",
    output       = "figure1.png",
    mersi_root   = "/data/Data_yuq/mersi",
)

# From pre-loaded arrays:
from figure_1 import make_figure1_from_arrays
make_figure1_from_arrays(
    rgb          = rgb_array,    # (H,W,3) float32 or None
    recal_clm    = recal_clm,    # (H,W) int32
    onboard_clm  = onboard_clm,
    lat          = lat,
    lon          = lon,
    output       = "figure1.png",
    date_str     = "2023-06-06  14:40 UTC",
)
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from plot_utils import (
    apply_nature_style, PANEL_WIDTH_IN, PANEL_HEIGHT_IN,
    make_geo_ax_with_caption, add_gridlines, panel_label, panel_title,
    panel_caption,
    plot_rgb, plot_rgb_placeholder, plot_clm, plot_diff,
    add_clm_colorbar, add_diff_colorbar,
    stats_caption_text, agreement_caption_text,
    get_extent, save_figure,
)
from io_mersi import (
    load_clm_hdf5, load_mersi_l1b, find_l1b_for_clm,
    parse_mersi_datetime, print_clm_distribution,
)


# ─────────────────────────────────────────────────────────────────────────────
# Core figure builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_figure1(
    lat:         np.ndarray,
    lon:         np.ndarray,
    rgb:         np.ndarray | None,   # (H,W,3) float32 [0,1] or None
    recal_clm:   np.ndarray,          # (H,W) int32
    onboard_clm: np.ndarray,          # (H,W) int32
    date_str:    str,
    output:      str,
    step:        int = 4,
) -> None:
    """Render and save the 4-panel Figure 1."""
    apply_nature_style()

    extent = get_extent(lat, lon, recal_clm, step=step)

    fig_w = PANEL_WIDTH_IN * 2 + 1.2
    fig_h = PANEL_HEIGHT_IN * 2 + 0.7
    fig   = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    gs = GridSpec(2, 2, figure=fig,
                  left=0.06, right=0.88,
                  top=0.91,  bottom=0.05,
                  wspace=0.28, hspace=0.36)

    titles  = [
        "RGB true colour",
        "Recalibration CLM",
        "Onboard CLM",
        "Class agreement",
    ]
    letters = ["a", "b", "c", "d"]
    specs   = [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]
    axs     = [make_geo_ax_with_caption(fig, sp) for sp in specs]

    # ── (a) RGB ─────────────────────────────────────────────────────
    ax = axs[0]
    if rgb is not None:
        plot_rgb(ax, lat, lon, rgb, step=step)
    else:
        plot_rgb_placeholder(ax, lat, lon, recal_clm, step=step)
    if extent:
        ax.set_extent(extent)
    add_gridlines(ax)

    # ── (b) Recalibration CLM ───────────────────────────────────────
    ax = axs[1]
    plot_clm(ax, lat, lon, recal_clm, step=step)
    if extent:
        ax.set_extent(extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    panel_caption(ax, stats_caption_text(recal_clm))

    # ── (c) Onboard CLM ─────────────────────────────────────────────
    ax = axs[2]
    plot_clm(ax, lat, lon, onboard_clm, step=step)
    if extent:
        ax.set_extent(extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    panel_caption(ax, stats_caption_text(onboard_clm))

    # ── (d) Agreement map  (Recal − Onboard) ────────────────────────
    ax = axs[3]
    sm, mask, diff = plot_diff(ax, lat, lon, recal_clm, onboard_clm, step=step)
    if extent:
        ax.set_extent(extent)
    add_gridlines(ax)
    add_diff_colorbar(fig, ax, sm, label="Recal − Onboard (class)")
    panel_caption(ax, agreement_caption_text(recal_clm, onboard_clm))

    # ── Decorations ──────────────────────────────────────────────────
    for ax, ltr, ttl in zip(axs, letters, titles):
        panel_title(ax, ttl)
        panel_label(ax, ltr)

    lat_c = float(np.nanmedian(lat))
    lon_c = float(np.nanmedian(lon))
    fig.suptitle(
        f"FY-3D MERSI-II cloud mask comparison   {date_str}   "
        f"centre {lat_c:.1f}°N  {lon_c:.1f}°E",
        fontsize=10.5, fontweight="normal", color="#333333", y=0.98)

    save_figure(fig, output)


# ─────────────────────────────────────────────────────────────────────────────
# Public entry points
# ─────────────────────────────────────────────────────────────────────────────

def make_figure1_from_arrays(
    recal_clm:   np.ndarray,
    onboard_clm: np.ndarray,
    lat:         np.ndarray,
    lon:         np.ndarray,
    rgb:         np.ndarray | None = None,
    output:      str               = "figure1.png",
    date_str:    str               = "",
    step:        int               = 4,
) -> None:
    """
    Build Figure 1 from pre-loaded NumPy arrays.
    Call this from backend_compare_and_viz.py.
    """
    print(f"\n[FIG1] Orbit: {output}")
    print_clm_distribution("Recalibration CLM", recal_clm)
    print_clm_distribution("Onboard CLM",        onboard_clm)
    _build_figure1(lat, lon, rgb, recal_clm, onboard_clm,
                   date_str, output, step=step)


def make_figure1_from_files(
    recal_path:   str,
    onboard_path: str,
    output:       str  = "figure1.png",
    mersi_root:   str  = "/data/Data_yuq/mersi",
    step:         int  = 4,
) -> None:
    """
    Build Figure 1 by reading HDF5 CLM files from disk.
    Mirrors the original visualize_clm_nature.py::visualize() API.
    """
    recal_data   = load_clm_hdf5(recal_path)
    onboard_data = load_clm_hdf5(onboard_path)
    if recal_data is None or onboard_data is None:
        print("[ERROR] Could not load CLM data — Figure 1 skipped.")
        return

    l1b_path  = find_l1b_for_clm(recal_path, mersi_root)
    l1b_data  = load_mersi_l1b(l1b_path) if l1b_path else None
    rgb       = l1b_data["rgb"] if l1b_data else None
    lat       = recal_data["lat"]
    lon       = recal_data["lon"]

    import re
    import os
    m = re.search(r'(\d{8})_(\d{4})', os.path.basename(onboard_path))
    date_str = ""
    if m:
        d, t = m.group(1), m.group(2)
        date_str = f"{d[:4]}-{d[4:6]}-{d[6:8]}  {t[:2]}:{t[2:]} UTC"

    make_figure1_from_arrays(
        recal_clm   = recal_data["clm"],
        onboard_clm = onboard_data["clm"],
        lat         = lat,
        lon         = lon,
        rgb         = rgb,
        output      = output,
        date_str    = date_str,
        step        = step,
    )
