"""
figure_2.py — MYD35 (truth) vs MERSI CLM validation figure  (6-panel)
=======================================================================
Plots ONLY the spatially overlapping region between the MERSI swath
and the matched MYD35 granule.

Panel layout  (2 rows × 3 cols):
  (a) MERSI RGB true colour             — overlap footprint
  (b) MYD35 cloud mask (reference)      — native MYD35 grid
  (c) MERSI recalibration CLM           — overlap region
  (d) MERSI recal − MYD35  (diff)       — overlap region
  (e) MERSI onboard − MYD35 (diff)      — overlap region
  (f) Confusion matrices (recal & onboard vs MYD35 truth)

v2 bug-fix notes
----------------
- The confusion-matrix panel (f) previously hand-placed two 4×4 grids
  with raw axes-coordinate math; long rotated row/column labels bled
  outside their allotted column and the two tables visually overlapped.
  It is now built with matplotlib's `ax.table()`, which guarantees
  fixed, non-overlapping cell geometry regardless of label length.
- Per-class % stats are no longer drawn as floating boxes ON TOP of
  the maps (this caused the text to overlap colorbars on narrow polar
  swaths). They are now one-line captions below each panel via
  plot_utils.panel_caption().
- Figure canvas is sized per-panel rather than a fixed total width, so
  every map gets enough room regardless of how elongated the swath is.

Usage
-----
from figure_2 import make_figure2

make_figure2(
    mersi_lat       = lat,
    mersi_lon       = lon,
    mersi_rgb       = rgb,             # (H,W,3) float32 or None
    recal_clm       = recal_clm,       # (H,W) int32 on MERSI grid
    onboard_clm     = onboard_clm,     # (H,W) int32 on MERSI grid
    myd35_data      = myd35_dict,      # output of io_myd35.load_best_myd35_for_mersi()
    output          = "figure2.png",
    mersi_date_str  = "2023-06-06  14:40 UTC",
    step            = 4,
)
"""

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

from plot_utils import (
    apply_nature_style, PANEL_WIDTH_IN, PANEL_HEIGHT_IN,
    make_geo_ax_with_caption, add_gridlines, panel_label, panel_title,
    panel_caption,
    plot_rgb, plot_rgb_placeholder, plot_clm, plot_diff,
    add_clm_colorbar, add_diff_colorbar,
    stats_caption_text, agreement_caption_text,
    get_extent, subsample, save_figure,
    CLM_LABEL_SHORT,
)
from io_mersi import print_clm_distribution


# ─────────────────────────────────────────────────────────────────────────────
# Overlap mask helpers
# ─────────────────────────────────────────────────────────────────────────────

def compute_overlap_extent(
    mersi_lat:  np.ndarray,
    mersi_lon:  np.ndarray,
    mersi_clm:  np.ndarray,
    myd35_lat:  np.ndarray,
    myd35_lon:  np.ndarray,
    myd35_clm:  np.ndarray,
    step:       int   = 4,
    pad:        float = 0.5,
) -> list | None:
    """Bounding-box intersection of two swath footprints."""
    def valid_bbox(lat, lon, clm, step):
        la = subsample(lat, step)
        lo = subsample(lon, step)
        cl = subsample(clm, step)
        m  = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
        if not m.any():
            return None
        return lo[m].min(), lo[m].max(), la[m].min(), la[m].max()

    bb1 = valid_bbox(mersi_lat, mersi_lon, mersi_clm, step)
    bb2 = valid_bbox(myd35_lat, myd35_lon, myd35_clm, step)
    if bb1 is None or bb2 is None:
        return None

    lon_min = max(bb1[0], bb2[0]) - pad
    lon_max = min(bb1[1], bb2[1]) + pad
    lat_min = max(bb1[2], bb2[2]) - pad
    lat_max = min(bb1[3], bb2[3]) + pad

    if lon_max <= lon_min or lat_max <= lat_min:
        return None
    return [float(lon_min), float(lon_max), float(lat_min), float(lat_max)]


def overlap_mask_on_mersi_grid(
    mersi_lat:           np.ndarray,
    mersi_lon:           np.ndarray,
    myd35_resampled_clm: np.ndarray,
) -> np.ndarray:
    """Boolean mask: MERSI pixels with a valid co-located MYD35 value."""
    return (
        np.isfinite(mersi_lat) &
        np.isfinite(mersi_lon) &
        (myd35_resampled_clm >= 0)
    )


# ─────────────────────────────────────────────────────────────────────────────
# Validation statistics
# ─────────────────────────────────────────────────────────────────────────────

def compute_validation_stats(
    mersi_clm: np.ndarray,   # (H,W) on MERSI grid
    myd35_clm: np.ndarray,   # (H,W) resampled to MERSI grid
    label:     str = "MERSI",
) -> dict:
    """
    Compute agreement metrics between MERSI CLM and MYD35 truth.
    Returns dict with agree_pct, pod, far, csi, hss, conf_mat, etc.
    """
    mask = (mersi_clm >= 0) & (myd35_clm >= 0)
    if not mask.any():
        return {}

    a = mersi_clm[mask]   # predicted
    b = myd35_clm[mask]   # reference (truth)
    n = len(a)

    agree_pct = 100 * np.mean(a == b)

    a_cloud = a <= 1
    b_cloud = b <= 1
    TP = int(( a_cloud &  b_cloud).sum())
    FP = int(( a_cloud & ~b_cloud).sum())
    FN = int((~a_cloud &  b_cloud).sum())
    TN = int((~a_cloud & ~b_cloud).sum())

    pod = 100 * TP / (TP + FN + 1e-9)
    far = 100 * FP / (TP + FP + 1e-9)
    csi = 100 * TP / (TP + FP + FN + 1e-9)

    expected = ((TP + FP) * (TP + FN) + (TN + FP) * (TN + FN)) / (n + 1e-9)
    hss = (TP + TN - expected) / (n - expected + 1e-9)

    conf_mat = np.zeros((4, 4), dtype=np.int64)
    for i in range(4):
        for j in range(4):
            conf_mat[i, j] = int(((a == i) & (b == j)).sum())

    stats = dict(
        n_pixels=int(n), agree_pct=agree_pct, pod=pod, far=far,
        csi=csi, hss=hss, TP=TP, FP=FP, FN=FN, TN=TN, conf_mat=conf_mat,
    )

    sep = "─" * 55
    print(sep)
    print(f"  Validation stats — {label} vs MYD35")
    print(f"  Overlap pixels : {n:>10,}")
    print(f"  Agreement      : {agree_pct:>7.2f}%")
    print(f"  POD (cloud)    : {pod:>7.2f}%")
    print(f"  FAR (cloud)    : {far:>7.2f}%")
    print(f"  CSI            : {csi:>7.2f}%")
    print(f"  HSS            : {hss:>7.4f}")
    print(sep)

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Panel renderers specific to Figure 2
# ─────────────────────────────────────────────────────────────────────────────

def _plot_myd35_clm(ax, myd35_data: dict, step: int = 4) -> None:
    plot_clm(ax, myd35_data["lat"], myd35_data["lon"],
             myd35_data["clm_native"], step=step)


def _validation_diff_panel(
    ax, fig,
    mersi_lat: np.ndarray,
    mersi_lon: np.ndarray,
    mersi_clm: np.ndarray,
    myd35_resampled: np.ndarray,
    step: int = 4,
    title_suffix: str = "",
) -> dict:
    """Plot difference panel + colorbar + caption. Returns validation stats."""
    ov_mask    = overlap_mask_on_mersi_grid(mersi_lat, mersi_lon, myd35_resampled)
    clm_masked = np.where(ov_mask, mersi_clm,       -1)
    myd_masked = np.where(ov_mask, myd35_resampled, -1)

    sm, _, _ = plot_diff(ax, mersi_lat, mersi_lon, clm_masked, myd_masked, step=step)
    add_diff_colorbar(fig, ax, sm, label=f"MERSI{title_suffix} − MYD35 (class)")

    stats = compute_validation_stats(mersi_clm, myd35_resampled,
                                     label=f"MERSI {title_suffix}")
    if stats:
        cap = (f"Agreement {stats['agree_pct']:.1f}%  ·  "
               f"POD {stats['pod']:.1f}%  ·  FAR {stats['far']:.1f}%  ·  "
               f"HSS {stats['hss']:.3f}")
        panel_caption(ax, cap)
    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Confusion-matrix panel  (built with ax.table — no manual coordinate math)
# ─────────────────────────────────────────────────────────────────────────────

def _draw_confusion_panel(
    ax,
    stats_recal:   dict,
    stats_onboard: dict,
) -> None:
    """
    Render two 4×4 confusion-matrix tables (row-normalised %, rows =
    MYD35 truth class, columns = MERSI predicted class) stacked
    vertically using matplotlib's table API, which guarantees cells
    never overlap regardless of label length.
    """
    ax.set_axis_off()

    labels = [CLM_LABEL_SHORT[v] for v in range(4)]
    cmap_conf = plt.get_cmap("Blues")

    def _cell_colors(conf_mat: np.ndarray) -> np.ndarray:
        cm = conf_mat.astype(float)
        row_tot = cm.sum(axis=1, keepdims=True)
        row_tot = np.where(row_tot == 0, 1, row_tot)
        pct = cm / row_tot
        return pct, cmap_conf(0.15 + 0.7 * pct)

    panels = [(stats_recal, "MERSI Recal  vs  MYD35"),
              (stats_onboard, "MERSI Onboard  vs  MYD35")]

    n_valid = sum(1 for s, _ in panels if s and "conf_mat" in s)
    if n_valid == 0:
        ax.text(0.5, 0.5, "No validation statistics available",
                ha="center", va="center", fontsize=9, color="#666666",
                transform=ax.transAxes)
        return

    # Stack the (up to 2) tables vertically, each in its own sub-axes
    # band, so they can never collide horizontally OR vertically.
    band_h = 1.0 / max(n_valid, 1)
    band_i = 0

    for stats, label in panels:
        if not stats or "conf_mat" not in stats:
            continue

        pct, colors = _cell_colors(stats["conf_mat"])

        y_top    = 1.0 - band_i * band_h
        y_bottom = y_top - band_h
        # Reserve top margin in the band for the title text
        table_top    = y_top - band_h * 0.18
        table_bottom = y_bottom + band_h * 0.06

        ax.text(0.5, y_top - band_h * 0.04, label,
                transform=ax.transAxes, ha="center", va="top",
                fontsize=9, fontweight="bold", color="#222222")

        cell_text = [[f"{pct[i, j]*100:.0f}%" for j in range(4)]
                     for i in range(4)]

        tbl = ax.table(
            cellText=cell_text,
            cellColours=colors,
            rowLabels=[f"truth: {l}" for l in labels],
            colLabels=[f"pred: {l}" for l in labels],
            bbox=[0.16, table_bottom, 0.80, table_top - table_bottom],
        )
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7.5)
        for (row, col), cell in tbl.get_celld().items():
            cell.set_edgecolor("#BBBBBB")
            cell.set_linewidth(0.4)
            if row == 0 or col == -1:
                cell.set_text_props(fontsize=7, color="#444444")
                cell.set_facecolor("#F5F5F5")
            else:
                # Darker fill → white text for contrast
                v = pct[row - 1, col]
                cell.set_text_props(
                    color="white" if v > 0.55 else "#1A1A1A")

        agree = stats.get("agree_pct", 0)
        hss   = stats.get("hss", 0)
        ax.text(0.5, table_bottom - band_h * 0.02,
                f"Overall agreement {agree:.1f}%   ·   HSS {hss:.3f}",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7.5, color="#444444")

        band_i += 1


# ─────────────────────────────────────────────────────────────────────────────
# Core figure builder
# ─────────────────────────────────────────────────────────────────────────────

def _build_figure2(
    mersi_lat:      np.ndarray,
    mersi_lon:      np.ndarray,
    mersi_rgb:      np.ndarray | None,
    recal_clm:      np.ndarray,
    onboard_clm:    np.ndarray,
    myd35_data:     dict,
    mersi_date_str: str,
    output:         str,
    step:           int = 4,
) -> dict:
    """Build and save Figure 2.  Returns dict of validation stats."""
    apply_nature_style()

    myd35_resampled = myd35_data["clm_resampled"]
    myd_dt = myd35_data.get("dt")
    dt_min = myd35_data.get("dt_diff_min", 0)

    overlap_extent = compute_overlap_extent(
        mersi_lat, mersi_lon, recal_clm,
        myd35_data["lat"], myd35_data["lon"], myd35_data["clm_native"],
        step=step)
    if overlap_extent is None:
        overlap_extent = get_extent(mersi_lat, mersi_lon, recal_clm, step=step)

    fig_w = PANEL_WIDTH_IN * 3 + 1.6
    fig_h = PANEL_HEIGHT_IN * 2 + 0.7
    fig   = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    gs = GridSpec(2, 3, figure=fig,
                  left=0.04, right=0.92,
                  top=0.91,  bottom=0.05,
                  wspace=0.40, hspace=0.36)

    panel_specs = [gs[0, 0], gs[0, 1], gs[0, 2],
                   gs[1, 0], gs[1, 1], gs[1, 2]]
    letters = ["a", "b", "c", "d", "e", "f"]
    titles  = [
        "MERSI RGB true colour",
        "MYD35 truth",
        "MERSI recal CLM",
        "Recal − MYD35 diff",
        "Onboard − MYD35 diff",
        "Confusion matrices",
    ]

    axs_geo = [make_geo_ax_with_caption(fig, panel_specs[i]) for i in range(5)]
    ax_conf = fig.add_subplot(panel_specs[5])

    # ── (a) MERSI RGB ───────────────────────────────────────────────
    ax = axs_geo[0]
    if mersi_rgb is not None:
        plot_rgb(ax, mersi_lat, mersi_lon, mersi_rgb, step=step)
    else:
        plot_rgb_placeholder(ax, mersi_lat, mersi_lon, recal_clm, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (b) MYD35 CLM on native grid ────────────────────────────────
    ax = axs_geo[1]
    _plot_myd35_clm(ax, myd35_data, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    myd_dt_str = f"Δt = {dt_min:.1f} min" if myd_dt is not None else ""
    panel_caption(ax, stats_caption_text(myd35_data["clm_native"]) +
                  (f"   ({myd_dt_str})" if myd_dt_str else ""))

    # ── (c) MERSI recalibration CLM (overlap region) ────────────────
    ax = axs_geo[2]
    ov_mask  = overlap_mask_on_mersi_grid(mersi_lat, mersi_lon, myd35_resampled)
    recal_ov = np.where(ov_mask, recal_clm, -1)
    plot_clm(ax, mersi_lat, mersi_lon, recal_ov, step=step)
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)
    add_clm_colorbar(fig, ax)
    panel_caption(ax, stats_caption_text(recal_ov))

    # ── (d) Recal − MYD35 diff ──────────────────────────────────────
    ax = axs_geo[3]
    stats_recal = _validation_diff_panel(
        ax, fig, mersi_lat, mersi_lon, recal_clm, myd35_resampled,
        step=step, title_suffix=" recal")
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (e) Onboard − MYD35 diff ─────────────────────────────────────
    ax = axs_geo[4]
    stats_onboard = _validation_diff_panel(
        ax, fig, mersi_lat, mersi_lon, onboard_clm, myd35_resampled,
        step=step, title_suffix=" onboard")
    if overlap_extent:
        ax.set_extent(overlap_extent)
    add_gridlines(ax)

    # ── (f) Confusion matrices ───────────────────────────────────────
    _draw_confusion_panel(ax_conf, stats_recal, stats_onboard)

    # ── Titles, labels, suptitle ─────────────────────────────────────
    all_axes = axs_geo + [ax_conf]
    for ax, ltr, ttl in zip(all_axes, letters, titles):
        panel_title(ax, ttl)
        panel_label(ax, ltr)

    lat_c = float(np.nanmedian(mersi_lat))
    lon_c = float(np.nanmedian(mersi_lon))
    fig.suptitle(
        f"FY-3D MERSI-II  vs  MYD35 (truth)   "
        f"MERSI: {mersi_date_str}   "
        f"centre {lat_c:.1f}°N {lon_c:.1f}°E",
        fontsize=10.5, fontweight="normal", color="#333333", y=0.98)

    save_figure(fig, output)
    return {"recal": stats_recal, "onboard": stats_onboard}


# ─────────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────────

def make_figure2(
    mersi_lat:      np.ndarray,
    mersi_lon:      np.ndarray,
    recal_clm:      np.ndarray,
    onboard_clm:    np.ndarray,
    myd35_data:     dict,
    mersi_rgb:      np.ndarray | None = None,
    output:         str               = "figure2.png",
    mersi_date_str: str               = "",
    step:           int               = 4,
) -> dict:
    """
    Build Figure 2: MYD35 validation comparison.

    Parameters
    ----------
    myd35_data : dict
        Output of io_myd35.load_best_myd35_for_mersi().
        Must contain keys: clm_native, clm_resampled, lat, lon, dt, dt_diff_min.

    Returns
    -------
    dict with 'recal' and 'onboard' validation stats.
    """
    print(f"\n[FIG2] Building MYD35 validation figure → {output}")
    print_clm_distribution("MERSI recal CLM",   recal_clm)
    print_clm_distribution("MERSI onboard CLM", onboard_clm)
    print_clm_distribution("MYD35 CLM (native)", myd35_data["clm_native"])
    print_clm_distribution("MYD35 CLM (resampled to MERSI grid)",
                            myd35_data["clm_resampled"])

    return _build_figure2(
        mersi_lat, mersi_lon, mersi_rgb,
        recal_clm, onboard_clm, myd35_data,
        mersi_date_str, output, step=step,
    )
