"""
FY-3D MERSI-II Cloud Mask Visualization  —  Nature-journal style
=================================================================
Four-panel figure: (a) RGB true colour | (b) Recalibration CLM
                   (c) Onboard CLM     | (d) Δ confidence heatmap

Design targets:
  • Nature / Science figure standards (300 dpi, 174 mm wide for full-page)
  • Helvetica Neue / Arial (Nature house font), 6–8 pt body, 7 pt labels
  • Muted, publication-quality colour palette
  • Shared cartopy basemap (coastlines, borders, gridlines)
  • Panel labels (a)(b)(c)(d) bold, top-left, outside axes
  • Stats printed to stdout BEFORE any rendering (unchanged)

Usage — identical to visualize_clm.py:
    python visualize_clm_nature.py \\
        --onboard /path/to/onboard_CLM.h5 \\
        --recal   /path/to/recal_CLM.h5 \\
        --output  figure1.png

    python visualize_clm_nature.py \\
        --data_dir /data/Data_yuq/fy3_cloud/recal_test/ \\
        --output_dir ./clm_plots_nature/

Integration with backend_compare_and_viz.py:
    Pass rgb=make_rgb(pxldat), recal_clm, onboard_clm arrays directly via
    the Python API:

        from visualize_clm_nature import visualize_from_arrays
        visualize_from_arrays(
            rgb=rgb_array,           # (n_elem, n_line, 3) float32 [0,1]
            recal_clm=py_cm,         # (n_elem, n_line) int, values 0-3
            onboard_clm=native_cm,   # (n_elem, n_line) int, values 0-3
            lat=lat, lon=lon,
            output='figure1.png',
            date_str='2023-06-06  14:40 UTC',
        )
"""

from __future__ import annotations
import argparse
import os
import re
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# ── L1b root ────────────────────────────────────────────────────
MERSI_ROOT = Path('/data/Data_yuq/mersi')


# ════════════════════════════════════════════════════════════════
# Nature figure global style
# ════════════════════════════════════════════════════════════════

# Full-page width in Nature = 174 mm; half-page = 85 mm
# We output 174 mm × ~120 mm at 300 dpi
NATURE_FULLWIDTH_IN = 174 / 25.4   # ~6.85 in
NATURE_DPI = 300

# Nature house font stack
FONT_FAMILY = ["Helvetica Neue", "Arial", "DejaVu Sans", "sans-serif"]

# Colour palette  ── muted, print-safe
# Cloud mask: 4-class discrete (CMYK-friendly)
CLM_HEX = {
    0: "#2166AC",   # Cloudy          – deep blue
    1: "#74ADD1",   # Prob. Cloudy    – light blue
    2: "#FEE090",   # Prob. Clear     – pale amber
    3: "#D73027",   # Confident Clear – brick red
}
CLM_LABEL = {
    0: "Cloudy",
    1: "Prob. Cloudy",
    2: "Prob. Clear",
    3: "Conf. Clear",
}
CLM_CMAP = mcolors.ListedColormap([CLM_HEX[v] for v in range(4)])
CLM_NORM = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], CLM_CMAP.N)

# Diverging colourmap for Δ confidence panel
DELTA_CMAP = "RdBu_r"

# Basemap style
LAND_COLOR   = "#F7F4EF"
OCEAN_COLOR  = "#EBF2F8"
LAKE_COLOR   = "#DDE8F2"
COAST_COLOR  = "#666666"
BORDER_COLOR = "#AAAAAA"
GRID_COLOR   = "#BBBBBB"


def _apply_nature_rcparams():
    """Set matplotlib rcParams to Nature style."""
    plt.rcParams.update({
        "font.family":         "sans-serif",
        "font.sans-serif":     FONT_FAMILY,
        "font.size":           7,
        "axes.titlesize":      7,
        "axes.labelsize":      6,
        "xtick.labelsize":     5,
        "ytick.labelsize":     5,
        "legend.fontsize":     6,
        "legend.framealpha":   0.9,
        "legend.edgecolor":    "#CCCCCC",
        "legend.handlelength": 1.2,
        "legend.handletextpad":0.4,
        "axes.linewidth":      0.5,
        "xtick.major.width":   0.5,
        "ytick.major.width":   0.5,
        "xtick.major.size":    2.0,
        "ytick.major.size":    2.0,
        "figure.dpi":          NATURE_DPI,
        "savefig.dpi":         NATURE_DPI,
        "savefig.bbox":        "tight",
        "savefig.facecolor":   "white",
        "pdf.fonttype":        42,   # embed fonts (Nature requirement)
        "ps.fonttype":         42,
        "lines.linewidth":     0.8,
        "patch.linewidth":     0.5,
    })


# ════════════════════════════════════════════════════════════════
# CLM color scheme  (same semantics as visualize_clm.py)
#
#  class 0 → Cloudy            deep blue
#  class 1 → Prob. Cloudy      light blue
#  class 2 → Prob. Clear       pale amber
#  class 3 → Confident Clear   brick red
#
#  NOTE: class 2/3 encode confidence only; no land/water split
#  is performed in the decode layer.  If Cloud_Mask_Value in the
#  HDF5 already encodes land/water, relabel accordingly.
# ════════════════════════════════════════════════════════════════


# ════════════════════════════════════════════════════════════════
# Distribution printing  (unchanged from visualize_clm.py)
# ════════════════════════════════════════════════════════════════

def print_clm_distribution(label: str, clm: np.ndarray) -> None:
    """Print pixel category distribution to stdout BEFORE any plotting."""
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
    print(sep)
    print()


# ════════════════════════════════════════════════════════════════
# Data I/O  (unchanged logic, adapted from visualize_clm.py)
# ════════════════════════════════════════════════════════════════

def decode_cloud_mask_from_bitmask(cm_raw: np.ndarray) -> np.ndarray:
    byte0 = cm_raw[:, :, 0]
    b0 = (byte0 >> 0) & 1
    b1 = (byte0 >> 1) & 1
    b2 = (byte0 >> 2) & 1
    result = np.full(byte0.shape, -1, dtype=np.int32)
    processed = b0 == 1
    result[processed & (b2 == 0) & (b1 == 0)] = 0
    result[processed & (b2 == 0) & (b1 == 1)] = 1
    result[processed & (b2 == 1) & (b1 == 0)] = 2
    result[processed & (b2 == 1) & (b1 == 1)] = 3
    return result


def derive_clm_from_confidence(conf: np.ndarray) -> np.ndarray:
    clm = np.full(conf.shape, -1, dtype=np.int32)
    valid = np.isfinite(conf)
    clm[valid & (conf < 0.33)]                         = 0
    clm[valid & (conf >= 0.33) & (conf < 0.66)]        = 1
    clm[valid & (conf >= 0.66) & (conf < 0.90)]        = 2
    clm[valid & (conf >= 0.90)]                         = 3
    return clm


def load_clm_data(path: str) -> dict:
    with h5py.File(path, "r") as f:
        grp = "Cloud_Mask_1km"
        lat = f[grp]["Latitude"][:]
        lon = f[grp]["Longitude"][:]
        clm = f[grp]["Cloud_Mask_Value"][:]
        if np.all(clm == 0) and "Cloud_Mask" in f[grp]:
            decoded = decode_cloud_mask_from_bitmask(f[grp]["Cloud_Mask"][:])
            if np.any(decoded != -1):
                clm = decoded
        if np.all(clm == 0) and "Confidence" in f[grp]:
            derived = derive_clm_from_confidence(f[grp]["Confidence"][:])
            if np.any(derived >= 0):
                clm = derived
    lat = np.where((lat < -90) | (lat > 90),   np.nan, lat)
    lon = np.where((lon < -180) | (lon > 180), np.nan, lon)
    clm = np.where((clm < 0) | (clm > 3),      -1,    clm)
    return {"lat": lat, "lon": lon, "clm": clm}


def load_rgb_from_l1b(l1b_path: str) -> np.ndarray | None:
    if not os.path.exists(l1b_path):
        return None
    try:
        with h5py.File(l1b_path, "r") as f:
            vis_250 = f["Data/EV_250_Aggr.1KM_RefSB"][:].astype(np.float64)
            vis_cal = f["Calibration/VIS_Cal_Coeff"][:]
            esd     = f.attrs.get("EarthSun Distance Ratio", 1.0)
        esd2   = float(np.squeeze(esd)) ** 2
        n_elem = vis_250.shape[2]
        n_line = vis_250.shape[1]
        rgb = np.zeros((n_elem, n_line, 3), dtype=np.float32)
        for i, b in enumerate([0, 1, 2]):
            c0, c1, c2 = vis_cal[b]
            dn   = vis_250[b]
            refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
            rgb[:, :, i] = refl.T.astype(np.float32)
        for ch in range(3):
            band  = rgb[:, :, ch]
            valid = band > 0
            if valid.any():
                p2, p98 = np.percentile(band[valid], [2, 98])
                rgb[:, :, ch] = np.clip((band - p2) / (p98 - p2 + 1e-10), 0, 1) if p98 > p2 else np.zeros_like(band)
            else:
                rgb[:, :, ch] = 0
        return rgb
    except Exception:
        return None


def find_l1b_path(clm_path: str, mersi_root: Path = MERSI_ROOT) -> str | None:
    m = re.search(r"(\d{8})_(\d{4})", os.path.basename(clm_path))
    if not m:
        return None
    date_str, time_tag = m.group(1), m.group(2)
    p = mersi_root / date_str / f"FY3D_MERSI_GBAL_L1_{date_str}_{time_tag}_1000M_MS.HDF"
    return str(p) if p.exists() else None


def _subsample(arr, step):
    return arr[::step, ::step]


# ════════════════════════════════════════════════════════════════
# Basemap helpers
# ════════════════════════════════════════════════════════════════

def _make_geo_ax(fig, spec):
    """Cartopy PlateCarree axis with Nature-style basemap."""
    ax = fig.add_subplot(spec, projection=ccrs.PlateCarree())
    ax.set_facecolor(OCEAN_COLOR)
    ax.add_feature(cfeature.OCEAN,     facecolor=OCEAN_COLOR,  zorder=0)
    ax.add_feature(cfeature.LAND,      facecolor=LAND_COLOR,   zorder=1)
    ax.add_feature(cfeature.LAKES,     facecolor=LAKE_COLOR,   zorder=1)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor=COAST_COLOR,  zorder=4)
    ax.add_feature(cfeature.BORDERS,   linewidth=0.25, edgecolor=BORDER_COLOR,
                   linestyle="--", zorder=4)
    return ax


def _add_gridlines(ax):
    gl = ax.gridlines(draw_labels=True, linewidth=0.25,
                      color=GRID_COLOR, alpha=0.8, linestyle=":")
    gl.top_labels   = False
    gl.right_labels = False
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {"size": 5, "color": "#444444"}
    gl.ylabel_style = {"size": 5, "color": "#444444"}


def _get_extent(lat, lon, clm, step=4, pad=1.5):
    la = _subsample(lat, step)
    lo = _subsample(lon, step)
    cl = _subsample(clm, step) if clm is not None else np.ones_like(la)
    mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0 if clm is not None else True)
    if not mask.any():
        return None
    return [lo[mask].min() - pad, lo[mask].max() + pad,
            la[mask].min() - pad, la[mask].max() + pad]


def _panel_label(ax, letter, fontsize=8):
    """Bold panel label (a), (b), … in top-left outside the axes."""
    ax.text(-0.04, 1.04, f"({letter})",
            transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold",
            va="bottom", ha="right",
            fontfamily="sans-serif")


# ════════════════════════════════════════════════════════════════
# Individual panel renderers
# ════════════════════════════════════════════════════════════════

def _plot_rgb(ax, lat, lon, rgb, step=4):
    """Scatter RGB true-colour pixels onto geo axis."""
    la = _subsample(lat, step)
    lo = _subsample(lon, step)
    rg = _subsample(rgb, step)
    mask = np.isfinite(la) & np.isfinite(lo)
    if not mask.any():
        ax.text(0.5, 0.5, "No L1b\navailable",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=6, color="#666666")
        return
    laf, lof, rgf = la[mask], lo[mask], rg[mask]
    if len(laf) > 180_000:
        rng = np.random.default_rng(0)
        idx = rng.choice(len(laf), 180_000, replace=False)
        laf, lof, rgf = laf[idx], lof[idx], rgf[idx]
    ax.scatter(lof, laf, c=rgf, s=0.12, linewidths=0,
               transform=ccrs.PlateCarree(), zorder=3, edgecolors="none",
               rasterized=True)


def _plot_rgb_placeholder(ax, lat, lon, clm, step=4):
    """Grey footprint when no L1b is available."""
    la = _subsample(lat, step)
    lo = _subsample(lon, step)
    cl = _subsample(clm, step)
    mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
    ax.scatter(lo[mask], la[mask], c="#CCCCCC", s=0.15, linewidths=0,
               transform=ccrs.PlateCarree(), zorder=3, alpha=0.6,
               rasterized=True)
    ax.text(0.5, 0.5, "No L1b file\n(CLM footprint)",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=6, color="#555555",
            bbox=dict(fc="white", ec="#CCCCCC", pad=3, lw=0.4, alpha=0.85))


def _plot_clm(ax, lat, lon, clm, step=4):
    """Scatter CLM-coloured pixels onto geo axis."""
    la = _subsample(lat, step)
    lo = _subsample(lon, step)
    cl = _subsample(clm, step)
    mask = np.isfinite(la) & np.isfinite(lo) & (cl >= 0)
    la, lo, cl = la[mask], lo[mask], cl[mask]

    # Draw clear first, cloudy on top for visual clarity
    for v in [2, 3, 1, 0]:
        m = cl == v
        if not m.any():
            continue
        ax.scatter(lo[m], la[m], c=CLM_HEX[v], s=0.15, linewidths=0,
                   transform=ccrs.PlateCarree(), zorder=3, edgecolors="none",
                   rasterized=True)


def _add_clm_colorbar(fig, ax, orientation="vertical", shrink=0.85, pad=0.03):
    """Attach a compact CLM colorbar to an axis."""
    sm = plt.cm.ScalarMappable(cmap=CLM_CMAP, norm=CLM_NORM)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation=orientation,
                        shrink=shrink, pad=pad, aspect=30,
                        ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(
        [CLM_LABEL[v] for v in range(4)],
        fontsize=5.5)
    cbar.ax.tick_params(length=2, width=0.4)
    cbar.outline.set_linewidth(0.4)
    return cbar


def _stats_annotation(ax, clm, fontsize=5.5):
    """Inset text box with per-category pixel percentages."""
    total = int((clm >= 0).sum())
    if total == 0:
        return
    lines = []
    for v in range(4):
        cnt = int((clm == v).sum())
        pct = 100.0 * cnt / total
        lines.append(f"{CLM_LABEL[v]}: {pct:.1f}%")
    txt = "\n".join(lines)
    ax.text(0.015, 0.015, txt,
            transform=ax.transAxes, va="bottom", ha="left",
            fontsize=fontsize, linespacing=1.5,
            bbox=dict(fc="white", ec="#AAAAAA", pad=3,
                      lw=0.4, alpha=0.88, boxstyle="round,pad=0.3"))


# ════════════════════════════════════════════════════════════════
# Main figure composer
# ════════════════════════════════════════════════════════════════

def _build_figure(recal_data: dict, onboard_data: dict,
                  rgb: np.ndarray | None,
                  date_str: str,
                  output: str,
                  step: int = 4):
    """
    Compose the 4-panel Nature-style figure and save.

    Panels:
      (a) RGB true colour
      (b) Recalibration CLM
      (c) Onboard CLM
      (d) Confidence difference  (recal − onboard, float, if both available)
    """
    _apply_nature_rcparams()

    lat = recal_data["lat"]
    lon = recal_data["lon"]
    extent = _get_extent(lat, lon, recal_data["clm"], step=step)

    # ── Figure size: full Nature width, 2×2 grid ──────────────────
    fig_w = NATURE_FULLWIDTH_IN
    fig_h = fig_w * 0.62          # aspect tuned for swath data
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    gs = GridSpec(2, 2,
                  figure=fig,
                  left=0.06, right=0.94,
                  top=0.90,  bottom=0.06,
                  wspace=0.18, hspace=0.32)

    axes_specs = [gs[0, 0], gs[0, 1], gs[1, 0], gs[1, 1]]
    letters    = ["a", "b", "c", "d"]
    titles     = [
        "RGB true colour (MERSI-II)",
        "Cloud mask — recalibration",
        "Cloud mask — onboard calibration",
        "Class agreement map",
    ]

    axs = [_make_geo_ax(fig, sp) for sp in axes_specs]

    # ── (a) RGB ───────────────────────────────────────────────────
    ax = axs[0]
    if rgb is not None:
        _plot_rgb(ax, lat, lon, rgb, step=step)
    else:
        _plot_rgb_placeholder(ax, lat, lon, recal_data["clm"], step=step)
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    _add_gridlines(ax)

    # ── (b) Recalibration CLM ─────────────────────────────────────
    ax = axs[1]
    _plot_clm(ax, lat, lon, recal_data["clm"], step=step)
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    _add_gridlines(ax)
    _stats_annotation(ax, recal_data["clm"])
    _add_clm_colorbar(fig, ax, shrink=0.82)

    # ── (c) Onboard CLM ───────────────────────────────────────────
    ax = axs[2]
    _plot_clm(ax, onboard_data["lat"], onboard_data["lon"],
              onboard_data["clm"], step=step)
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    _add_gridlines(ax)
    _stats_annotation(ax, onboard_data["clm"])
    _add_clm_colorbar(fig, ax, shrink=0.82)

    # ── (d) Agreement / difference map ───────────────────────────
    ax = axs[3]
    # Show agreement: same class = 1, differ = 0; use recal lat/lon
    la = _subsample(lat, step)
    lo = _subsample(lon, step)
    cl_r = _subsample(recal_data["clm"], step)
    cl_o = _subsample(onboard_data["clm"], step)
    mask = (np.isfinite(la) & np.isfinite(lo) & (cl_r >= 0) & (cl_o >= 0))

    diff = cl_r.astype(np.float32) - cl_o.astype(np.float32)   # −3 … +3
    diff_plot = np.where(mask, diff, np.nan)

    # Scatter with diverging cmap centred on 0
    vmax = 3.0
    d_cmap = plt.get_cmap(DELTA_CMAP)
    d_norm = mcolors.Normalize(vmin=-vmax, vmax=vmax)

    colors_d = d_cmap(d_norm(diff_plot[mask]))
    ax.scatter(lo[mask], la[mask], c=colors_d, s=0.15, linewidths=0,
               transform=ccrs.PlateCarree(), zorder=3, edgecolors="none",
               rasterized=True)

    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    _add_gridlines(ax)

    # Colorbar for difference
    sm_d = plt.cm.ScalarMappable(cmap=d_cmap, norm=d_norm)
    sm_d.set_array([])
    cbar_d = fig.colorbar(sm_d, ax=ax, orientation="vertical",
                          shrink=0.82, pad=0.03, aspect=30)
    cbar_d.set_label("Recal − Onboard (class)", fontsize=5.5, labelpad=3)
    cbar_d.ax.tick_params(labelsize=5, length=2, width=0.4)
    cbar_d.outline.set_linewidth(0.4)

    # Agreement statistics inset
    agree_pct = 100 * np.mean(cl_r[mask] == cl_o[mask])
    ax.text(0.015, 0.015,
            f"Agreement: {agree_pct:.1f}%",
            transform=ax.transAxes, va="bottom", ha="left",
            fontsize=5.5,
            bbox=dict(fc="white", ec="#AAAAAA", pad=3,
                      lw=0.4, alpha=0.88, boxstyle="round,pad=0.3"))

    # ── Shared titles, panel labels, suptitle ─────────────────────
    for ax, ltr, ttl in zip(axs, letters, titles):
        ax.set_title(ttl, fontsize=6.5, fontweight="bold",
                     color="#222222", pad=4, loc="left")
        _panel_label(ax, ltr, fontsize=8)

    # Suptitle: satellite / date metadata
    lat_c = float(np.nanmedian(lat))
    lon_c = float(np.nanmedian(lon))
    fig.suptitle(
        f"FY-3D MERSI-II cloud mask comparison   {date_str}   "
        f"centre {lat_c:.1f}°N  {lon_c:.1f}°E",
        fontsize=7.5, fontweight="normal", color="#333333",
        y=0.96)

    # ── Save ──────────────────────────────────────────────────────
    fig.savefig(output, dpi=NATURE_DPI, bbox_inches="tight",
                facecolor="white", pil_kwargs={"compression": 6})
    plt.close(fig)
    print(f"[SAVE] {output}  ({os.path.getsize(output)/1024:.0f} KB)")


# ════════════════════════════════════════════════════════════════
# Public Python API  (for backend_compare_and_viz.py integration)
# ════════════════════════════════════════════════════════════════

def visualize_from_arrays(
    rgb: np.ndarray | None,
    recal_clm: np.ndarray,
    onboard_clm: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    output: str = "figure1.png",
    date_str: str = "",
    step: int = 4,
):
    """
    Direct-array entry point — no HDF5 files needed.

    Parameters
    ----------
    rgb         : (n_elem, n_line, 3) float32 [0, 1]  or None
    recal_clm   : (n_elem, n_line) int  values 0-3 / -1
    onboard_clm : (n_elem, n_line) int  values 0-3 / -1
    lat, lon    : (n_elem, n_line) float64 geolocation
    output      : output PNG path
    date_str    : e.g. "2023-06-06  14:40 UTC"
    step        : subsampling step for scatter (4 = every 4th pixel)
    """
    # Print stats BEFORE any plotting
    print(f"\n[INFO] Orbit: {os.path.basename(output)}")
    print_clm_distribution("② Recalibration CLM", recal_clm)
    print_clm_distribution("③ Onboard CLM",        onboard_clm)

    recal_data   = {"lat": lat, "lon": lon, "clm": recal_clm}
    onboard_data = {"lat": lat, "lon": lon, "clm": onboard_clm}

    _build_figure(recal_data, onboard_data, rgb, date_str, output, step=step)


# ════════════════════════════════════════════════════════════════
# HDF5-file entry point  (mirror of visualize_clm.py)
# ════════════════════════════════════════════════════════════════

def visualize(onboard_path, recal_path, output="figure1.png",
              step=4, mersi_root=MERSI_ROOT):
    recal_data   = load_clm_data(recal_path)
    onboard_data = load_clm_data(onboard_path)

    # Print stats BEFORE any plotting
    print(f"\n[INFO] Orbit: {os.path.basename(onboard_path)}")
    print_clm_distribution("② Recalibration CLM", recal_data["clm"])
    print_clm_distribution("③ Onboard CLM",        onboard_data["clm"])

    l1b_path = find_l1b_path(recal_path, mersi_root)
    rgb = load_rgb_from_l1b(l1b_path) if l1b_path else None

    m = re.search(r"(\d{8})_(\d{4})", os.path.basename(onboard_path))
    date_str = ""
    if m:
        d, t = m.group(1), m.group(2)
        date_str = f"{d[:4]}-{d[4:6]}-{d[6:8]}  {t[:2]}:{t[2:]} UTC"

    _build_figure(recal_data, onboard_data, rgb, date_str, output, step=step)


# ════════════════════════════════════════════════════════════════
# Batch
# ════════════════════════════════════════════════════════════════

def find_orbit_pairs(data_dir: str) -> list:
    pairs = []
    for onboard in sorted(Path(data_dir).rglob("*_CLM_CLA.h5")):
        if "_recal" in onboard.name:
            continue
        recal = onboard.with_name(onboard.name.replace("_CLM_CLA.h5", "_CLM_CLA_recal.h5"))
        if not recal.exists():
            print(f"[WARN] Missing recal for {onboard.name}, skipping")
            continue
        m = re.search(r"(\d{8})_(\d{4})", onboard.name)
        if not m:
            continue
        pairs.append((str(onboard), str(recal), m.group(1), m.group(2)))
    return pairs


def batch_visualize(data_dir: str, output_dir: str, step=4, mersi_root=MERSI_ROOT):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs = find_orbit_pairs(data_dir)
    if not pairs:
        print(f"[ERROR] No CLM pairs found in {data_dir}")
        return
    print(f"[INFO] Found {len(pairs)} orbit pair(s)\n")
    for onboard, recal, date_str, time_str in pairs:
        outpath = output_dir / f"CLM_nature_{date_str}_{time_str}.png"
        if outpath.exists():
            print(f"[SKIP] {outpath.name}")
            continue
        print(f"[PLOT] {date_str} {time_str} ...")
        try:
            visualize(onboard, recal, str(outpath), step=step, mersi_root=mersi_root)
        except Exception as e:
            print(f"  [ERROR] {e}")
    print(f"\n[DONE] Output in {output_dir}/")


# ════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="FY-3D MERSI-II CLM — Nature-journal style visualization")
    parser.add_argument("--onboard")
    parser.add_argument("--recal")
    parser.add_argument("--output",     default="figure1.png")
    parser.add_argument("--data_dir")
    parser.add_argument("--output_dir", default="./clm_plots_nature")
    parser.add_argument("--step",       type=int, default=4)
    parser.add_argument("--mersi_root", default=str(MERSI_ROOT))
    args = parser.parse_args()

    mersi_root = Path(args.mersi_root)
    if args.data_dir:
        batch_visualize(args.data_dir, args.output_dir,
                        step=args.step, mersi_root=mersi_root)
    elif args.onboard and args.recal:
        visualize(args.onboard, args.recal, args.output,
                  step=args.step, mersi_root=mersi_root)
    else:
        parser.print_help()
        print("\nError: provide --data_dir, or --onboard + --recal.")