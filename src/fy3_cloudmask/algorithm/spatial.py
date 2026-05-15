"""Spatial analysis utilities: tview lookup table, 3x3 neighborhood stats, uniformity check.

Port of tview.f, get_regdif.f, get_regstd.f, spatial_var.f, check_reg_uniformity().
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ..constants import BAD_DATA


# ---------------------------------------------------------------------------
# APOLLO 11-12um BTD lookup table (from tview.f)
# ---------------------------------------------------------------------------
# Secant of viewing angle
TVIEW_UTAB = np.array([2.00, 1.75, 1.50, 1.25, 1.00], dtype=np.float64)

# 11um brightness temperature (K)
TVIEW_TTAB = np.array([
    190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0,
    260.0, 270.0, 280.0, 290.0, 300.0, 310.0,
], dtype=np.float64)

# Threshold table (11-12um BTD in K) - 5 rows (angle) x 13 cols (temperature)
TVIEW_TAB = np.array([
    [18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0],
    [16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0],
    [14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0],
    [12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
    [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0],
], dtype=np.float64)


@njit(cache=True)
def tview(secant_vza: float, bt_11um: float) -> float:
    """APOLLO 11-12um BTD lookup table interpolation.

    Bilinear interpolation over a 5x13 table indexed by secant viewing angle
    and 11um brightness temperature.

    Port of tview.f with key=1 (linear interpolation).

    Args:
        secant_vza: Secant of satellite viewing angle (1/cos(VZA)).
        bt_11um: 11um brightness temperature (K).

    Returns:
        Interpolated 11-12um BTD threshold (K). Returns 99.0 if out of range.
    """
    # Clamp inputs to table bounds
    u = secant_vza
    t = bt_11um

    if u < TVIEW_UTAB[4] or u > TVIEW_UTAB[0]:
        return 99.0
    if t < TVIEW_TTAB[0] or t > TVIEW_TTAB[12]:
        return 99.0

    # Find bounding indices for secant angle (descending order)
    iu = 0
    for i in range(4):
        if u >= TVIEW_UTAB[i + 1]:
            iu = i
            break
    if u >= TVIEW_UTAB[0]:
        iu = 3

    # Find bounding indices for temperature (ascending order)
    it = 0
    for i in range(12):
        if t >= TVIEW_TTAB[i] and t < TVIEW_TTAB[i + 1]:
            it = i
            break
    if t >= TVIEW_TTAB[12]:
        it = 11

    # Bilinear interpolation
    u1 = TVIEW_UTAB[iu]
    u2 = TVIEW_UTAB[iu + 1]
    t1 = TVIEW_TTAB[it]
    t2 = TVIEW_TTAB[it + 1]

    du = u2 - u1
    dt = t2 - t1

    if du == 0.0 or dt == 0.0:
        return TVIEW_TAB[iu, it]

    fu = (u - u1) / du
    ft = (t - t1) / dt

    # Interpolate
    result = (TVIEW_TAB[iu, it] * (1.0 - fu) * (1.0 - ft) +
              TVIEW_TAB[iu + 1, it] * fu * (1.0 - ft) +
              TVIEW_TAB[iu, it + 1] * (1.0 - fu) * ft +
              TVIEW_TAB[iu + 1, it + 1] * fu * ft)

    return result


@njit(cache=True)
def get_regional_mean(data_3x3: np.ndarray) -> float:
    """Compute mean of 3x3 neighborhood (excluding center pixel).

    Args:
        data_3x3: 3x3 array of values.

    Returns:
        Mean of the 8 surrounding pixels.
    """
    total = 0.0
    count = 0
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue  # skip center
            val = data_3x3[i, j]
            if val > BAD_DATA + 1.0:
                total += val
                count += 1
    if count == 0:
        return BAD_DATA
    return total / count


@njit(cache=True)
def get_regional_std(data_3x3: np.ndarray) -> float:
    """Compute standard deviation of 3x3 neighborhood (excluding center).

    Args:
        data_3x3: 3x3 array of values.

    Returns:
        Standard deviation of the 8 surrounding pixels.
    """
    mean = get_regional_mean(data_3x3)
    if mean < BAD_DATA + 1.0:
        return BAD_DATA

    total = 0.0
    count = 0
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue
            val = data_3x3[i, j]
            if val > BAD_DATA + 1.0:
                diff = val - mean
                total += diff * diff
                count += 1

    if count < 2:
        return BAD_DATA
    return np.sqrt(total / (count - 1))


@njit(cache=True)
def get_regional_diff(data_3x3: np.ndarray) -> float:
    """Compute max difference between center and 3x3 neighborhood.

    Args:
        data_3x3: 3x3 array of values.

    Returns:
        Maximum absolute difference between center and surrounding pixels.
    """
    center = data_3x3[1, 1]
    if center < BAD_DATA + 1.0:
        return BAD_DATA

    max_diff = 0.0
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue
            val = data_3x3[i, j]
            if val > BAD_DATA + 1.0:
                diff = abs(val - center)
                if diff > max_diff:
                    max_diff = diff

    return max_diff


@njit(cache=True)
def check_reg_uniformity(
    eco_center: int,
    eco_neighbors: np.ndarray,
    snow_center: int,
    snow_neighbors: np.ndarray,
    lsf_center: int,
    lsf_neighbors: np.ndarray,
    is_edge: bool,
    is_snow: bool,
    is_ice: bool,
) -> tuple[bool, bool, bool, bool, int, int, int]:
    """Check 3x3 neighborhood consistency.

    Port of check_reg_uniformity() in fylat_fy3mersi_cloud_mask.f90.

    Args:
        eco_center: Ecosystem type of center pixel.
        eco_neighbors: 3x3 array of ecosystem types.
        snow_center: Snow mask value of center pixel.
        snow_neighbors: 3x3 array of snow mask values.
        lsf_center: Land/sea flag of center pixel.
        lsf_neighbors: 3x3 array of land/sea flags.
        is_edge: Whether pixel is on image edge.
        is_snow: Whether center pixel is snow.
        is_ice: Whether center pixel is ice.

    Returns:
        Tuple of (uniform, is_coast, is_land, is_water, n_land, n_water, n_coast).
    """
    uniform = True
    is_coast = False
    is_land = False
    is_water = False

    # Edge pixels are not uniform
    if is_edge:
        uniform = False

    # Snow/ice pixels are not uniform
    if is_snow or is_ice:
        uniform = False

    # Count surface types in 3x3 box
    n_land = 0
    n_water = 0
    n_coast = 0
    n_other = 0

    for i in range(3):
        for j in range(3):
            lsf = lsf_neighbors[i, j]
            if lsf == 1 or lsf == 4:
                n_land += 1
            elif lsf == 2:
                n_coast += 1
            elif lsf == 0:
                n_water += 1
            elif lsf == 3:
                n_land += 1  # lake is treated as land
            else:
                n_other += 1

    # Check for mixed surface types
    if not (n_land == 9 or n_water == 9 or n_coast == 9):
        if n_other == 0:
            uniform = False

    # Check ecosystem consistency
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue
            if eco_neighbors[i, j] != eco_center:
                uniform = False
                break

    # Check snow mask consistency
    for i in range(3):
        for j in range(3):
            if i == 1 and j == 1:
                continue
            if snow_neighbors[i, j] != snow_center:
                uniform = False
                break

    # Coastline handling
    if n_water + n_coast == 9 and n_water != 9:
        is_coast = True
        is_land = True
        is_water = False
        uniform = False
    elif n_land == 9:
        is_land = True
    elif n_water == 9:
        is_water = True

    return uniform, is_coast, is_land, is_water, n_land, n_water, n_coast


@njit(cache=True)
def spatial_var_test(
    indat_11um: np.ndarray,
    threshold: float,
) -> bool:
    """Spatial variability test for ocean pixels.

    Computes the standard deviation of 11um BT in a 3x3 neighborhood.
    If std > threshold, the pixel is spatially variable (likely cloud).

    Args:
        indat_11um: 3x3 array of 11um brightness temperatures.
        threshold: Variability threshold (e.g., 0.6 K).

    Returns:
        True if spatially variable (std > threshold).
    """
    std = get_regional_std(indat_11um)
    if std < BAD_DATA + 1.0:
        return False
    return std > threshold
