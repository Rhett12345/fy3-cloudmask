"""Cloud mask restoral and post-processing tests.

Port of chk_land.f90, chk_land_nite.f90, chk_coast.f90, chk_sunglint.f90,
chk_spatial_var.f90, chk_shallow_water.f, shadows.f90, noncld_obs_chk.f90,
thin_ci_chk_ir.f90.

These tests can override the initial cloud mask decision based on
additional criteria (spatial context, adjacency, shadows, etc.).
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from ..confidence import conf_test_thresholds, encode_confidence
from ..bitops import set_bit, clear_bit, check_bit
from ..spatial import get_regional_mean, get_regional_std, get_regional_diff
from ...constants import (
    BAD_DATA, BIT_LAND_RESTORAL, BIT_SHADOW, BIT_NCO, BIT_THIN_CIRRUS_IR,
    BIT_CLOUD_ADJ, BIT_TEMPORAL,
    IR_11, IR_12, IR_38, IR_85,
    BAND_064, BAND_086, BAND_138,
    VIS_VALID_MIN, VIS_VALID_MAX,
)


def chk_land_restoral(
    confdnc: float,
    pxldat: np.ndarray,
    sfctmp: float,
    nwp_pmsl: float,
    nwp_precip_water: float,
    vza: float,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Land cloud restoral check.

    Port of chk_land.f90. Can restore clear-sky confidence for
    pixels initially flagged as cloudy.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        sfctmp: Surface temperature (K).
        nwp_pmsl: NWP mean sea level pressure (hPa).
        nwp_precip_water: NWP precipitable water (mm).
        vza: Satellite viewing angle (degrees).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('land_restoral', {})

    # Only apply if confidence is in cloudy range
    if confdnc > 0.95:
        return confdnc

    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]
    masir4 = pxldat[IR_38]

    if masir11 < BAD_DATA + 1.0:
        return confdnc

    # Temperature difference test
    if sfctmp > 0.0 and sfctmp < 350.0:
        sfcdif = sfctmp - masir11
        restoral_thr = thr.get('temp_diff', 5.0)

        if sfcdif < restoral_thr:
            # Restore to clear
            new_conf = thr.get('restored_confidence', 0.97)
            confdnc = max(confdnc, new_conf)
            set_bit(testbits, BIT_LAND_RESTORAL)
            set_bit(qa_bits, 26)

    # Precipitable water test
    if nwp_precip_water > 0:
        pw_thr = thr.get('precip_water', 6.0)
        if nwp_precip_water > pw_thr:
            # High moisture - be more conservative
            confdnc = min(confdnc, 0.95)

    return confdnc


def chk_land_nite_restoral(
    confdnc: float,
    pxldat: np.ndarray,
    sfctmp: float,
    nwp_pmsl: float,
    nwp_precip_water: float,
    vza: float,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Nighttime land cloud restoral check.

    Port of chk_land_nite.f90.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        sfctmp: Surface temperature (K).
        nwp_pmsl: NWP mean sea level pressure (hPa).
        nwp_precip_water: NWP precipitable water (mm).
        vza: Satellite viewing angle (degrees).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('land_nite_restoral', {})

    if confdnc > 0.95:
        return confdnc

    masir11 = pxldat[IR_11]
    if masir11 < BAD_DATA + 1.0:
        return confdnc

    if sfctmp > 0.0 and sfctmp < 350.0:
        sfcdif = sfctmp - masir11
        restoral_thr = thr.get('temp_diff', 5.0)

        if sfcdif < restoral_thr:
            new_conf = thr.get('restored_confidence', 0.97)
            confdnc = max(confdnc, new_conf)
            set_bit(testbits, BIT_LAND_RESTORAL)
            set_bit(qa_bits, 26)

    return confdnc


def chk_coast_restoral(
    confdnc: float,
    pxldat: np.ndarray,
    sfctmp: float,
    sst: float,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Coastal cloud restoral check.

    Port of chk_coast.f90.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        sfctmp: Surface temperature (K).
        sst: Sea surface temperature (K).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('coast_restoral', {})

    if confdnc > 0.95:
        return confdnc

    masir11 = pxldat[IR_11]
    if masir11 < BAD_DATA + 1.0:
        return confdnc

    # Use SST if available, otherwise use NWP surface temp
    temp = sst if sst > 100.0 else sfctmp
    if temp > 0.0 and temp < 350.0:
        sfcdif = temp - masir11
        restoral_thr = thr.get('temp_diff', 5.0)

        if sfcdif < restoral_thr:
            new_conf = thr.get('restored_confidence', 0.97)
            confdnc = max(confdnc, new_conf)

    return confdnc


def chk_sunglint_restoral(
    confdnc: float,
    pxldat: np.ndarray,
    refang: float,
    snglnt: bool,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Sun glint restoral check.

    Port of chk_sunglint.f90.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        refang: Reflectance angle.
        snglnt: Whether pixel is in sun glint region.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    if not snglnt:
        return confdnc

    thr = thresholds.get('sunglint_restoral', {})

    # In sun glint region, reduce confidence
    max_conf = thr.get('max_confidence', 0.95)
    confdnc = min(confdnc, max_conf)

    return confdnc


def chk_shallow_water(
    confdnc: float,
    pxldat: np.ndarray,
    sh_ocean: bool,
    sh_lake: bool,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Shallow water restoral check.

    Port of chk_shallow_water.f.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        sh_ocean: Whether shallow ocean.
        sh_lake: Whether shallow lake.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    if not sh_ocean and not sh_lake:
        return confdnc

    thr = thresholds.get('shallow_water', {})

    # Shallow water reduces confidence due to bottom reflectance
    max_conf = thr.get('max_confidence', 0.95)
    confdnc = min(confdnc, max_conf)

    return confdnc


def chk_spatial_var(
    confdnc: float,
    indat_3x3: np.ndarray,
    uniform: bool,
    is_edge: bool,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Spatial variability check.

    Port of chk_spatial_var.f90.

    Args:
        confdnc: Current confidence value.
        indat_3x3: 3x3 array of 11um BT.
        uniform: Whether pixel neighborhood is uniform.
        is_edge: Whether pixel is on image edge.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    if is_edge:
        return confdnc

    thr = thresholds.get('spatial_var', {})

    std = get_regional_std(indat_3x3)
    if std < BAD_DATA + 1.0:
        return confdnc

    # High spatial variability suggests cloud
    var_threshold = thr.get('std_threshold', 0.6)
    if std > var_threshold and confdnc > 0.95:
        # Reduce confidence if spatially variable
        confdnc = min(confdnc, 0.95)

    return confdnc


def chk_cloud_adj(
    confdnc: float,
    cm_array: np.ndarray,
    row: int,
    col: int,
    n_rows: int,
    n_cols: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Cloud adjacency check.

    Port of noncld_obs_chk.f90. If neighboring pixels are cloudy,
    reduce confidence of current pixel.

    Args:
        confdnc: Current confidence value.
        cm_array: 2D cloud mask array (already processed).
        row: Current row index.
        col: Current column index.
        n_rows: Total number of rows.
        n_cols: Total number of columns.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('cloud_adjacency', {})

    # Count cloudy neighbors in 3x3 box
    n_cloudy = 0
    n_total = 0
    for i in range(max(0, row - 1), min(n_rows, row + 2)):
        for j in range(max(0, col - 1), min(n_cols, col + 2)):
            if i == row and j == col:
                continue
            n_total += 1
            if cm_array[i, j] < 2:  # cloudy or probably cloudy
                n_cloudy += 1

    if n_total == 0:
        return confdnc

    # If most neighbors are cloudy, reduce confidence
    cloudy_frac = n_cloudy / n_total
    adj_threshold = thr.get('cloudy_fraction', 0.75)

    if cloudy_frac > adj_threshold:
        max_conf = thr.get('max_confidence', 0.95)
        confdnc = min(confdnc, max_conf)
        clear_bit(testbits, BIT_CLOUD_ADJ)

    return confdnc


def chk_thin_cirrus_ir(
    confdnc: float,
    pxldat: np.ndarray,
    vza: float,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Thin cirrus IR check.

    Port of thin_ci_chk_ir.f90.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('thin_cirrus_ir', {})

    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    if masir11 < BAD_DATA + 1.0 or masir12 < BAD_DATA + 1.0:
        return confdnc

    # 11-12um BTD for thin cirrus
    btd = masir11 - masir12
    cirrus_threshold = thr.get('btd_threshold', 1.0)

    if btd < cirrus_threshold and confdnc > 0.95:
        # Thin cirrus detected, reduce confidence
        clear_bit(testbits, BIT_THIN_CIRRUS_IR)
        confdnc = min(confdnc, 0.95)

    return confdnc


def chk_shadow(
    confdnc: float,
    pxldat: np.ndarray,
    indat_3x3_vis: np.ndarray,
    indat_3x3_11um: np.ndarray,
    sza: float,
    vza: float,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> float:
    """Cloud shadow detection.

    Port of shadows.f90.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        indat_3x3_vis: 3x3 array of visible reflectance.
        indat_3x3_11um: 3x3 array of 11um BT.
        sza: Solar zenith angle (degrees).
        vza: Satellite viewing angle (degrees).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value.
    """
    thr = thresholds.get('shadow', {})

    masv66 = pxldat[BAND_064]
    masir11 = pxldat[IR_11]

    if masv66 < VIS_VALID_MIN or masir11 < BAD_DATA + 1.0:
        return confdnc

    # Check for anomalously low visible reflectance
    vis_mean = get_regional_mean(indat_3x3_vis)
    if vis_mean < BAD_DATA + 1.0:
        return confdnc

    vis_diff = vis_mean - masv66
    shadow_vis_threshold = thr.get('vis_diff_threshold', 0.1)

    if vis_diff > shadow_vis_threshold:
        # Possible shadow - check 11um uniformity
        ir_std = get_regional_std(indat_3x3_11um)
        if ir_std < 1.0:  # Uniform IR field
            clear_bit(testbits, BIT_SHADOW)

    return confdnc
