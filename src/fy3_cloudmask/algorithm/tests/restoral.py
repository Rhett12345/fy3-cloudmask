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
from ..bitops import set_bit, clear_bit, check_bit, set_bit_qa
from ..spatial import get_regional_mean, get_regional_std, get_regional_diff, tview
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
    """Sun glint clear-sky restoral check.

    Port of chk_sunglint.f90. The Fortran version performs clear-sky
    restoral tests in sun-glint conditions — it can RAISE confidence
    to 0.96 when clear-sky conditions are verified. It does NOT cap
    confidence.

    Full Fortran logic requires spatial variability test results and
    0.895/0.935um band data not yet available in the pipeline. This
    simplified version removes the incorrect cap.

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

    # Fortran chk_sunglint.f90 does NOT cap confidence.
    # It performs clear-sky restoral tests and can raise confidence to 0.96.
    # The previous min(conf, 0.95) was incorrect.
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

    Port of chk_spatial_var.f90. If spatially uniform, boosts confidence.
    Does NOT reduce confidence for variable scenes.

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

    var_threshold = thr.get('std_threshold', 0.6)

    # Set QA bit indicating test was applied
    set_bit_qa(qa_bits, 25)

    if std <= var_threshold:
        # Spatially uniform - boost confidence
        set_bit(testbits, 25)
        if confdnc > 0.66:
            confdnc = max(confdnc, 0.96)
        else:
            confdnc = max(confdnc, 0.67)

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
    """Non-cloud obstruction check (dust/smoke).

    Port of noncld_obs_chk.f90. The Fortran version checks for dust/smoke
    obstructions using IR BTD tests. It does NOT modify confidence — it
    only clears test bit 28 if the 11-12um BTD is below the dust threshold.

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
        Updated confidence value (unchanged per Fortran behavior).
    """
    # Fortran noncld_obs_chk.f90 does NOT modify confidence.
    # It only checks for dust/smoke using IR BTD and clears bit 28.
    # The cloud adjacency confidence cap was a Python-only addition
    # that created a feedback loop (91% cloudy → all pixels capped at 0.95).
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

    Port of thin_ci_chk_ir.f90. Only sets flag bits; does NOT modify confidence.

    Args:
        confdnc: Current confidence value.
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Updated confidence value (unchanged).
    """
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    if masir11 < BAD_DATA + 1.0 or masir12 < BAD_DATA + 1.0:
        return confdnc
    if vza <= 0.0:
        return confdnc

    masdf1 = masir11 - masir12
    cosvza = math.cos(vza * math.pi / 180.0)
    schi = 1.0 / cosvza if abs(cosvza) > 1e-6 else 99.0

    # APOLLO lookup table threshold
    diftemp = tview(schi, masir11)

    if diftemp < 0.1 or abs(schi - 99.0) < 0.0001:
        return confdnc

    ci1 = diftemp
    ci2 = diftemp + 0.3 * diftemp

    # Set QA bit indicating thin cirrus test applied
    set_bit_qa(qa_bits, 11)

    # Only flag if BTD falls in thin cirrus range
    if masdf1 > ci1 and masdf1 <= ci2:
        clear_bit(testbits, BIT_THIN_CIRRUS_IR)

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
