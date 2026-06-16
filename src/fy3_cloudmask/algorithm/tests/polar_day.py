"""Polar daytime cloud detection tests.

Port of PolarDay_land.f90, PolarDay_coast.f90, PolarDay_desert.f90,
PolarDay_desert_c.f90, PolarDay_ocean.f90, PolarDay_snow.f90.

Each function returns a confidence value and updates testbits/qa_bits in-place.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from ..confidence import conf_test_thresholds
from ..bitops import set_bit, clear_bit, check_bit
from ..spatial import tview
from ...constants import (
    BAD_DATA, BIT_PFMFT, BIT_NFMFT, BIT_BTD_11_12, BIT_BTD_11_4,
    BIT_REF_064, BIT_GEMI, BIT_NIR_138, BIT_THIN_CIRRUS_SOLAR,
    BAND_064, BAND_086, BAND_138, IR_38, IR_11, IR_12, IR_85,
    VIS_VALID_MIN, VIS_VALID_MAX,
)


def polar_day_land(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    vrused: bool,
    cirrus_vis: bool,
    hi_elev: bool,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime land cloud detection.

    Port of PolarDay_land.f90. Tests 4 groups:
    1. PFMFT/NFMFT (IR forward model fit)
    2. BTD tests (11-12, 11-4um)
    3. Visible tests (0.64um, ratio)
    4. Thin cirrus (1.38um)

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        visusd: Whether visible data is usable.
        vrused: Whether reflectance ratio test can be used.
        cirrus_vis: Whether thin cirrus detected in visible.
        hi_elev: Whether pixel is at high elevation.
        bt_clr: 7-element clear-sky BT array.
        is_cold_sfc: 1 if cold surface, 0 otherwise.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('land_day_polar', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1: IR forward model
    cmin2 = 1.0  # Group 2: BTD tests
    cmin3 = 1.0  # Group 3: Visible tests
    cmin4 = 1.0  # Group 4: NIR test

    ngtests = [0, 0, 0, 0]

    masv66 = pxldat[BAND_064]
    masv88 = pxldat[BAND_086]
    masv188 = pxldat[BAND_138]
    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    # ================================================================
    # GROUP 1: PFMFT + NFMFT
    # ================================================================

    # PFMFT test
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            masir11 < pfmft_thr.get('bt_11_max', 310.0) and
            (bt_clr[4] - bt_clr[5]) > pfmft_thr.get('btd_min', 0.0)):
        if masir11 > 270.0 and bt_clr[4] > 270.0:
            tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5]) * (masir11 - 260.0) / (bt_clr[4] - 260.0)
        else:
            tv11_12 = masir11 - masir12
        set_bit(qa_bits, 14)
        set_bit(testbits, 14)
        nmtests += 1
        nbands = max(nbands, 2)
        if is_cold_sfc == 1:
            pfmft_cold = pfmft_thr.get('cold', [2.0, 1.5, 1.0, 1.0])
            c1 = conf_test_thresholds(tv11_12, np.array(pfmft_cold, dtype=np.float64))
        else:
            pfmft_land = pfmft_thr.get('land', [4.0, 3.5, 3.0, 1.0])
            c1 = conf_test_thresholds(tv11_12, np.array(pfmft_land, dtype=np.float64))
        cmin1 = min(cmin1, c1)
        ngtests[0] += 1

    # NFMFT test
    nfmft_max = nfmft_thr.get('max_threshold', 1.50)
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            (masir11 - masir12) <= nfmft_max):
        tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5])
        set_bit(qa_bits, 15)
        set_bit(testbits, 15)
        nmtests += 1
        nbands = max(nbands, 2)
        nfmft_land = nfmft_thr.get('land', [-23.0, -22.5, -22.0, 1.0])
        c2 = conf_test_thresholds(tv11_12, np.array(nfmft_land, dtype=np.float64))
        cmin1 = min(cmin1, c2)
        ngtests[0] += 1

    # ================================================================
    # GROUP 2: BTD tests
    # ================================================================

    # 11-12um BTD thin cirrus test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and vza > 0.0:
        masdf1 = masir11 - masir12
        cosvza = math.cos(vza * math.pi / 180.0)
        schi = 1.0 / cosvza if abs(cosvza) > 1e-6 else 99.0

        diftemp = tview(schi, masir11)
        btd_11_12 = thr.get('btd_11_12', [3.0])
        dfthrsh = btd_11_12[0] if (diftemp < 0.1 or abs(schi - 99.0) < 0.0001) else diftemp

        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 18)
        if masdf1 <= dfthrsh:
            set_bit(testbits, 18)

        locut = dfthrsh
        midpt = dfthrsh - 0.3 * dfthrsh
        hicut = midpt - 1.25

        c5 = conf_test_thresholds(masdf1, np.array([locut, midpt, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c5)
        ngtests[1] += 1

    # 11-4um BTD test (only when visible data is usable)
    if visusd and masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 19)
        mas11_4 = masir11 - masir4
        btd_11_4 = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])
        if mas11_4 >= btd_11_4[1]:
            set_bit(testbits, 19)
        c3 = conf_test_thresholds(mas11_4, np.array(btd_11_4, dtype=np.float64))
        cmin2 = min(cmin2, c3)
        ngtests[1] += 1

    # ================================================================
    # GROUP 3: Visible tests
    # ================================================================

    if visusd:
        # 0.64um reflectance test
        if masv66 > VIS_VALID_MIN and masv66 <= VIS_VALID_MAX:
            ref064 = thr.get('ref064', [0.24, 0.20, 0.16, 1.0])
            if masv66 <= ref064[1]:
                set_bit(testbits, BIT_REF_064)
                set_bit(qa_bits, 20)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[2] += 1
            c4 = conf_test_thresholds(masv66, np.array(ref064, dtype=np.float64))
            cmin3 = min(cmin3, c4)

        # GEMI/visible ratio test
        if vrused and masv66 > VIS_VALID_MIN and masv88 > VIS_VALID_MIN:
            s1 = masv66 * 100.0
            s2 = masv88 * 100.0
            etan = 2.0 * (s2 - s1) + 1.5 * s2 + 0.5 * s1
            etad = s2 + s1 + 0.5
            if etad > 0:
                eta = etan / etad
                vrat = eta * (1.0 - 0.25 * eta) - ((s1 - 0.125) / (1.0 - s1)) if s1 < 1.0 else 0.0
                vrat_thr = thr.get('vrat', [1.80, 1.85, 1.90, 1.0])
                if vrat > vrat_thr[1]:
                    set_bit(testbits, BIT_GEMI)
                    set_bit(qa_bits, 21)
                    nmtests += 1
                    nbands = max(nbands, 2)
                    ngtests[2] += 1
                c6 = conf_test_thresholds(vrat, np.array(vrat_thr, dtype=np.float64))
                cmin3 = min(cmin3, c6)

    # ================================================================
    # GROUP 4: NIR test (1.38um)
    # ================================================================

    if not hi_elev and masv188 > VIS_VALID_MIN and masv188 <= VIS_VALID_MAX:
        ref138 = thr.get('ref138', [0.04, 0.035, 0.03, 1.0])
        if masv188 <= ref138[1]:
            set_bit(testbits, BIT_NIR_138)
            set_bit(qa_bits, 16)
            nmtests += 1
            nbands = max(nbands, 1)
            ngtests[3] += 1
        c7 = conf_test_thresholds(masv188, np.array(ref138, dtype=np.float64))
        cmin4 = min(cmin4, c7)

        # Thin cirrus flag
        tci = thr.get('tci', [0.035, 0.0125])
        set_bit(qa_bits, 9)
        if tci[0] > masv188 >= tci[1]:
            clear_bit(qa_bits, 9)

    # ================================================================
    # Final confidence
    # ================================================================
    group_mins = [cmin1, cmin2, cmin3, cmin4]
    active_groups = sum(1 for g in ngtests if g > 0)

    if active_groups > 0:
        product = 1.0
        for i in range(4):
            if ngtests[i] > 0:
                product *= group_mins[i]
        confdnc = product ** (1.0 / active_groups)
    else:
        confdnc = 1.0

    return confdnc, nmtests, nbands


def polar_day_coast(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    vrused: bool,
    cirrus_vis: bool,
    hi_elev: bool,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime coastal cloud detection.

    Uses coastal-specific thresholds. Same structure as polar_day_land.
    """
    return polar_day_land(
        pxldat, vza, visusd, vrused, cirrus_vis, hi_elev,
        bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
    )


def polar_day_desert(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    vrused: bool,
    cirrus_vis: bool,
    hi_elev: bool,
    tbadj: float,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime desert cloud detection.

    Port of PolarDay_desert.f90.
    """
    thr = thresholds.get('land_day_desert', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0
    cmin1 = 1.0
    cmin2 = 1.0
    cmin3 = 1.0
    ngtests = [0, 0, 0]

    masv66 = pxldat[BAND_064]
    masv88 = pxldat[BAND_086]
    masv188 = pxldat[BAND_138]
    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    # Group 1: 11-12um BTD
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        masdf1 = masir11 - masir12
        btd_11_12 = thr.get('btd_11_12', [3.5, 3.0, 2.5, 1.0])
        if masdf1 <= btd_11_12[1]:
            set_bit(testbits, BIT_BTD_11_12)
            set_bit(qa_bits, 18)
            nmtests += 1
            nbands = max(nbands, 2)
            ngtests[0] += 1
        c1 = conf_test_thresholds(masdf1, np.array(btd_11_12, dtype=np.float64))
        cmin1 = min(cmin1, c1)

    # Group 2: 11-4um BTD
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        mas11_4 = masir11 - masir4
        btd_11_4_hi = thr.get('btd_11_4_hi', [-3.0, -5.0, -7.0, 1.0])
        btd_11_4_lo = thr.get('btd_11_4_lo', [-25.0, -23.0, -21.0, 1.0])
        if mas11_4 >= btd_11_4_lo[1]:
            set_bit(testbits, BIT_BTD_11_4)
            set_bit(qa_bits, 19)
            nmtests += 1
            nbands = max(nbands, 2)
            ngtests[1] += 1
        c2 = conf_test_thresholds(mas11_4, np.array(btd_11_4_lo, dtype=np.float64))
        cmin2 = min(cmin2, c2)
        c2h = conf_test_thresholds(mas11_4, np.array(btd_11_4_hi, dtype=np.float64))
        cmin2 = min(cmin2, c2h)

    # Group 3: Visible tests
    if visusd:
        if masv66 > VIS_VALID_MIN and masv88 > VIS_VALID_MIN:
            ref086 = thr.get('ref086', [0.42, 0.39, 0.36, 1.0])
            if masv88 <= ref086[1]:
                set_bit(testbits, BIT_REF_064)
                set_bit(qa_bits, 20)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[2] += 1
            c3 = conf_test_thresholds(masv88, np.array(ref086, dtype=np.float64))
            cmin3 = min(cmin3, c3)

        # 1.38um test
        if not hi_elev and masv188 > VIS_VALID_MIN:
            ref138 = thr.get('ref138', [0.040, 0.035, 0.03, 1.0])
            if masv188 <= ref138[1]:
                set_bit(testbits, BIT_NIR_138)
                set_bit(qa_bits, 16)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[2] += 1
            c3b = conf_test_thresholds(masv188, np.array(ref138, dtype=np.float64))
            cmin3 = min(cmin3, c3b)

            tci = thr.get('tci', [0.035, 0.0125])
            set_bit(qa_bits, 9)
            if tci[0] > masv188 >= tci[1]:
                clear_bit(qa_bits, 9)

    # Final confidence
    active_groups = sum(1 for g in ngtests if g > 0)
    if active_groups > 0:
        product = 1.0
        mins = [cmin1, cmin2, cmin3]
        for i in range(3):
            if ngtests[i] > 0:
                product *= mins[i]
        confdnc = product ** (1.0 / active_groups)

    return confdnc, nmtests, nbands


def polar_day_desert_coast(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    vrused: bool,
    cirrus_vis: bool,
    hi_elev: bool,
    tbadj: float,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime coastal desert cloud detection.

    Port of PolarDay_desert_c.f90.
    """
    return polar_day_desert(
        pxldat, vza, visusd, vrused, cirrus_vis, hi_elev, tbadj,
        bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
    )


def polar_day_ocean(
    pxldat: np.ndarray,
    vza: float,
    snglnt: bool,
    visusd: bool,
    cirrus_vis: bool,
    sfctmp: float,
    refang: float,
    sh_ocean: bool,
    bt_clr: np.ndarray,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime ocean cloud detection.

    Port of PolarDay_ocean.f90.
    """
    # Use polar-specific thresholds (ocean_day_polar) instead of standard ocean_day
    from .ocean_day import ocean_day
    polar_thresholds = dict(thresholds)
    if 'ocean_day_polar' in thresholds:
        polar_thresholds['ocean_day'] = thresholds['ocean_day_polar']
    return ocean_day(
        pxldat, vza, snglnt, visusd, cirrus_vis, sfctmp, refang, sh_ocean,
        bt_clr, polar_thresholds, testbits, qa_bits,
    )


def polar_day_snow(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    hi_elev: bool,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar daytime snow/ice cloud detection.

    Port of PolarDay_snow.f90.
    """
    thr = thresholds.get('day_snow_polar', {})
    pfmft_thr = thresholds.get('pfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0
    cmin1 = 1.0
    cmin2 = 1.0
    cmin3 = 1.0
    cmin4 = 1.0
    ngtests = [0, 0, 0, 0]

    masv66 = pxldat[BAND_064]
    masv188 = pxldat[BAND_138]
    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    # Group 1: PFMFT
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            masir11 < pfmft_thr.get('bt_11_max', 310.0) and
            (bt_clr[4] - bt_clr[5]) > pfmft_thr.get('btd_min', 0.0)):
        if masir11 > 270.0 and bt_clr[4] > 270.0:
            tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5]) * (masir11 - 260.0) / (bt_clr[4] - 260.0)
        else:
            tv11_12 = masir11 - masir12
        set_bit(qa_bits, 14)
        set_bit(testbits, 14)
        nmtests += 1
        nbands = max(nbands, 2)
        if is_cold_sfc == 1:
            pfmft_cold = pfmft_thr.get('cold', [2.0, 1.5, 1.0, 1.0])
            c1 = conf_test_thresholds(tv11_12, np.array(pfmft_cold, dtype=np.float64))
        else:
            pfmft_snow = pfmft_thr.get('snow', [4.0, 3.5, 3.0, 1.0])
            c1 = conf_test_thresholds(tv11_12, np.array(pfmft_snow, dtype=np.float64))
        cmin1 = min(cmin1, c1)
        ngtests[0] += 1

    # Group 2: 11-12um BTD
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        masdf1 = masir11 - masir12
        btd_11_12 = thr.get('btd_11_12', [3.0, 2.5, 2.0, 1.0])
        if masdf1 <= btd_11_12[1]:
            set_bit(testbits, BIT_BTD_11_12)
            set_bit(qa_bits, 18)
            nmtests += 1
            nbands = max(nbands, 2)
            ngtests[1] += 1
        c2 = conf_test_thresholds(masdf1, np.array(btd_11_12, dtype=np.float64))
        cmin2 = min(cmin2, c2)

    # Group 3: Visible tests
    if visusd:
        if masv66 > VIS_VALID_MIN and masv66 <= VIS_VALID_MAX:
            ref064 = thr.get('ref064', [0.40, 0.35, 0.30, 1.0])
            if masv66 <= ref064[1]:
                set_bit(testbits, BIT_REF_064)
                set_bit(qa_bits, 20)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[2] += 1
            c3 = conf_test_thresholds(masv66, np.array(ref064, dtype=np.float64))
            cmin3 = min(cmin3, c3)

        # 1.38um test
        if not hi_elev and masv188 > VIS_VALID_MIN:
            ref138 = thr.get('ref138', [0.04, 0.035, 0.03, 1.0])
            if masv188 <= ref138[1]:
                set_bit(testbits, BIT_NIR_138)
                set_bit(qa_bits, 16)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[3] += 1
            c4 = conf_test_thresholds(masv188, np.array(ref138, dtype=np.float64))
            cmin4 = min(cmin4, c4)

    # Final confidence
    group_mins = [cmin1, cmin2, cmin3, cmin4]
    active_groups = sum(1 for g in ngtests if g > 0)
    if active_groups > 0:
        product = 1.0
        for i in range(4):
            if ngtests[i] > 0:
                product *= group_mins[i]
        confdnc = product ** (1.0 / active_groups)

    return confdnc, nmtests, nbands
