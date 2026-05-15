"""Polar nighttime cloud detection tests.

Port of PolarNite_land.f90, PolarNite_ocean.f90, PolarNite_snow.f90.

Each function returns a confidence value and updates testbits/qa_bits in-place.
"""

from __future__ import annotations

import math

import numpy as np
from numba import njit

from ..confidence import conf_test_thresholds
from ..bitops import set_bit, clear_bit, check_bit
from ..spatial import tview, get_regional_std
from ...constants import (
    BAD_DATA, BIT_PFMFT, BIT_NFMFT, BIT_BTD_11_12, BIT_BTD_11_4,
    IR_38, IR_73, IR_85, IR_11, IR_12,
)


def polar_nite_land(
    pxldat: np.ndarray,
    vza: float,
    sfctmp: float,
    hi_elev: bool,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar nighttime land cloud detection.

    Port of PolarNite_land.f90. Tests 2 groups:
    1. High thick cloud (PFMFT/NFMFT, surface temperature)
    2. BTD tests (11-12, 11-4um)

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        sfctmp: Surface temperature from NWP (K).
        hi_elev: Whether pixel is at high elevation.
        bt_clr: 7-element clear-sky BT array.
        is_cold_sfc: 1 if cold surface, 0 otherwise.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('polar_nite_land', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1: high thick cloud
    cmin2 = 1.0  # Group 2: BTD tests

    ngtests = [0, 0]

    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    max_vza = 65.49

    # ================================================================
    # GROUP 1: High thick cloud tests
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

    # Surface temperature test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and not hi_elev:
        if 0.0 < sfctmp < 350.0:
            masdf1 = masir11 - masir12
            masdf2 = masir11 - (masir4 - 1.5)

            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 27)

            delta_t = 0.0
            lst_thrsh = 12.0 + delta_t

            if masdf1 >= 1.0:
                midpt = lst_thrsh + 2.0 * int(masdf1)
            else:
                midpt = lst_thrsh

            a = vza / max_vza
            corr = (a ** 4) * 3.0
            midpt = midpt + corr
            locut = midpt + 2.0
            hicut = midpt - 2.0

            sfcdif = sfctmp - masir11
            if sfcdif < midpt:
                set_bit(testbits, 27)

            c9 = conf_test_thresholds(sfcdif, np.array([locut, midpt, hicut, 1.0], dtype=np.float64))
            cmin1 = min(cmin1, c9)
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

    # 11-4um BTD test
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 19)
        mas11_4 = masir11 - masir4
        btd_11_4 = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])
        if mas11_4 <= btd_11_4[1]:
            set_bit(testbits, 19)
        c3 = conf_test_thresholds(mas11_4, np.array(btd_11_4, dtype=np.float64))
        cmin2 = min(cmin2, c3)
        ngtests[1] += 1

    # ================================================================
    # Final confidence
    # ================================================================
    group_mins = [cmin1, cmin2]
    active_groups = sum(1 for g in ngtests if g > 0)

    if active_groups > 0:
        product = 1.0
        for i in range(2):
            if ngtests[i] > 0:
                product *= group_mins[i]
        confdnc = product ** (1.0 / active_groups)
    else:
        confdnc = 1.0

    return confdnc, nmtests, nbands


def polar_nite_ocean(
    indat_3x3: np.ndarray,
    pxldat: np.ndarray,
    vza: float,
    sfctmp: float,
    sh_ocean: bool,
    uniform: bool,
    bt_clr: np.ndarray,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar nighttime ocean cloud detection.

    Port of PolarNite_ocean.f90.
    """
    from .ocean_nite import ocean_nite
    return ocean_nite(
        indat_3x3, pxldat, vza, sfctmp, sh_ocean, uniform,
        bt_clr, thresholds, testbits, qa_bits,
    )


def polar_nite_snow(
    pxldat: np.ndarray,
    vza: float,
    sfctmp: float,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Polar nighttime snow/ice cloud detection.

    Port of PolarNite_snow.f90.
    """
    thr = thresholds.get('polar_nite_snow', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0
    cmin1 = 1.0
    cmin2 = 1.0
    ngtests = [0, 0]

    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    max_vza = 65.49

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
            pfmft_snow = pfmft_thr.get('snow', [4.0, 3.5, 3.0, 1.0])
            c1 = conf_test_thresholds(tv11_12, np.array(pfmft_snow, dtype=np.float64))
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
        nfmft_snow = nfmft_thr.get('snow', [-23.0, -22.5, -22.0, 1.0])
        c2 = conf_test_thresholds(tv11_12, np.array(nfmft_snow, dtype=np.float64))
        cmin1 = min(cmin1, c2)
        ngtests[0] += 1

    # 11-12um BTD
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        masdf1 = masir11 - masir12
        btd_11_12 = thr.get('btd_11_12', [3.0, 2.5, 2.0, 1.0])
        if masdf1 <= btd_11_12[1]:
            set_bit(testbits, BIT_BTD_11_12)
            set_bit(qa_bits, 18)
            nmtests += 1
            nbands = max(nbands, 2)
            ngtests[1] += 1
        c3 = conf_test_thresholds(masdf1, np.array(btd_11_12, dtype=np.float64))
        cmin2 = min(cmin2, c3)

    # Final confidence
    group_mins = [cmin1, cmin2]
    active_groups = sum(1 for g in ngtests if g > 0)
    if active_groups > 0:
        product = 1.0
        for i in range(2):
            if ngtests[i] > 0:
                product *= group_mins[i]
        confdnc = product ** (1.0 / active_groups)

    return confdnc, nmtests, nbands
