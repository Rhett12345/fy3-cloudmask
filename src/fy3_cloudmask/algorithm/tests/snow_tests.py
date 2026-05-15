"""Snow/ice surface cloud detection tests.

Port of Day_snow.f90, Nite_snow.f90, Antarctic_day.f90.

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
    BIT_REF_064, BIT_NIR_138, BIT_THIN_CIRRUS_SOLAR,
    BAND_064, BAND_086, BAND_138, IR_38, IR_11, IR_12,
    VIS_VALID_MIN, VIS_VALID_MAX,
)


def day_snow(
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
    """Daytime snow/ice surface cloud detection.

    Port of Day_snow.f90. Tests 4 groups:
    1. PFMFT/NFMFT (IR forward model fit)
    2. BTD tests (11-12, 11-4um)
    3. Visible tests (0.64um)
    4. Thin cirrus (1.38um)

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        visusd: Whether visible data is usable.
        hi_elev: Whether pixel is at high elevation.
        bt_clr: 7-element clear-sky BT array.
        is_cold_sfc: 1 if cold surface, 0 otherwise.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('day_snow', {})
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

        locut = dfthrsh + 0.3 * dfthrsh
        hicut = dfthrsh - 1.25
        c5 = conf_test_thresholds(masdf1, np.array([locut, dfthrsh, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c5)
        ngtests[1] += 1

    # 11-4um BTD test
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
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
            ref064 = thr.get('ref064', [0.40, 0.35, 0.30, 1.0])
            if masv66 <= ref064[1]:
                set_bit(testbits, BIT_REF_064)
                set_bit(qa_bits, 20)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[2] += 1
            c4 = conf_test_thresholds(masv66, np.array(ref064, dtype=np.float64))
            cmin3 = min(cmin3, c4)

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


def nite_snow(
    pxldat: np.ndarray,
    vza: float,
    sfctmp: float,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Nighttime snow/ice surface cloud detection.

    Port of Nite_snow.f90. Tests 2 groups:
    1. PFMFT/NFMFT, surface temperature
    2. BTD tests (11-12, 11-4, 4-12um)

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        sfctmp: Surface temperature from NWP (K).
        bt_clr: 7-element clear-sky BT array.
        is_cold_sfc: 1 if cold surface, 0 otherwise.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('nite_snow', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1
    cmin2 = 1.0  # Group 2

    ngtests = [0, 0]

    masir4 = pxldat[IR_38]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    max_vza = 65.49

    # ================================================================
    # GROUP 1: PFMFT/NFMFT + SST
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

    # Surface temperature test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        if 0.0 < sfctmp < 350.0:
            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 27)

            lst_thrsh = 12.0
            masdf1 = masir11 - masir12
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

        locut = dfthrsh + 0.3 * dfthrsh
        hicut = dfthrsh - 1.25
        c5 = conf_test_thresholds(masdf1, np.array([locut, dfthrsh, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c5)
        ngtests[1] += 1

    # 11-4um BTD test
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
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

    # 4-12um BTD test
    if masir12 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        mas4_12 = masir4 - masir12
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 17)
        btd_4_12 = thr.get('btd_4_12', [4.5, 4.0, 3.5, 1.0])
        if mas4_12 <= btd_4_12[1]:
            set_bit(testbits, 17)
        c4 = conf_test_thresholds(mas4_12, np.array(btd_4_12, dtype=np.float64))
        cmin2 = min(cmin2, c4)
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


def antarctic_day(
    pxldat: np.ndarray,
    vza: float,
    visusd: bool,
    bt_clr: np.ndarray,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Antarctic daytime cloud detection.

    Port of Antarctic_day.f90. Simplified tests for Antarctic ice sheet.

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        visusd: Whether visible data is usable.
        bt_clr: 7-element clear-sky BT array.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('antarctic_day', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0
    cmin1 = 1.0
    cmin2 = 1.0
    ngtests = [0, 0]

    masv66 = pxldat[BAND_064]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    # Group 1: 11-12um BTD
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        masdf1 = masir11 - masir12
        btd_11_12 = thr.get('btd_11_12', [3.0, 2.5, 2.0, 1.0])
        if masdf1 <= btd_11_12[1]:
            set_bit(testbits, BIT_BTD_11_12)
            set_bit(qa_bits, 18)
            nmtests += 1
            nbands = max(nbands, 2)
            ngtests[0] += 1
        c1 = conf_test_thresholds(masdf1, np.array(btd_11_12, dtype=np.float64))
        cmin1 = min(cmin1, c1)

    # Group 2: Visible test
    if visusd:
        if masv66 > VIS_VALID_MIN and masv66 <= VIS_VALID_MAX:
            ref064 = thr.get('ref064', [0.50, 0.45, 0.40, 1.0])
            if masv66 <= ref064[1]:
                set_bit(testbits, BIT_REF_064)
                set_bit(qa_bits, 20)
                nmtests += 1
                nbands = max(nbands, 1)
                ngtests[1] += 1
            c2 = conf_test_thresholds(masv66, np.array(ref064, dtype=np.float64))
            cmin2 = min(cmin2, c2)

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
