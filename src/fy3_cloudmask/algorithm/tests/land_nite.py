"""Nighttime land cloud detection tests.

Port of LandNite.f90.

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
    IR_38, IR_73, IR_85, IR_11, IR_12,
)


def land_nite(
    pxldat: np.ndarray,
    plat: float,
    vza: float,
    coast: bool,
    desert: bool,
    hi_elev: bool,
    sh_lake: bool,
    sfctmp: float,
    eco_type: int,
    ptwp: float,
    bt_clr: np.ndarray,
    is_cold_sfc: int,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> tuple[float, int, int]:
    """Nighttime land cloud detection.

    Port of LandNite.f90. Tests 3 groups:
    1. High thick cloud (PFMFT/NFMFT, surface temperature)
    2. Low cloud BTD (11-12, 11-4, 7.3-11um)
    3. High thin cloud (4-12um)

    Args:
        pxldat: 25-element pixel data array.
        plat: Latitude (degrees).
        vza: Satellite viewing angle (degrees).
        coast: Whether pixel is coastal.
        desert: Whether pixel is desert.
        hi_elev: Whether pixel is at high elevation.
        sh_lake: Whether pixel is shallow lake.
        sfctmp: Surface temperature from NWP (K).
        eco_type: Ecosystem type index.
        ptwp: Precipitable water (mm).
        bt_clr: 7-element clear-sky BT array.
        is_cold_sfc: 1 if cold surface, 0 otherwise.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('land_nite', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1: high thick cloud
    cmin2 = 1.0  # Group 2: BTD tests
    cmin5 = 1.0  # Group 5: high thin cloud

    ngtests = [0, 0, 0]

    masir4 = pxldat[IR_38]
    masir73 = pxldat[IR_73]
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
        # Note: In Fortran, cmin1 update is commented out
        # cmin1 = min(cmin1, c1)
        # ngtests[0] += 1

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
        # Note: In Fortran, cmin1 update is commented out
        # cmin1 = min(cmin1, c2)
        # ngtests[0] += 1

    # Surface temperature test
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            not hi_elev and eco_type != 8):
        if 0.0 < sfctmp < 350.0:
            masdf1 = masir11 - masir12
            masdf2 = masir11 - (masir4 - 1.5)  # corrected masir4

            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 27)

            delta_t = 0.0  # NWP-dependent adjustment
            if desert:
                lst_thrsh = 20.0 + delta_t
            elif masdf1 >= -0.2 or (masdf1 < -0.2 and (masdf2 <= -0.5 or masdf2 >= 1.0)):
                lst_thrsh = 12.0 + delta_t
            else:
                lst_thrsh = 20.0 + delta_t

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
        nl11_12hi = thr.get('btd_11_12', [3.0])
        dfthrsh = nl11_12hi[0] if (diftemp < 0.1 or abs(schi - 99.0) < 0.0001) else diftemp

        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 18)
        if masdf1 <= dfthrsh:
            set_bit(testbits, 18)

        locut = dfthrsh
        midpt = dfthrsh - 0.3 * dfthrsh
        if masir11 < 270.0:
            if abs(plat) <= 30.0:
                hicut = midpt - 1.25
            else:
                a = (90.0 - abs(plat)) / 60.0
                hicut = -0.1 - ((a ** 4) * 1.15)
        else:
            hicut = midpt - 1.25

        c5 = conf_test_thresholds(masdf1, np.array([locut, midpt, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c5)
        ngtests[1] += 1

    # 11-4um BTD fog/low cloud test
    if (masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0 and
            masir12 > BAD_DATA + 1.0):
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 19)
        mas11_4 = masir11 - masir4
        masdf1 = masir11 - masir12

        # Get thresholds based on 11-12um BTD
        nl_thresholds = _get_nl_thresholds(masdf1, thr)
        locut = nl_thresholds[0]
        midpt = nl_thresholds[1]
        hicut = nl_thresholds[2]
        power = nl_thresholds[3]

        if sh_lake:
            locut += 2.0
            midpt += 2.0
            hicut += 2.0

        if mas11_4 <= midpt:
            set_bit(testbits, 19)

        c3 = conf_test_thresholds(mas11_4, np.array([locut, midpt, hicut, power], dtype=np.float64))
        cmin2 = min(cmin2, c3)
        ngtests[1] += 1

    # 7.3-11um BTD thick mid-level cloud test
    if (masir11 > BAD_DATA + 1.0 and masir73 > BAD_DATA + 1.0 and
            masir4 > BAD_DATA + 1.0):
        mas11_4 = masir11 - masir4
        if mas11_4 <= -2.0:
            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 23)
            mas7_11 = masir73 - masir11
            nl7_11s = thr.get('btd_73_11', [-5.0, -3.0, -1.0, 1.0])
            if mas7_11 <= nl7_11s[1]:
                set_bit(testbits, 23)
            c6 = conf_test_thresholds(mas7_11, np.array(nl7_11s, dtype=np.float64))
            # Note: In Fortran, cmin2 update is commented out
            # cmin2 = min(cmin2, c6)
            # ngtests[1] += 1

    # ================================================================
    # GROUP 5: High thin cloud test (4-12um)
    # ================================================================

    if masir12 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        mas4_12 = masir4 - masir12
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 17)
        nl4_12hi = thr.get('btd_4_12', [4.5, 4.0, 3.5, 1.0])
        if mas4_12 <= nl4_12hi[1]:
            set_bit(testbits, 17)
        c4 = conf_test_thresholds(mas4_12, np.array(nl4_12hi, dtype=np.float64))
        cmin5 = min(cmin5, c4)
        ngtests[2] += 1

    # ================================================================
    # Final confidence
    # ================================================================
    group_mins = [cmin1, cmin2, cmin5]
    active_groups = sum(1 for g in ngtests if g > 0)

    if active_groups > 0:
        product = 1.0
        for i in range(3):
            if ngtests[i] > 0:
                product *= group_mins[i]
        confdnc = product ** (1.0 / active_groups)
    else:
        confdnc = 1.0

    return confdnc, nmtests, nbands


def _get_nl_thresholds(btd_11_12: float, thr: dict) -> list[float]:
    """Get nighttime land 11-4um thresholds based on 11-12um BTD.

    Port of get_nl_thresholds() function.

    Args:
        btd_11_12: 11-12um brightness temperature difference (K).
        thr: Threshold dictionary.

    Returns:
        List of [locut, midpt, hicut, power].
    """
    # Default thresholds
    nl11_4 = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])

    # Adjust based on 11-12um BTD
    if btd_11_12 >= 1.0:
        offset = 2.0 * int(btd_11_12)
        return [nl11_4[0] + offset, nl11_4[1] + offset, nl11_4[2] + offset, nl11_4[3]]
    else:
        return nl11_4
