"""Nighttime ocean cloud detection tests.

Port of ocean_nite.f90.

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
    BAD_DATA, BIT_BTD_11_12, BIT_BTD_11_4, IR_38, IR_73, IR_85, IR_11, IR_12,
)


def ocean_nite(
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
    """Nighttime ocean cloud detection.

    Port of ocean_nite.f90. Tests 2 groups:
    1. High thick cloud (11um BT, PFMFT/NFMFT, SST)
    2. BTD tests (8-11, 11-12, 11-4, 8.6-7.3, variability)

    Args:
        indat_3x3: 3x3 array of 11um BT for spatial variability test.
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        sfctmp: Sea surface temperature (K).
        sh_ocean: Whether shallow ocean.
        uniform: Whether pixel neighborhood is uniform.
        bt_clr: 7-element clear-sky BT array.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('ocean_nite', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1: high thick cloud
    cmin2 = 1.0  # Group 2: BTD tests

    ngtests = [0, 0]

    masir4 = pxldat[IR_38]
    masir73 = pxldat[IR_73]
    masir8 = pxldat[IR_85]
    masir11 = pxldat[IR_11]
    masir12 = pxldat[IR_12]

    max_vza = 65.49

    # ================================================================
    # GROUP 1: High thick cloud tests
    # ================================================================

    # 11um BT threshold test
    if masir11 > BAD_DATA + 1.0:
        nobt11 = thr.get('bt11', [235.0, 270.0, 265.0, 260.0])
        nmtests += 1
        nbands = max(nbands, 1)
        set_bit(qa_bits, 13)
        if masir11 >= nobt11[1]:
            set_bit(testbits, 13)
        c1 = conf_test_thresholds(masir11, np.array(nobt11, dtype=np.float64))
        cmin1 = min(cmin1, c1)
        ngtests[0] += 1

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
        pfmft_ocean = pfmft_thr.get('ocean', [4.0, 3.5, 3.0, 1.0])
        c2 = conf_test_thresholds(tv11_12, np.array(pfmft_ocean, dtype=np.float64))
        # Note: In Fortran, cmin1 update is commented out
        # cmin1 = min(cmin1, c2)
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
        nfmft_ocean = nfmft_thr.get('ocean', [-23.0, -22.5, -22.0, 1.0])
        c3 = conf_test_thresholds(tv11_12, np.array(nfmft_ocean, dtype=np.float64))
        # Note: In Fortran, cmin1 update is commented out
        # cmin1 = min(cmin1, c3)
        # ngtests[0] += 1

    # SST test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        if 0.0 < sfctmp < 350.0:
            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 27)

            sst_thrsh = 10.0 if sh_ocean else 6.0
            masdf1 = masir11 - masir12
            if masdf1 >= 1.0:
                midpt = sst_thrsh + 2.0 * int(masdf1)
            else:
                midpt = sst_thrsh

            a = vza / max_vza
            corr = (a ** 4) * 3.0
            midpt = midpt + corr
            locut = midpt + 1.0
            hicut = midpt - 2.0

            sfcdif = sfctmp - masir11
            if sfcdif < midpt:
                set_bit(testbits, 27)

            c10 = conf_test_thresholds(sfcdif, np.array([locut, midpt, hicut, 1.0], dtype=np.float64))
            cmin1 = min(cmin1, c10)
            ngtests[0] += 1

    # ================================================================
    # GROUP 2: BTD tests
    # ================================================================

    # 8-11um tri-spectral test
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            masir8 > BAD_DATA + 1.0):
        masdf2 = masir8 - masir11
        masdf1 = masir11 - masir12
        tri_thres = _trispc(masdf1)
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 18)
        if masdf2 < tri_thres:
            set_bit(testbits, 18)
        locut = tri_thres + 0.5
        hicut = tri_thres - 0.5
        c4 = conf_test_thresholds(masdf2, np.array([locut, tri_thres, hicut, 1.0], dtype=np.float64))
        # Note: In Fortran, this test is commented out for cmin2
        # cmin2 = min(cmin2, c4)
        # ngtests[1] += 1

    # 11-12um BTD thin cirrus test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and vza > 0.0:
        masdf1 = masir11 - masir12
        cosvza = math.cos(vza * math.pi / 180.0)
        schi = 1.0 / cosvza if abs(cosvza) > 1e-6 else 99.0

        diftemp = tview(schi, masir11)
        no11_12hi = thr.get('btd_11_12', [3.0])
        dfthrsh = no11_12hi[0] if (diftemp < 0.1 or abs(schi - 99.0) < 0.0001) else diftemp

        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 18)
        if masdf1 <= dfthrsh:
            pass  # nptests
        else:
            if check_bit(testbits, 18):
                clear_bit(testbits, 18)

        locut = dfthrsh + 0.3 * dfthrsh
        hicut = dfthrsh - 1.25
        c7 = conf_test_thresholds(masdf1, np.array([locut, dfthrsh, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c7)
        ngtests[1] += 1

    # 11-4um BTD fog/low cloud test
    if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 19)
        mas11_4 = masir11 - masir4
        no11_4lo = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])
        if mas11_4 <= no11_4lo[1]:
            set_bit(testbits, 19)
        c6 = conf_test_thresholds(mas11_4, np.array(no11_4lo, dtype=np.float64))
        # Note: In Fortran, this test is commented out for cmin2
        # cmin2 = min(cmin2, c6)
        # ngtests[1] += 1

    # Water vapor cloud test (8.6-7.3um)
    if masir73 > BAD_DATA + 1.0 and masir8 > BAD_DATA + 1.0:
        dwvs = masir8 - masir73
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 29)
        no86_73 = thr.get('btd_86_73', [-5.0, -3.0, -1.0, 1.0])
        if dwvs > no86_73[1]:
            set_bit(testbits, 29)
        c9 = conf_test_thresholds(dwvs, np.array(no86_73, dtype=np.float64))
        cmin2 = min(cmin2, c9)
        ngtests[1] += 1

    # Spatial variability test
    if uniform:
        nmtests += 1
        nbands = max(nbands, 1)
        set_bit(qa_bits, 30)

        no_11var = thr.get('variability', [0.3, 0.6, 0.9, 1.0])
        std = get_regional_std(indat_3x3)
        np_val = std if std > BAD_DATA + 1.0 else 0.0

        if np_val > no_11var[1]:
            set_bit(testbits, 30)

        c11 = conf_test_thresholds(np_val, np.array(no_11var, dtype=np.float64))
        cmin2 = min(cmin2, c11)
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


def _trispc(btd_11_12: float) -> float:
    """Tri-spectral test threshold based on 11-12um BTD.

    Port of trispc() function. Returns 8-11um BTD threshold.

    Args:
        btd_11_12: 11-12um brightness temperature difference (K).

    Returns:
        Threshold for 8-11um BTD test.
    """
    return 0.558 * btd_11_12 - 3.36
