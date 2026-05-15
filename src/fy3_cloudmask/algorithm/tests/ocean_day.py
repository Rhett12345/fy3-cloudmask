"""Daytime ocean cloud detection tests.

Port of ocean_day.f90.

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
    BAD_DATA, BIT_BTD_11_12, BIT_BTD_11_4, BIT_NIR_138,
    BAND_064, BAND_086, BAND_138, IR_38, IR_11, IR_12, IR_85, IR_73,
    VIS_VALID_MIN, VIS_VALID_MAX,
)


def ocean_day(
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
    """Daytime ocean cloud detection.

    Port of ocean_day.f90. Tests 4 groups:
    1. High thick cloud (11um BT, PFMFT/NFMFT, SST)
    2. Low cloud BTD (8-11, 11-12, 11-4um)
    3. Visible (0.86um refl, 0.86/0.66 ratio)
    4. Thin cirrus (1.38um)

    Args:
        pxldat: 25-element pixel data array.
        vza: Satellite viewing angle (degrees).
        snglnt: Whether pixel is in sun glint region.
        visusd: Whether visible data is usable.
        cirrus_vis: Whether thin cirrus detected in visible.
        sfctmp: Sea surface temperature (K).
        refang: Reflectance angle.
        sh_ocean: Whether shallow ocean (< 50m depth or near shore).
        bt_clr: 7-element clear-sky BT array.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    thr = thresholds.get('ocean_day', {})
    pfmft_thr = thresholds.get('pfmft', {})
    nfmft_thr = thresholds.get('nfmft', {})
    snglnt_thr = thresholds.get('sunglint', {})

    nmtests = 0
    nbands = 0
    confdnc = 1.0

    cmin1 = 1.0  # Group 1: high thick cloud
    cmin2 = 1.0  # Group 2: BTD tests
    cmin3 = 1.0  # Group 3: visible tests
    cmin4 = 1.0  # Group 4: thin cirrus

    ngtests = [0, 0, 0, 0]

    masv66 = pxldat[BAND_064]
    masv88 = pxldat[BAND_086]
    masv188 = pxldat[BAND_138]
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
        dobt11 = thr.get('bt11', [235.0, 270.0, 265.0, 260.0])
        nmtests += 1
        nbands = max(nbands, 1)
        set_bit(qa_bits, 13)
        if masir11 >= dobt11[1]:
            set_bit(testbits, 13)
        c1 = conf_test_thresholds(masir11, np.array(dobt11, dtype=np.float64))
        cmin1 = min(cmin1, c1)
        ngtests[0] += 1

    # PFMFT test
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            masir11 < pfmft_thr.get('bt_11_max', 310.0) and
            (masir11 - masir12) < pfmft_thr.get('btd_min', 0.0)):
        if masir11 > 270.0 and bt_clr[4] > 270.0:
            tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5]) * (masir11 - 260.0) / (bt_clr[4] - 260.0)
        else:
            tv11_12 = masir11 - masir12
        set_bit(qa_bits, 14)
        set_bit(testbits, 14)
        nmtests += 1
        nbands = max(nbands, 2)

    # NFMFT test
    nfmft_max = nfmft_thr.get('max_threshold', 1.50)
    if (masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and
            (masir11 - masir12) <= nfmft_max):
        tv11_12 = (masir11 - masir12) - (bt_clr[4] - bt_clr[5])
        set_bit(qa_bits, 15)
        set_bit(testbits, 15)
        nmtests += 1
        nbands = max(nbands, 2)

    # SST test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0:
        if 0.0 < sfctmp < 350.0:
            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 27)

            sst_thrsh = 10.0 if sh_ocean else 6.0
            r24_25 = masir11 - masir12
            if r24_25 >= 1.0:
                midpt = sst_thrsh + 2.0 * int(r24_25)
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
        r23_24 = masir8 - masir11
        r24_25 = masir11 - masir12
        tri_thres = _trispc(r24_25)
        nmtests += 1
        nbands = max(nbands, 2)
        set_bit(qa_bits, 18)
        if r23_24 < tri_thres:
            set_bit(testbits, 18)
        locut = tri_thres + 0.5
        hicut = tri_thres - 0.5
        c4 = conf_test_thresholds(r23_24, np.array([locut, tri_thres, hicut, 1.0], dtype=np.float64))
        # Note: In Fortran, this test is commented out for cmin2
        # cmin2 = min(cmin2, c4)
        # ngtests[1] += 1

    # 11-12um BTD thin cirrus test
    if masir11 > BAD_DATA + 1.0 and masir12 > BAD_DATA + 1.0 and vza > 0.0:
        masdf1 = masir11 - masir12
        cosvza = math.cos(vza * math.pi / 180.0)
        schi = 1.0 / cosvza if abs(cosvza) > 1e-6 else 99.0

        diftemp = tview(schi, masir11)
        do11_12hi = thr.get('btd_11_12', [3.0])
        dfthrsh = do11_12hi[0] if (diftemp < 0.1 or abs(schi - 99.0) < 0.0001) else diftemp

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
        c6 = conf_test_thresholds(masdf1, np.array([locut, dfthrsh, hicut, 1.0], dtype=np.float64))
        cmin2 = min(cmin2, c6)
        ngtests[1] += 1

    # 11-4um BTD fog/low cloud test (only when vis usable and not sunglint)
    if visusd and not snglnt:
        if masir11 > BAD_DATA + 1.0 and masir4 > BAD_DATA + 1.0:
            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 19)
            r24_21 = masir11 - masir4
            do11_4lo = thr.get('btd_11_4', [-14.0, -12.0, -10.0, 1.0])
            if r24_21 >= do11_4lo[1]:
                set_bit(testbits, 19)
            c7 = conf_test_thresholds(r24_21, np.array(do11_4lo, dtype=np.float64))
            cmin2 = min(cmin2, c7)
            ngtests[1] += 1

    # ================================================================
    # GROUP 3: Visible tests
    # ================================================================

    if visusd:
        # 0.86um NIR reflectance test
        if masv88 > VIS_VALID_MIN:
            if snglnt:
                sg_thr = snglnt_thr.get('ref086', [0.03, 0.09, 0.15, 1.0])
                locut = sg_thr[0]
                midpt = sg_thr[1]
                hicut = sg_thr[2]
                power = sg_thr[3]
            else:
                doref2 = thr.get('ref086', [0.03, 0.09, 0.15, 1.0])
                locut = doref2[0]
                midpt = doref2[1]
                hicut = doref2[2]
                power = doref2[3]

            nmtests += 1
            nbands = max(nbands, 1)
            set_bit(qa_bits, 20)
            if masv88 <= midpt:
                set_bit(testbits, 20)
            c8 = conf_test_thresholds(masv88, np.array([locut, midpt, hicut, power], dtype=np.float64))
            cmin3 = min(cmin3, c8)
            ngtests[2] += 1

        # Visible channel ratio test (0.86/0.66)
        if masv66 > VIS_VALID_MIN and masv88 > VIS_VALID_MIN:
            if snglnt:
                sg_thr = snglnt_thr.get('vrat', [0.7, 0.9, 1.1, 1.0])
                locuta = [sg_thr[0], sg_thr[0]]
                hicuta = [sg_thr[2], sg_thr[2]]
                midpta = [sg_thr[1], sg_thr[1]]
            else:
                dovratlo = thr.get('vrat_lo', [0.7, 0.9, 1.1, 1.0])
                dovrathi = thr.get('vrat_hi', [1.5, 1.8, 2.1, 1.0])
                locuta = [dovratlo[0], dovrathi[0]]
                hicuta = [dovratlo[2], dovrathi[2]]
                midpta = [dovratlo[1], dovrathi[1]]

            nmtests += 1
            nbands = max(nbands, 2)
            set_bit(qa_bits, 21)
            vrat = masv88 / masv66 if masv66 > 0 else 0.0
            if vrat < midpta[0] or vrat > midpta[1]:
                set_bit(testbits, 21)
            c9 = _conf_test_2val(vrat, locuta, hicuta, 1.0, midpta)
            cmin3 = min(cmin3, c9)
            ngtests[2] += 1

    # ================================================================
    # GROUP 4: Thin cirrus test (1.38um)
    # ================================================================

    if visusd:
        if masv188 > VIS_VALID_MIN:
            nmtests += 1
            nbands = max(nbands, 1)
            set_bit(qa_bits, 16)
            doref3 = thr.get('ref138', [0.04, 0.035, 0.03, 1.0])
            if masv188 <= doref3[1]:
                set_bit(testbits, 16)
            c11 = conf_test_thresholds(masv188, np.array(doref3, dtype=np.float64))
            cmin4 = min(cmin4, c11)
            ngtests[3] += 1

        # Thin cirrus flag
        dotci = thr.get('tci', [0.035, 0.0125])
        set_bit(qa_bits, 9)
        if dotci[0] > masv188 >= dotci[1]:
            clear_bit(qa_bits, 9)
            cirrus_vis = True

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


def _trispc(btd_11_12: float) -> float:
    """Tri-spectral test threshold based on 11-12um BTD.

    Port of trispc() function. Returns 8-11um BTD threshold.

    Args:
        btd_11_12: 11-12um brightness temperature difference (K).

    Returns:
        Threshold for 8-11um BTD test.
    """
    # Linear regression from global HIRS data
    # threshold = 0.558 * btd_11_12 - 3.36
    return 0.558 * btd_11_12 - 3.36


def _conf_test_2val(
    val: float,
    locut: list[float],
    hicut: list[float],
    power: float,
    midpt: list[float],
) -> float:
    """Two-value confidence test for ratio tests.

    Used when clear-sky values fall between two thresholds.

    Args:
        val: Test value.
        locut: Two-element list of low cutoffs.
        hicut: Two-element list of high cutoffs.
        power: S-curve power.
        midpt: Two-element list of midpoints.

    Returns:
        Confidence value in [0.0, 1.0].
    """
    # For ratio tests, we have two ranges: below midpt[0] and above midpt[1]
    # Compute confidence for each range and take the maximum
    c1 = conf_test_thresholds(val, np.array([locut[0], midpt[0], hicut[0], power], dtype=np.float64))
    c2 = conf_test_thresholds(val, np.array([locut[1], midpt[1], hicut[1], power], dtype=np.float64))
    return max(c1, c2)
