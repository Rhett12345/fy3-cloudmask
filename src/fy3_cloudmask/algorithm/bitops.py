"""Bit manipulation utilities for cloud mask products.

Exact port of Fortran set_bit.f, clear_bit.f, check_bits.f,
fill_bit_pixel.f90, set_quality_A.f, set_unused_bits.f, proc_path.f.

The cloud mask uses 6-byte (48-bit) testbits and 10-byte (80-bit) qa_bits arrays.
Bit numbering: bits 0-7 in byte 0, bits 8-15 in byte 1, etc.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ..constants import (
    BIT_PROCESSED, BIT_CONF_LSB, BIT_CONF_MSB,
    BIT_DAY, BIT_NO_SUNGLINT, BIT_NO_SNOW_ICE,
    BIT_COAST, BIT_DESERT, BIT_NCO, BIT_THIN_CIRRUS_SOLAR,
    BIT_SHADOW, BIT_THIN_CIRRUS_IR, BIT_CLOUD_ADJ,
    BIT_TEMPORAL, BIT_SUSPENDED_DUST, INITIAL_BITS,
)


@njit(cache=True)
def set_bit(testbits: np.ndarray, bit_num: int) -> None:
    """Set a bit in the testbits array.

    Args:
        testbits: 6-byte array (uint8).
        bit_num: Bit number (0-47).
    """
    iword = bit_num // 8
    ipos = bit_num % 8
    testbits[iword] = testbits[iword] | (1 << ipos)


@njit(cache=True)
def clear_bit(testbits: np.ndarray, bit_num: int) -> None:
    """Clear a bit in the testbits array.

    Args:
        testbits: 6-byte array (uint8).
        bit_num: Bit number (0-47).
    """
    iword = bit_num // 8
    ipos = bit_num % 8
    testbits[iword] = testbits[iword] & ~(1 << ipos)


@njit(cache=True)
def check_bit(testbits: np.ndarray, bit_num: int) -> bool:
    """Check if a bit is set in the testbits array.

    Args:
        testbits: 6-byte array (uint8).
        bit_num: Bit number (0-47).

    Returns:
        True if bit is set.
    """
    iword = bit_num // 8
    ipos = bit_num % 8
    return (testbits[iword] & (1 << ipos)) != 0


@njit(cache=True)
def init_testbits() -> np.ndarray:
    """Initialize testbits with default values (matching pxinit.f).

    Sets initial bits: NCO(8), thin_cirrus_solar(9), shadow(10),
    thin_cirrus_ir(11), cloud_adj(12), temporal(24), suspended_dust(28),
    spare(31) all to 1.

    Returns:
        Initialized 6-byte array.
    """
    testbits = np.zeros(6, dtype=np.uint8)
    for bit_num in INITIAL_BITS:
        set_bit(testbits, bit_num)
    return testbits


@njit(cache=True)
def init_qa_bits() -> np.ndarray:
    """Initialize QA bits array.

    Returns:
        Initialized 10-byte array.
    """
    return np.zeros(10, dtype=np.uint8)


@njit(cache=True)
def fill_bit_pixel(
    nmtests: int,
    nbands: int,
    bad_geo: bool,
    snglnt: bool,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
) -> None:
    """Final QA quality assembly for one pixel.

    Port of fill_bit_pixel.f90. Modifies testbits and qa_bits in-place.

    Quality encoding based on nmtests and nbands:
    - nmtests==0 or nbands==0 or bad_geo: all testbits zeroed, qa_bits[0]=0
    - nmtests < 3: quality=4 (bits 0,3 set)
    - nmtests < 7: quality=6 (bits 0,2,3 set)
    - else: quality=7 (bits 0,1,2,3 set)

    Sun glint reduction: if snglnt and qa_bits[0]==15, reduce to 13.

    Args:
        nmtests: Number of tests run.
        nbands: Number of bands used.
        bad_geo: Whether geolocation is bad.
        snglnt: Whether pixel is in sun glint region.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).
    """
    if nmtests == 0 or nbands == 0 or bad_geo:
        # No valid tests or bad geometry - zero everything
        for i in range(6):
            testbits[i] = 0
        qa_bits[0] = 0
        return

    # Set processed bit
    set_bit(testbits, BIT_PROCESSED)

    # Quality encoding
    if nmtests < 3:
        # Low quality
        qa_bits[0] = 0
        set_bit_qa(qa_bits, 0)  # bit 0
        set_bit_qa(qa_bits, 3)  # bit 3
    elif nmtests < 7:
        # Medium quality
        qa_bits[0] = 0
        set_bit_qa(qa_bits, 0)  # bit 0
        set_bit_qa(qa_bits, 2)  # bit 2
        set_bit_qa(qa_bits, 3)  # bit 3
    else:
        # High quality
        qa_bits[0] = 0
        set_bit_qa(qa_bits, 0)  # bit 0
        set_bit_qa(qa_bits, 1)  # bit 1
        set_bit_qa(qa_bits, 2)  # bit 2
        set_bit_qa(qa_bits, 3)  # bit 3

    # Sun glint quality reduction
    if snglnt and qa_bits[0] == 15:
        qa_bits[0] = 13


@njit(cache=True)
def set_bit_qa(qa_bits: np.ndarray, bit_num: int) -> None:
    """Set a bit in the QA bits array.

    Args:
        qa_bits: 10-byte array (uint8).
        bit_num: Bit number (0-79).
    """
    iword = bit_num // 8
    ipos = bit_num % 8
    qa_bits[iword] = qa_bits[iword] | (1 << ipos)


@njit(cache=True)
def proc_path(
    testbits: np.ndarray,
    day: bool,
    snglnt: bool,
    water: bool,
    snow: bool,
    ice: bool,
    coast: bool,
    desert: bool,
    land: bool,
    shadow: bool,
    smoke: bool,
) -> None:
    """Set processing path bits in testbits.

    Port of proc_path.f:
    - bit 5: no snow/ice
    - bit 3: day
    - bit 4: no sunglint (or not water)
    - bit 6: coast
    - bit 7: desert
    - bits 6+7: land (both set)
    - bit 10: cleared if shadow
    - bit 8: cleared if smoke

    Args:
        testbits: 6-byte array (modified in-place).
        day: Is daytime.
        snglnt: Is in sun glint.
        water: Is water surface.
        snow: Is snow.
        ice: Is ice.
        coast: Is coast.
        desert: Is desert.
        land: Is land.
        shadow: Has cloud shadow.
        smoke: Has smoke/non-cloud obstruction.
    """
    # Bit 5: no snow/ice
    if not snow and not ice:
        set_bit(testbits, BIT_NO_SNOW_ICE)

    # Bit 3: day
    if day:
        set_bit(testbits, BIT_DAY)

    # Bit 4: no sunglint
    if not snglnt or not water:
        set_bit(testbits, BIT_NO_SUNGLINT)

    # Bits 6,7: surface type
    if land:
        set_bit(testbits, BIT_COAST)   # bit 6
        set_bit(testbits, BIT_DESERT)  # bit 7
    elif coast:
        set_bit(testbits, BIT_COAST)   # bit 6
    elif desert:
        set_bit(testbits, BIT_DESERT)  # bit 7
    # Default: water (bits 6,7 both 0)

    # Bit 10: shadow (cleared if shadow detected)
    if shadow:
        clear_bit(testbits, BIT_SHADOW)

    # Bit 8: NCO (cleared if smoke detected)
    if smoke:
        clear_bit(testbits, BIT_NCO)


@njit(cache=True)
def set_unused_bits(testbits: np.ndarray) -> None:
    """Set unused/spare bits (12, 24, 31).

    Args:
        testbits: 6-byte array (modified in-place).
    """
    set_bit(testbits, 12)  # cloud adjacency (unused, always set)
    set_bit(testbits, 24)  # temporal consistency (unused, always set)
    set_bit(testbits, 31)  # spare bit


@njit(cache=True)
def convert_cloud_mask(testbits: np.ndarray) -> tuple[int, int]:
    """Convert 6-byte cloud mask to scalar value (0-3).

    Port of convert_cloud_mask() in fylat_fy3mersi_cloud_mask.f90.

    Reads bits from byte 0:
    - b0 = bit 0 (processed)
    - b1 = bit 1 (confidence LSB)
    - b2 = bit 2 (confidence MSB)

    Output:
    - b0==0: 0 (not processed, treated as cloudy)
    - b0==1, b2==0, b1==0: 0 (cloudy)
    - b0==1, b2==0, b1==1: 1 (probably cloudy)
    - b0==1, b2==1, b1==0: 2 (probably clear)
    - b0==1, b2==1, b1==1: 3 (confident clear)

    Args:
        testbits: 6-byte array.

    Returns:
        Tuple of (cloud_mask_value, processed_flag).
    """
    byte = testbits[0]
    b0 = (byte >> 0) & 1  # processed
    b1 = (byte >> 1) & 1  # confidence LSB
    b2 = (byte >> 2) & 1  # confidence MSB

    if b0 == 0:
        return (0, 0)
    elif b2 == 0 and b1 == 0:
        return (0, 1)
    elif b2 == 0 and b1 == 1:
        return (1, 1)
    elif b2 == 1 and b1 == 0:
        return (2, 1)
    else:  # b2==1 and b1==1
        return (3, 1)
