"""S-curve confidence test function and confidence encoding.

Exact port of Fortran conf_test.f and set_confdnc.f.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ..constants import CONF_CLEAR, CONF_PROB_CLEAR, CONF_CLOUDY


@njit(cache=True)
def conf_test(val: float, locut: float, midpt: float, hicut: float, power: float) -> float:
    """S-curve confidence interpolation.

    Computes a confidence value in [0, 1] based on the test value's position
    relative to an S-curve defined by the threshold parameters.

    The S-curve maps a test value to a confidence:
    - val <= locut: confidence = 0.0 (certainly cloudy)
    - val >= hicut: confidence = 1.0 (certainly clear)
    - In between: S-shaped interpolation

    When hicut < locut, the curve is flipped (higher values = more cloudy).

    Args:
        val: Test value to evaluate.
        locut: Lower cutoff (0% confidence boundary).
        midpt: Midpoint (50% confidence point).
        hicut: Upper cutoff (100% confidence boundary).
        power: S-curve shape parameter (typically 1.0 for linear).

    Returns:
        Confidence value in [0.0, 1.0].
    """
    coeff = 2.0 ** (power - 1.0)

    if hicut > locut:
        # Normal direction: higher value = more clear
        gamma = hicut
        alpha = locut
        flipped = False
    else:
        # Flipped direction: lower value = more clear
        gamma = locut
        alpha = hicut
        flipped = True

    beta = midpt

    # Out-of-range clamping
    if not flipped:
        if val > gamma:
            return 1.0
        if val < alpha:
            return 0.0
    else:
        if val > gamma:
            return 0.0
        if val < alpha:
            return 1.0

    # In-range S-curve interpolation
    if val <= beta:
        range_val = 2.0 * (beta - alpha)
        if range_val == 0.0:
            return 0.5
        s1 = (val - alpha) / range_val
        if not flipped:
            c = coeff * (s1 ** power)
        else:
            c = 1.0 - (coeff * (s1 ** power))
    else:
        range_val = 2.0 * (beta - gamma)
        if range_val == 0.0:
            return 0.5
        s1 = (val - gamma) / range_val
        if not flipped:
            c = 1.0 - (coeff * (s1 ** power))
        else:
            c = coeff * (s1 ** power)

    # Clamp to [0, 1]
    if c < 0.0:
        c = 0.0
    elif c > 1.0:
        c = 1.0

    return c


@njit(cache=True)
def conf_test_thresholds(val: float, thresholds: np.ndarray) -> float:
    """Convenience wrapper: conf_test with a 4-element threshold array.

    Args:
        val: Test value.
        thresholds: Array of [locut, midpt, hicut, power].

    Returns:
        Confidence value in [0.0, 1.0].
    """
    return conf_test(val, thresholds[0], thresholds[1], thresholds[2], thresholds[3])


def encode_confidence(confdnc: float) -> tuple[int, int]:
    """Encode confidence value into two bit flags.

    Port of set_confdnc.f:
    - confdnc > 0.99: bits (1, 1) = confident clear
    - confdnc > 0.95: bits (0, 1) = probably clear
    - confdnc > 0.66: bits (1, 0) = probably cloudy
    - else:           bits (0, 0) = confident cloudy

    Args:
        confdnc: Confidence value in [0.0, 1.0].

    Returns:
        Tuple of (bit1, bit2) for encoding into testbits.
    """
    if confdnc > CONF_CLEAR:
        return (1, 1)  # confident clear
    elif confdnc > CONF_PROB_CLEAR:
        return (0, 1)  # probably clear
    elif confdnc > CONF_CLOUDY:
        return (1, 0)  # probably cloudy
    else:
        return (0, 0)  # confident cloudy


def compute_group_confidence(group_confidences: list[float]) -> float:
    """Compute final confidence as geometric mean of group minimums.

    In the Fortran code, each test group maintains a minimum confidence (cmin).
    The final confidence is the geometric mean of all non-trivial group minimums.

    Args:
        group_confidences: List of minimum confidence values from each test group.

    Returns:
        Final confidence value.
    """
    # Count groups that actually ran (had tests)
    active = [c for c in group_confidences if c < 1.0 or True]  # all groups count
    n_groups = len(active)
    if n_groups == 0:
        return 1.0

    product = 1.0
    for c in active:
        product *= c

    return product ** (1.0 / n_groups)
