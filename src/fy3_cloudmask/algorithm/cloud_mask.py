"""Cloud mask result types.

The cloud mask algorithm is implemented in Fortran (src/fortran/cloudmask/).
This module provides Python-side types and the native backend bridge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class CloudMaskResult:
    """Result of cloud mask processing for a single pixel."""
    cloud_mask: int = 0          # 0=cloudy, 1=prob cloudy, 2=prob clear, 3=confident clear
    confidence: float = 1.0      # Raw confidence value
    n_tests: int = 0             # Number of tests applied
    n_bands: int = 0             # Number of bands used
    testbits: np.ndarray = None  # 6-byte test bits
    qa_bits: np.ndarray = None   # 10-byte QA bits
