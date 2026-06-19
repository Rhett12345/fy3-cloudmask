"""Tests for the cloud mask algorithm (Fortran native backend)."""

import numpy as np
import pytest

from fy3_cloudmask.algorithm import CloudMaskResult, is_native_available


def test_cloud_mask_result_dataclass():
    """Test CloudMaskResult dataclass construction."""
    result = CloudMaskResult(cloud_mask=3, confidence=0.99, n_tests=5, n_bands=6)
    assert result.cloud_mask == 3
    assert result.confidence == 0.99
    assert result.n_tests == 5
    assert result.n_bands == 6


def test_native_backend_available():
    """Test that native backend is available."""
    assert is_native_available(), "Native Fortran backend must be available"
