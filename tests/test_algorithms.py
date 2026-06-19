"""Tests for cloud mask algorithm (Fortran native backend bridge).

The algorithm is implemented in Fortran (src/fortran/cloudmask/).
These tests verify the Python-to-Fortran bridge is working correctly.
"""

import numpy as np
import pytest

from fy3_cloudmask.algorithm import is_native_available, get_backend_info, CloudMaskResult


def test_native_backend_available():
    """Test that the Fortran native backend is available."""
    assert is_native_available(), "Native Fortran backend must be available"


def test_backend_info():
    """Test backend info returns valid data."""
    info = get_backend_info()
    assert 'backend' in info
    assert 'C++/Fortran' in info['backend']


def test_cloud_mask_result():
    """Test CloudMaskResult dataclass."""
    result = CloudMaskResult(cloud_mask=3, confidence=0.99, n_tests=5, n_bands=6)
    assert result.cloud_mask == 3
    assert result.confidence == 0.99


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
