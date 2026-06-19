"""End-to-end pipeline test with real FY-3D MERSI-II data (Fortran native backend).

Tests the complete cloud mask workflow using the native C++/Fortran backend.
"""

import logging
import time
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml

from fy3_cloudmask.algorithm import process_swath_native, is_native_available
from fy3_cloudmask.output import (
    compute_cloud_amount,
    write_cloud_mask,
    write_cloud_amount,
    write_combined_product,
)
from fy3_cloudmask.constants import BAD_DATA, SZA_NIGHT

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# Test data paths
L1B_FILE = '/home/liusy2020/yuq/fydata/mersi_20230606/FY3D_MERSI_GBAL_L1_20230606_0500_1000M_MS.HDF'
GEO_FILE = '/home/liusy2020/yuq/fydata/mersi_20230606/FY3D_MERSI_GBAL_L1_20230606_0500_GEO1K_MS.HDF'
THRESHOLDS_FILE = '/home/liusy2020/yuq/cloudmask/fy3_cloudmask/config/thresholds/mersi_ii3d_v8.yaml'
OUTPUT_DIR = '/home/liusy2020/yuq/data-yuq/cloudmask_test/output'


@pytest.mark.skip(reason="Requires real FY-3D data and native backend build")
def test_native_backend_basic():
    """Basic smoke test for native backend."""
    assert is_native_available(), "Native backend must be available"
