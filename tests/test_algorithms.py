"""Comprehensive tests for cloud mask algorithm components."""

import numpy as np
import pytest

from fy3_cloudmask.algorithm.confidence import (
    conf_test, conf_test_thresholds, encode_confidence, compute_group_confidence,
)
from fy3_cloudmask.algorithm.bitops import (
    set_bit, clear_bit, check_bit, init_testbits, init_qa_bits,
    fill_bit_pixel, proc_path, convert_cloud_mask, set_unused_bits,
)
from fy3_cloudmask.algorithm.spatial import (
    tview, get_regional_mean, get_regional_std, get_regional_diff,
    check_reg_uniformity, spatial_var_test,
)
from fy3_cloudmask.algorithm.surface_classifier import (
    PixelFlags, classify_pixel_surface, detect_snow_ndsi,
    is_desert_ecosystem, is_vrat_disabled, is_greenland, is_new_zealand,
)
from fy3_cloudmask.constants import BAD_DATA


class TestConfidenceFunctions:
    """Test confidence calculation functions."""

    def test_conf_test_linear(self):
        """Test linear S-curve confidence."""
        # Value at lower cutoff
        c = conf_test(0.0, 0.0, 0.5, 1.0, 1.0)
        assert c == 0.0

        # Value at upper cutoff
        c = conf_test(1.0, 0.0, 0.5, 1.0, 1.0)
        assert c == 1.0

        # Value at midpoint
        c = conf_test(0.5, 0.0, 0.5, 1.0, 1.0)
        assert abs(c - 0.5) < 0.01

    def test_conf_test_flipped(self):
        """Test flipped S-curve (lower values = more clear)."""
        # When hicut < locut, curve is flipped
        c = conf_test(0.0, 1.0, 0.5, 0.0, 1.0)
        assert c == 1.0

        c = conf_test(1.0, 1.0, 0.5, 0.0, 1.0)
        assert c == 0.0

    def test_conf_test_out_of_range(self):
        """Test out-of-range values."""
        # Above range
        c = conf_test(2.0, 0.0, 0.5, 1.0, 1.0)
        assert c == 1.0

        # Below range
        c = conf_test(-1.0, 0.0, 0.5, 1.0, 1.0)
        assert c == 0.0

    def test_conf_test_thresholds(self):
        """Test threshold-based confidence."""
        thresholds = np.array([0.0, 0.5, 1.0, 1.0], dtype=np.float64)
        c = conf_test_thresholds(0.5, thresholds)
        assert abs(c - 0.5) < 0.01

    def test_encode_confidence_levels(self):
        """Test confidence encoding at all levels."""
        # Confident clear
        assert encode_confidence(0.999) == (1, 1)
        assert encode_confidence(0.995) == (1, 1)

        # Probably clear
        assert encode_confidence(0.97) == (0, 1)
        assert encode_confidence(0.96) == (0, 1)

        # Probably cloudy
        assert encode_confidence(0.80) == (1, 0)
        assert encode_confidence(0.67) == (1, 0)

        # Cloudy
        assert encode_confidence(0.50) == (0, 0)
        assert encode_confidence(0.00) == (0, 0)

    def test_compute_group_confidence(self):
        """Test group confidence computation."""
        # Single group
        c = compute_group_confidence([0.8])
        assert abs(c - 0.8) < 0.01

        # Multiple groups - geometric mean
        c = compute_group_confidence([0.8, 0.9])
        expected = (0.8 * 0.9) ** 0.5
        assert abs(c - expected) < 0.01

        # Three groups
        c = compute_group_confidence([0.8, 0.9, 0.7])
        expected = (0.8 * 0.9 * 0.7) ** (1/3)
        assert abs(c - expected) < 0.01


class TestBitOperations:
    """Test bit manipulation functions."""

    def test_set_bit(self):
        """Test setting bits."""
        tb = np.zeros(6, dtype=np.uint8)
        set_bit(tb, 0)
        assert tb[0] == 1

        set_bit(tb, 7)
        assert tb[0] == 129  # 1 + 128

        set_bit(tb, 8)
        assert tb[1] == 1

    def test_clear_bit(self):
        """Test clearing bits."""
        tb = np.array([255, 255, 0, 0, 0, 0], dtype=np.uint8)
        clear_bit(tb, 0)
        assert tb[0] == 254

        clear_bit(tb, 7)
        assert tb[0] == 126

    def test_check_bit(self):
        """Test checking bits."""
        tb = np.array([5, 0, 0, 0, 0, 0], dtype=np.uint8)  # 0b00000101
        assert check_bit(tb, 0) == True
        assert check_bit(tb, 1) == False
        assert check_bit(tb, 2) == True
        assert check_bit(tb, 3) == False

    def test_init_testbits(self):
        """Test testbits initialization."""
        tb = init_testbits()
        assert len(tb) == 6
        assert tb.dtype == np.uint8
        # Check initial bits are set
        assert check_bit(tb, 8)   # NCO
        assert check_bit(tb, 9)   # thin_cirrus_solar
        assert check_bit(tb, 10)  # shadow
        assert check_bit(tb, 11)  # thin_cirrus_ir
        assert check_bit(tb, 12)  # cloud_adj
        assert check_bit(tb, 24)  # temporal
        assert check_bit(tb, 28)  # suspended_dust
        assert check_bit(tb, 31)  # spare

    def test_fill_bit_pixel(self):
        """Test QA quality assembly."""
        tb = init_testbits()
        qa = init_qa_bits()

        # High quality (7+ tests)
        fill_bit_pixel(10, 2, False, False, tb, qa)
        assert qa[0] == 15  # 0b00001111

        # Medium quality (3-6 tests)
        qa = init_qa_bits()
        fill_bit_pixel(5, 2, False, False, tb, qa)
        assert qa[0] == 13  # 0b00001101

        # Low quality (1-2 tests)
        qa = init_qa_bits()
        fill_bit_pixel(2, 2, False, False, tb, qa)
        assert qa[0] == 9   # 0b00001001

    def test_proc_path(self):
        """Test processing path bits."""
        tb = init_testbits()

        # Land, day, no snow
        proc_path(tb, day=True, snglnt=False, water=False,
                  snow=False, ice=False, coast=False, desert=False,
                  land=True, shadow=False, smoke=False)
        assert check_bit(tb, 3)  # day
        assert check_bit(tb, 4)  # no sunglint
        assert check_bit(tb, 5)  # no snow/ice
        assert check_bit(tb, 6)  # coast (set for land)
        assert check_bit(tb, 7)  # desert (set for land)

    def test_convert_cloud_mask(self):
        """Test cloud mask conversion."""
        # Not processed
        tb = np.zeros(6, dtype=np.uint8)
        mask, processed = convert_cloud_mask(tb)
        assert mask == 0
        assert processed == 0

        # Confident clear (bits 0, 1, 2)
        tb = np.zeros(6, dtype=np.uint8)
        set_bit(tb, 0)  # processed
        set_bit(tb, 1)  # conf_lsb
        set_bit(tb, 2)  # conf_msb
        mask, processed = convert_cloud_mask(tb)
        assert mask == 3
        assert processed == 1

        # Cloudy (processed, conf=0)
        tb = np.zeros(6, dtype=np.uint8)
        set_bit(tb, 0)  # processed
        mask, processed = convert_cloud_mask(tb)
        assert mask == 0
        assert processed == 1


class TestSpatialFunctions:
    """Test spatial analysis functions."""

    def test_tview_interpolation(self):
        """Test APOLLO lookup table interpolation."""
        # Center of table
        val = tview(1.0, 250.0)
        assert val > 0
        assert val < 20

        # Out of range
        val = tview(3.0, 250.0)
        assert val == 99.0

    def test_get_regional_mean(self):
        """Test 3x3 mean calculation."""
        data = np.array([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ], dtype=np.float64)
        mean = get_regional_mean(data)
        # Excludes center (5)
        expected = (1+2+3+4+6+7+8+9) / 8
        assert abs(mean - expected) < 0.01

    def test_get_regional_std(self):
        """Test 3x3 standard deviation calculation."""
        data = np.array([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ], dtype=np.float64)
        std = get_regional_std(data)
        assert std > 0

    def test_get_regional_diff(self):
        """Test 3x3 max difference calculation."""
        data = np.array([
            [1, 2, 3],
            [4, 10, 6],
            [7, 8, 9],
        ], dtype=np.float64)
        diff = get_regional_diff(data)
        # Center is 10, neighbors range from 1 to 9
        # Max diff = |10 - 1| = 9
        assert diff == 9

    def test_check_reg_uniformity(self):
        """Test 3x3 uniformity check."""
        eco_center = 1
        eco_neighbors = np.ones((3, 3), dtype=np.int32)
        snow_center = 0
        snow_neighbors = np.zeros((3, 3), dtype=np.int32)
        lsf_center = 1
        lsf_neighbors = np.ones((3, 3), dtype=np.int32)

        uniform, is_coast, is_land, is_water, n_land, n_water, n_coast = check_reg_uniformity(
            eco_center, eco_neighbors, snow_center, snow_neighbors,
            lsf_center, lsf_neighbors, False, False, False,
        )
        assert uniform == True
        assert is_land == True

    def test_spatial_var_test(self):
        """Test spatial variability test."""
        # Uniform field
        data = np.full((3, 3), 280.0, dtype=np.float64)
        assert spatial_var_test(data, 0.6) == False

        # Variable field
        data = np.array([
            [280, 285, 280],
            [285, 290, 285],
            [280, 285, 280],
        ], dtype=np.float64)
        assert spatial_var_test(data, 0.6) == True


class TestSurfaceClassifier:
    """Test surface classification functions."""

    def test_is_desert_ecosystem(self):
        """Test desert ecosystem detection."""
        assert is_desert_ecosystem(8) == True
        assert is_desert_ecosystem(46) == True
        assert is_desert_ecosystem(1) == False

    def test_is_vrat_disabled(self):
        """Test VRAT disabled ecosystems."""
        assert is_vrat_disabled(2) == True
        assert is_vrat_disabled(8) == True
        assert is_vrat_disabled(1) == False

    def test_is_greenland(self):
        """Test Greenland detection."""
        assert is_greenland(70.0, -40.0, 0) == True
        assert is_greenland(0.0, 0.0, 0) == False

    def test_is_new_zealand(self):
        """Test New Zealand detection."""
        assert is_new_zealand(-40.0, 175.0) == True
        assert is_new_zealand(0.0, 0.0) == False

    def test_pixel_flags_defaults(self):
        """Test default pixel flags."""
        flags = PixelFlags()
        assert flags.land == False
        assert flags.water == False
        assert flags.coast == False
        assert flags.desert == False
        assert flags.day == False
        assert flags.night == False
        assert flags.snow == False
        assert flags.ice == False
        assert flags.snglnt == False
        assert flags.visusd == True
        assert flags.vrused == True
        assert flags.hi_elev == False
        assert flags.antarctic == False
        assert flags.uniform == True
        assert flags.bad_value == False
        assert flags.process == True


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
