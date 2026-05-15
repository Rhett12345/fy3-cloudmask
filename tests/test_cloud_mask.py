"""Tests for the cloud mask algorithm."""

import numpy as np
import pytest

from fy3_cloudmask.algorithm import run_cloud_mask_pixel, CloudMaskResult
from fy3_cloudmask.algorithm.bitops import init_testbits, init_qa_bits, convert_cloud_mask
from fy3_cloudmask.algorithm.confidence import conf_test, encode_confidence
from fy3_cloudmask.algorithm.surface_classifier import PixelFlags, classify_pixel_surface
from fy3_cloudmask.constants import BAD_DATA


class TestConfidence:
    """Test confidence functions."""

    def test_conf_test_clear(self):
        """Test confidence for clear-sky values."""
        # Value above hicut should return 1.0
        c = conf_test(0.5, 0.0, 0.2, 0.4, 1.0)
        assert c == 1.0

    def test_conf_test_cloudy(self):
        """Test confidence for cloudy values."""
        # Value below locut should return 0.0
        c = conf_test(-0.1, 0.0, 0.2, 0.4, 1.0)
        assert c == 0.0

    def test_conf_test_midrange(self):
        """Test confidence for mid-range values."""
        # Value at midpoint should return 0.5
        c = conf_test(0.2, 0.0, 0.2, 0.4, 1.0)
        assert abs(c - 0.5) < 0.01

    def test_encode_confidence(self):
        """Test confidence encoding."""
        # Confident clear
        bit1, bit2 = encode_confidence(0.995)
        assert bit1 == 1 and bit2 == 1

        # Probably clear
        bit1, bit2 = encode_confidence(0.97)
        assert bit1 == 0 and bit2 == 1

        # Probably cloudy
        bit1, bit2 = encode_confidence(0.80)
        assert bit1 == 1 and bit2 == 0

        # Cloudy
        bit1, bit2 = encode_confidence(0.50)
        assert bit1 == 0 and bit2 == 0


class TestBitops:
    """Test bit operations."""

    def test_init_testbits(self):
        """Test testbits initialization."""
        tb = init_testbits()
        assert len(tb) == 6
        assert tb.dtype == np.uint8

    def test_init_qa_bits(self):
        """Test QA bits initialization."""
        qa = init_qa_bits()
        assert len(qa) == 10
        assert qa.dtype == np.uint8

    def test_convert_cloud_mask(self):
        """Test cloud mask conversion."""
        # Not processed (all zeros)
        tb = np.zeros(6, dtype=np.uint8)
        mask, processed = convert_cloud_mask(tb)
        assert mask == 0
        assert processed == 0

        # Processed, confident clear (bits 0, 1, 2 set)
        tb = np.zeros(6, dtype=np.uint8)
        tb[0] = 0b00000111  # bits 0, 1, 2
        mask, processed = convert_cloud_mask(tb)
        assert mask == 3  # Confident clear
        assert processed == 1


class TestPixelFlags:
    """Test pixel flags."""

    def test_default_flags(self):
        """Test default flag values."""
        flags = PixelFlags()
        assert flags.land == False
        assert flags.water == False
        assert flags.coast == False
        assert flags.desert == False
        assert flags.day == False
        assert flags.night == False
        assert flags.snow == False
        assert flags.ice == False


class TestCloudMaskPixel:
    """Test cloud mask pixel processing."""

    def create_test_pixel(self):
        """Create a test pixel with reasonable values."""
        pxldat = np.zeros(25, dtype=np.float64)

        # Visible channels (0-18)
        pxldat[0] = 0.10   # 0.47um
        pxldat[1] = 0.12   # 0.55um
        pxldat[2] = 0.15   # 0.65um
        pxldat[3] = 0.18   # 0.86um
        pxldat[4] = 0.02   # 1.38um
        pxldat[5] = 0.10   # 1.64um
        pxldat[6] = 0.12   # 2.13um

        # IR channels (19-24)
        pxldat[19] = 280.0  # 3.8um
        pxldat[20] = 282.0  # 4.05um
        pxldat[21] = 260.0  # 7.3um
        pxldat[22] = 275.0  # 8.5um
        pxldat[23] = 290.0  # 11.0um
        pxldat[24] = 288.0  # 12.0um

        return pxldat

    def test_daytime_land_pixel(self):
        """Test processing of daytime land pixel."""
        pxldat = self.create_test_pixel()

        # Clear-sky BT from RTM (simplified)
        bt_clr = np.array([280, 282, 260, 275, 290, 288, 0], dtype=np.float64)

        # Simple thresholds
        thresholds = {
            'land_day': {
                'ref064': [0.24, 0.20, 0.16, 1.0],
                'ref138': [0.04, 0.035, 0.03, 1.0],
            },
            'pfmft': {
                'bt_11_max': 310.0,
                'btd_min': 0.0,
                'land': [4.0, 3.5, 3.0, 1.0],
            },
            'nfmft': {
                'max_threshold': 1.50,
                'land': [-23.0, -22.5, -22.0, 1.0],
            },
        }

        result = run_cloud_mask_pixel(
            pxldat=pxldat,
            lat=35.0,
            lon=115.0,
            elevation=100.0,
            lsf=1,  # Land
            sza=30.0,  # Daytime
            vza=10.0,
            glint_angle=60.0,
            eco_type=1,  # Forest
            snow_mask_val=0,
            sst=BAD_DATA,
            nwp_sfctmp=290.0,
            nwp_pmsl=1013.0,
            nwp_u_wind=5.0,
            nwp_v_wind=3.0,
            nwp_precip_water=20.0,
            sensor_id=21,
            bt_clr=bt_clr,
            thresholds=thresholds,
        )

        assert isinstance(result, CloudMaskResult)
        assert result.cloud_mask in [0, 1, 2, 3]
        assert 0.0 <= result.confidence <= 1.0
        assert result.n_tests >= 0
        assert len(result.testbits) == 6
        assert len(result.qa_bits) == 10

    def test_nighttime_ocean_pixel(self):
        """Test processing of nighttime ocean pixel."""
        pxldat = self.create_test_pixel()

        bt_clr = np.array([280, 282, 260, 275, 290, 288, 0], dtype=np.float64)

        thresholds = {
            'ocean_nite': {
                'bt11': [235.0, 270.0, 265.0, 260.0],
                'btd_11_12': [3.0],
                'btd_11_4': [-14.0, -12.0, -10.0, 1.0],
                'btd_86_73': [-5.0, -3.0, -1.0, 1.0],
                'variability': [0.3, 0.6, 0.9, 1.0],
            },
            'pfmft': {
                'bt_11_max': 310.0,
                'btd_min': 0.0,
                'ocean': [4.0, 3.5, 3.0, 1.0],
            },
            'nfmft': {
                'max_threshold': 1.50,
                'ocean': [-23.0, -22.5, -22.0, 1.0],
            },
        }

        result = run_cloud_mask_pixel(
            pxldat=pxldat,
            lat=0.0,
            lon=170.0,
            elevation=-100.0,
            lsf=0,  # Ocean
            sza=100.0,  # Nighttime
            vza=10.0,
            glint_angle=60.0,
            eco_type=0,
            snow_mask_val=0,
            sst=295.0,
            nwp_sfctmp=295.0,
            nwp_pmsl=1013.0,
            nwp_u_wind=5.0,
            nwp_v_wind=3.0,
            nwp_precip_water=30.0,
            sensor_id=21,
            bt_clr=bt_clr,
            thresholds=thresholds,
        )

        assert isinstance(result, CloudMaskResult)
        assert result.cloud_mask in [0, 1, 2, 3]
        assert 0.0 <= result.confidence <= 1.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
