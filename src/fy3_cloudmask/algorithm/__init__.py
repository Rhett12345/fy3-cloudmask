"""Cloud mask algorithm modules."""

from .cloud_mask import run_cloud_mask_pixel, run_cloud_mask_swath, CloudMaskResult
from .confidence import conf_test, conf_test_thresholds, encode_confidence
from .bitops import set_bit, clear_bit, check_bit, init_testbits, init_qa_bits, convert_cloud_mask
from .spatial import tview, get_regional_mean, get_regional_std, get_regional_diff
from .surface_classifier import classify_pixel_surface, PixelFlags, detect_snow_ndsi
from .native_backend import is_native_available, process_swath_native, get_backend_info

__all__ = [
    # Main algorithm
    'run_cloud_mask_pixel', 'run_cloud_mask_swath', 'CloudMaskResult',
    # Native backend
    'is_native_available', 'process_swath_native', 'get_backend_info',
    # Confidence
    'conf_test', 'conf_test_thresholds', 'encode_confidence',
    # Bit operations
    'set_bit', 'clear_bit', 'check_bit', 'init_testbits', 'init_qa_bits', 'convert_cloud_mask',
    # Spatial
    'tview', 'get_regional_mean', 'get_regional_std', 'get_regional_diff',
    # Surface classification
    'classify_pixel_surface', 'PixelFlags', 'detect_snow_ndsi',
]
