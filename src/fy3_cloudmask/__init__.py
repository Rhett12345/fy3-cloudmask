"""FY-3D MERSI-II Cloud Mask Retrieval System - Python Implementation."""

__version__ = "3.2.0"

from .algorithm import run_cloud_mask_pixel, run_cloud_mask_swath, CloudMaskResult
from .algorithm.native_backend import is_native_available, get_backend_info
from .config import load_config, FY3Config
from .output import (
    write_cloud_mask, write_cloud_amount, write_combined_product,
    compute_cloud_amount, compute_cloud_amount_with_coords,
)

__all__ = [
    # Algorithm
    'run_cloud_mask_pixel',
    'run_cloud_mask_swath',
    'CloudMaskResult',
    # Native backend
    'is_native_available',
    'get_backend_info',
    # Config
    'load_config',
    'FY3Config',
    # Output
    'write_cloud_mask',
    'write_cloud_amount',
    'write_combined_product',
    'compute_cloud_amount',
    'compute_cloud_amount_with_coords',
]
