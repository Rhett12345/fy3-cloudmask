"""FY-3D MERSI-II Cloud Mask Retrieval System — Fortran native backend."""

__version__ = "3.4.0"

from .algorithm import CloudMaskResult, is_native_available, get_backend_info, process_swath_native
from .config import load_config, FY3Config
from .output import (
    write_cloud_mask, write_cloud_amount, write_combined_product,
    compute_cloud_amount, compute_cloud_amount_with_coords,
)

__all__ = [
    # Algorithm
    'CloudMaskResult',
    # Native backend
    'is_native_available',
    'get_backend_info',
    'process_swath_native',
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
