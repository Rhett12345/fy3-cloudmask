"""Cloud mask algorithm — Fortran native backend only.

The algorithm is implemented in Fortran (src/fortran/cloudmask/).
This package provides the Python-to-Fortran bridge via pybind11.
"""

from .cloud_mask import CloudMaskResult
from .native_backend import is_native_available, get_backend_info, process_swath_native

__all__ = [
    'CloudMaskResult',
    'is_native_available',
    'get_backend_info',
    'process_swath_native',
]
