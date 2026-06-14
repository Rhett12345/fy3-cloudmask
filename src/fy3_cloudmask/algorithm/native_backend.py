"""Native C++/Fortran backend for the cloud mask algorithm.

This module provides the high-performance cloud mask engine implemented in
Fortran with C++ OpenMP parallelization. When available, it replaces the
pure Python/Numba implementation for a 30-100x speedup.

Architecture:
    Python (config/CLI/IO) -> C++ (OpenMP pixel loop) -> Fortran (core algorithm)

Usage:
    from fy3_cloudmask.algorithm.native_backend import process_swath_native

    result = process_swath_native(
        ref_vis, tbb_ir, lat, lon, satzen, solzen, relaz, glint,
        sfctmp, pmsl, uwind, vwind, tpw, elev, eco, snow_mask, btclr,
        n_elem, n_line
    )
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

# Try to import the native module
try:
    from fy3_cloudmask import _cloudmask_native
    NATIVE_AVAILABLE = True
    logger.info("Native C++/Fortran backend available (version %s)",
                getattr(_cloudmask_native, '__version__', 'unknown'))
except ImportError:
    NATIVE_AVAILABLE = False
    logger.info("Native backend not available, using pure Python/Numba")


def is_native_available() -> bool:
    """Check if the native C++/Fortran backend is available."""
    return NATIVE_AVAILABLE


def get_backend_info() -> dict:
    """Get information about the available backend."""
    if NATIVE_AVAILABLE:
        return {
            'backend': 'C++/Fortran (OpenMP)',
            'version': getattr(_cloudmask_native, '__version__', 'unknown'),
            'module': str(_cloudmask_native),
        }
    else:
        return {
            'backend': 'Python/Numba',
            'version': 'N/A',
            'module': 'fy3_cloudmask.algorithm.cloud_mask',
        }


def process_swath_native(
    ref_vis: np.ndarray,      # (nElem, nLine, 19) float32
    tbb_ir: np.ndarray,       # (nElem, nLine, 6) float32
    lat: np.ndarray,          # (nElem, nLine) float32
    lon: np.ndarray,          # (nElem, nLine) float32
    satzen: np.ndarray,       # (nElem, nLine) float32
    solzen: np.ndarray,       # (nElem, nLine) float32
    relaz: np.ndarray,        # (nElem, nLine) float32
    glint: np.ndarray,        # (nElem, nLine) float32
    sfctmp: np.ndarray,       # (nElem, nLine) float32
    pmsl: np.ndarray,         # (nElem, nLine) float32
    uwind: np.ndarray,        # (nElem, nLine) float32
    vwind: np.ndarray,        # (nElem, nLine) float32
    tpw: np.ndarray,          # (nElem, nLine) float32
    elev: np.ndarray,         # (nElem, nLine) float32
    eco: np.ndarray,          # (nElem, nLine) int8
    lsf: np.ndarray,          # (nElem, nLine) int8
    snow_mask: np.ndarray,    # (nElem, nLine) int8
    btclr: np.ndarray,        # (nElem, nLine, 7) float32
    n_elem: int,
    n_line: int,
) -> dict:
    """Process a full swath through the native cloud mask engine.

    Parameters
    ----------
    ref_vis : ndarray, shape (nElem, nLine, 19), float32
        Visible/NIR reflectance for 19 channels.
    tbb_ir : ndarray, shape (nElem, nLine, 6), float32
        IR brightness temperature for 6 channels.
    lat, lon : ndarray, shape (nElem, nLine), float32
        Geolocation.
    satzen, solzen, relaz, glint : ndarray, shape (nElem, nLine), float32
        Geometry angles (degrees).
    sfctmp, pmsl : ndarray, shape (nElem, nLine), float32
        NWP surface temperature (K) and pressure (hPa).
    uwind, vwind : ndarray, shape (nElem, nLine), float32
        NWP wind components (m/s).
    tpw : ndarray, shape (nElem, nLine), float32
        Total precipitable water (cm).
    elev : ndarray, shape (nElem, nLine), float32
        Elevation (m).
    eco : ndarray, shape (nElem, nLine), int8
        IGBP ecosystem type.
    lsf : ndarray, shape (nElem, nLine), int8
        Land-sea flag (0=water, 1=land, 2=coast, 3=shallow_lake, 4=land).
    snow_mask : ndarray, shape (nElem, nLine), int8
        NISE snow/ice mask.
    btclr : ndarray, shape (nElem, nLine, 7), float32
        Clear-sky brightness temperatures from RTM.
    n_elem, n_line : int
        Swath dimensions.

    Returns
    -------
    dict with keys:
        cm_bitarray : ndarray, shape (nElem, nLine, 6), int8
        qa_bitarray : ndarray, shape (nElem, nLine, 10), int8
        cloud_mask : ndarray, shape (nElem, nLine), int32
        confidence : ndarray, shape (nElem, nLine), float32
        nmtests, nbands, shadow, smoke : ndarray, shape (nElem, nLine), int32
    """
    if not NATIVE_AVAILABLE:
        raise RuntimeError(
            "Native backend not available. Build with: cd ext && make install"
        )

    # Ensure arrays are contiguous and correct dtype
    ref_vis = np.ascontiguousarray(ref_vis, dtype=np.float32)
    tbb_ir = np.ascontiguousarray(tbb_ir, dtype=np.float32)
    lat = np.ascontiguousarray(lat, dtype=np.float32)
    lon = np.ascontiguousarray(lon, dtype=np.float32)
    satzen = np.ascontiguousarray(satzen, dtype=np.float32)
    solzen = np.ascontiguousarray(solzen, dtype=np.float32)
    relaz = np.ascontiguousarray(relaz, dtype=np.float32)
    glint = np.ascontiguousarray(glint, dtype=np.float32)
    sfctmp = np.ascontiguousarray(sfctmp, dtype=np.float32)
    pmsl = np.ascontiguousarray(pmsl, dtype=np.float32)
    uwind = np.ascontiguousarray(uwind, dtype=np.float32)
    vwind = np.ascontiguousarray(vwind, dtype=np.float32)
    tpw = np.ascontiguousarray(tpw, dtype=np.float32)
    elev = np.ascontiguousarray(elev, dtype=np.float32)
    eco = np.ascontiguousarray(eco, dtype=np.int8)
    lsf = np.ascontiguousarray(lsf, dtype=np.int8)
    snow_mask = np.ascontiguousarray(snow_mask, dtype=np.int8)
    btclr = np.ascontiguousarray(btclr, dtype=np.float32)

    # Get code root path for threshold file lookup
    import os
    code_root = os.environ.get('FY3_CODE_ROOT', '')
    if not code_root:
        code_root = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'coeff')
        code_root = os.path.normpath(code_root)
    if not code_root.endswith('/'):
        code_root += '/'

    # Call native engine
    return _cloudmask_native.process_swath(
        ref_vis, tbb_ir, lat, lon, satzen, solzen, relaz, glint,
        sfctmp, pmsl, uwind, vwind, tpw, elev, eco, lsf, snow_mask, btclr,
        n_elem, n_line, code_root
    )
