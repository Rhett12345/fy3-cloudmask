"""Cloud amount computation from 5x5 pixel boxes.

Port of fylat_fy3mersi_cloud_amount.f90.

Computes 5km cloud cover from 5x5 pixel boxes of the 1km cloud mask.
"""

from __future__ import annotations

import numpy as np
from numba import njit

from ..constants import (
    CLOUD_AMOUNT_BOX_SIZE,
    CLOUD_AMOUNT_MIN_VALID,
    CLOUD_AMOUNT_MAX_VALID,
)


def compute_cloud_amount(
    cm_tmp: np.ndarray,
    qa_processed: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute 5km cloud amount from 1km cloud mask.

    Port of fy3mersi_cloud_amount() in fylat_fy3mersi_cloud_mask.f90.

    The 1km cloud mask (2048x2000) is divided into 5x5 pixel boxes.
    For each box, cloud cover is computed as the percentage of cloudy pixels
    among valid pixels.

    Args:
        cm_tmp: (n_elem, n_line) int32 cloud mask values (0-3).
        qa_processed: (n_elem, n_line) int32 processed flag (1=processed, 0=not).

    Returns:
        Tuple of (cloud_amount, cloud_amount_qa):
        - cloud_amount: (n_elem_5km, n_line_5km) uint8 cloud cover percentage (0-100).
        - cloud_amount_qa: (n_elem_5km, n_line_5km) uint8 quality flag (0=bad, 1=low, 2=high).
    """
    n_elem, n_line = cm_tmp.shape
    box_size = CLOUD_AMOUNT_BOX_SIZE

    # Output dimensions (5km grid)
    n_elem_5km = n_elem // box_size
    n_line_5km = n_line // box_size

    cloud_amount = np.full((n_elem_5km, n_line_5km), 255, dtype=np.uint8)
    cloud_amount_qa = np.zeros((n_elem_5km, n_line_5km), dtype=np.uint8)

    # Process each 5x5 box
    for j in range(n_line_5km):
        for i in range(n_elem_5km):
            # Extract 5x5 box
            row_start = j * box_size
            row_end = row_start + box_size
            col_start = i * box_size
            col_end = col_start + box_size

            cm_box = cm_tmp[col_start:col_end, row_start:row_end]
            qa_box = qa_processed[col_start:col_end, row_start:row_end]

            # Compute cloud cover
            cc, qflag = _calculate_cloud_cover_5x5(cm_box, qa_box)

            cloud_amount[i, j] = cc
            cloud_amount_qa[i, j] = qflag

    return cloud_amount, cloud_amount_qa


def _calculate_cloud_cover_5x5(
    cm_box: np.ndarray,
    qa_box: np.ndarray,
) -> tuple[int, int]:
    """Calculate cloud cover for a 5x5 pixel box.

    Port of calculate_cloud_cover_5x5() subroutine.

    Args:
        cm_box: (5, 5) cloud mask values (0-3).
        qa_box: (5, 5) processed flags (1=processed, 0=not).

    Returns:
        Tuple of (cloud_cover_percent, quality_flag).
        - cloud_cover_percent: 0-100, or 255 for invalid.
        - quality_flag: 0=bad, 1=low quality, 2=high quality.
    """
    num_valid = 0
    num_cloudy = 0

    for ii in range(5):
        for jj in range(5):
            if qa_box[ii, jj] == 1:
                num_valid += 1
                if cm_box[ii, jj] < 2:  # cloudy (0) or probably cloudy (1)
                    num_cloudy += 1

    # Quality determination
    if num_valid == CLOUD_AMOUNT_MAX_VALID:
        # High quality: all 25 pixels valid
        qflag = 2
        cc = int((num_cloudy / num_valid) * 100.0)
    elif num_valid > CLOUD_AMOUNT_MIN_VALID:
        # Low quality: 16-24 valid pixels
        qflag = 1
        cc = int((num_cloudy / num_valid) * 100.0)
    else:
        # Bad quality: 15 or fewer valid pixels
        qflag = 0
        cc = 255  # Invalid

    return cc, qflag


def compute_cloud_amount_with_coords(
    cm_tmp: np.ndarray,
    qa_processed: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute 5km cloud amount with geolocation.

    Args:
        cm_tmp: (n_elem, n_line) int32 cloud mask values (0-3).
        qa_processed: (n_elem, n_line) int32 processed flag.
        lon: (n_elem, n_line) float64 longitude array.
        lat: (n_elem, n_line) float64 latitude array.

    Returns:
        Tuple of (cloud_amount, cloud_amount_qa, lon_5km, lat_5km).
    """
    n_elem, n_line = cm_tmp.shape
    box_size = CLOUD_AMOUNT_BOX_SIZE

    n_elem_5km = n_elem // box_size
    n_line_5km = n_line // box_size

    cloud_amount = np.full((n_elem_5km, n_line_5km), 255, dtype=np.uint8)
    cloud_amount_qa = np.zeros((n_elem_5km, n_line_5km), dtype=np.uint8)
    lon_5km = np.zeros((n_elem_5km, n_line_5km), dtype=np.float64)
    lat_5km = np.zeros((n_elem_5km, n_line_5km), dtype=np.float64)

    for j in range(n_line_5km):
        for i in range(n_elem_5km):
            row_start = j * box_size
            row_end = row_start + box_size
            col_start = i * box_size
            col_end = col_start + box_size

            cm_box = cm_tmp[col_start:col_end, row_start:row_end]
            qa_box = qa_processed[col_start:col_end, row_start:row_end]

            cc, qflag = _calculate_cloud_cover_5x5(cm_box, qa_box)

            cloud_amount[i, j] = cc
            cloud_amount_qa[i, j] = qflag

            # Use center pixel of 5x5 box for geolocation
            center_col = col_start + 2
            center_row = row_start + 2
            lon_5km[i, j] = lon[center_col, center_row]
            lat_5km[i, j] = lat[center_col, center_row]

    return cloud_amount, cloud_amount_qa, lon_5km, lat_5km
