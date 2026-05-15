"""HDF5 output writer for cloud mask products.

Port of io_module.f90 write functions.

Writes CLM (cloud mask) and CLA (cloud amount) products in HDF5 format.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import h5py
    HAS_H5PY = True
except ImportError:
    HAS_H5PY = False


def write_cloud_mask(
    filepath: str,
    cm_bitarray: np.ndarray,
    cm_qa_bitarray: np.ndarray,
    cm_tmp: np.ndarray,
    confidence: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    attributes: Optional[dict] = None,
) -> None:
    """Write CLM (cloud mask) product as HDF5.

    Args:
        filepath: Output file path.
        cm_bitarray: (n_elem, n_line, 6) uint8 testbits array.
        cm_qa_bitarray: (n_elem, n_line, 10) uint8 QA bits array.
        cm_tmp: (n_elem, n_line) int32 cloud mask values (0-3).
        confidence: (n_elem, n_line) float64 confidence values.
        lon: (n_elem, n_line) float64 longitude array.
        lat: (n_elem, n_line) float64 latitude array.
        attributes: Optional global attributes dictionary.

    Raises:
        ImportError: If h5py is not installed.
    """
    if not HAS_H5PY:
        raise ImportError("h5py is required for writing HDF5 output")

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, 'w') as f:
        # Global attributes
        f.attrs['title'] = 'FY-3D MERSI-II Cloud Mask Product'
        f.attrs['institution'] = 'National Satellite Meteorological Center, CMA'
        f.attrs['source'] = 'FYLAT Cloud Mask Retrieval System V3.2 (Python)'
        f.attrs['history'] = f'Created {datetime.now().isoformat()}'
        f.attrs['Conventions'] = 'CF-1.8'

        if attributes:
            for key, value in attributes.items():
                f.attrs[key] = value

        # Cloud Mask dataset
        ds_cm = f.create_dataset(
            'Cloud_Mask',
            data=cm_bitarray,
            dtype=np.uint8,
            compression='gzip',
            compression_opts=4,
        )
        ds_cm.attrs['long_name'] = 'Cloud Mask Test Bits'
        ds_cm.attrs['units'] = 'bit flags'
        ds_cm.attrs['_FillValue'] = 0
        ds_cm.attrs['valid_range'] = [0, 255]

        # Quality Assurance dataset
        ds_qa = f.create_dataset(
            'Quality_Assurance',
            data=cm_qa_bitarray,
            dtype=np.uint8,
            compression='gzip',
            compression_opts=4,
        )
        ds_qa.attrs['long_name'] = 'Quality Assurance Flags'
        ds_qa.attrs['units'] = 'bit flags'
        ds_qa.attrs['_FillValue'] = 0

        # Cloud Mask Values (0-3)
        ds_tmp = f.create_dataset(
            'Cloud_Mask_Value',
            data=cm_tmp,
            dtype=np.int32,
            compression='gzip',
            compression_opts=4,
        )
        ds_tmp.attrs['long_name'] = 'Cloud Mask Value'
        ds_tmp.attrs['units'] = 'category'
        ds_tmp.attrs['_FillValue'] = -999
        ds_tmp.attrs['flag_values'] = [0, 1, 2, 3]
        ds_tmp.attrs['flag_meanings'] = 'cloudy probably_cloudy probably_clear confident_clear'

        # Confidence
        ds_conf = f.create_dataset(
            'Confidence',
            data=confidence,
            dtype=np.float64,
            compression='gzip',
            compression_opts=4,
        )
        ds_conf.attrs['long_name'] = 'Cloud Detection Confidence'
        ds_conf.attrs['units'] = '1'
        ds_conf.attrs['_FillValue'] = -999.0
        ds_conf.attrs['valid_range'] = [0.0, 1.0]

        # Geolocation
        ds_lon = f.create_dataset(
            'Longitude',
            data=lon,
            dtype=np.float64,
        )
        ds_lon.attrs['long_name'] = 'Longitude'
        ds_lon.attrs['units'] = 'degrees_east'

        ds_lat = f.create_dataset(
            'Latitude',
            data=lat,
            dtype=np.float64,
        )
        ds_lat.attrs['long_name'] = 'Latitude'
        ds_lat.attrs['units'] = 'degrees_north'


def write_cloud_amount(
    filepath: str,
    cloud_amount: np.ndarray,
    cloud_amount_qa: np.ndarray,
    lon_5km: np.ndarray,
    lat_5km: np.ndarray,
    attributes: Optional[dict] = None,
) -> None:
    """Write CLA (cloud amount) product as HDF5.

    Args:
        filepath: Output file path.
        cloud_amount: (n_elem_5km, n_line_5km) uint8 cloud cover percentage.
        cloud_amount_qa: (n_elem_5km, n_line_5km) uint8 quality flag.
        lon_5km: (n_elem_5km, n_line_5km) float64 longitude array.
        lat_5km: (n_elem_5km, n_line_5km) float64 latitude array.
        attributes: Optional global attributes dictionary.

    Raises:
        ImportError: If h5py is not installed.
    """
    if not HAS_H5PY:
        raise ImportError("h5py is required for writing HDF5 output")

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, 'w') as f:
        # Global attributes
        f.attrs['title'] = 'FY-3D MERSI-II Cloud Amount Product'
        f.attrs['institution'] = 'National Satellite Meteorological Center, CMA'
        f.attrs['source'] = 'FYLAT Cloud Mask Retrieval System V3.2 (Python)'
        f.attrs['history'] = f'Created {datetime.now().isoformat()}'
        f.attrs['Conventions'] = 'CF-1.8'
        f.attrs['grid_resolution'] = '5km'

        if attributes:
            for key, value in attributes.items():
                f.attrs[key] = value

        # Cloud Amount dataset
        ds_ca = f.create_dataset(
            'Cloud_Amount',
            data=cloud_amount,
            dtype=np.uint8,
            compression='gzip',
            compression_opts=4,
        )
        ds_ca.attrs['long_name'] = 'Cloud Amount'
        ds_ca.attrs['units'] = 'percent'
        ds_ca.attrs['_FillValue'] = 255
        ds_ca.attrs['valid_range'] = [0, 100]

        # Quality dataset
        ds_qa = f.create_dataset(
            'Cloud_Amount_QA',
            data=cloud_amount_qa,
            dtype=np.uint8,
            compression='gzip',
            compression_opts=4,
        )
        ds_qa.attrs['long_name'] = 'Cloud Amount Quality Flag'
        ds_qa.attrs['units'] = 'category'
        ds_qa.attrs['flag_values'] = [0, 1, 2]
        ds_qa.attrs['flag_meanings'] = 'bad_quality low_quality high_quality'

        # Geolocation
        ds_lon = f.create_dataset(
            'Longitude',
            data=lon_5km,
            dtype=np.float64,
        )
        ds_lon.attrs['long_name'] = 'Longitude'
        ds_lon.attrs['units'] = 'degrees_east'

        ds_lat = f.create_dataset(
            'Latitude',
            data=lat_5km,
            dtype=np.float64,
        )
        ds_lat.attrs['long_name'] = 'Latitude'
        ds_lat.attrs['units'] = 'degrees_north'


def write_combined_product(
    filepath: str,
    cm_bitarray: np.ndarray,
    cm_qa_bitarray: np.ndarray,
    cm_tmp: np.ndarray,
    confidence: np.ndarray,
    cloud_amount: np.ndarray,
    cloud_amount_qa: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    lon_5km: np.ndarray,
    lat_5km: np.ndarray,
    attributes: Optional[dict] = None,
) -> None:
    """Write combined CLM+CLA product as HDF5.

    Args:
        filepath: Output file path.
        cm_bitarray: (n_elem, n_line, 6) uint8 testbits array.
        cm_qa_bitarray: (n_elem, n_line, 10) uint8 QA bits array.
        cm_tmp: (n_elem, n_line) int32 cloud mask values.
        confidence: (n_elem, n_line) float64 confidence values.
        cloud_amount: (n_elem_5km, n_line_5km) uint8 cloud cover.
        cloud_amount_qa: (n_elem_5km, n_line_5km) uint8 quality flag.
        lon: (n_elem, n_line) float64 longitude array.
        lat: (n_elem, n_line) float64 latitude array.
        lon_5km: (n_elem_5km, n_line_5km) float64 longitude array.
        lat_5km: (n_elem_5km, n_line_5km) float64 latitude array.
        attributes: Optional global attributes dictionary.
    """
    if not HAS_H5PY:
        raise ImportError("h5py is required for writing HDF5 output")

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(filepath, 'w') as f:
        # Global attributes
        f.attrs['title'] = 'FY-3D MERSI-II Cloud Mask and Amount Product'
        f.attrs['institution'] = 'National Satellite Meteorological Center, CMA'
        f.attrs['source'] = 'FYLAT Cloud Mask Retrieval System V3.2 (Python)'
        f.attrs['history'] = f'Created {datetime.now().isoformat()}'
        f.attrs['Conventions'] = 'CF-1.8'

        if attributes:
            for key, value in attributes.items():
                f.attrs[key] = value

        # 1km Cloud Mask group
        grp_1km = f.create_group('Cloud_Mask_1km')
        grp_1km.attrs['grid_resolution'] = '1km'

        grp_1km.create_dataset('Cloud_Mask', data=cm_bitarray, dtype=np.uint8,
                               compression='gzip', compression_opts=4)
        grp_1km.create_dataset('Quality_Assurance', data=cm_qa_bitarray, dtype=np.uint8,
                               compression='gzip', compression_opts=4)
        grp_1km.create_dataset('Cloud_Mask_Value', data=cm_tmp, dtype=np.int32,
                               compression='gzip', compression_opts=4)
        grp_1km.create_dataset('Confidence', data=confidence, dtype=np.float64,
                               compression='gzip', compression_opts=4)
        grp_1km.create_dataset('Longitude', data=lon, dtype=np.float64)
        grp_1km.create_dataset('Latitude', data=lat, dtype=np.float64)

        # 5km Cloud Amount group
        grp_5km = f.create_group('Cloud_Amount_5km')
        grp_5km.attrs['grid_resolution'] = '5km'

        grp_5km.create_dataset('Cloud_Amount', data=cloud_amount, dtype=np.uint8,
                               compression='gzip', compression_opts=4)
        grp_5km.create_dataset('Cloud_Amount_QA', data=cloud_amount_qa, dtype=np.uint8,
                               compression='gzip', compression_opts=4)
        grp_5km.create_dataset('Longitude', data=lon_5km, dtype=np.float64)
        grp_5km.create_dataset('Latitude', data=lat_5km, dtype=np.float64)
