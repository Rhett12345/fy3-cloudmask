"""End-to-end pipeline for cloud mask processing.

Orchestrates the complete cloud mask workflow:
1. Read satellite data (GEO + L1b)
2. Read and interpolate NWP data
3. Read ancillary data
4. Run cloud mask algorithm
5. Compute cloud amount
6. Write outputs
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np

from .config import FY3Config, load_config
from .algorithm import run_cloud_mask_swath, CloudMaskResult
from .algorithm.native_backend import is_native_available, process_swath_native, get_backend_info
from .output import (
    compute_cloud_amount, compute_cloud_amount_with_coords,
    write_cloud_mask, write_cloud_amount, write_combined_product,
)

logger = logging.getLogger(__name__)


@dataclass
class ProcessingResult:
    """Result of processing a single orbit."""
    success: bool
    l1b_path: str
    geo_path: str
    output_path: str
    processing_time: float  # seconds
    n_pixels_processed: int
    n_cloudy: int
    n_clear: int
    error_message: Optional[str] = None


class CloudMaskPipeline:
    """Cloud mask processing pipeline.

    This class orchestrates the complete cloud mask workflow for
    processing FY-3D MERSI-II satellite data.

    Example:
        >>> pipeline = CloudMaskPipeline('config/default.yaml')
        >>> result = pipeline.process_orbit(
        ...     l1b_path='data/FY3D_MERSI_20230101_0000_1000M_MS.HDF',
        ...     geo_path='data/FY3D_MERSI_GEO_20230101_0000_1000M.HDF',
        ...     output_dir='output/',
        ... )
    """

    def __init__(self, config_path: str):
        """Initialize pipeline with configuration.

        Args:
            config_path: Path to YAML configuration file.
        """
        self.config = load_config(config_path)
        self.thresholds = self._load_thresholds()
        self._backend_info = get_backend_info()
        logger.info(f"Pipeline initialized with config: {config_path}")
        logger.info(f"Backend: {self._backend_info['backend']}")

    def _load_thresholds(self) -> dict:
        """Load threshold configuration.

        Returns:
            Threshold dictionary.
        """
        import yaml
        thresholds_path = Path(self.config.paths.coeff_dir) / 'thresholds' / 'mersi_ii3d_v8.yaml'
        if thresholds_path.exists():
            with open(thresholds_path) as f:
                return yaml.safe_load(f)
        else:
            logger.warning(f"Threshold file not found: {thresholds_path}")
            return {}

    def process_orbit(
        self,
        l1b_path: str,
        geo_path: str,
        output_dir: str,
        nwp_path1: Optional[str] = None,
        nwp_path2: Optional[str] = None,
        oisst_path: Optional[str] = None,
    ) -> ProcessingResult:
        """Process a single orbit.

        Args:
            l1b_path: Path to L1b HDF5 file.
            geo_path: Path to GEO HDF5 file.
            output_dir: Output directory.
            nwp_path1: Path to first NWP GRIB file (optional).
            nwp_path2: Path to second NWP GRIB file (optional).
            oisst_path: Path to OISST file (optional).

        Returns:
            ProcessingResult with processing statistics.
        """
        import time
        start_time = time.time()

        logger.info(f"Processing orbit: {l1b_path}")
        logger.info(f"  GEO: {geo_path}")
        logger.info(f"  Output: {output_dir}")

        try:
            # Step 1: Read satellite data
            logger.info("Step 1: Reading satellite data...")
            sat_data = self._read_satellite_data(l1b_path, geo_path)

            # Step 2: Read NWP data (if available)
            logger.info("Step 2: Reading NWP data...")
            nwp_data = self._read_nwp_data(nwp_path1, nwp_path2)

            # Step 3: Read ancillary data
            logger.info("Step 3: Reading ancillary data...")
            ancillary = self._read_ancillary_data(oisst_path)

            # Step 4: Run cloud mask algorithm
            logger.info("Step 4: Running cloud mask algorithm...")
            if is_native_available():
                logger.info("  Using native C++/Fortran backend (OpenMP)")
                cm_result = self._run_native_backend(sat_data, nwp_data)
                cm_bitarray = cm_result['cm_bitarray']
                cm_qa_bitarray = cm_result['qa_bitarray']
                cm_tmp = cm_result['cloud_mask']
                confidence = cm_result['confidence']
            else:
                logger.info("  Using Python/Numba backend")
                cm_bitarray, cm_qa_bitarray, cm_tmp, confidence = run_cloud_mask_swath(
                    pxldat_swath=sat_data['pxldat'],
                    lat_swath=sat_data['lat'],
                    lon_swath=sat_data['lon'],
                    elevation_swath=sat_data['elevation'],
                    lsf_swath=sat_data['lsf'],
                    sza_swath=sat_data['sza'],
                    vza_swath=sat_data['vza'],
                    glint_angle_swath=sat_data['glint_angle'],
                    eco_type_swath=sat_data['eco_type'],
                    snow_mask_swath=sat_data['snow_mask'],
                    sst_swath=sat_data.get('sst', np.zeros_like(sat_data['lat'])),
                    nwp_sfctmp_swath=nwp_data.get('sfctmp', np.zeros_like(sat_data['lat'])),
                    nwp_pmsl_swath=nwp_data.get('pmsl', np.zeros_like(sat_data['lat'])),
                    nwp_u_wind_swath=nwp_data.get('u_wind', np.zeros_like(sat_data['lat'])),
                    nwp_v_wind_swath=nwp_data.get('v_wind', np.zeros_like(sat_data['lat'])),
                    nwp_precip_water_swath=nwp_data.get('precip_water', np.zeros_like(sat_data['lat'])),
                    bt_clr_swath=nwp_data.get('bt_clr', np.zeros((*sat_data['lat'].shape, 7))),
                    sensor_id=self.config.sensor.sensor_id,
                    thresholds=self.thresholds,
                )

            # Step 5: Compute cloud amount
            logger.info("Step 5: Computing cloud amount...")
            qa_processed = (cm_tmp >= 0).astype(np.int32)
            cloud_amount, cloud_amount_qa = compute_cloud_amount(cm_tmp, qa_processed)

            # Step 6: Write outputs
            logger.info("Step 6: Writing outputs...")
            output_path = self._write_outputs(
                output_dir, l1b_path,
                cm_bitarray, cm_qa_bitarray, cm_tmp, confidence,
                cloud_amount, cloud_amount_qa,
                sat_data['lon'], sat_data['lat'],
            )

            # Compute statistics
            n_pixels = cm_tmp.size
            n_cloudy = np.sum(cm_tmp < 2)
            n_clear = np.sum(cm_tmp >= 2)

            elapsed = time.time() - start_time
            logger.info(f"Processing complete in {elapsed:.1f}s")
            logger.info(f"  Total pixels: {n_pixels}")
            logger.info(f"  Cloudy: {n_cloudy} ({100*n_cloudy/n_pixels:.1f}%)")
            logger.info(f"  Clear: {n_clear} ({100*n_clear/n_pixels:.1f}%)")

            return ProcessingResult(
                success=True,
                l1b_path=l1b_path,
                geo_path=geo_path,
                output_path=output_path,
                processing_time=elapsed,
                n_pixels_processed=n_pixels,
                n_cloudy=n_cloudy,
                n_clear=n_clear,
            )

        except Exception as e:
            elapsed = time.time() - start_time
            logger.error(f"Processing failed: {e}")
            return ProcessingResult(
                success=False,
                l1b_path=l1b_path,
                geo_path=geo_path,
                output_path='',
                processing_time=elapsed,
                n_pixels_processed=0,
                n_cloudy=0,
                n_clear=0,
                error_message=str(e),
            )

    def _run_native_backend(self, sat_data: dict, nwp_data: dict) -> dict:
        """Run cloud mask using the native C++/Fortran backend.

        Args:
            sat_data: Satellite data dictionary.
            nwp_data: NWP data dictionary.

        Returns:
            Dictionary with cloud mask results.
        """
        n_elem, n_line = sat_data['lat'].shape

        # Extract VIS and IR bands from pxldat
        pxldat = sat_data['pxldat']
        if pxldat.ndim == 3 and pxldat.shape[2] == 25:
            ref_vis = pxldat[:, :, :19].astype(np.float32)
            tbb_ir = pxldat[:, :, 19:25].astype(np.float32)
        else:
            raise ValueError(f"Unexpected pxldat shape: {pxldat.shape}")

        # Prepare geometry arrays
        lat = sat_data['lat'].astype(np.float32)
        lon = sat_data['lon'].astype(np.float32)
        satzen = sat_data['vza'].astype(np.float32)
        solzen = sat_data['sza'].astype(np.float32)
        relaz = np.zeros_like(lat, dtype=np.float32)  # Placeholder
        glint = sat_data['glint_angle'].astype(np.float32)

        # NWP arrays
        sfctmp = nwp_data.get('sfctmp', np.full_like(lat, 300.0)).astype(np.float32)
        pmsl = nwp_data.get('pmsl', np.full_like(lat, 1013.0)).astype(np.float32)
        uwind = nwp_data.get('u_wind', np.zeros_like(lat)).astype(np.float32)
        vwind = nwp_data.get('v_wind', np.zeros_like(lat)).astype(np.float32)
        tpw = nwp_data.get('precip_water', np.full_like(lat, 3.0)).astype(np.float32)
        btclr = nwp_data.get('bt_clr', np.full((*lat.shape, 7), 280.0, dtype=np.float32))

        # Ancillary arrays
        elev = sat_data['elevation'].astype(np.float32)
        eco = sat_data['eco_type'].astype(np.int8)
        lsf = sat_data['lsf'].astype(np.int8)
        snow_mask = sat_data['snow_mask'].astype(np.int8)

        # Call native engine
        return process_swath_native(
            ref_vis, tbb_ir, lat, lon, satzen, solzen, relaz, glint,
            sfctmp, pmsl, uwind, vwind, tpw, elev, eco, lsf, snow_mask, btclr,
            n_elem, n_line
        )

    def _read_satellite_data(self, l1b_path: str, geo_path: str) -> dict:
        """Read satellite data from HDF5 files.

        Args:
            l1b_path: L1b file path.
            geo_path: GEO file path.

        Returns:
            Dictionary with satellite data arrays.
        """
        # Placeholder - actual implementation depends on data format
        # This will be implemented in Phase 2 (Data I/O)
        logger.warning("Using placeholder satellite data reader")
        n_elem = self.config.sensor.n_elem
        n_line = self.config.sensor.n_line

        return {
            'pxldat': np.zeros((n_elem, n_line, 25), dtype=np.float64),
            'lat': np.zeros((n_elem, n_line), dtype=np.float64),
            'lon': np.zeros((n_elem, n_line), dtype=np.float64),
            'elevation': np.zeros((n_elem, n_line), dtype=np.float64),
            'lsf': np.ones((n_elem, n_line), dtype=np.int32),
            'sza': np.full((n_elem, n_line), 30.0, dtype=np.float64),
            'vza': np.full((n_elem, n_line), 10.0, dtype=np.float64),
            'glint_angle': np.full((n_elem, n_line), 60.0, dtype=np.float64),
            'eco_type': np.ones((n_elem, n_line), dtype=np.int32),
            'snow_mask': np.zeros((n_elem, n_line), dtype=np.int32),
        }

    def _read_nwp_data(self, path1: Optional[str], path2: Optional[str]) -> dict:
        """Read NWP data from GRIB files.

        Args:
            path1: First NWP file path.
            path2: Second NWP file path.

        Returns:
            Dictionary with NWP data arrays.
        """
        # Placeholder - actual implementation depends on NWP format
        logger.warning("Using placeholder NWP data reader")
        return {}

    def _read_ancillary_data(self, oisst_path: Optional[str]) -> dict:
        """Read ancillary data.

        Args:
            oisst_path: OISST file path.

        Returns:
            Dictionary with ancillary data.
        """
        # Placeholder - actual implementation depends on ancillary data format
        logger.warning("Using placeholder ancillary data reader")
        return {}

    def _write_outputs(
        self,
        output_dir: str,
        l1b_path: str,
        cm_bitarray: np.ndarray,
        cm_qa_bitarray: np.ndarray,
        cm_tmp: np.ndarray,
        confidence: np.ndarray,
        cloud_amount: np.ndarray,
        cloud_amount_qa: np.ndarray,
        lon: np.ndarray,
        lat: np.ndarray,
    ) -> str:
        """Write output files.

        Args:
            output_dir: Output directory.
            l1b_path: Input L1b path (for naming).
            cm_bitarray: Cloud mask testbits.
            cm_qa_bitarray: QA bits.
            cm_tmp: Cloud mask values.
            confidence: Confidence values.
            cloud_amount: Cloud amount.
            cloud_amount_qa: Cloud amount QA.
            lon: Longitude array.
            lat: Latitude array.

        Returns:
            Output file path.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate output filename from input
        l1b_name = Path(l1b_path).stem
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        # Write combined product
        output_path = output_dir / f'{l1b_name}_CLM_CLA_{timestamp}.h5'

        # Compute 5km geolocation
        n_elem, n_line = lon.shape
        box_size = 5
        n_elem_5km = n_elem // box_size
        n_line_5km = n_line // box_size
        lon_5km = lon[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]
        lat_5km = lat[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]

        write_combined_product(
            str(output_path),
            cm_bitarray, cm_qa_bitarray, cm_tmp, confidence,
            cloud_amount, cloud_amount_qa,
            lon, lat, lon_5km, lat_5km,
            attributes={
                'input_l1b': l1b_path,
                'processing_time': datetime.now().isoformat(),
            },
        )

        logger.info(f"Output written to: {output_path}")
        return str(output_path)

    def process_batch(
        self,
        start_date: str,
        end_date: str,
        data_root: str,
        output_root: str,
        n_workers: int = 1,
    ) -> list[ProcessingResult]:
        """Process multiple orbits.

        Args:
            start_date: Start date (YYYY-MM-DD).
            end_date: End date (YYYY-MM-DD).
            data_root: Root directory for input data.
            output_root: Root directory for output.
            n_workers: Number of parallel workers.

        Returns:
            List of ProcessingResult for each orbit.
        """
        from datetime import datetime, timedelta

        start = datetime.strptime(start_date, '%Y-%m-%d')
        end = datetime.strptime(end_date, '%Y-%m-%d')

        results = []
        current = start

        while current <= end:
            # Find orbits for this date
            date_str = current.strftime('%Y%m%d')
            orbits = self._find_orbits(data_root, date_str)

            for l1b_path, geo_path in orbits:
                result = self.process_orbit(
                    l1b_path, geo_path,
                    output_dir=str(Path(output_root) / date_str),
                )
                results.append(result)

            current += timedelta(days=1)

        # Summary
        n_success = sum(1 for r in results if r.success)
        n_failed = sum(1 for r in results if not r.success)
        logger.info(f"Batch processing complete: {n_success} succeeded, {n_failed} failed")

        return results

    def _find_orbits(self, data_root: str, date_str: str) -> list[tuple[str, str]]:
        """Find orbit files for a given date.

        Args:
            data_root: Root directory.
            date_str: Date string (YYYYMMDD).

        Returns:
            List of (l1b_path, geo_path) tuples.
        """
        data_path = Path(data_root)
        orbits = []

        # Look for L1b files matching date pattern
        pattern = f'*{date_str}*1000M_MS.HDF'
        for l1b_path in data_path.rglob(pattern):
            # Find matching GEO file
            geo_name = l1b_path.name.replace('1000M_MS', 'GEO')
            geo_path = l1b_path.parent / geo_name
            if geo_path.exists():
                orbits.append((str(l1b_path), str(geo_path)))

        return orbits
