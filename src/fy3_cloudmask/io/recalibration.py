"""Recalibration coefficient manager for FY-3D MERSI-II.

Loads daily recalibration coefficients for 7 solar reflectance bands from
CSV files organized as: {base_dir}/YYYYMM/RAD_YYYYMMDD.csv

Formula: reflectance = cal0 + cal1 * DN
(cal2 is always 0; coefficients already include 0.01/esd² normalization)
"""

from __future__ import annotations

import csv
import logging
import os
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)

# Number of solar reflectance bands that use recalibration
N_SOLAR_BANDS = 7


class RecalibrationManager:
    """Manages daily recalibration coefficients for 7 solar reflectance bands.

    Parameters
    ----------
    base_dir : str
        Path to the recalibration data directory (e.g., '../fy3d_recali').
        Expected structure: {base_dir}/YYYYMM/RAD_YYYYMMDD.csv
    """

    def __init__(self, base_dir: str):
        self.base_dir = Path(base_dir).resolve()
        if not self.base_dir.is_dir():
            raise FileNotFoundError(f"Recalibration directory not found: {self.base_dir}")

    def _csv_path(self, date_str: str) -> Path:
        """Get the CSV file path for a given date (YYYYMMDD)."""
        yyyymm = date_str[:6]
        return self.base_dir / yyyymm / f"RAD_{date_str}.csv"

    def load_coefficients(self, date_str: str) -> tuple[np.ndarray, np.ndarray]:
        """Load recalibration coefficients for a specific date.

        Parameters
        ----------
        date_str : str
            Date in YYYYMMDD format.

        Returns
        -------
        cal0 : ndarray, shape (7,), float64
            Intercept coefficients for channels 1-7.
        cal1 : ndarray, shape (7,), float64
            Slope coefficients for channels 1-7.
        """
        csv_path = self._csv_path(date_str)
        if not csv_path.exists():
            # Try nearest available date
            csv_path = self._find_nearest(date_str)
            logger.warning("Exact date %s not found, using nearest: %s", date_str, csv_path.stem)

        return self._parse_csv(csv_path)

    def _parse_csv(self, csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
        """Parse a recalibration CSV file.

        Expected format:
            ,cal0,cal1,cal2
            ch01,-3.264,0.0273,0.0
            ch02,-4.324,0.0259,0.0
            ...
            ch07,-2.602,0.0208,0.0

        Note: Coefficients are in percentage scale (reflectance 0-100).
        We divide by 100 to convert to decimal scale (0-1) matching onboard calibration.
        """
        cal0 = np.zeros(N_SOLAR_BANDS, dtype=np.float64)
        cal1 = np.zeros(N_SOLAR_BANDS, dtype=np.float64)

        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)  # skip header
            for i, row in enumerate(reader):
                if i >= N_SOLAR_BANDS:
                    break
                # row[0] = channel name (e.g., 'ch01'), row[1] = cal0, row[2] = cal1
                # Divide by 100 to convert from percentage to decimal
                cal0[i] = float(row[1]) / 100.0
                cal1[i] = float(row[2]) / 100.0

        logger.debug("Loaded recalibration from %s: cal0=%s, cal1=%s",
                      csv_path.name, cal0[:3], cal1[:3])
        return cal0, cal1

    def _find_nearest(self, date_str: str) -> Path:
        """Find the nearest available recalibration file to the requested date.

        Searches forward and backward from the target date, up to 30 days.
        """
        target = datetime.strptime(date_str, '%Y%m%d')
        yyyymm = date_str[:6]

        # Check same month first
        month_dir = self.base_dir / yyyymm
        if month_dir.is_dir():
            # Look for closest file in same month
            candidates = []
            for f in month_dir.glob('RAD_*.csv'):
                fdate = f.stem.replace('RAD_', '')
                try:
                    fd = datetime.strptime(fdate, '%Y%m%d')
                    candidates.append((abs((fd - target).days), f))
                except ValueError:
                    continue
            if candidates:
                candidates.sort(key=lambda x: x[0])
                return candidates[0][1]

        # Search nearby months
        for delta in range(1, 31):
            for sign in [1, -1]:
                check_date = target + timedelta(days=sign * delta)
                check_path = self._csv_path(check_date.strftime('%Y%m%d'))
                if check_path.exists():
                    return check_path

        raise FileNotFoundError(
            f"No recalibration file found within 30 days of {date_str} "
            f"(searched in {self.base_dir})"
        )

    def apply_to_pxldat(
        self,
        pxldat: np.ndarray,
        dn_vis: np.ndarray,
        date_str: str,
    ) -> None:
        """Apply recalibration to the first 7 bands of pxldat in-place.

        Parameters
        ----------
        pxldat : ndarray, shape (nElem, nLine, 25)
            Full pixel data array. Bands 0-6 will be overwritten.
        dn_vis : ndarray, shape (7, nElem, nLine) or (nElem, nLine, 7)
            Raw DN values for the 7 solar reflectance bands.
            If 3D with shape (7, nElem, nLine), will be transposed.
        date_str : str
            Date in YYYYMMDD format for coefficient lookup.
        """
        cal0, cal1 = self.load_coefficients(date_str)

        # Handle both input shapes: (7, nElem, nLine) from HDF or (nElem, nLine, 7)
        if dn_vis.ndim == 3 and dn_vis.shape[0] == N_SOLAR_BANDS:
            dn_vis = dn_vis.transpose(1, 2, 0)  # -> (nElem, nLine, 7)

        for b in range(N_SOLAR_BANDS):
            dn = dn_vis[:, :, b].astype(np.float64)
            refl = cal0[b] + cal1[b] * dn
            pxldat[:, :, b] = np.clip(refl, -0.1, 2.0)
