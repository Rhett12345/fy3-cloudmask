"""Surface type classification for each pixel.

Port of get_pxldat() in fylat_fy3mersi_cloud_mask.f90.
Determines surface type flags (land, water, coast, desert, snow, ice, etc.)
based on land/sea mask, ecosystem type, geographic location, and NDSI snow detection.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from numba import njit

from ..constants import (
    BAD_DATA, SZA_NIGHT, POLAR_LAT, ANTARCTIC_LAT,
    GLINT_ANGLE_MAX, HIGH_ELEVATION, GREENLAND_ELEVATION,
    DESERT_ECOSYSTEM_TYPES, VRAT_DISABLED_ECOSYSTEMS,
    BAND_055, BAND_064, BAND_086, BAND_138, BAND_213,
    BAND_103, IR_38, IR_85, IR_11, IR_12, BAND_NDSI_NIR,
    VIS_VALID_MIN, VIS_VALID_MAX, IR_VALID_MIN, IR_VALID_MAX,
)


@dataclass
class PixelFlags:
    """All surface type and condition flags for one pixel."""
    polar: bool = False
    land: bool = False
    water: bool = False
    coast: bool = False
    desert: bool = False
    day: bool = False
    night: bool = False
    snow: bool = False
    ice: bool = False
    snglnt: bool = False
    visusd: bool = True
    vrused: bool = True
    hi_elev: bool = False
    antarctic: bool = False
    uniform: bool = True
    bad_value: bool = False
    bad_geo: bool = False
    process: bool = True
    map_ice: bool = False
    map_snow: bool = False
    ndsi_snow: bool = False
    sh_ocean: bool = False
    sh_lake: bool = False
    sg_bad_data: bool = False
    cirrus_ir: bool = False
    cirrus_vis: bool = False
    shadow: bool = False
    smoke: bool = False
    no_250: bool = False
    New_Zealand: bool = False
    Greenland: bool = False


def is_desert_ecosystem(eco_type: int) -> bool:
    """Check if ecosystem type indicates desert."""
    return eco_type in DESERT_ECOSYSTEM_TYPES


def is_vrat_disabled(eco_type: int) -> bool:
    """Check if VRAT test should be disabled for this ecosystem type."""
    return eco_type in VRAT_DISABLED_ECOSYSTEMS


def is_greenland(lat: float, lon: float, eco_type: int) -> bool:
    """Check if pixel is in Greenland region."""
    # Approximate Greenland bounds
    if lat > 60.0 and lat < 84.0:
        if lon > -75.0 and lon < -10.0:
            return True
    # Also check ecosystem type for Greenland
    return False


def is_new_zealand(lat: float, lon: float) -> bool:
    """Check if pixel is in New Zealand region."""
    if lat > -48.0 and lat < -34.0:
        if lon > 165.0 and lon < 180.0:
            return True
    return False


def classify_pixel_surface(
    pxldat: np.ndarray,
    lat: float,
    lon: float,
    elevation: float,
    lsf: int,
    sza: float,
    vza: float,
    glint_angle: float,
    eco_type: int,
    snow_mask_val: int,
    sst: float,
    nwp_sfctmp: float,
    nwp_pmsl: float,
    nwp_u_wind: float,
    nwp_v_wind: float,
    nwp_precip_water: float,
    sensor_id: int,
    thresholds: dict,
) -> tuple[PixelFlags, np.ndarray]:
    """Classify surface type and conditions for a pixel.

    Port of get_pxldat() in fylat_fy3mersi_cloud_mask.f90.

    Args:
        pxldat: 25-element array of pixel data (19 VIS reflectances + 6 IR BTs).
        lat: Latitude (degrees).
        lon: Longitude (degrees).
        elevation: Surface elevation (m).
        lsf: Land/sea flag from GEO data.
        sza: Solar zenith angle (degrees).
        vza: Satellite viewing angle (degrees).
        glint_angle: Sun glint angle (degrees).
        eco_type: IGBP ecosystem type.
        snow_mask_val: Snow/ice mask value from ancillary data.
        sst: Sea surface temperature (K).
        nwp_sfctmp: NWP surface temperature (K).
        nwp_pmsl: NWP mean sea level pressure (hPa).
        nwp_u_wind: NWP U-wind component (m/s).
        nwp_v_wind: NWP V-wind component (m/s).
        nwp_precip_water: NWP precipitable water (mm).
        sensor_id: Sensor ID (21=FY-3D, 22=FY-3E).
        thresholds: Threshold dictionary.

    Returns:
        Tuple of (PixelFlags, modified pxldat array).
    """
    flags = PixelFlags()
    pxldat = pxldat.copy()

    # Note: Band swap removed - using correct band indices directly
    # Band 5 (1.38um) is at index 4, Band 19 (1.03um) is at index 18

    # Initialize NWP-derived values
    sfctmp = nwp_sfctmp if nwp_sfctmp > 0 else BAD_DATA
    pmsl = nwp_pmsl if nwp_pmsl > 0 else BAD_DATA

    # SST override for water/coast pixels
    if (lsf == 0 or lsf == 2 or lsf == 3) and sst > 100.0:
        sfctmp = sst

    # Elevation adjustment for desert 11um threshold
    tbadj = (elevation / 1000.0) * 5.0 if elevation > 0 else 0.0

    # Geographic classification
    flags.antarctic = lat < ANTARCTIC_LAT

    # Day/night classification
    if sza > SZA_NIGHT:
        flags.night = True
        flags.day = False
        flags.visusd = False
    else:
        flags.day = True
        flags.night = False
        flags.visusd = True

    # Polar classification
    flags.polar = abs(lat) > POLAR_LAT

    # Sun glint
    flags.snglnt = glint_angle <= GLINT_ANGLE_MAX

    # High elevation
    flags.hi_elev = elevation > HIGH_ELEVATION
    flags.Greenland = is_greenland(lat, lon, eco_type)
    if flags.Greenland and elevation > GREENLAND_ELEVATION:
        flags.hi_elev = True
    # Ellesmere Island special case
    if lat > 76.0 and lat < 84.0 and lon > -95.0 and lon < -60.0:
        if elevation > GREENLAND_ELEVATION:
            flags.hi_elev = True

    # Desert classification
    flags.desert = is_desert_ecosystem(eco_type)
    # Additional regional desert rules
    if not flags.desert:
        # Africa
        if lat > -35.0 and lat < 37.0 and lon > -20.0 and lon < 55.0:
            if eco_type == 42 and elevation > 1000.0:
                flags.desert = True
        # Eurasia
        elif lat > 15.0 and lat < 55.0 and lon > 40.0 and lon < 140.0:
            if eco_type == 42 and elevation > 1000.0:
                flags.desert = True
        # Australia
        elif lat > -40.0 and lat < -10.0 and lon > 112.0 and lon < 155.0:
            if eco_type == 42:
                flags.desert = True

    # VRAT disabled for certain ecosystems
    flags.vrused = not is_vrat_disabled(eco_type)

    # Land/sea classification from land-sea mask
    if lsf == 1 or lsf == 4:
        flags.land = True
        flags.water = False
        flags.coast = False
        flags.sh_ocean = False
    elif lsf == 2:
        flags.land = True
        flags.coast = True
        flags.water = False
    elif lsf == 3:
        flags.land = True
        flags.sh_lake = True
        if flags.day:
            flags.coast = True
        else:
            flags.coast = False
        flags.water = False
    elif lsf == 0:
        flags.water = True
        flags.land = False
        flags.coast = False
        flags.sh_ocean = True
    else:
        # Fallback: use VRAT to determine land/water
        flags.water = True
        flags.land = False
        if pxldat[BAND_064] > VIS_VALID_MIN and pxldat[BAND_086] > VIS_VALID_MIN:
            if pxldat[BAND_064] > 0.001:
                vrat = pxldat[BAND_086] / pxldat[BAND_064]
                if vrat > 0.9:
                    flags.land = True
                    flags.water = False

    # New Zealand check
    flags.New_Zealand = is_new_zealand(lat, lon)

    # Snow/ice from ancillary map
    if flags.water:
        if 25 <= snow_mask_val <= 100:
            flags.map_ice = True
        elif snow_mask_val in (101, 103, 104):
            flags.map_snow = True
    else:
        if 25 <= snow_mask_val <= 100:
            flags.map_snow = True
        elif snow_mask_val in (101, 103, 104):
            flags.map_snow = True
        elif snow_mask_val == 200 and abs(lat) > POLAR_LAT:
            flags.map_snow = True

    # Daytime NDSI snow detection
    if flags.day:
        ndsi_snow = detect_snow_ndsi(
            pxldat, flags, thresholds, flags.Greenland, flags.hi_elev, lat,
        )
        flags.ndsi_snow = ndsi_snow

        # Combine map and NDSI snow/ice
        if flags.water:
            flags.ice = flags.map_ice
            if flags.ndsi_snow and not flags.map_ice:
                flags.snow = True
            else:
                flags.snow = flags.map_snow
        else:
            flags.snow = flags.map_snow or flags.ndsi_snow
            flags.ice = False
    else:
        # Nighttime: use map data directly
        flags.ice = flags.map_ice
        flags.snow = flags.map_snow

    # Check for bad data values
    if pxldat[IR_11] < IR_VALID_MIN or pxldat[IR_11] > IR_VALID_MAX:
        flags.bad_value = True

    return flags, pxldat


def detect_snow_ndsi(
    pxldat: np.ndarray,
    flags: PixelFlags,
    thresholds: dict,
    is_greenland: bool,
    is_hi_elev: bool,
    lat: float,
) -> bool:
    """NDSI-based snow detection with false snow removal.

    Port of snow_mask() in fylat_fy3mersi_cloud_mask.f90.

    Args:
        pxldat: 25-element pixel data array.
        flags: Pixel flags.
        thresholds: Snow mask thresholds.
        is_greenland: Whether pixel is in Greenland.
        is_hi_elev: Whether pixel is at high elevation.
        lat: Latitude.

    Returns:
        True if snow detected.
    """
    sm = thresholds.get('snow_mask', {})

    # Gate: 11um BT must be below threshold
    bt11_thresh = sm.get('bt11_threshold', 280.0)
    if pxldat[IR_11] > bt11_thresh:
        return False

    # Compute NDSI
    masv55 = pxldat[BAND_055]
    masnir = pxldat[BAND_NDSI_NIR]
    if masv55 + masnir <= 0:
        return False
    ndsi = (masv55 - masnir) / (masv55 + masnir)

    # Snow detection criteria
    ndsi_thresh = sm.get('ndsi_threshold', 0.4)
    ref086_thresh = sm.get('ref086_threshold', 0.11)

    if ndsi <= ndsi_thresh:
        return False
    if pxldat[BAND_086] <= ref086_thresh:
        return False

    # --- False snow removal filters ---

    # 1. Thin cirrus filter (not Greenland)
    ref138_thresh = sm.get('ref138_threshold', 0.0525)
    if not is_greenland:
        if pxldat[BAND_138] > ref138_thresh:
            return False

    # 2. Sun glint filter
    if flags.water and flags.snglnt:
        if -60.0 < lat < 50.0:
            return False

    # 3. Ice cloud filter
    btd_85_11_thresh = sm.get('btd_85_11_threshold', 0.5)
    gl_offset = 1.5 if is_greenland else 0.0
    if pxldat[IR_85] - pxldat[IR_11] >= btd_85_11_thresh + gl_offset:
        return False

    # 4. Water cloud filter (3.8um available)
    btd_37_11_thresh = sm.get('btd_37_11_threshold', 9.0)
    btd_37_11_hel = sm.get('btd_37_11_hel_threshold', 10.0)
    if not is_greenland:
        thresh = btd_37_11_hel if is_hi_elev else btd_37_11_thresh
        if pxldat[IR_38] - pxldat[IR_11] >= thresh:
            return False

    # 5. NIR brightness filter (not high elevation, not Greenland)
    nir_thresh = sm.get('nir_threshold', 0.17)
    if not is_hi_elev and not is_greenland:
        if masnir > nir_thresh:
            return False

    return True
