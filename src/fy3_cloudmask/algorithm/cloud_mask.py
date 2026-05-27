"""Main cloud mask algorithm driver.

Port of fylat_fy3mersi_cloud_mask.f90.

This module orchestrates the cloud mask algorithm:
1. Classify pixel surface type
2. Dispatch to appropriate surface/lighting test
3. Apply restoral tests
4. Assemble final cloud mask bits
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from ..constants import (
    BAD_DATA, SZA_NIGHT, POLAR_LAT, ANTARCTIC_LAT,
    COLD_SURFACE_TEMP, HIGH_ELEVATION,
    BIT_PROCESSED, BIT_CONF_LSB, BIT_CONF_MSB, BIT_DAY,
    BIT_NO_SUNGLINT, BIT_NO_SNOW_ICE, BIT_COAST, BIT_DESERT,
    CONF_CLEAR, CONF_PROB_CLEAR, CONF_CLOUDY,
    CM_CLOUDY, CM_PROB_CLOUDY, CM_PROB_CLEAR, CM_CLEAR,
)
from .confidence import encode_confidence
from .bitops import (
    set_bit, clear_bit, check_bit, init_testbits, init_qa_bits,
    fill_bit_pixel, proc_path, set_unused_bits, convert_cloud_mask,
)
from .surface_classifier import classify_pixel_surface, PixelFlags
from .spatial import check_reg_uniformity, get_regional_std
from .tests import (
    land_day_standard, land_day_coast, land_day_desert, land_day_desert_coast,
    land_nite,
    ocean_day, ocean_nite,
    polar_day_land, polar_day_coast, polar_day_desert, polar_day_desert_coast,
    polar_day_ocean, polar_day_snow,
    polar_nite_land, polar_nite_ocean, polar_nite_snow,
    day_snow, nite_snow, antarctic_day,
    chk_land_restoral, chk_land_nite_restoral, chk_coast_restoral,
    chk_sunglint_restoral, chk_shallow_water, chk_spatial_var,
    chk_cloud_adj, chk_thin_cirrus_ir, chk_shadow,
)


@dataclass
class CloudMaskResult:
    """Result of cloud mask processing for a single pixel."""
    cloud_mask: int = 0          # 0=cloudy, 1=prob cloudy, 2=prob clear, 3=confident clear
    confidence: float = 1.0      # Raw confidence value
    n_tests: int = 0             # Number of tests applied
    n_bands: int = 0             # Number of bands used
    testbits: np.ndarray = None  # 6-byte test bits
    qa_bits: np.ndarray = None   # 10-byte QA bits


def run_cloud_mask_pixel(
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
    bt_clr: np.ndarray,
    thresholds: dict,
    indat_3x3_11um: np.ndarray = None,
    indat_3x3_vis: np.ndarray = None,
    cm_array: np.ndarray = None,
    row: int = 0,
    col: int = 0,
    n_rows: int = 2000,
    n_cols: int = 2048,
) -> CloudMaskResult:
    """Run cloud mask algorithm for a single pixel.

    Port of the main pixel loop in fylat_fy3mersi_cloud_mask.f90.

    Args:
        pxldat: 25-element pixel data array.
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
        bt_clr: 7-element clear-sky BT array (from RTM or NWP).
        thresholds: Threshold dictionary.
        indat_3x3_11um: 3x3 array of 11um BT for spatial tests.
        indat_3x3_vis: 3x3 array of visible reflectance for shadow test.
        cm_array: 2D cloud mask array (for adjacency test).
        row: Current row index.
        col: Current column index.
        n_rows: Total number of rows.
        n_cols: Total number of columns.

    Returns:
        CloudMaskResult with all outputs.
    """
    result = CloudMaskResult()
    result.testbits = init_testbits()
    result.qa_bits = init_qa_bits()

    # ================================================================
    # Step 1: Classify pixel surface type
    # ================================================================
    flags, pxldat = classify_pixel_surface(
        pxldat, lat, lon, elevation, lsf, sza, vza, glint_angle,
        eco_type, snow_mask_val, sst, nwp_sfctmp, nwp_pmsl,
        nwp_u_wind, nwp_v_wind, nwp_precip_water, sensor_id, thresholds,
    )

    # Skip processing for bad data
    if flags.bad_value or flags.bad_geo:
        result.cloud_mask = CM_CLOUDY
        result.confidence = 0.0
        return result

    # ================================================================
    # Step 2: Set processing path bits
    # ================================================================
    proc_path(
        result.testbits, flags.day, flags.snglnt, flags.water,
        flags.snow, flags.ice, flags.coast, flags.desert, flags.land,
        flags.shadow, flags.smoke,
    )

    # ================================================================
    # Step 3: Dispatch to appropriate test function
    # ================================================================
    confdnc, nmtests, nbands = _dispatch_test(
        pxldat, flags, vza, sza, lat, lon, sfctmp=nwp_sfctmp,
        sst=sst, eco_type=eco_type, bt_clr=bt_clr,
        thresholds=thresholds, testbits=result.testbits,
        qa_bits=result.qa_bits, indat_3x3_11um=indat_3x3_11um,
    )

    # ================================================================
    # Step 4: Apply restoral tests
    # ================================================================
    if flags.day:
        # Land restoral
        if flags.land and not (flags.snow or flags.ice):
            if confdnc <= 0.95:
                confdnc = chk_land_restoral(
                    confdnc, pxldat, nwp_sfctmp, nwp_pmsl,
                    nwp_precip_water, vza, thresholds,
                    result.testbits, result.qa_bits,
                )

        # Night land restoral (should not happen in day, but for completeness)
    else:
        # Night land restoral
        if flags.land and not (flags.snow or flags.ice):
            if confdnc <= 0.95:
                confdnc = chk_land_nite_restoral(
                    confdnc, pxldat, nwp_sfctmp, nwp_pmsl,
                    nwp_precip_water, vza, thresholds,
                    result.testbits, result.qa_bits,
                )

    # Coastal restoral
    if flags.coast and not (flags.snow or flags.ice):
        confdnc = chk_coast_restoral(
            confdnc, pxldat, nwp_sfctmp, sst, thresholds,
            result.testbits, result.qa_bits,
        )

    # Sun glint restoral
    if flags.snglnt:
        refang = 0.0  # TODO: compute from geometry
        confdnc = chk_sunglint_restoral(
            confdnc, pxldat, refang, flags.snglnt, thresholds,
            result.testbits, result.qa_bits,
        )

    # Shallow water restoral
    if flags.sh_ocean or flags.sh_lake:
        confdnc = chk_shallow_water(
            confdnc, pxldat, flags.sh_ocean, flags.sh_lake, thresholds,
            result.testbits, result.qa_bits,
        )

    # Spatial variability check
    if indat_3x3_11um is not None and flags.water:
        confdnc = chk_spatial_var(
            confdnc, indat_3x3_11um, flags.uniform, False, thresholds,
            result.testbits, result.qa_bits,
        )

    # Thin cirrus IR check
    confdnc = chk_thin_cirrus_ir(
        confdnc, pxldat, vza, thresholds,
        result.testbits, result.qa_bits,
    )

    # Shadow check
    if indat_3x3_vis is not None and indat_3x3_11um is not None:
        confdnc = chk_shadow(
            confdnc, pxldat, indat_3x3_vis, indat_3x3_11um,
            sza, vza, thresholds, result.testbits, result.qa_bits,
        )

    # Cloud adjacency check
    if cm_array is not None:
        confdnc = chk_cloud_adj(
            confdnc, cm_array, row, col, n_rows, n_cols, thresholds,
            result.testbits, result.qa_bits,
        )

    # ================================================================
    # Step 5: Encode confidence and assemble final bits
    # ================================================================
    result.confidence = confdnc
    result.n_tests = nmtests
    result.n_bands = nbands

    # Encode confidence into testbits
    bit1, bit2 = encode_confidence(confdnc)
    if bit1:
        set_bit(result.testbits, BIT_CONF_LSB)
    if bit2:
        set_bit(result.testbits, BIT_CONF_MSB)

    # Set day bit
    if flags.day:
        set_bit(result.testbits, BIT_DAY)

    # Set surface type bits
    if flags.land:
        set_bit(result.testbits, BIT_COAST)
        set_bit(result.testbits, BIT_DESERT)
    elif flags.coast:
        set_bit(result.testbits, BIT_COAST)
    elif flags.desert:
        set_bit(result.testbits, BIT_DESERT)

    # Fill bit pixel (quality assembly)
    fill_bit_pixel(
        nmtests, nbands, flags.bad_geo, flags.snglnt,
        result.testbits, result.qa_bits,
    )

    # Set unused bits
    set_unused_bits(result.testbits)

    # Convert to scalar cloud mask
    result.cloud_mask, _ = convert_cloud_mask(result.testbits)

    return result


def _dispatch_test(
    pxldat: np.ndarray,
    flags: PixelFlags,
    vza: float,
    sza: float,
    lat: float,
    lon: float,
    sfctmp: float,
    sst: float,
    eco_type: int,
    bt_clr: np.ndarray,
    thresholds: dict,
    testbits: np.ndarray,
    qa_bits: np.ndarray,
    indat_3x3_11um: np.ndarray = None,
    refang: float = 0.0,
) -> tuple[float, int, int]:
    """Dispatch to appropriate test function based on surface type and lighting.

    Args:
        pxldat: 25-element pixel data array.
        flags: Pixel surface type flags.
        vza: Satellite viewing angle (degrees).
        sza: Solar zenith angle (degrees).
        lat: Latitude (degrees).
        lon: Longitude (degrees).
        sfctmp: Surface temperature (K).
        sst: Sea surface temperature (K).
        eco_type: Ecosystem type.
        bt_clr: 7-element clear-sky BT array.
        thresholds: Threshold dictionary.
        testbits: 6-byte array (modified in-place).
        qa_bits: 10-byte array (modified in-place).
        indat_3x3_11um: 3x3 array of 11um BT for spatial tests.
        refang: Reflectance angle.

    Returns:
        Tuple of (confidence, n_tests, n_bands).
    """
    confdnc = 1.0
    nmtests = 0
    nbands = 0

    # Determine if cold surface
    is_cold_sfc = 1 if sfctmp < COLD_SURFACE_TEMP and sfctmp > 0 else 0

    # Antarctic special case
    if flags.antarctic and flags.day:
        confdnc, nmtests, nbands = antarctic_day(
            pxldat, vza, flags.visusd, bt_clr, thresholds, testbits, qa_bits,
        )
        return confdnc, nmtests, nbands

    # Polar region
    if flags.polar:
        if flags.day:
            # Polar daytime
            if flags.snow or flags.ice:
                confdnc, nmtests, nbands = polar_day_snow(
                    pxldat, vza, flags.visusd, flags.hi_elev,
                    bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
                )
            elif flags.desert and flags.coast:
                confdnc, nmtests, nbands = polar_day_desert_coast(
                    pxldat, vza, flags.visusd, flags.vrused, False,
                    flags.hi_elev, 0.0, bt_clr, is_cold_sfc,
                    thresholds, testbits, qa_bits,
                )
            elif flags.coast:
                confdnc, nmtests, nbands = polar_day_coast(
                    pxldat, vza, flags.visusd, flags.vrused, False,
                    flags.hi_elev, bt_clr, is_cold_sfc,
                    thresholds, testbits, qa_bits,
                )
            elif flags.desert:
                confdnc, nmtests, nbands = polar_day_desert(
                    pxldat, vza, flags.visusd, flags.vrused, False,
                    flags.hi_elev, 0.0, bt_clr, is_cold_sfc,
                    thresholds, testbits, qa_bits,
                )
            elif flags.water:
                confdnc, nmtests, nbands = polar_day_ocean(
                    pxldat, vza, flags.snglnt, flags.visusd, False,
                    sst, refang, flags.sh_ocean, bt_clr,
                    thresholds, testbits, qa_bits,
                )
            else:
                # Polar land (default)
                confdnc, nmtests, nbands = polar_day_land(
                    pxldat, vza, flags.visusd, flags.vrused, False,
                    flags.hi_elev, bt_clr, is_cold_sfc,
                    thresholds, testbits, qa_bits,
                )
        else:
            # Polar nighttime
            if flags.snow or flags.ice:
                confdnc, nmtests, nbands = polar_nite_snow(
                    pxldat, vza, sfctmp, bt_clr, is_cold_sfc,
                    thresholds, testbits, qa_bits,
                )
            elif flags.water:
                confdnc, nmtests, nbands = polar_nite_ocean(
                    indat_3x3_11um if indat_3x3_11um is not None else np.zeros((3, 3)),
                    pxldat, vza, sfctmp, flags.sh_ocean, flags.uniform,
                    bt_clr, thresholds, testbits, qa_bits,
                )
            else:
                # Polar land (default)
                confdnc, nmtests, nbands = polar_nite_land(
                    pxldat, vza, sfctmp, flags.hi_elev,
                    bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
                )
        return confdnc, nmtests, nbands

    # Non-polar region
    if flags.day:
        # Daytime
        if flags.snow or flags.ice:
            confdnc, nmtests, nbands = day_snow(
                pxldat, vza, flags.visusd, flags.hi_elev,
                bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
            )
        elif flags.water:
            confdnc, nmtests, nbands = ocean_day(
                pxldat, vza, flags.snglnt, flags.visusd, False,
                sst, refang, flags.sh_ocean, bt_clr,
                thresholds, testbits, qa_bits,
            )
        else:
            # Land variants
            if flags.desert and flags.coast:
                confdnc, nmtests, nbands = land_day_desert_coast(
                    pxldat, bt_clr, vza, is_cold_sfc, flags.hi_elev,
                    0.0, thresholds, testbits, qa_bits,
                )
            elif flags.coast:
                confdnc, nmtests, nbands = land_day_coast(
                    pxldat, bt_clr, vza, is_cold_sfc, flags.hi_elev,
                    thresholds, testbits, qa_bits,
                )
            elif flags.desert:
                confdnc, nmtests, nbands = land_day_desert(
                    pxldat, bt_clr, vza, is_cold_sfc, flags.hi_elev,
                    0.0, thresholds, testbits, qa_bits,
                )
            else:
                confdnc, nmtests, nbands = land_day_standard(
                    pxldat, bt_clr, vza, is_cold_sfc, flags.hi_elev,
                    thresholds, testbits, qa_bits,
                )
    else:
        # Nighttime
        if flags.snow or flags.ice:
            confdnc, nmtests, nbands = nite_snow(
                pxldat, vza, sfctmp, bt_clr, is_cold_sfc,
                thresholds, testbits, qa_bits,
            )
        elif flags.water:
            confdnc, nmtests, nbands = ocean_nite(
                indat_3x3_11um if indat_3x3_11um is not None else np.zeros((3, 3)),
                pxldat, vza, sfctmp, flags.sh_ocean, flags.uniform,
                bt_clr, thresholds, testbits, qa_bits,
            )
        else:
            # Land (night)
            confdnc, nmtests, nbands = land_nite(
                pxldat, lat, vza, flags.coast, flags.desert,
                flags.hi_elev, flags.sh_lake, sfctmp, eco_type,
                0.0, bt_clr, is_cold_sfc, thresholds, testbits, qa_bits,
            )

    return confdnc, nmtests, nbands


def run_cloud_mask_swath(
    pxldat_swath: np.ndarray,
    lat_swath: np.ndarray,
    lon_swath: np.ndarray,
    elevation_swath: np.ndarray,
    lsf_swath: np.ndarray,
    sza_swath: np.ndarray,
    vza_swath: np.ndarray,
    glint_angle_swath: np.ndarray,
    eco_type_swath: np.ndarray,
    snow_mask_swath: np.ndarray,
    sst_swath: np.ndarray,
    nwp_sfctmp_swath: np.ndarray,
    nwp_pmsl_swath: np.ndarray,
    nwp_u_wind_swath: np.ndarray,
    nwp_v_wind_swath: np.ndarray,
    nwp_precip_water_swath: np.ndarray,
    bt_clr_swath: np.ndarray,
    sensor_id: int,
    thresholds: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run cloud mask algorithm for an entire swath.

    This is the main entry point for processing a full satellite swath.

    Args:
        pxldat_swath: (n_elem, n_line, 25) pixel data array.
        lat_swath: (n_elem, n_line) latitude array.
        lon_swath: (n_elem, n_line) longitude array.
        elevation_swath: (n_elem, n_line) elevation array.
        lsf_swath: (n_elem, n_line) land/sea flag array.
        sza_swath: (n_elem, n_line) solar zenith angle array.
        vza_swath: (n_elem, n_line) viewing zenith angle array.
        glint_angle_swath: (n_elem, n_line) sun glint angle array.
        eco_type_swath: (n_elem, n_line) ecosystem type array.
        snow_mask_swath: (n_elem, n_line) snow mask array.
        sst_swath: (n_elem, n_line) SST array.
        nwp_sfctmp_swath: (n_elem, n_line) NWP surface temperature array.
        nwp_pmsl_swath: (n_elem, n_line) NWP MSL pressure array.
        nwp_u_wind_swath: (n_elem, n_line) NWP U-wind array.
        nwp_v_wind_swath: (n_elem, n_line) NWP V-wind array.
        nwp_precip_water_swath: (n_elem, n_line) NWP precipitable water array.
        bt_clr_swath: (n_elem, n_line, 7) clear-sky BT array.
        sensor_id: Sensor ID (21=FY-3D, 22=FY-3E).
        thresholds: Threshold dictionary.

    Returns:
        Tuple of (cm_bitarray, cm_qa_bitarray, cm_tmp, confidence_array):
        - cm_bitarray: (n_elem, n_line, 6) uint8 testbits array
        - cm_qa_bitarray: (n_elem, n_line, 10) uint8 QA bits array
        - cm_tmp: (n_elem, n_line) int32 cloud mask values (0-3)
        - confidence_array: (n_elem, n_line) float64 confidence values
    """
    n_elem, n_line = lat_swath.shape

    # Output arrays
    cm_bitarray = np.zeros((n_elem, n_line, 6), dtype=np.uint8)
    cm_qa_bitarray = np.zeros((n_elem, n_line, 10), dtype=np.uint8)
    cm_tmp = np.zeros((n_elem, n_line), dtype=np.int32)
    confidence_array = np.zeros((n_elem, n_line), dtype=np.float64)

    # Process each pixel
    for j in range(n_line):
        for i in range(n_elem):
            # Extract 3x3 neighborhood for spatial tests
            indat_3x3_11um = _extract_3x3(pxldat_swath[:, :, 23], i, j, n_elem, n_line)
            indat_3x3_vis = _extract_3x3(pxldat_swath[:, :, 2], i, j, n_elem, n_line)

            result = run_cloud_mask_pixel(
                pxldat=pxldat_swath[i, j, :],
                lat=lat_swath[i, j],
                lon=lon_swath[i, j],
                elevation=elevation_swath[i, j],
                lsf=int(lsf_swath[i, j]),
                sza=sza_swath[i, j],
                vza=vza_swath[i, j],
                glint_angle=glint_angle_swath[i, j],
                eco_type=int(eco_type_swath[i, j]),
                snow_mask_val=int(snow_mask_swath[i, j]),
                sst=sst_swath[i, j],
                nwp_sfctmp=nwp_sfctmp_swath[i, j],
                nwp_pmsl=nwp_pmsl_swath[i, j],
                nwp_u_wind=nwp_u_wind_swath[i, j],
                nwp_v_wind=nwp_v_wind_swath[i, j],
                nwp_precip_water=nwp_precip_water_swath[i, j],
                sensor_id=sensor_id,
                bt_clr=bt_clr_swath[i, j, :],
                thresholds=thresholds,
                indat_3x3_11um=indat_3x3_11um,
                indat_3x3_vis=indat_3x3_vis,
                cm_array=cm_tmp,
                row=i,
                col=j,
                n_rows=n_elem,
                n_cols=n_line,
            )

            cm_bitarray[i, j, :] = result.testbits
            cm_qa_bitarray[i, j, :] = result.qa_bits
            cm_tmp[i, j] = result.cloud_mask
            confidence_array[i, j] = result.confidence

    return cm_bitarray, cm_qa_bitarray, cm_tmp, confidence_array


def _extract_3x3(data_2d: np.ndarray, col: int, row: int, n_cols: int, n_rows: int) -> np.ndarray:
    """Extract 3x3 neighborhood from 2D array.

    Args:
        data_2d: 2D data array.
        col: Center column index.
        row: Center row index.
        n_cols: Total columns.
        n_rows: Total rows.

    Returns:
        3x3 array with BAD_DATA for out-of-bounds pixels.
    """
    result = np.full((3, 3), BAD_DATA, dtype=np.float64)
    for di in range(-1, 2):
        for dj in range(-1, 2):
            ni = col + di
            nj = row + dj
            if 0 <= ni < n_cols and 0 <= nj < n_rows:
                result[di + 1, dj + 1] = data_2d[ni, nj]
    return result
