"""Constants and magic numbers for the FY-3D MERSI-II cloud mask system.

All values extracted from the Fortran source code to ensure algorithmic fidelity.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Sensor dimensions
# ---------------------------------------------------------------------------
SENSOR_CONFIGS = {
    1: {"name": "MODIS-to-MERSI", "n_elem": 1354, "n_line": 2030},
    2: {"name": "MODIS-to-MERSI (MERSI fmt)", "n_elem": 2048, "n_line": 2000},
    3: {"name": "VIIRS-to-MERSI", "n_elem": 3354, "n_line": 3030},
    21: {"name": "FY-3D MERSI-II", "n_elem": 2048, "n_line": 2000},
    22: {"name": "FY-3E MERSI-II", "n_elem": 2048, "n_line": 2000},
}

# ---------------------------------------------------------------------------
# Channel configuration (25 channels: 19 VIS + 6 IR)
# ---------------------------------------------------------------------------
N_CHANNELS = 25
N_VIS = 19
N_IR = 6

# Wavelengths in micrometers (0-indexed, bands 0-24)
VIS_WAVELENGTHS_UM = np.array([
    0.470, 0.550, 0.650, 0.860, 1.380, 1.640, 2.130,
    0.410, 0.440, 0.490, 0.555, 0.670, 0.709, 0.746,
    0.865, 0.905, 0.936, 0.940, 1.030,
], dtype=np.float64)

IR_WAVELENGTHS_UM = np.array([
    3.800, 4.050, 7.230, 8.560, 10.710, 11.950,
], dtype=np.float64)

IR_WAVENUMBERS_CM = np.array([
    2643.0, 2472.0, 1383.0, 1168.0, 933.0, 837.0,
], dtype=np.float64)

# Channel flags: 0 = VIS, 1 = IR (30 elements, first 19=VIS, next 6=IR, rest unused)
CHAN_FLAG = np.array([0]*19 + [1]*6 + [0]*5, dtype=np.int32)

# ---------------------------------------------------------------------------
# Special band indices (1-based, matching Fortran pxldat indexing)
# In Python arrays we use 0-based, so band_1 = index 0, band_2 = index 1, etc.
# These are 1-based for documentation; actual code uses 0-based.
# ---------------------------------------------------------------------------
# pxldat(1)=band1(0.47um), pxldat(2)=band2(0.55um), pxldat(3)=band3(0.65um),
# pxldat(4)=band4(0.86um), pxldat(5)=band5(1.38um), pxldat(7)=band7(2.13um),
# pxldat(19)=band19(1.03um), pxldat(20)=IR1(3.8um), pxldat(23)=IR4(8.5um),
# pxldat(24)=IR5(11.0um), pxldat(25)=IR6(12.0um)

# 0-based indices for direct numpy array access
BAND_047 = 0   # 0.47um
BAND_055 = 1   # 0.55um (used in NDSI)
BAND_064 = 2   # 0.64um (0.65um, main VIS test)
BAND_086 = 3   # 0.86um
BAND_138 = 4   # 1.38um (cirrus detection) -- but note: for sensor_id>20, swapped with band 19
BAND_164 = 5   # 1.64um
BAND_213 = 6   # 2.13um (NDSI for Aqua-like)
BAND_103 = 18  # 1.03um
IR_38 = 19     # 3.8um  (pxldat index 20, 0-based = 19)
IR_40 = 20     # 4.05um (pxldat index 21, 0-based = 20)
IR_73 = 21     # 7.3um  (pxldat index 22, 0-based = 21)
IR_85 = 22     # 8.5um  (pxldat index 23, 0-based = 22)
IR_11 = 23     # 11.0um (pxldat index 24, 0-based = 23)
IR_12 = 24     # 12.0um (pxldat index 25, 0-based = 24)

# For NDSI: 2.25um is band 7 for Aqua-like, but band 7 is 2.13um for MERSI-II
# The Fortran code uses pxldat(7) which is band 7 = 2.13um for MERSI-II
BAND_NDSI_NIR = 6  # 0-based index for 2.13um (used as NIR in NDSI)

# ---------------------------------------------------------------------------
# Quality / missing data values
# ---------------------------------------------------------------------------
BAD_DATA = -999.0
MISSING = -99.0
MISSING_INT = -99

# Valid data ranges
VIS_VALID_MIN = -99.0
VIS_VALID_MAX = 2.3
# IR_VALID_MIN = 0.0
# IR_VALID_MAX = 1000.0
IR_VALID_MIN = 150.5   # clip下限是150K，真实大气BT不可能达到150K
IR_VALID_MAX = 340.0   # 地球表面最高温度约330K，留少量余量

# ---------------------------------------------------------------------------
# Geometric / physical thresholds
# ---------------------------------------------------------------------------
SZA_NIGHT = 85.0           # Solar zenith angle threshold for day/night
POLAR_LAT = 60.0           # Latitude threshold for polar regions
ANTARCTIC_LAT = -60.0      # Latitude threshold for Antarctic
GLINT_ANGLE_MAX = 36.0     # Sun glint angle threshold
COLD_SURFACE_TEMP = 265.0  # Cold surface temperature threshold (K)
HIGH_ELEVATION = 2000.0    # High elevation threshold (m)
GREENLAND_ELEVATION = 200.0  # Greenland elevation threshold (m)

# Desert ecosystem types (IGBP classification)
DESERT_ECOSYSTEM_TYPES = [8, 46, 50, 51, 59, 71, 11, 9, 52]

# Ecosystem types where VRAT is disabled
VRAT_DISABLED_ECOSYSTEMS = [2, 8, 11, 40, 41, 46, 51, 52, 59, 71, 50]

# ---------------------------------------------------------------------------
# Confidence thresholds
# ---------------------------------------------------------------------------
CONF_CLOUDY = 0.66         # Below this = cloudy
CONF_PROB_CLEAR = 0.95     # Above this = probably clear
CONF_CLEAR = 0.99          # Above this = confident clear

# ---------------------------------------------------------------------------
# Cloud mask output encoding (cm_tmp values)
# ---------------------------------------------------------------------------
CM_CLOUDY = 0
CM_PROB_CLOUDY = 1
CM_PROB_CLEAR = 2
CM_CLEAR = 3

# ---------------------------------------------------------------------------
# Bit layout for testbits (6 bytes = 48 bits)
# Byte 0 (bits 0-7):  processed(0), conf_lsb(1), conf_msb(2), day(3), no_sunglint(4),
#                      no_snow_ice(5), coast(6), desert(7)
# Byte 1 (bits 8-15): nco(8), thin_cirrus_solar(9), shadow(10), thin_cirrus_ir(11),
#                      cloud_adj(12), unused(13), pfmft(14), nfmft(15)
# Byte 2 (bits 16-23): nir_138(16), unused(17), btd_11_12(18), btd_11_4(19),
#                       ref_064(20), gemi(21), unused(22-23)
# Byte 3 (bits 24-31): temporal(24), unused(25), land_restoral(26), unused(27),
#                       suspended_dust(28), unused(29-31)
# ---------------------------------------------------------------------------
BIT_PROCESSED = 0
BIT_CONF_LSB = 1
BIT_CONF_MSB = 2
BIT_DAY = 3
BIT_NO_SUNGLINT = 4
BIT_NO_SNOW_ICE = 5
BIT_COAST = 6
BIT_DESERT = 7
BIT_NCO = 8
BIT_THIN_CIRRUS_SOLAR = 9
BIT_SHADOW = 10
BIT_THIN_CIRRUS_IR = 11
BIT_CLOUD_ADJ = 12
BIT_PFMFT = 14
BIT_NFMFT = 15
BIT_NIR_138 = 16
BIT_BTD_11_12 = 18
BIT_BTD_11_4 = 19
BIT_REF_064 = 20
BIT_GEMI = 21
BIT_TEMPORAL = 24
BIT_LAND_RESTORAL = 26
BIT_SUSPENDED_DUST = 28

# Initial bits set in pxinit
INITIAL_BITS = np.array([BIT_NCO, BIT_THIN_CIRRUS_SOLAR, BIT_SHADOW,
                          BIT_THIN_CIRRUS_IR, BIT_CLOUD_ADJ, BIT_TEMPORAL,
                          BIT_SUSPENDED_DUST, 31], dtype=np.int32)

# ---------------------------------------------------------------------------
# NWP configuration
# ---------------------------------------------------------------------------
NWP_SOURCES = {
    1: {"name": "NCEP_Reanalysis", "resolution": 1.0, "levels": 26, "grib": 1},
    2: {"name": "GFS_1p00", "resolution": 1.0, "levels": 26, "grib": 2},
    3: {"name": "T639", "resolution": 0.125, "levels": 36, "grib": 2},
    4: {"name": "NCEP_Reanalysis_G2", "resolution": 1.0, "levels": 26, "grib": 2},
    5: {"name": "GFS_0p50", "resolution": 0.5, "levels": 26, "grib": 2},
    6: {"name": "GRAPES_GFS", "resolution": 0.25, "levels": 40, "grib": 2},
    7: {"name": "GDAS_0p25", "resolution": 0.25, "levels": 31, "grib": 2},
    8: {"name": "GFS_0p25", "resolution": 0.25, "levels": 31, "grib": 2},
    9: {"name": "GFS_0p50_41L", "resolution": 0.5, "levels": 41, "grib": 2},
    10: {"name": "GFS_0p25_41L", "resolution": 0.25, "levels": 41, "grib": 2},
}

# Standard pressure levels (hPa) for 26-level NWP
P_LEVELS_26 = np.array([
    10, 20, 30, 50, 70, 100, 150, 200, 250,
    300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
    800, 850, 900, 925, 950, 975, 1000,
], dtype=np.float64)

# 101-level pressure grid for RTM interpolation
NLEVELS_INTERP = 101
PTOP = 25.0   # hPa
PBOT = 300.0  # hPa
MIN_WATER_VAPOR = 0.0003  # g/kg minimum water vapor mixing ratio
WV_EXTRAP_POWER = 3.0     # Power law exponent for WV extrapolation

# CO2 mixing ratio (ppm)
CO2_PPM = 380.0

# ---------------------------------------------------------------------------
# Radiative transfer
# ---------------------------------------------------------------------------
RTM_VZA_BINSIZE = 0.01  # cosine VZA bin size
DEFAULT_EMISSIVITY = 0.99

# Planck table
PLANCK_T_MIN = 159.0  # K
PLANCK_T_MAX = 360.0  # K
PLANCK_T_STEP = 1.0   # K
PLANCK_N_ENTRIES = 201

# ---------------------------------------------------------------------------
# Cloud amount
# ---------------------------------------------------------------------------
CLOUD_AMOUNT_BOX_SIZE = 5  # 5x5 pixel box
CLOUD_AMOUNT_MIN_VALID = 15  # Minimum valid pixels for quality=1
CLOUD_AMOUNT_MAX_VALID = 25  # Maximum valid pixels for quality=2

# ---------------------------------------------------------------------------
# Surface data grids
# ---------------------------------------------------------------------------
NISE_GRID_SHAPE = (721, 721)
OISST_GRID_SHAPE = (720, 1440)  # (lat, lon) at 0.25 degree
OISST_C_TO_K_OFFSET = 273.15
EMISSIVITY_GRID_SHAPE = (3600, 7200)  # 0.05 degree
ALBEDO_GRID_SHAPE = (2700, 5400)      # 0.06666 degree

# Water surface emissivity for 16 channels (used for ocean pixels)
WATER_SURFACE_EMISSIVITY = np.array([
    0.990, 0.990, 0.990, 0.990, 0.990, 0.990,  # channels 1-6
    0.990, 0.990, 0.990, 0.990, 0.990, 0.990,  # channels 7-12
    0.990, 0.990, 0.990, 0.990,                  # channels 13-16
], dtype=np.float64)

# ---------------------------------------------------------------------------
# US Standard Atmosphere (101 levels, hardcoded from Fortran)
# These are reference profiles used as fallback when NWP data is unavailable.
# ---------------------------------------------------------------------------
# Pressure (hPa), Temperature (K), Water vapor (g/kg), Ozone (g/kg), Height (m)
# Loaded from Fortran data arrays in platform_module.f90
# For now, we define the structure; actual values loaded from file.
