"""End-to-end pipeline test with real FY-3D MERSI-II data.

Tests the complete cloud mask workflow:
1. Read L1b + GEO HDF5 data
2. Convert DN to physical units (reflectance/BT)
3. Run cloud mask algorithm
4. Write HDF5 output
5. Verify output validity
"""

import logging
import time
from pathlib import Path

import h5py
import numpy as np
import pytest
import yaml

from fy3_cloudmask.algorithm import run_cloud_mask_swath
from fy3_cloudmask.algorithm.cloud_mask import run_cloud_mask_pixel
from fy3_cloudmask.output import (
    compute_cloud_amount,
    write_cloud_mask,
    write_cloud_amount,
    write_combined_product,
)
from fy3_cloudmask.constants import BAD_DATA, SZA_NIGHT

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)

# Test data paths
L1B_FILE = '/home/liusy2020/yuq/fydata/mersi_20230606/FY3D_MERSI_GBAL_L1_20230606_0500_1000M_MS.HDF'
GEO_FILE = '/home/liusy2020/yuq/fydata/mersi_20230606/FY3D_MERSI_GBAL_L1_20230606_0500_GEO1K_MS.HDF'
THRESHOLDS_FILE = '/home/liusy2020/yuq/cloudmask/fy3_cloudmask/config/thresholds/mersi_ii3d_v8.yaml'
OUTPUT_DIR = '/home/liusy2020/yuq/data-yuq/cloudmask_test/output'
NWP_FILE = '/data/nwp/20230606/ORG/gfs0p25_41L_20230606_06_00'


def read_l1b_data(l1b_path: str) -> dict:
    """Read L1b HDF5 data and convert DN to physical units.

    VIS channels: DN -> reflectance (0-1+)
    IR channels: DN -> brightness temperature (K)

    For FY-3D MERSI-II (from Fortran source):
    - VIS: reflectance = (coef0 + coef1*DN + coef2*DN^2) * 0.01 / cos(sza) / esd^2
    - IR: DN IS radiance in mW/(m^2*sr*cm^-1). Use Planck inversion for BT.

    Returns:
        dict with pxldat (n_elem, n_line, 25) and calibration metadata.
    """
    # Planck constants (matching Fortran planck_module.f90 brite_m function)
    H_PLANCK = 6.62606876e-34   # Planck constant (J·s)
    C_LIGHT = 2.99792458e+08    # Speed of light (m/s)
    K_BOLTZMANN = 1.3806503e-23 # Boltzmann constant (J/K)
    C1_PLANCK = 2.0 * H_PLANCK * C_LIGHT * C_LIGHT  # 2hc² (W·m²/sr)
    C2_PLANCK = H_PLANCK * C_LIGHT / K_BOLTZMANN     # hc/k (m·K)

    # FY-3D MERSI-II central wavenumbers (cm⁻¹) from platform_module.f90
    IR_WAVENUMBERS = np.array([2643.4359, 2471.654, 1382.621, 1168.182, 933.364, 836.941])

    # FY-3D MERSI-II TBB temperature correction coefficients (fylat_sensor_id=21)
    # From planck_module.f90 cwn_fy3d/tcs_fy3d/tci_fy3d data statements
    # BT_final = (BT_raw - tci) / tcs  (rad2tbb direction)
    TCS_FY3D = np.array([0.9992917440, 0.9994814177, 0.9989956900, 0.9997135336, 0.9980397975, 0.9983777125])
    TCI_FY3D = np.array([0.50718071650, 0.3493280160, 0.40925130837, 0.1014073981, 0.57633464244, 0.4317181810])

    with h5py.File(l1b_path, 'r') as f:
        # Read raw DN data (band, line, pixel)
        vis_250 = f['Data/EV_250_Aggr.1KM_RefSB'][:].astype(np.float64)    # (4, 2000, 2048)
        vis_1km = f['Data/EV_1KM_RefSB'][:].astype(np.float64)              # (15, 2000, 2048)
        ir_250 = f['Data/EV_250_Aggr.1KM_Emissive'][:].astype(np.float64)   # (2, 2000, 2048)
        ir_1km = f['Data/EV_1KM_Emissive'][:].astype(np.float64)            # (4, 2000, 2048)

        # Calibration coefficients
        vis_cal = f['Calibration/VIS_Cal_Coeff'][:]       # (19, 3): intercept, slope, quadratic
        esd = f.attrs['EarthSun Distance Ratio']

        # Day/night flag from file
        day_flag = f.attrs.get('Day Or Night Flag', b'D')

    n_band_vis_250, n_line, n_pixel = vis_250.shape
    n_band_ir_250 = ir_250.shape[0]
    n_band_ir_1km = ir_1km.shape[0]

    # Assemble pxldat: (n_pixel, n_line, 25)
    pxldat = np.zeros((n_pixel, n_line, 25), dtype=np.float64)

    # VIS channels: convert DN to reflectance
    # For FY-3D: reflectance = (coef0 + coef1*DN + coef2*DN^2) * 0.01 / esd^2
    # (cos(sza) correction applied later per-pixel)
    esd2 = esd * esd

    # Bands 1-4 (250m aggregated to 1km)
    for b in range(4):
        c0, c1, c2 = vis_cal[b]
        dn = vis_250[b]  # (n_line, n_pixel)
        refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        refl = np.clip(refl, -0.1, 2.0)
        pxldat[:, :, b] = refl.T  # transpose to (n_pixel, n_line)

    # Bands 5-19 (1km)
    for b in range(15):
        c0, c1, c2 = vis_cal[b + 4]
        dn = vis_1km[b]
        refl = (c0 + c1 * dn + c2 * dn * dn) * 0.01 / esd2
        refl = np.clip(refl, -0.1, 2.0)
        pxldat[:, :, b + 4] = refl.T

    # IR channels: DN -> radiance -> BT using Planck formula
    # Step 1: radiance_mw = (DN + intercept) * slope [mW/(m²·sr·cm⁻¹)]
    # Step 2: BT = c2*vs / log(c1*vs^3 / (1e-5*R) + 1)  [Fortran brite_m]
    # Step 3: BT_final = (BT_raw - tci) / tcs

    # Read Slope/Intercept from dataset attributes
    with h5py.File(l1b_path, 'r') as f:
        slope_250 = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Slope']      # (2,)
        intercept_250 = f['Data/EV_250_Aggr.1KM_Emissive'].attrs['Intercept']  # (2,)
        slope_1km = f['Data/EV_1KM_Emissive'].attrs['Slope']                # (4,)
        intercept_1km = f['Data/EV_1KM_Emissive'].attrs['Intercept']        # (4,)

    # Band layout (0-based pxldat indices):
    #   band 20 (3.8um)  -> pxldat[:,:,19], IR_WAVENUMBERS[0], ir_1km[0]
    #   band 21 (4.0um)  -> pxldat[:,:,20], IR_WAVENUMBERS[1], ir_1km[1]
    #   band 22 (7.2um)  -> pxldat[:,:,21], IR_WAVENUMBERS[2], ir_1km[2]
    #   band 23 (8.6um)  -> pxldat[:,:,22], IR_WAVENUMBERS[3], ir_1km[3]
    #   band 24 (10.7um) -> pxldat[:,:,23], IR_WAVENUMBERS[4], ir_250[0]
    #   band 25 (12.0um) -> pxldat[:,:,24], IR_WAVENUMBERS[5], ir_250[1]
    # Fortran brite_m formula: BT = c2*vs / log(c1*vs^3 / (1e-5*R) + 1)
    for b in range(n_band_ir_1km):  # bands 20-23
        radiance_mw = (ir_1km[b] + intercept_1km[b]) * slope_1km[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b]
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs**3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b]) / TCS_FY3D[b]
        bt = np.clip(bt, 150.0, 350.0)
        pxldat[:, :, b + 19] = bt.T

    for b in range(n_band_ir_250):  # bands 24-25
        radiance_mw = (ir_250[b] + intercept_250[b]) * slope_250[b]
        radiance_mw = np.maximum(radiance_mw, 0.01)
        wvn = IR_WAVENUMBERS[b + 4]
        vs = 100.0 * wvn
        bt_raw = C2_PLANCK * vs / np.log(C1_PLANCK * vs**3 / (1e-5 * radiance_mw) + 1.0)
        bt = (bt_raw - TCI_FY3D[b + 4]) / TCS_FY3D[b + 4]
        bt = np.clip(bt, 150.0, 350.0)
        pxldat[:, :, b + 23] = bt.T

    return {
        'pxldat': pxldat,
        'day_flag': day_flag,
    }


def read_geo_data(geo_path: str) -> dict:
    """Read GEO HDF5 data.

    Returns:
        dict with geolocation and geometry arrays (n_pixel, n_line).
    """
    with h5py.File(geo_path, 'r') as f:
        lat = f['Geolocation/Latitude'][:].astype(np.float64)       # (2000, 2048)
        lon = f['Geolocation/Longitude'][:].astype(np.float64)      # (2000, 2048)
        dem = f['Geolocation/DEM'][:].astype(np.float64)            # (2000, 2048)
        lsf = f['Geolocation/LandSeaMask'][:].astype(np.int32)      # (2000, 2048)
        eco = f['Geolocation/LandCover'][:].astype(np.int32)        # (2000, 2048)
        sza_raw = f['Geolocation/SolarZenith'][:].astype(np.float64)
        vza_raw = f['Geolocation/SensorZenith'][:].astype(np.float64)
        saa_raw = f['Geolocation/SolarAzimuth'][:].astype(np.float64)
        vaa_raw = f['Geolocation/SensorAzimuth'][:].astype(np.float64)

    # Scale angles (stored as degrees * 100)
    sza = sza_raw / 100.0
    vza = vza_raw / 100.0
    saa = saa_raw / 100.0
    vaa = vaa_raw / 100.0

    # Compute sun glint angle (simplified)
    # glint_angle ≈ arccos(cos(sza)*cos(vza) + sin(sza)*sin(vza)*cos(saa-vaa))
    sza_rad = np.radians(sza)
    vza_rad = np.radians(vza)
    daa_rad = np.radians(saa - vaa)
    cos_glint = np.cos(sza_rad) * np.cos(vza_rad) + np.sin(sza_rad) * np.sin(vza_rad) * np.cos(daa_rad)
    cos_glint = np.clip(cos_glint, -1.0, 1.0)
    glint_angle = np.degrees(np.arccos(cos_glint))

    # Transpose from (line, pixel) to (pixel, line)
    return {
        'lat': lat.T,
        'lon': lon.T,
        'elevation': dem.T,
        'lsf': lsf.T,
        'eco_type': eco.T,
        'sza': sza.T,
        'vza': vza.T,
        'glint_angle': glint_angle.T,
    }


def load_thresholds(thresholds_path: str) -> dict:
    """Load algorithm thresholds from YAML file."""
    with open(thresholds_path) as f:
        return yaml.safe_load(f)


def read_nwp_binary(nwp_path: str, nvar: int = 283) -> dict:
    """Read GFS 0.25-degree 41-layer NWP binary file.

    Binary format: Fortran direct access, (nlon, nlat, nvar) float32.
    Fields in order: psfc(1), pmsl(2), tsfc(3), zsfc(4), albedo(5),
    t_sigma(6), rh_sigma(7), u_sigma(8), v_sigma(9), tpw(10), ...

    Args:
        nwp_path: Path to binary NWP file.
        nvar: Number of variables in file (default 283 for 41-layer).

    Returns:
        dict with NWP fields on 0.25-degree grid (lon: -180..180, lat: -90..90).
    """
    nlon, nlat = 1440, 721
    data = np.fromfile(nwp_path, dtype=np.float32)
    expected = nlon * nlat * nvar
    if data.size < expected:
        raise ValueError(f"NWP file too small: {data.size} < {expected}")
    # Reshape: raw data is (nlon, nlat, nvar) in Fortran column-major order
    # In Python C-order: first axis varies fastest -> (nvar, nlat, nlon) after transpose
    arr = data[:expected].reshape(nlon, nlat, nvar).transpose(2, 1, 0)  # (nvar, nlat, nlon)

    # NWP grid: lat from -90 to 90 (index 0 = -90), lon from 0 to 359.75
    # Shift longitude to -180..180
    # Fortran code: out(1:180) = in(181:360), out(181:360) = in(1:180)
    # In 0-indexed: out[0:720] = in[720:1440], out[720:1440] = in[0:720]
    arr_shifted = np.empty_like(arr)
    arr_shifted[:, :, :720] = arr[:, :, 720:]
    arr_shifted[:, :, 720:] = arr[:, :, :720]
    lon = np.concatenate([np.arange(720, 1440) * 0.25 - 180, np.arange(0, 720) * 0.25 - 180 + 360])
    # Actually simpler: just create the shifted lon grid
    lon = np.linspace(-180, 179.75, nlon)
    lat = np.linspace(-90, 90, nlat)

    # Extract key fields (0-indexed)
    result = {
        'lon': lon,
        'lat': lat,
        'psfc': arr_shifted[0],     # surface pressure (Pa)
        'pmsl': arr_shifted[1],     # mean sea level pressure (Pa)
        'tsfc': arr_shifted[2],     # surface temperature (K)
        'zsfc': arr_shifted[3],     # surface height (m)
        'u_wind': arr_shifted[7],   # U-wind at sigma 0.995 (m/s)
        'v_wind': arr_shifted[8],   # V-wind at sigma 0.995 (m/s)
        'tpw': arr_shifted[9],      # total precipitable water (mm)
    }
    return result


def interpolate_nwp_to_pixels(nwp: dict, lat_px: np.ndarray, lon_px: np.ndarray) -> dict:
    """Interpolate NWP fields to satellite pixel locations using nearest-neighbor.

    Args:
        nwp: NWP data dict from read_nwp_binary.
        lat_px: Satellite pixel latitudes (n_elem, n_line).
        lon_px: Satellite pixel longitudes (n_elem, n_line).

    Returns:
        dict with interpolated NWP fields at pixel locations.
    """
    nwp_lat = nwp['lat']  # (721,) -90 to 90
    nwp_lon = nwp['lon']  # (1440,) -180 to 179.75

    # Find nearest grid indices for each pixel
    # lat: -90 to 90, monotonically increasing
    lat_idx = np.searchsorted(nwp_lat, lat_px.ravel()).clip(0, len(nwp_lat) - 1)
    # lon: -180 to 179.75, monotonically increasing
    lon_idx = np.searchsorted(nwp_lon, lon_px.ravel()).clip(0, len(nwp_lon) - 1)

    shape = lat_px.shape
    fields = {}
    for key in ['tsfc', 'pmsl', 'u_wind', 'v_wind', 'tpw']:
        fields[key] = nwp[key][lat_idx, lon_idx].reshape(shape)

    return fields


class TestPipelineE2E:
    """End-to-end pipeline tests with real FY-3D data."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Load data and thresholds."""
        self.l1b_path = L1B_FILE
        self.geo_path = GEO_FILE
        self.thresholds_path = THRESHOLDS_FILE
        self.nwp_path = NWP_FILE
        self.output_dir = Path(OUTPUT_DIR)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check data files exist
        if not Path(self.l1b_path).exists():
            pytest.skip(f"L1b file not found: {self.l1b_path}")
        if not Path(self.geo_path).exists():
            pytest.skip(f"GEO file not found: {self.geo_path}")

        self.thresholds = load_thresholds(self.thresholds_path)

        # Load NWP data
        if Path(self.nwp_path).exists():
            logger.info(f"Loading NWP data from {self.nwp_path}...")
            self.nwp = read_nwp_binary(self.nwp_path)
            logger.info(f"  NWP tsfc range: {self.nwp['tsfc'].min():.1f} - {self.nwp['tsfc'].max():.1f} K")
            logger.info(f"  NWP pmsl range: {self.nwp['pmsl'].min()/100:.1f} - {self.nwp['pmsl'].max()/100:.1f} hPa")
            logger.info(f"  NWP tpw range: {self.nwp['tpw'].min():.1f} - {self.nwp['tpw'].max():.1f} mm")
        else:
            logger.warning(f"NWP file not found: {self.nwp_path}, using dummy values")
            self.nwp = None

    def test_read_geo_data(self):
        """Test reading GEO data from HDF5."""
        geo = read_geo_data(self.geo_path)

        assert geo['lat'].shape == (2048, 2000), f"Unexpected lat shape: {geo['lat'].shape}"
        assert geo['lon'].shape == (2048, 2000)
        assert geo['sza'].shape == (2048, 2000)
        assert geo['lsf'].shape == (2048, 2000)

        # Check value ranges
        assert -90 <= geo['lat'].min() <= geo['lat'].max() <= 90
        assert -180 <= geo['lon'].min() <= geo['lon'].max() <= 180
        assert 0 <= geo['sza'].min() <= geo['sza'].max() <= 180
        assert 0 <= geo['vza'].min() <= geo['vza'].max() <= 90

        logger.info(f"GEO data loaded: lat=[{geo['lat'].min():.1f}, {geo['lat'].max():.1f}], "
                     f"lon=[{geo['lon'].min():.1f}, {geo['lon'].max():.1f}]")
        logger.info(f"  SZA=[{geo['sza'].min():.1f}, {geo['sza'].max():.1f}], "
                     f"LandSeaMask values: {np.unique(geo['lsf'])}")

    def test_read_l1b_data(self):
        """Test reading L1b data from HDF5."""
        l1b = read_l1b_data(self.l1b_path)

        assert l1b['pxldat'].shape == (2048, 2000, 25), f"Unexpected shape: {l1b['pxldat'].shape}"

        # Check VIS channels are in reasonable reflectance range
        vis_data = l1b['pxldat'][:, :, :19]
        vis_valid = vis_data[(vis_data > -0.1) & (vis_data < 2.0)]
        vis_mean = np.mean(vis_valid) if len(vis_valid) > 0 else 0
        logger.info(f"VIS mean reflectance: {vis_mean:.4f} (valid pixels: {len(vis_valid)}/{vis_data.size})")

        # Check IR channels are in reasonable BT range
        ir_data = l1b['pxldat'][:, :, 19:]
        ir_valid = ir_data[(ir_data > 150) & (ir_data < 350)]
        if len(ir_valid) > 0:
            logger.info(f"IR mean BT: {np.mean(ir_valid):.1f}K (valid pixels: {len(ir_valid)}/{ir_data.size})")
            assert np.mean(ir_valid) > 200, f"IR BT too low: {np.mean(ir_valid)}"

    def test_pixel_level_cloud_mask(self):
        """Test cloud mask on individual pixels from real data."""
        l1b = read_l1b_data(self.l1b_path)
        geo = read_geo_data(self.geo_path)

        # Pick a few test pixels from different locations
        test_pixels = [
            (1024, 1000, "Center pixel"),
            (500, 500, "Upper-left region"),
            (1500, 1500, "Lower-right region"),
        ]

        for col, row, desc in test_pixels:
            pxldat = l1b['pxldat'][col, row, :]
            lat = float(geo['lat'][col, row])
            lon = float(geo['lon'][col, row])
            sza = float(geo['sza'][col, row])
            vza = float(geo['vza'][col, row])
            glint = float(geo['glint_angle'][col, row])
            lsf = int(geo['lsf'][col, row])
            eco = int(geo['eco_type'][col, row])
            elev = float(geo['elevation'][col, row])

            # Get NWP data at this pixel
            if self.nwp is not None:
                nwp_interp = interpolate_nwp_to_pixels(
                    self.nwp,
                    np.array([[lat]]),
                    np.array([[lon]]),
                )
                nwp_sfctmp = float(nwp_interp['tsfc'][0, 0])
                nwp_pmsl = float(nwp_interp['pmsl'][0, 0]) / 100.0  # Pa -> hPa
                nwp_u = float(nwp_interp['u_wind'][0, 0])
                nwp_v = float(nwp_interp['v_wind'][0, 0])
                nwp_pw = float(nwp_interp['tpw'][0, 0])
            else:
                bt_11 = pxldat[23]
                nwp_sfctmp = bt_11 if bt_11 > 200 else 290.0
                nwp_pmsl = 1013.0
                nwp_u, nwp_v = 5.0, 3.0
                nwp_pw = 20.0

            # Compute clear-sky BT from NWP surface temperature
            bt_clr = np.array([
                nwp_sfctmp - 10.0,  # 3.8um
                nwp_sfctmp - 20.0,  # 6.7um
                nwp_sfctmp - 10.0,  # 7.3um
                nwp_sfctmp - 5.0,   # 8.5um
                nwp_sfctmp,          # 11um
                nwp_sfctmp - 1.0,   # 12um
                0.0,                 # unused
            ], dtype=np.float64)

            result = run_cloud_mask_pixel(
                pxldat=pxldat,
                lat=lat,
                lon=lon,
                elevation=elev,
                lsf=lsf,
                sza=sza,
                vza=vza,
                glint_angle=glint,
                eco_type=eco,
                snow_mask_val=0,
                sst=BAD_DATA,
                nwp_sfctmp=nwp_sfctmp,
                nwp_pmsl=nwp_pmsl,
                nwp_u_wind=nwp_u,
                nwp_v_wind=nwp_v,
                nwp_precip_water=nwp_pw,
                sensor_id=21,
                bt_clr=bt_clr,
                thresholds=self.thresholds,
            )

            logger.info(f"{desc} ({col},{row}): lat={lat:.2f}, lon={lon:.2f}, "
                        f"sza={sza:.1f}, lsf={lsf}, cloud_mask={result.cloud_mask}, "
                        f"confidence={result.confidence:.3f}, n_tests={result.n_tests}")

            assert result.cloud_mask in [0, 1, 2, 3], f"Invalid cloud_mask: {result.cloud_mask}"
            assert 0.0 <= result.confidence <= 1.0, f"Invalid confidence: {result.confidence}"
            assert result.testbits is not None
            assert result.qa_bits is not None

    def test_swath_level_cloud_mask_small(self):
        """Test cloud mask on a small subset (100x100 pixels) of real data."""
        l1b = read_l1b_data(self.l1b_path)
        geo = read_geo_data(self.geo_path)

        # Extract small subset (100x100 pixels from center)
        col_start, col_end = 974, 1074
        row_start, row_end = 950, 1050

        pxldat_sub = l1b['pxldat'][col_start:col_end, row_start:row_end, :]
        lat_sub = geo['lat'][col_start:col_end, row_start:row_end]
        lon_sub = geo['lon'][col_start:col_end, row_start:row_end]
        elev_sub = geo['elevation'][col_start:col_end, row_start:row_end]
        lsf_sub = geo['lsf'][col_start:col_end, row_start:row_end]
        sza_sub = geo['sza'][col_start:col_end, row_start:row_end]
        vza_sub = geo['vza'][col_start:col_end, row_start:row_end]
        glint_sub = geo['glint_angle'][col_start:col_end, row_start:row_end]
        eco_sub = geo['eco_type'][col_start:col_end, row_start:row_end]

        n_elem, n_line = lat_sub.shape
        logger.info(f"Running swath cloud mask on {n_elem}x{n_line} subset...")

        # Interpolate NWP data to satellite pixels
        if self.nwp is not None:
            logger.info("Interpolating NWP data to satellite pixels...")
            nwp_interp = interpolate_nwp_to_pixels(self.nwp, lat_sub, lon_sub)
            nwp_sfctmp = nwp_interp['tsfc']
            nwp_pmsl = nwp_interp['pmsl'] / 100.0  # Pa -> hPa
            nwp_u = nwp_interp['u_wind']
            nwp_v = nwp_interp['v_wind']
            nwp_pw = nwp_interp['tpw']
            logger.info(f"  NWP sfctmp at center: {nwp_sfctmp[n_elem//2, n_line//2]:.1f} K")
            logger.info(f"  NWP pmsl at center: {nwp_pmsl[n_elem//2, n_line//2]:.1f} hPa")
            logger.info(f"  NWP tpw at center: {nwp_pw[n_elem//2, n_line//2]:.1f} mm")
        else:
            nwp_sfctmp = pxldat_sub[:, :, 23].copy()
            nwp_sfctmp[nwp_sfctmp < 200] = 290.0
            nwp_pmsl = np.full((n_elem, n_line), 1013.0)
            nwp_u = np.full((n_elem, n_line), 5.0)
            nwp_v = np.full((n_elem, n_line), 3.0)
            nwp_pw = np.full((n_elem, n_line), 20.0)

        # Compute clear-sky BT from NWP surface temperature
        # bt_clr[4] = clear-sky 11um BT ≈ NWP surface temperature
        # bt_clr[5] = clear-sky 12um BT ≈ 11um BT - 1K (water vapor absorption)
        bt_clr = np.zeros((n_elem, n_line, 7), dtype=np.float64)
        bt_clr[:, :, 4] = nwp_sfctmp  # clear-sky 11um ~ surface temp
        bt_clr[:, :, 5] = nwp_sfctmp - 1.0  # clear-sky 12um slightly cooler
        bt_clr[:, :, 0] = nwp_sfctmp - 10.0  # 3.8um
        bt_clr[:, :, 1] = nwp_sfctmp - 20.0  # 6.7um water vapor
        bt_clr[:, :, 2] = nwp_sfctmp - 10.0  # 7.3um
        bt_clr[:, :, 3] = nwp_sfctmp - 5.0   # 8.5um

        start_time = time.time()
        cm_bitarray, cm_qa_bitarray, cm_tmp, confidence = run_cloud_mask_swath(
            pxldat_swath=pxldat_sub,
            lat_swath=lat_sub,
            lon_swath=lon_sub,
            elevation_swath=elev_sub,
            lsf_swath=lsf_sub,
            sza_swath=sza_sub,
            vza_swath=vza_sub,
            glint_angle_swath=glint_sub,
            eco_type_swath=eco_sub,
            snow_mask_swath=np.zeros((n_elem, n_line), dtype=np.int32),
            sst_swath=np.zeros((n_elem, n_line), dtype=np.float64),
            nwp_sfctmp_swath=nwp_sfctmp,
            nwp_pmsl_swath=nwp_pmsl,
            nwp_u_wind_swath=nwp_u,
            nwp_v_wind_swath=nwp_v,
            nwp_precip_water_swath=nwp_pw,
            bt_clr_swath=bt_clr,
            sensor_id=21,
            thresholds=self.thresholds,
        )
        elapsed = time.time() - start_time

        # Verify output shapes
        assert cm_tmp.shape == (n_elem, n_line)
        assert confidence.shape == (n_elem, n_line)
        assert cm_bitarray.shape == (n_elem, n_line, 6)
        assert cm_qa_bitarray.shape == (n_elem, n_line, 10)

        # Compute statistics
        n_total = cm_tmp.size
        n_cloudy = np.sum(cm_tmp == 0)
        n_prob_cloudy = np.sum(cm_tmp == 1)
        n_prob_clear = np.sum(cm_tmp == 2)
        n_clear = np.sum(cm_tmp == 3)

        logger.info(f"Swath processing completed in {elapsed:.1f}s ({n_total/elapsed:.0f} pixels/s)")
        logger.info(f"  Total pixels: {n_total}")
        logger.info(f"  Cloudy:       {n_cloudy} ({100*n_cloudy/n_total:.1f}%)")
        logger.info(f"  Prob cloudy:  {n_prob_cloudy} ({100*n_prob_cloudy/n_total:.1f}%)")
        logger.info(f"  Prob clear:   {n_prob_clear} ({100*n_prob_clear/n_total:.1f}%)")
        logger.info(f"  Clear:        {n_clear} ({100*n_clear/n_total:.1f}%)")
        logger.info(f"  Mean confidence: {np.mean(confidence):.3f}")

        # Basic validity checks
        assert np.all(np.isin(cm_tmp, [0, 1, 2, 3])), "Invalid cloud mask values"
        assert np.all(confidence >= 0) and np.all(confidence <= 1), "Invalid confidence values"

    def test_write_hdf5_output(self):
        """Test writing cloud mask output to HDF5."""
        l1b = read_l1b_data(self.l1b_path)
        geo = read_geo_data(self.geo_path)

        # Use a very small subset for output test (50x50)
        col_start, col_end = 1000, 1050
        row_start, row_end = 975, 1025

        pxldat_sub = l1b['pxldat'][col_start:col_end, row_start:row_end, :]
        lat_sub = geo['lat'][col_start:col_end, row_start:row_end]
        lon_sub = geo['lon'][col_start:col_end, row_start:row_end]
        elev_sub = geo['elevation'][col_start:col_end, row_start:row_end]
        lsf_sub = geo['lsf'][col_start:col_end, row_start:row_end]
        sza_sub = geo['sza'][col_start:col_end, row_start:row_end]
        vza_sub = geo['vza'][col_start:col_end, row_start:row_end]
        glint_sub = geo['glint_angle'][col_start:col_end, row_start:row_end]
        eco_sub = geo['eco_type'][col_start:col_end, row_start:row_end]

        n_elem, n_line = lat_sub.shape

        if self.nwp is not None:
            nwp_interp = interpolate_nwp_to_pixels(self.nwp, lat_sub, lon_sub)
            nwp_sfctmp = nwp_interp['tsfc']
            nwp_pmsl = nwp_interp['pmsl'] / 100.0
            nwp_u = nwp_interp['u_wind']
            nwp_v = nwp_interp['v_wind']
            nwp_pw = nwp_interp['tpw']
        else:
            nwp_sfctmp = pxldat_sub[:, :, 23].copy()
            nwp_sfctmp[nwp_sfctmp < 200] = 290.0
            nwp_pmsl = np.full((n_elem, n_line), 1013.0)
            nwp_u = np.full((n_elem, n_line), 5.0)
            nwp_v = np.full((n_elem, n_line), 3.0)
            nwp_pw = np.full((n_elem, n_line), 20.0)

        bt_clr = np.zeros((n_elem, n_line, 7), dtype=np.float64)
        bt_clr[:, :, 4] = nwp_sfctmp
        bt_clr[:, :, 5] = nwp_sfctmp - 1.0
        bt_clr[:, :, 0] = nwp_sfctmp - 10.0
        bt_clr[:, :, 1] = nwp_sfctmp - 20.0
        bt_clr[:, :, 2] = nwp_sfctmp - 10.0
        bt_clr[:, :, 3] = nwp_sfctmp - 5.0

        cm_bitarray, cm_qa_bitarray, cm_tmp, confidence = run_cloud_mask_swath(
            pxldat_swath=pxldat_sub,
            lat_swath=lat_sub,
            lon_swath=lon_sub,
            elevation_swath=elev_sub,
            lsf_swath=lsf_sub,
            sza_swath=sza_sub,
            vza_swath=vza_sub,
            glint_angle_swath=glint_sub,
            eco_type_swath=eco_sub,
            snow_mask_swath=np.zeros((n_elem, n_line), dtype=np.int32),
            sst_swath=np.zeros((n_elem, n_line), dtype=np.float64),
            nwp_sfctmp_swath=nwp_sfctmp,
            nwp_pmsl_swath=nwp_pmsl,
            nwp_u_wind_swath=nwp_u,
            nwp_v_wind_swath=nwp_v,
            nwp_precip_water_swath=nwp_pw,
            bt_clr_swath=bt_clr,
            sensor_id=21,
            thresholds=self.thresholds,
        )

        # Compute cloud amount
        qa_processed = (cm_tmp >= 0).astype(np.int32)
        cloud_amount, cloud_amount_qa = compute_cloud_amount(cm_tmp, qa_processed)

        # Compute 5km geolocation
        box_size = 5
        n_elem_5km = n_elem // box_size
        n_line_5km = n_line // box_size
        lon_5km = lon_sub[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]
        lat_5km = lat_sub[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]

        # Write output
        output_path = self.output_dir / 'test_cloudmask_output.h5'
        write_combined_product(
            str(output_path),
            cm_bitarray, cm_qa_bitarray, cm_tmp, confidence,
            cloud_amount, cloud_amount_qa,
            lon_sub, lat_sub, lon_5km, lat_5km,
            attributes={
                'input_l1b': self.l1b_path,
                'input_geo': self.geo_path,
                'test_description': 'E2E pipeline test output',
            },
        )

        assert output_path.exists(), f"Output file not created: {output_path}"

        # Verify output file
        with h5py.File(str(output_path), 'r') as f:
            logger.info(f"Output file: {output_path}")
            logger.info(f"  Size: {output_path.stat().st_size / 1024:.1f} KB")
            logger.info(f"  Datasets:")
            for name, obj in f.items():
                if isinstance(obj, h5py.Group):
                    for ds_name, ds in obj.items():
                        if isinstance(ds, h5py.Dataset):
                            logger.info(f"    {name}/{ds_name}: shape={ds.shape}, dtype={ds.dtype}")
                elif isinstance(obj, h5py.Dataset):
                    logger.info(f"    {name}: shape={obj.shape}, dtype={obj.dtype}")

            assert 'Cloud_Mask_1km' in f, "Cloud mask group missing"
            assert 'Cloud_Mask_1km/Cloud_Mask_Value' in f, "Cloud mask values missing"
            assert 'Cloud_Amount_5km' in f, "Cloud amount group missing"

        logger.info(f"Output written successfully: {output_path}")

    def test_full_orbit_with_nwp(self):
        """Test full orbit (2048x2000) cloud mask with real NWP data."""
        l1b = read_l1b_data(self.l1b_path)
        geo = read_geo_data(self.geo_path)

        pxldat = l1b['pxldat']
        lat = geo['lat']
        lon = geo['lon']
        elev = geo['elevation']
        lsf = geo['lsf']
        sza = geo['sza']
        vza = geo['vza']
        glint = geo['glint_angle']
        eco = geo['eco_type']

        n_elem, n_line = lat.shape
        logger.info(f"Running full orbit cloud mask on {n_elem}x{n_line} = {n_elem*n_line} pixels...")

        # Interpolate NWP data
        if self.nwp is not None:
            logger.info("Interpolating NWP data to full orbit...")
            t0 = time.time()
            nwp_interp = interpolate_nwp_to_pixels(self.nwp, lat, lon)
            nwp_sfctmp = nwp_interp['tsfc']
            nwp_pmsl = nwp_interp['pmsl'] / 100.0
            nwp_u = nwp_interp['u_wind']
            nwp_v = nwp_interp['v_wind']
            nwp_pw = nwp_interp['tpw']
            logger.info(f"  NWP interpolation took {time.time()-t0:.1f}s")
            logger.info(f"  NWP sfctmp range: {nwp_sfctmp.min():.1f} - {nwp_sfctmp.max():.1f} K")
            logger.info(f"  NWP pmsl range: {nwp_pmsl.min():.1f} - {nwp_pmsl.max():.1f} hPa")
            logger.info(f"  NWP tpw range: {nwp_pw.min():.1f} - {nwp_pw.max():.1f} mm")
        else:
            logger.warning("No NWP data, using dummy values")
            nwp_sfctmp = pxldat[:, :, 23].copy()
            nwp_sfctmp[nwp_sfctmp < 200] = 290.0
            nwp_pmsl = np.full((n_elem, n_line), 1013.0)
            nwp_u = np.full((n_elem, n_line), 5.0)
            nwp_v = np.full((n_elem, n_line), 3.0)
            nwp_pw = np.full((n_elem, n_line), 20.0)

        bt_clr = np.zeros((n_elem, n_line, 7), dtype=np.float64)
        bt_clr[:, :, 4] = nwp_sfctmp
        bt_clr[:, :, 5] = nwp_sfctmp - 1.0
        bt_clr[:, :, 0] = nwp_sfctmp - 10.0
        bt_clr[:, :, 1] = nwp_sfctmp - 20.0
        bt_clr[:, :, 2] = nwp_sfctmp - 10.0
        bt_clr[:, :, 3] = nwp_sfctmp - 5.0

        start_time = time.time()
        cm_bitarray, cm_qa_bitarray, cm_tmp, confidence = run_cloud_mask_swath(
            pxldat_swath=pxldat,
            lat_swath=lat,
            lon_swath=lon,
            elevation_swath=elev,
            lsf_swath=lsf,
            sza_swath=sza,
            vza_swath=vza,
            glint_angle_swath=glint,
            eco_type_swath=eco,
            snow_mask_swath=np.zeros((n_elem, n_line), dtype=np.int32),
            sst_swath=np.zeros((n_elem, n_line), dtype=np.float64),
            nwp_sfctmp_swath=nwp_sfctmp,
            nwp_pmsl_swath=nwp_pmsl,
            nwp_u_wind_swath=nwp_u,
            nwp_v_wind_swath=nwp_v,
            nwp_precip_water_swath=nwp_pw,
            bt_clr_swath=bt_clr,
            sensor_id=21,
            thresholds=self.thresholds,
        )
        elapsed = time.time() - start_time

        # Statistics
        n_total = cm_tmp.size
        n_cloudy = np.sum(cm_tmp == 0)
        n_prob_cloudy = np.sum(cm_tmp == 1)
        n_prob_clear = np.sum(cm_tmp == 2)
        n_clear = np.sum(cm_tmp == 3)

        logger.info(f"Full orbit processing completed in {elapsed:.1f}s ({n_total/elapsed:.0f} pixels/s)")
        logger.info(f"  Total pixels: {n_total}")
        logger.info(f"  Cloudy:       {n_cloudy} ({100*n_cloudy/n_total:.1f}%)")
        logger.info(f"  Prob cloudy:  {n_prob_cloudy} ({100*n_prob_cloudy/n_total:.1f}%)")
        logger.info(f"  Prob clear:   {n_prob_clear} ({100*n_prob_clear/n_total:.1f}%)")
        logger.info(f"  Clear:        {n_clear} ({100*n_clear/n_total:.1f}%)")
        logger.info(f"  Mean confidence: {np.mean(confidence):.3f}")

        # Verify output
        assert np.all(np.isin(cm_tmp, [0, 1, 2, 3])), "Invalid cloud mask values"
        assert np.all(confidence >= 0) and np.all(confidence <= 1), "Invalid confidence values"

        # Write full orbit output
        qa_processed = (cm_tmp >= 0).astype(np.int32)
        cloud_amount, cloud_amount_qa = compute_cloud_amount(cm_tmp, qa_processed)

        box_size = 5
        n_elem_5km = n_elem // box_size
        n_line_5km = n_line // box_size
        lon_5km = lon[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]
        lat_5km = lat[2::box_size, 2::box_size][:n_elem_5km, :n_line_5km]

        output_path = self.output_dir / 'FY3D_MERSI_20230606_0500_CLM_CLA_nwp.h5'
        write_combined_product(
            str(output_path),
            cm_bitarray, cm_qa_bitarray, cm_tmp, confidence,
            cloud_amount, cloud_amount_qa,
            lon, lat, lon_5km, lat_5km,
            attributes={
                'input_l1b': self.l1b_path,
                'input_geo': self.geo_path,
                'input_nwp': self.nwp_path,
                'test_description': 'Full orbit with real NWP data',
            },
        )
        logger.info(f"Full orbit output: {output_path} ({output_path.stat().st_size/1024/1024:.1f} MB)")


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
