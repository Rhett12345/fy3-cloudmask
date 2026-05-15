# FY-3D MERSI-II Cloud Mask Retrieval System

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A modern Python implementation of the FYLAT FY-3D MERSI-II Cloud Mask Retrieval System (V3.2), ported from the original Fortran codebase.

## Overview

This system generates cloud mask products from FY-3D MERSI-II satellite data. It implements the MODIS MOD35-derived cloud detection algorithm with:

- **25-channel processing** (19 VIS + 6 IR)
- **Surface type classification** (land, water, coast, desert, snow, ice, polar)
- **Day/night processing paths** with different spectral tests
- **S-curve confidence interpolation** for cloud probability
- **48-bit testbits** and **80-bit QA bits** cloud mask encoding
- **4-level cloud confidence**: cloudy(0), probably cloudy(1), probably clear(2), confident clear(3)

## Key Features

- **Exact algorithm fidelity** - Preserves original Fortran algorithm logic
- **NumPy/Numba acceleration** - Performance-sensitive parts use JIT compilation
- **Modular design** - Clean separation of concerns
- **YAML configuration** - Human-readable config files
- **HDF5 output** - Standard satellite data format
- **CLI interface** - Command-line tool for batch processing
- **Comprehensive tests** - 34 unit tests covering all components

## Project Structure

```
fy3_cloudmask/
├── pyproject.toml                    # Package configuration
├── README.md                         # This file
├── config/
│   ├── default.yaml                  # Default configuration
│   ├── sensors/
│   │   └── fy3d_mersi_ii.yaml       # Sensor parameters
│   └── thresholds/
│       └── mersi_ii3d_v8.yaml       # Algorithm thresholds
├── src/
│   └── fy3_cloudmask/
│       ├── __init__.py
│       ├── __main__.py              # Package entry point
│       ├── cli.py                   # Command-line interface
│       ├── config.py                # Configuration loader
│       ├── constants.py             # All magic numbers and constants
│       ├── pipeline.py              # End-to-end processing pipeline
│       ├── algorithm/
│       │   ├── __init__.py
│       │   ├── cloud_mask.py        # Main cloud mask driver
│       │   ├── confidence.py        # S-curve confidence functions
│       │   ├── bitops.py            # Bit manipulation utilities
│       │   ├── spatial.py           # Spatial analysis (tview, 3x3 stats)
│       │   ├── surface_classifier.py # Surface type classification
│       │   └── tests/
│       │       ├── land_day.py      # Daytime land tests
│       │       ├── land_nite.py     # Nighttime land tests
│       │       ├── ocean_day.py     # Daytime ocean tests
│       │       ├── ocean_nite.py    # Nighttime ocean tests
│       │       ├── polar_day.py     # Polar daytime tests
│       │       ├── polar_nite.py    # Polar nighttime tests
│       │       ├── snow_tests.py    # Snow/ice surface tests
│       │       └── restoral.py      # Post-processing restoral tests
│       ├── data/                    # Data I/O (future)
│       └── output/
│           ├── writer.py            # HDF5 output writer
│           └── cloud_amount.py      # Cloud amount computation
├── tests/
│   ├── test_cloud_mask.py           # Cloud mask integration tests
│   └── test_algorithms.py           # Algorithm component tests
└── scripts/
    └── convert_thresholds.py        # Threshold file converter
```

## Installation

### Prerequisites

- Python 3.9 or later
- conda (recommended)

### Install Dependencies

```bash
# Create conda environment
conda create -n cloudmask python=3.10
conda activate cloudmask

# Install dependencies
pip install numpy numba h5py pyyaml click pytest

# Install package in development mode
pip install -e .
```

### Quick Install

```bash
pip install -e .
```

## Usage

### Command Line Interface

```bash
# Process single orbit
python -m fy3_cloudmask process \
    --config config/default.yaml \
    --l1b data/FY3D_MERSI_20230101_0000_1000M_MS.HDF \
    --geo data/FY3D_MERSI_GEO_20230101_0000_1000M.HDF \
    --output output/

# Process batch (multiple orbits)
python -m fy3_cloudmask batch \
    --config config/default.yaml \
    --start 2023-01-01 \
    --end 2023-01-31 \
    --data data/ \
    --output output/ \
    --workers 4

# Convert threshold file
python -m fy3_cloudmask convert-thresholds \
    --input coeff/fylat_thresholds.mersi.ii3d.v8 \
    --output config/thresholds/mersi_ii3d_v8.yaml
```

### Python API

```python
from fy3_cloudmask import run_cloud_mask_pixel, run_cloud_mask_swath
from fy3_cloudmask import load_config, write_cloud_mask

# Load configuration
config = load_config('config/default.yaml')

# Process single pixel
result = run_cloud_mask_pixel(
    pxldat=pxldat,          # 25-element array
    lat=35.0,
    lon=115.0,
    elevation=100.0,
    lsf=1,                  # Land
    sza=30.0,               # Solar zenith angle
    vza=10.0,               # Viewing zenith angle
    glint_angle=60.0,
    eco_type=1,             # Forest
    snow_mask_val=0,
    sst=290.0,              # Sea surface temperature
    nwp_sfctmp=290.0,       # NWP surface temperature
    nwp_pmsl=1013.0,        # NWP mean sea level pressure
    nwp_u_wind=5.0,
    nwp_v_wind=3.0,
    nwp_precip_water=20.0,
    sensor_id=21,           # FY-3D
    bt_clr=bt_clr,          # Clear-sky BT from RTM
    thresholds=thresholds,
)

print(f"Cloud mask: {result.cloud_mask}")  # 0-3
print(f"Confidence: {result.confidence}")  # 0.0-1.0

# Process entire swath
cm_bitarray, cm_qa, cm_tmp, confidence = run_cloud_mask_swath(
    pxldat_swath=pxldat_swath,  # (2048, 2000, 25)
    lat_swath=lat_swath,
    lon_swath=lon_swath,
    # ... other parameters
    sensor_id=21,
    thresholds=thresholds,
)

# Write output
write_cloud_mask('output/CLM.h5', cm_bitarray, cm_qa, cm_tmp, confidence, lon, lat)
```

## Algorithm Description

### Cloud Mask Decision Tree

```
1. Classify pixel surface type
   ├── Land / Water / Coast / Desert
   ├── Snow / Ice
   ├── Polar (>60° latitude)
   └── Day / Night (SZA > 85°)

2. Dispatch to appropriate test function
   ├── Daytime Land: land_day_standard, land_day_coast, land_day_desert
   ├── Daytime Ocean: ocean_day
   ├── Nighttime Land: land_nite
   ├── Nighttime Ocean: ocean_nite
   ├── Polar variants: polar_day_*, polar_nite_*
   └── Snow/Ice: day_snow, nite_snow, antarctic_day

3. Apply spectral tests (per group)
   Group 1: IR forward model fit (PFMFT/NFMFT)
   Group 2: BTD tests (11-12μm, 11-4μm, 8-11μm)
   Group 3: Visible tests (0.64μm, GEMI ratio)
   Group 4: NIR test (1.38μm cirrus)

4. Compute confidence
   - Each test returns confidence [0, 1]
   - Group confidence = minimum of test confidences
   - Final confidence = geometric mean of group confidences

5. Post-processing (restoral tests)
   - Land/coast restoral
   - Sun glint adjustment
   - Shallow water correction
   - Spatial variability check
   - Cloud adjacency check
   - Thin cirrus IR check
   - Shadow detection

6. Encode output
   - Confidence → 2-bit encoding (cloudy/prob_cloudy/prob_clear/clear)
   - Assemble 48-bit testbits
   - Assemble 80-bit QA bits
   - Compute 5km cloud amount
```

### Confidence Encoding

| Confidence Range | Cloud Mask Value | Meaning |
|-----------------|------------------|---------|
| > 0.99 | 3 | Confident clear |
| > 0.95 | 2 | Probably clear |
| > 0.66 | 1 | Probably cloudy |
| ≤ 0.66 | 0 | Cloudy |

### Bit Layout

**Testbits (6 bytes = 48 bits):**
- Byte 0 (bits 0-7): processed(0), conf_lsb(1), conf_msb(2), day(3), no_sunglint(4), no_snow_ice(5), coast(6), desert(7)
- Byte 1 (bits 8-15): nco(8), thin_cirrus_solar(9), shadow(10), thin_cirrus_ir(11), cloud_adj(12), pfmft(14), nfmft(15)
- Byte 2 (bits 16-23): nir_138(16), btd_11_12(18), btd_11_4(19), ref_064(20), gemi(21)
- Byte 3 (bits 24-31): temporal(24), land_restoral(26), suspended_dust(28)

## Configuration

### Main Configuration (config/default.yaml)

```yaml
sensor:
  sensor_id: 21          # 21=FY-3D, 22=FY-3E
  n_elem: 2048           # Pixels per line
  n_line: 2000           # Lines per swath

paths:
  coeff_dir: ./coeff     # Coefficient files
  output_dir: ./output   # Output directory

algorithm:
  cloudmask_id: 1        # Enable cloud mask
  cloudamount_id: 1      # Enable cloud amount
```

### Threshold Configuration (config/thresholds/mersi_ii3d_v8.yaml)

```yaml
snow_mask:
  bt11_threshold: 261.0
  ndsi_threshold: 0.4
  ref086_threshold: 0.11

land_day:
  ref064: [0.24, 0.20, 0.16, 1.0]    # [locut, midpt, hicut, power]
  ref138: [0.04, 0.035, 0.03, 1.0]
  vrat: [1.80, 1.85, 1.90, 1.0]

pfmft:
  bt_11_max: 310.0
  btd_min: 0.0
  land: [4.0, 3.5, 3.0, 1.0]
  cold: [2.0, 1.5, 1.0, 1.0]
```

## Testing

```bash
# Run all tests
PYTHONPATH=src python -m pytest tests/ -v

# Run specific test file
PYTHONPATH=src python -m pytest tests/test_algorithms.py -v

# Run with coverage
PYTHONPATH=src python -m pytest tests/ --cov=fy3_cloudmask
```

### Test Coverage

- **34 unit tests** covering:
  - Confidence functions (6 tests)
  - Bit operations (7 tests)
  - Spatial analysis (5 tests)
  - Surface classification (5 tests)
  - Cloud mask pixel processing (2 tests)
  - Additional algorithm tests (9 tests)

## Output Products

### CLM (Cloud Mask)

- **File format**: HDF5
- **Datasets**:
  - `Cloud_Mask`: (n_elem, n_line, 6) uint8 - Test bits
  - `Quality_Assurance`: (n_elem, n_line, 10) uint8 - QA bits
  - `Cloud_Mask_Value`: (n_elem, n_line) int32 - Cloud mask (0-3)
  - `Confidence`: (n_elem, n_line) float64 - Confidence (0-1)
  - `Longitude`, `Latitude`: Geolocation

### CLA (Cloud Amount)

- **File format**: HDF5
- **Resolution**: 5km (5x5 pixel boxes)
- **Datasets**:
  - `Cloud_Amount`: (n_elem_5km, n_line_5km) uint8 - Cloud cover percentage (0-100)
  - `Cloud_Amount_QA`: Quality flag (0=bad, 1=low, 2=high)

## Performance

- **Single orbit processing**: ~2-5 minutes (depending on hardware)
- **Memory usage**: ~4-8 GB for single orbit
- **Acceleration**: Numba JIT compilation for hot loops

## Migration from Fortran

This Python implementation is a complete rewrite of the original Fortran system:

| Component | Original Fortran | Python Implementation |
|-----------|------------------|----------------------|
| Lines of code | ~37,000 | ~5,000 |
| Configuration | Namelist (.nml) | YAML |
| Data I/O | HDF5 Fortran | h5py |
| NWP processing | wgrib2 shell scripts | xarray + cfgrib |
| Ancillary data | HDF4 Fortran | pyhdf |
| RTM (PFAAST/PLoD) | Included | Deferred (future) |

### What's Preserved

- All cloud detection algorithms
- All threshold values
- Bit layout and encoding
- Decision tree logic
- Surface classification rules

### What's Improved

- Modular, maintainable code structure
- Comprehensive unit tests
- YAML configuration
- Python ecosystem integration
- Better error handling
- Documentation

## Dependencies

### Required

- `numpy` >= 1.20
- `numba` >= 0.50
- `h5py` >= 3.0
- `pyyaml` >= 5.0

### Optional

- `click` >= 8.0 (for CLI)
- `pytest` >= 7.0 (for testing)
- `xarray` + `cfgrib` (for NWP GRIB reading)
- `pyhdf` (for HDF4 ancillary data)

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Original Fortran code by Min Min (minmin@cma.gov.cn), National Satellite Meteorological Center, CMA
- MODIS MOD35 Cloud Mask Algorithm (NASA)
- APOLLO 11-12μm BTD lookup table

## Contact

For questions or issues, please open an issue on GitHub or contact the maintainer.

## Changelog

### V3.2.0 (2026-05-15)

- Initial Python release
- Complete rewrite from Fortran
- All cloud mask algorithms implemented
- Comprehensive test suite
- YAML configuration
- CLI interface
- HDF5 output

---

**Note**: This is a research codebase. For production use, please validate against reference datasets.
