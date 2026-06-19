# FY-3D MERSI-II Cloud Mask Retrieval System

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Fortran-native cloud mask retrieval system for FY-3D MERSI-II satellite data, implementing the MODIS MOD35-derived cloud detection algorithm. Python is used only for data I/O and orchestration; all algorithm computation runs in Fortran with OpenMP parallelization.

## Overview

This system generates cloud mask products from FY-3D MERSI-II satellite data:

- **25-channel processing** (19 VIS + 6 IR)
- **Surface type classification** (land, water, coast, desert, snow, ice, polar)
- **Day/night processing paths** with different spectral tests
- **S-curve confidence interpolation** for cloud probability
- **48-bit testbits** and **80-bit QA bits** cloud mask encoding
- **4-level cloud confidence**: cloudy(0), probably cloudy(1), probably clear(2), confident clear(3)
- **NWP integration** (FNL GRIB2 data for surface temperature, pressure, winds, precipitable water)
- **Recalibration support** for daily solar reflectance band coefficients

## Quick Start

```bash
# Build native Fortran backend
cd ext/ && ./build.sh --install

# Process single orbit
python -m fy3_cloudmask process \
    --config config/default.yaml \
    --l1b data/L1b.HDF --geo data/GEO.HDF --output output/

# Run tests
PYTHONPATH=src python -m pytest tests/ -v
```

## Architecture

### Fortran-Native Design

All cloud mask algorithm computation runs in Fortran with OpenMP parallelization. Python handles data reading (L1B, GEO, NWP), calibration, and HDF5 output writing.

```
Python: data I/O → calibration → NWP interpolation
    ↓
Fortran (pybind11): cloud mask algorithm (~37,000 lines)
    ↓                 18 spectral test paths + 8 restoral tests
    ↓                 S-curve confidence + bit encoding
Python: HDF5 output + cloud amount computation
```

Performance: ~10-30s per orbit (2048 × 2000 pixels) with OpenMP.

### Algorithm Flow (per pixel)

```
L1b HDF5 + GEO HDF5
    ↓
SurfaceClassifier → land/water/coast/desert/snow/ice/polar + day/night
    ↓
18 spectral test paths (LandDay, ocean_day, LandNite, polar_day, etc.)
    ↓
S-curve confidence → [0, 1]
    ↓
8 restoral tests (sunglint, cirrus, shadow, etc.)
    ↓
Bit encoding → 48-bit testbits + 80-bit QA → 4-level cloud mask
    ↓
HDF5 output (cm + conf datasets)
```

### Key Modules

| Module | Description |
|--------|-------------|
| `config.py` | `FY3Config` dataclass, YAML loading |
| `constants.py` | Band indices, bit positions, physical constants |
| `pipeline.py` | `CloudMaskPipeline` orchestrates full processing |
| `algorithm/native_backend.py` | pybind11 bridge to Fortran native engine |
| `algorithm/cloud_mask.py` | `CloudMaskResult` dataclass |
| `io/recalibration.py` | Daily recalibration coefficient manager |
| `output/writer.py` | HDF5 output writer |
| `output/cloud_amount.py` | 5km cloud amount from 5×5 pixel boxes |

### Fortran Source (`src/fortran/`)

| Directory | Description |
|-----------|-------------|
| `core/` | Foundation modules (names, constant, planck, numerical) |
| `cloudmask/` | Cloud mask algorithm (~37,000 lines): 18 test paths, 8 restoral tests, spatial analysis, threshold reader |
| `utils/` | String utilities + C sources (median filter) |
| `c_api/` | ISO_C_BINDING wrappers with OpenMP parallelization |

## Native Backend Build

```bash
cd ext/

# Build and install
./build.sh --install

# Debug build
./build.sh --debug --install

# Clean
./build.sh --clean
```

Requirements: `gfortran`, `g++`, `gcc`, `hdf5`, `pybind11`.

## Configuration

### Hierarchy

```
config/default.yaml              → runtime paths, processing options
config/sensors/fy3d_mersi_ii.yaml → band wavelengths, sensor geometry
config/thresholds/mersi_ii3d_v8.yaml → 790+ algorithm thresholds
```

## Output Products

### HDF5 Structure

```
cm          (nElem, nLine) int32    — 0=cloudy, 1=prob_cloudy, 2=prob_clear, 3=conf_clear
conf        (nElem, nLine) float32  — confidence [0, 1]
lat         (nElem, nLine) float32  — latitude
lon         (nElem, nLine) float32  — longitude
```

## Recalibration

Supports daily recalibration coefficients for 7 solar reflectance bands (channels 1-7).

```python
from fy3_cloudmask.io.recalibration import RecalibrationManager

mgr = RecalibrationManager('../fy3d_recali')
cal0, cal1 = mgr.load_coefficients('20200308')
pxldat = read_l1b_data(l1b_path, recal_cal0=cal0, recal_cal1=cal1)
```

## Project Structure

```
fy3_cloudmask/
├── README.md
├── config/
│   ├── default.yaml
│   ├── sensors/fy3d_mersi_ii.yaml
│   └── thresholds/mersi_ii3d_v8.yaml
├── src/fy3_cloudmask/
│   ├── algorithm/
│   │   ├── cloud_mask.py          # CloudMaskResult dataclass
│   │   ├── native_backend.py      # pybind11 Fortran bridge
│   │   └── __init__.py
│   ├── io/
│   │   └── recalibration.py       # Daily recalibration manager
│   ├── output/
│   │   ├── writer.py              # HDF5 writer
│   │   └── cloud_amount.py        # 5km cloud amount
│   ├── config.py
│   ├── constants.py
│   └── pipeline.py
├── src/fortran/
│   ├── core/                      # Foundation modules
│   ├── cloudmask/                 # Cloud mask algorithm (18 test paths + 8 restoral)
│   ├── utils/                     # String/C utilities
│   └── c_api/                     # ISO_C_BINDING + OpenMP wrappers
├── ext/
│   ├── build.sh                   # Fortran backend build script
│   ├── cloudmask_pybind.cpp       # pybind11 wrapper
│   └── include/
├── scripts/
│   ├── process_20200308_f90.py    # Fortran processing pipeline
│   ├── run_fortran_only.py        # Batch Fortran-native processing
│   ├── validate_cloudmask.py      # MYD35 validation
│   ├── spatial_analysis.py        # Spatial noise analysis
│   └── convert_thresholds.py      # Fortran→YAML threshold converter
└── tests/
    ├── test_algorithms.py
    ├── test_cloud_mask.py
    └── test_pipeline_e2e.py
```

## Testing

```bash
# Unit tests
PYTHONPATH=src python -m pytest tests/ -v

# Native backend validation
PYTHONPATH=src python scripts/run_fortran_only.py
```

## Dependencies

**Required**: `numpy`, `h5py`, `pyyaml`, `cfgrib`, `scipy`

**Build**: `gfortran`, `g++`, `hdf5`, `pybind11`

## License

MIT License. Original Fortran code by Min Min (minmin@cma.gov.cn), National Satellite Meteorological Center, CMA.
