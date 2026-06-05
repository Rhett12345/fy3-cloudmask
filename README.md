# FY-3D MERSI-II Cloud Mask Retrieval System

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python implementation of the FYLAT FY-3D MERSI-II Cloud Mask Retrieval System (V3.2), ported from the original ~37,000-line Fortran codebase. Supports dual backend architecture with optional C++/Fortran native acceleration.

## Overview

This system generates cloud mask products from FY-3D MERSI-II satellite data. It implements the MODIS MOD35-derived cloud detection algorithm with:

- **25-channel processing** (19 VIS + 6 IR)
- **Surface type classification** (land, water, coast, desert, snow, ice, polar)
- **Day/night processing paths** with different spectral tests
- **S-curve confidence interpolation** for cloud probability
- **48-bit testbits** and **80-bit QA bits** cloud mask encoding
- **4-level cloud confidence**: cloudy(0), probably cloudy(1), probably clear(2), confident clear(3)
- **NWP integration** (GFS 0.25° data for surface temperature, pressure, precipitable water)
- **Recalibration support** for daily solar reflectance band coefficients

## Quick Start

```bash
# Install
pip install -e .

# Build native backend (optional, 30-100x speedup)
cd ext/ && ./build.sh --install

# Process single orbit
python -m fy3_cloudmask process \
    --config config/default.yaml \
    --l1b data/L1b.HDF --geo data/GEO.HDF --output output/

# Run tests
PYTHONPATH=src python -m pytest tests/test_algorithms.py tests/test_cloud_mask.py -v
```

## Architecture

### Dual Backend Design

The system has two mathematically-equivalent backends, auto-selected at runtime:

| Backend | Speed | Use Case |
|---------|-------|----------|
| **Python/Numba** (default) | ~12 min/orbit | Development, debugging, no build required |
| **C++/Fortran + OpenMP** | ~10-30s/orbit | Production, batch processing |

```python
# Auto-selects best available backend
from fy3_cloudmask.algorithm.native_backend import is_native_available, process_swath_native
```

### Algorithm Flow (per pixel)

```
L1b HDF5 + GEO HDF5 + NWP Binary
    ↓
SurfaceClassifier → land/water/coast/desert/snow/ice/polar + day/night
    ↓
18 spectral test paths (land_day, ocean_nite, polar_day, etc.)
    ↓
S-curve confidence → [0, 1]
    ↓
8 restoral tests (sunglint, cirrus, shadow, etc.)
    ↓
Bit encoding → 48-bit testbits + 80-bit QA → 4-level cloud mask
    ↓
HDF5 output (Cloud_Mask_1km/ + Cloud_Amount_5km/)
```

### Key Modules

| Module | Description |
|--------|-------------|
| `config.py` | `FY3Config` dataclass, YAML loading with recursive merge |
| `constants.py` | Band indices, bit positions, physical constants |
| `pipeline.py` | `CloudMaskPipeline` orchestrates full processing |
| `algorithm/cloud_mask.py` | Main driver: `run_cloud_mask_pixel()` and `run_cloud_mask_swath()` |
| `algorithm/tests/` | 18 test path modules (e.g., `land_day.py`, `ocean_nite.py`) |
| `algorithm/native_backend.py` | Dual backend auto-selection and pybind11 interface |
| `io/recalibration.py` | Daily recalibration coefficient manager |
| `output/writer.py` | HDF5 output writer (CLM + CLA products) |
| `output/cloud_amount.py` | 5km cloud amount from 5×5 pixel boxes |

## Recalibration

Supports daily recalibration coefficients for 7 solar reflectance bands (channels 1-7), replacing onboard calibration coefficients.

### Data Structure

```
fy3d_recali/
├── 202208/
│   ├── RAD_20220803.csv
│   ├── RAD_20220808.csv
│   └── ...
├── 202209/
└── ...
```

### CSV Format

```csv
,cal0,cal1,cal2
ch01,-3.264,0.0273,0.0
ch02,-4.324,0.0259,0.0
...
ch07,-2.602,0.0208,0.0
```

### Usage

```python
from fy3_cloudmask.io.recalibration import RecalibrationManager

mgr = RecalibrationManager('../fy3d_recali')
cal0, cal1 = mgr.load_coefficients('20220803')

# Apply to L1b data (replaces onboard VIS_Cal_Coeff)
pxldat = read_l1b_data(l1b_path, recal_cal0=cal0, recal_cal1=cal1)
```

### Recalibration Test

```bash
# Compare onboard vs recalibration (single date)
python scripts/test_recalibration.py --date 20220803 --output /tmp/recal_test

# Batch test with HDF5 output
OMP_NUM_THREADS=1 python scripts/test_recalibration.py --date 20220803 --max-orbits 2
```

Output: HDF5 files with `_onboard.h5` and `_recal.h5` suffixes, plus `comparison_summary.json`.

## Native Backend Build

The optional C++/Fortran backend provides 30-100x speedup via OpenMP parallelization:

```bash
cd ext/

# Build and install
./build.sh --install

# Debug build
./build.sh --debug --install

# Clean
./build.sh --clean
```

Requirements: `gfortran`, `g++`, `gcc`, `hdf5`, `pybind11`. The build script auto-detects conda or system compilers.

## Configuration

### Hierarchy

```
config/default.yaml              → runtime paths, processing options
config/sensors/fy3d_mersi_ii.yaml → band wavelengths, sensor geometry
config/thresholds/mersi_ii3d_v8.yaml → 790+ algorithm thresholds
```

### Main Config (config/default.yaml)

```yaml
sensor:
  sensor_id: 21          # 21=FY-3D, 22=FY-3E
  n_elem: 2048
  n_line: 2000

paths:
  coeff_dir: ./coeff
  output_dir: ./output
```

## Output Products

### HDF5 Structure

```
Cloud_Mask_1km/
  ├── Cloud_Mask          (nElem, nLine, 6) uint8   — test bits
  ├── Quality_Assurance   (nElem, nLine, 10) uint8  — QA bits
  ├── Cloud_Mask_Value    (nElem, nLine) int32       — 0-3
  ├── Confidence          (nElem, nLine) float64     — 0-1
  ├── Longitude           (nElem, nLine) float64
  └── Latitude            (nElem, nLine) float64
Cloud_Amount_5km/
  ├── Cloud_Amount        (nElem/5, nLine/5) uint8   — 0-100%
  ├── Cloud_Amount_QA     (nElem/5, nLine/5) uint8
  ├── Longitude           (nElem/5, nLine/5) float64
  └── Latitude            (nElem/5, nLine/5) float64
```

## Project Structure

```
fy3_cloudmask/
├── pyproject.toml
├── README.md
├── config/
│   ├── default.yaml
│   ├── sensors/fy3d_mersi_ii.yaml
│   └── thresholds/mersi_ii3d_v8.yaml
├── src/fy3_cloudmask/
│   ├── algorithm/
│   │   ├── cloud_mask.py          # Main driver
│   │   ├── confidence.py          # S-curve confidence
│   │   ├── bitops.py              # Bit manipulation
│   │   ├── spatial.py             # Spatial analysis
│   │   ├── surface_classifier.py  # Surface classification
│   │   ├── native_backend.py      # Dual backend interface
│   │   └── tests/                 # 18 spectral test paths
│   ├── io/
│   │   └── recalibration.py       # Daily recalibration manager
│   ├── output/
│   │   ├── writer.py              # HDF5 writer
│   │   └── cloud_amount.py        # 5km cloud amount
│   ├── config.py
│   ├── constants.py
│   └── pipeline.py
├── ext/
│   ├── build.sh                   # Native backend build script
│   ├── cloudmask_pybind.cpp       # pybind11 wrapper
│   └── include/cloudmask_engine.hpp
├── src/fortran/
│   ├── core/                      # Foundation modules
│   ├── cloudmask/                 # Cloud mask algorithm
│   ├── utils/                     # String/C utilities
│   └── c_api/                     # ISO_C_BINDING + OpenMP wrappers
├── scripts/
│   ├── run_fortran_only.py        # Batch Fortran-native processing
│   ├── test_recalibration.py      # Recalibration comparison test
│   ├── find_matched_files.py      # L1b/GEO/NWP file matching
│   └── convert_thresholds.py      # Fortran→YAML threshold converter
└── tests/
    ├── test_algorithms.py
    ├── test_cloud_mask.py
    └── test_pipeline_e2e.py
```

## Testing

```bash
# Unit tests
PYTHONPATH=src python -m pytest tests/test_algorithms.py tests/test_cloud_mask.py -v

# E2E with real FY-3D data (~12 min)
PYTHONPATH=src python -m pytest tests/test_pipeline_e2e.py -v

# Native backend validation
PYTHONPATH=src python scripts/run_fortran_only.py --output /tmp/validation
```

## Dependencies

**Required**: `numpy`, `numba`, `h5py`, `pyyaml`

**Optional**: `click` (CLI), `pytest` (testing), `pybind11` (native backend build)

## License

MIT License. Original Fortran code by Min Min (minmin@cma.gov.cn), National Satellite Meteorological Center, CMA.
