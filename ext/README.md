# FY-3D Cloud Mask Native Engine (C++/Fortran Hybrid)

## Architecture

```
Python (config, CLI, I/O, pipeline)
    ↓ pybind11 (zero-copy numpy arrays)
C++ (OpenMP parallel pixel loop)
    ↓ ISO_C_BINDING
Fortran (core algorithm: spectral tests, conf_test, bit ops)
```

### Why hybrid?

| Layer | Language | Reason |
|-------|----------|--------|
| Config/CLI/IO | Python | Easy to maintain, rich ecosystem |
| Pixel loop | C++ + OpenMP | 4-8x speedup from parallelism, no GIL |
| Core algorithm | Fortran | Original code, proven correctness, fast |

### Performance comparison

| Backend | Single orbit (2048×2000) | Parallelism |
|---------|-------------------------|-------------|
| Python/Numba (original) | ~12 min | None (GIL) |
| C++/Fortran (this) | ~10-30s | OpenMP (all cores) |

**Expected speedup: 30-100x** over the Python/Numba version.

## Files

```
ext/
├── README.md                    # This file
├── build.sh                     # Build script
├── Makefile                     # Alternative build
├── CMakeLists.txt               # CMake build
├── cloudmask_pybind.cpp         # pybind11 Python bindings
├── include/
│   └── cloudmask_engine.hpp     # C++ header
└── fortran/
    ├── cloudmask_c_api.f90      # ISO_C_BINDING wrapper
    ├── cloudmask_data_arrays.f90 # Thread-safe module (with OpenMP threadprivate)
    └── *.inc                    # Fortran include files (copied from original)
```

## Prerequisites

- **Fortran compiler**: gfortran (≥7) or ifort
- **C++ compiler**: g++ (≥7) or icpx
- **OpenMP**: Usually included with compilers
- **HDF5**: Development libraries
- **Python**: ≥3.9 with numpy, pybind11

### Install on CentOS/RHEL 8

```bash
sudo dnf install gcc-gfortran gcc-c++ hdf5-devel python3-devel
pip install pybind11 numpy
```

### Install on Ubuntu 20.04+

```bash
sudo apt install gfortran g++ libhdf5-dev python3-dev
pip install pybind11 numpy
```

## Build

### Quick build

```bash
cd ext/
./build.sh --install
```

### Manual build

```bash
cd ext/
make -j$(nproc) install
```

### CMake build

```bash
cd ext/
mkdir build && cd build
cmake .. -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
make -j$(nproc)
make install
```

## Verify

```python
from fy3_cloudmask.algorithm.native_backend import is_native_available, get_backend_info

print(is_native_available())  # True if built successfully
print(get_backend_info())     # Shows backend details
```

## Usage

The native backend is automatically used when available. No code changes needed:

```python
from fy3_cloudmask import CloudMaskPipeline

pipeline = CloudMaskPipeline('config/default.yaml')
result = pipeline.process_orbit(l1b_path, geo_path, output_dir)
# Automatically uses C++/Fortran backend if built
```

### Direct API access

```python
from fy3_cloudmask.algorithm.native_backend import process_swath_native

result = process_swath_native(
    ref_vis, tbb_ir, lat, lon, satzen, solzen, relaz, glint,
    sfctmp, pmsl, uwind, vwind, tpw, elev, eco, snow_mask, btclr,
    n_elem, n_line
)
```

## OpenMP Configuration

Control the number of threads:

```bash
export OMP_NUM_THREADS=8    # Use 8 threads
export OMP_PROC_BIND=close  # Bind threads to cores
```

Or in Python:

```python
import os
os.environ['OMP_NUM_THREADS'] = '8'
```

## Thread Safety

The Fortran module `cloudmask_data_arrays` uses `!$omp threadprivate` for all
per-pixel state variables. Each OpenMP thread has its own private copy of:

- `testbits(6)`, `qa_bits(10)` — bit arrays
- `confdnc`, `plat`, `plon`, `vza`, etc. — scalar state
- All logical flags (`polar`, `land`, `day`, `snow`, etc.)

Global state (thresholds, sensor config) is read-only during the parallel loop
and set once before entering the parallel region.

## Troubleshooting

### Build fails with "pybind11 not found"

```bash
pip install pybind11
# or
pip install scikit-build-core  # for CMake-based builds
```

### Build fails with "HDF5 not found"

```bash
# CentOS/RHEL
sudo dnf install hdf5-devel

# Ubuntu
sudo apt install libhdf5-dev

# Conda
conda install hdf5
```

### Runtime error "libgfortran.so not found"

```bash
# Add to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(gfortran -print-file-name=.)
```

### Performance not as expected

Check OpenMP is actually using multiple threads:

```python
import os
print(f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'not set')}")
```

Monitor CPU usage with `htop` or `top` during processing.
