# FY-3D Cloud Mask Native Engine — Deployment & Performance Guide

**Date**: 2026-06-05
**Version**: v3.2.0
**Scope**: C++/Fortran native backend engineering optimization and deployment

---

## 1. Engineering Optimizations Summary

### 1.1 Debug Print Removal

Removed 6 debug print statements from `src/fortran/c_api/cloudmask_c_api.f90`:
- `write(*,*) 'DBG: ...'` and `flush(6)` calls scattered throughout the C API entry point
- These caused unnecessary I/O overhead in production and polluted stdout

### 1.2 Build System Rewrite (`ext/build.sh`)

Complete rewrite for portable deployment:

| Change | Before | After |
|--------|--------|-------|
| Compiler detection | Hardcoded `gfortran`/`g++`/`gcc` | Auto-detect: conda prefix → system PATH → fallback |
| RPATH | Hardcoded conda lib path | `$ORIGIN/../lib` (portable) |
| Optimization | `-O2` only | `-O3` for Fortran, `-O2` for C++, `-DNDEBUG` in release |
| HDF5 discovery | Hardcoded path | Dynamic search: `$CONDA_PREFIX/lib`, `/usr/lib64`, etc. |
| Build artifacts | Stale `.mod` files accumulate | Clean stale `.mod` before compile |

Key deployment features:
- **`--debug`** flag: enables `-g -O0 -fcheck=all -Wall -Wextra`
- **`--clean`** flag: removes entire build directory
- **`--install`** flag: copies `.so` to `src/fy3_cloudmask/` and runs quick verification

### 1.3 Cache-Blocked C++ Transpose (`ext/cloudmask_pybind.cpp`)

Added cache-blocked tiling (BLOCK=32) for C-order ↔ Fortran-order memory transposition:

```cpp
static constexpr int BLOCK = 32;

// 2D transpose: (nElem, nLine) C-order → Fortran column-major
for (int jj = 0; jj < nLine; jj += BLOCK) {
    for (int ii = 0; ii < nElem; ii += BLOCK) {
        for (int j = jj; j < jEnd; j++) {
            for (int i = ii; i < iEnd; i++) {
                dst[i * nLine + j] = src[j * nElem + i];
            }
        }
    }
}
```

3D transpose uses the same tiling, processing all K slices in the inner loop for cache locality.

### 1.4 Dead Code Cleanup

Deleted files:
- `ext/Makefile` — had `-fdefault-integer-8` ABI mismatch bug (Fortran default integers are 4-byte, not 8-byte)
- `ext/fort.41` — leftover Fortran runtime artifact
- 6 one-off scripts: `backend_compare_and_viz.py`, `batch_backend_compare.py`, `detailed_compare.py`, `select_test_dates.py`, `threshold_audit.py`, `convert_yaml_to_fortran.py`
- `scripts/test_dates_manifest.json` — unused manifest

### 1.5 Directory Rename

`retrieval_system_V3.1_cldmask/` → `coeff/`

The directory only contains threshold coefficient files (`.dat`, `.inc`). Updated all references in 6+ files including `native_backend.py`, `cloudmask_c_api.f90`, `build.sh`, and test scripts.

### 1.6 Environment Specification

Created `environment.yml` for reproducible conda environments:
- Python 3.11, NumPy <2.0, HDF5 1.14.*
- Compiler toolchain: `gfortran_linux-64>=12.0`, `gxx_linux-64>=12.0`, `gcc_linux-64>=12.0`
- pybind11 >=2.10

---

## 2. Performance Test Results

### 2.1 Test Configuration

- **Test data**: 10 dates, 20 orbits (2 orbits per date)
- **Data source**: FY-3D MERSI-II L1b + GEO HDF files
- **Platform**: Linux 4.18.0-147.el8.x86_64, 80-core machine
- **Build**: Release mode (`-O3 -DNDEBUG`)

### 2.2 Single-Thread Baseline (OMP_NUM_THREADS=1)

| Metric | Value |
|--------|-------|
| Mean time per orbit | ~6.6 seconds |
| Throughput | ~0.63 Mpix/s |
| Memory per orbit | ~2.1 GB peak |

### 2.3 OpenMP Thread Scaling

Tested with `OMP_NUM_THREADS` = 1, 2, 4, 8, 16, 32, 64, 80:

| Threads | Time (s) | Speedup vs 1T | Efficiency |
|---------|----------|---------------|------------|
| 1 | 6.60 | 1.00x | 100% |
| 2 | 4.12 | 1.60x | 80% |
| 4 | 2.85 | 2.32x | 58% |
| 8 | 2.10 | 3.14x | 39% |
| 16 | 1.85 | 3.57x | 22% |
| 32 | 1.72 | 3.84x | 12% |
| 64 | 1.68 | 3.93x | 6% |
| 80 | 1.68 | 3.93x | 5% |

**Key findings**:
- Diminishing returns beyond 8 threads due to per-pixel Fortran algorithm overhead and memory bandwidth saturation
- The Fortran core algorithm processes pixels sequentially within each OpenMP task — the parallelism is at the swath-strip level, not pixel level
- **Recommendation**: Use `OMP_NUM_THREADS=4~8` for optimal throughput/core ratio on production machines

### 2.4 Critical: OMP_NUM_THREADS Behavior

**Must set `OMP_NUM_THREADS` BEFORE the Python process starts.**

Setting `os.environ['OMP_NUM_THREADS']` after the Fortran library loads has no effect — OpenMP reads the environment variable during runtime initialization, which happens at first `import _cloudmask_native`.

```bash
# Correct: set before Python starts
OMP_NUM_THREADS=4 python -m fy3_cloudmask process ...

# Wrong: has no effect
python -c "import os; os.environ['OMP_NUM_THREADS']='4'; import _cloudmask_native"
```

---

## 3. Deployment Checklist

### 3.1 Target Machine Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| OS | Linux (kernel >= 3.10) | CentOS 7+ / Ubuntu 18.04+ |
| CPU | 4 cores | 8+ cores |
| RAM | 4 GB | 8+ GB (for concurrent orbits) |
| Disk | 1 GB (code + coeff) | 10+ GB (output storage) |
| GCC/GFortran | >= 12.0 | >= 12.0 |
| HDF5 | 1.14.x | 1.14.x (with Fortran bindings) |
| Python | 3.10+ | 3.11 |

### 3.2 Build Steps on Target Machine

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate cloudmask

# 2. Build native engine
cd ext/
./build.sh --install

# 3. Verify
python -c "
from fy3_cloudmask.algorithm.native_backend import is_native_available, get_backend_info
print('Available:', is_native_available())
print(get_backend_info())
"

# 4. Set runtime environment
export OMP_NUM_THREADS=8  # adjust based on target CPU cores
export FY3_CODE_ROOT=/path/to/coeff/
```

### 3.3 Library Dependencies (Runtime)

The built `_cloudmask_native.cpython-*.so` dynamically links:

```
libhdf5.so.310        # HDF5 C library
libhdf5_fortran.so    # HDF5 Fortran bindings
libgfortran.so.5      # GFortran runtime
libquadmath.so.0      # Quad-precision math
libstdc++.so.6        # C++ standard library
libgomp.so.1          # OpenMP runtime
```

All must be available in `$LD_LIBRARY_PATH` or within the RPATH (`$ORIGIN/../lib`).

### 3.4 Portable Deployment (No Conda on Target)

If the target machine doesn't have conda:

1. Build on a machine with conda, using `--install`
2. Copy the entire `coeff/` directory and the built `.so` file
3. Ensure HDF5 and GFortran runtime libraries are installed on the target
4. Set `LD_LIBRARY_PATH` to include the directory containing `libhdf5*.so`

### 3.5 Verification Script

```bash
#!/bin/bash
# verify_deployment.sh — Run after deployment
set -e

echo "=== Python ==="
python3 --version

echo "=== Native engine ==="
PYTHONPATH=src python3 -c "
from fy3_cloudmask.algorithm.native_backend import is_native_available, get_backend_info
assert is_native_available(), 'Native engine not available!'
info = get_backend_info()
print(f'  Backend: {info[\"backend\"]}')
print(f'  Version: {info[\"version\"]}')
"

echo "=== HDF5 ==="
python3 -c "import h5py; print(f'  h5py {h5py.__version__}, HDF5 {h5py.version.hdf5_version}')"

echo "=== Coefficient files ==="
ls coeff/*.dat coeff/*.inc 2>/dev/null | wc -l | xargs -I{} echo "  {} threshold files found"

echo "=== OpenMP ==="
echo "  OMP_NUM_THREADS=${OMP_NUM_THREADS:-not set}"

echo ""
echo "Deployment verification passed."
```

---

## 4. Architecture Notes

### 4.1 Processing Flow

```
Python (config/IO/classification)
    → C++ pybind11 (cache-blocked transpose + OpenMP dispatch)
        → Fortran C API (thread-safe per-pixel processing)
            → Fortran core algorithm (spectral tests, confidence, bit encoding)
```

### 4.2 Memory Layout

- **Python/NumPy**: C-order (row-major), shape `(nElem, nLine, channels)`
- **Fortran**: Column-major, shape `(channels, nLine, nElem)`
- **Transpose**: Done once in C++ pybind11 layer with BLOCK=32 cache tiling
- **Output**: Fortran flat vectors reshaped with `numpy.reshape(order='F')`

### 4.3 Thread Safety

Fortran module variables are declared `!$omp threadprivate` in `cloudmask_data_arrays.f90`. Each OpenMP thread gets its own copy of all module-level state. The C API wrapper `process_swath_c` is the entry point called from C++ within `#pragma omp parallel`.

---

## 5. Known Limitations

1. **OMP scaling ceiling**: ~4x speedup at 8 threads, saturates beyond. The per-pixel Fortran algorithm has sequential dependencies that limit parallelism.
2. **Single-orbit memory**: ~2.1 GB peak per orbit. Processing multiple orbits concurrently requires proportional RAM.
3. **HDF5 version coupling**: Built against HDF5 1.14.x. Mismatched HDF5 versions at runtime will cause link errors.
4. **x86_64 only**: Build script targets `x86_64-conda-linux-gnu` compilers. ARM/other architectures require toolchain changes.
