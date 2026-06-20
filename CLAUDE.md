# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FY-3D MERSI-II Cloud Mask Retrieval System (v3.2.0) — a Python port of a ~37,000-line Fortran codebase that generates cloud mask products from FY-3D MERSI-II satellite data. Implements the MODIS MOD35-derived cloud detection algorithm with 25-channel processing, surface type classification, day/night paths, S-curve confidence, and dual backend architecture.

## Common Commands

```bash
# Install in dev mode
pip install -e .

# Run unit/integration tests
PYTHONPATH=src python -m pytest tests/test_algorithms.py tests/test_cloud_mask.py -v

# Run E2E tests (requires real FY-3D data, ~12 min)
PYTHONPATH=src python -m pytest tests/test_pipeline_e2e.py -v

# CLI: single orbit
python -m fy3_cloudmask process --config config/default.yaml --l1b data/L1b.HDF --geo data/GEO.HDF --output output/

# CLI: batch
python -m fy3_cloudmask batch --config config/default.yaml --start 2023-01-01 --end 2023-01-31 --data data/ --output output/ --workers 4

# Build native C++/Fortran backend (optional, 30-100x speedup)
cd ext/ && ./build.sh --install
```

Note: Tests require `PYTHONPATH=src` since the package uses a `src/` layout. pytest is configured in `pyproject.toml` with `testpaths = ["tests"]` and `addopts = "-v --tb=short"`.

## Architecture

### Dual Backend Design

The system has two mathematically-equivalent backends, auto-selected at runtime (`algorithm/native_backend.py`):

1. **Python/Numba** (default): Pure Python with `@njit` JIT for confidence, bitops, spatial analysis. ~12 min per orbit.
2. **C++/Fortran**: OpenMP-parallelized via pybind11. ~10-30s per orbit. Built from `ext/`.

### Algorithm Flow (per pixel)

```
L1b HDF5 + GEO HDF5
    ↓
CloudMaskPipeline (pipeline.py)
    ↓
SurfaceClassifier (surface_classifier.py)
  → classifies: land/water/coast/desert/snow/ice/polar + day/night
    ↓
_dispatch_test (cloud_mask.py)
  → routes to one of 18 spectral test paths in algorithm/tests/
    ↓
S-curve confidence (confidence.py)
  → maps test values → [0,1] confidence via threshold interpolation
    ↓
8 restoral tests (tests/restoral.py)
  → post-processing corrections (sunglint, cirrus, shadow, etc.)
    ↓
Bit encoding (bitops.py)
  → 48-bit testbits + 80-bit QA bits → 4-level cloud mask
    ↓
HDF5 output (output/writer.py) + cloud amount on 5km grid (output/cloud_amount.py)
```

### Key Modules

- `config.py` — `FY3Config` dataclass, loads YAML, supports recursive merge overrides
- `constants.py` — All band indices, bit positions, physical constants, sensor dimensions
- `pipeline.py` — `CloudMaskPipeline` orchestrates: data read → NWP interpolation → algorithm → output
- `algorithm/cloud_mask.py` — Main driver: `run_cloud_mask_pixel()` and `run_cloud_mask_swath()`
- `algorithm/tests/` — 18 test path modules organized by surface type + lighting (e.g., `land_day.py`, `ocean_nite.py`)

### Configuration Hierarchy

```
config/default.yaml              → runtime paths, processing options
config/sensors/fy3d_mersi_ii.yaml → band wavelengths, sensor geometry
config/thresholds/mersi_ii3d_v8.yaml → 790+ algorithm thresholds (by surface type / lighting)
```

Thresholds were converted from Fortran namelist format via `scripts/convert_thresholds.py`.

## Important Caveats

- **Known migration bugs**: See `docs/code_review/migration_bug_report.md` — 12 bugs identified in the Fortran-to-Python migration, including 4 critical (wrong coastal/desert thresholds, shadow detection missing preconditions). Fix these before trusting production output.
- **Source layout**: Package source is under `src/fy3_cloudmask/`, not directly in the root. Always set `PYTHONPATH=src` when running without `pip install -e .`.
- **Numba cold start**: First invocation of JIT-compiled functions is slow (~30s). Subsequent calls are fast.
- **E2E tests need real data**: `test_pipeline_e2e.py` expects FY-3D L1b/GEO HDF files in `tests/reference_data/`. Unit and integration tests use synthetic data and run without external files.
- **Constants are authoritative**: All band mappings, bit positions, and magic numbers live in `constants.py`. Do not hardcode numeric indices — import from constants.

代码的主要功能依赖是fortran，而不是python，我之后拿去部署的工程化代码也是使用fortran。
一定要保障fortran的可用性

## 行为规则

- **每次输出前必须以"打报告"开头**，然后再输出正文内容。

## 调参数据

- **小规模调参的测试数据日期固定为 2022-08-03（2022年8月3日）**。该日期的 FY-3D MERSI-II 数据作为参数调整的标准验证集。