# FY-3D MERSI-II Cloud Mask Retrieval System

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A modern Python implementation of the FYLAT FY-3D MERSI-II Cloud Mask Retrieval System (V3.1), ported from the original Fortran codebase.

## Overview

This system generates cloud mask products from FY-3D MERSI-II satellite data. It implements the MODIS MOD35-derived cloud detection algorithm with:

- **25-channel processing** (19 VIS + 6 IR)
- **Surface type classification** (land, water, coast, desert, snow, ice, polar)
- **Day/night processing paths** with different spectral tests
- **S-curve confidence interpolation** for cloud probability
- **48-bit testbits** and **80-bit QA bits** cloud mask encoding
- **4-level cloud confidence**: cloudy(0), probably cloudy(1), probably clear(2), confident clear(3)
- **NWP integration** (GFS 0.25° data for surface temperature, pressure, precipitable water)

## Key Features

- **Exact algorithm fidelity** - Preserves original Fortran algorithm logic
- **NumPy/Numba acceleration** - Performance-sensitive parts use JIT compilation
- **Modular design** - Clean separation of concerns
- **YAML configuration** - Human-readable config files
- **HDF5 output** - Standard satellite data format
- **CLI interface** - Command-line tool for batch processing
- **Comprehensive tests** - 34 unit tests covering all components

## Implementation Architecture

### Multi-Language Hybrid Design

```
┌─────────────────────────────────────────────────────┐
│                  Python 编排层                        │
│  CLI (click) → Pipeline → Config → Output Writer     │
│  负责：命令行交互 / 流程编排 / 配置管理 / HDF5 写入      │
└──────────────┬──────────────┬───────────────────────┘
               │              │
     ┌─────────▼────┐  ┌──────▼──────────────┐
     │  Numba 后端   │  │  C++/Fortran 后端    │
     │  (纯 Python)  │  │  (原生加速 30-100x)  │
     │              │  │                      │
     │ @njit JIT   │  │ pybind11 绑定        │
     │ 置信度/位操作 │  │ OpenMP 并行像素循环    │
     │ 空间分析     │  │ Fortran 核心算法      │
     └──────┬───────┘  └──────┬───────────────┘
            │                 │
            └────────┬────────┘
                     │
            ┌────────▼────────┐
            │  阈值配置 (YAML)  │
            │  790+ 物理参数    │
            └─────────────────┘
```

代码运行时会自动检测原生后端是否可用（`native_backend.py:30-37`），优先使用 C++/Fortran + OpenMP 版本，不可用时回退到纯 Python/Numba 版本。两条路径产出的结果**数学上等价**（同一套 Fortran 算法的两种翻译）。

### 核心部件详解

#### 部件 1：配置系统 (`config.py` + `config/default.yaml` + 阈值 YAML)

```
config/default.yaml          → FY3Config dataclass  → 所有模块读取
config/thresholds/*.yaml     → dict                 → 光谱检验阈值参数
```

YAML 配置文件替代原始 Fortran 的 namelist (`.nml`) 格式。`config.py:147` 的 `load_config()` 支持递归合并覆盖，CLI 参数可以覆盖配置文件中的任意字段。阈值文件 `mersi_ii3d_v8.yaml` 包含所有物理阈值参数，按地表类型/光照条件分层组织。

#### 部件 2：常量定义 (`constants.py`)

```
传感器维度 / 波段映射 / 波长-波数对照 / 位布局定义 / 物理常数
```

所有"魔法数字"集中管理。将 Fortran 中散落的数值提取为命名常量，保证算法精度（如 `BAND_064=2` 代表 0.64μm 通道在 pxldat 数组中的位置）。

#### 部件 3：地表分类器 (`surface_classifier.py`)

```
输入: lat, lon, elevation, lsf, sza, eco_type, snow_mask
  ↓
输出: PixelFlags (24个布尔标志)
  polar / land / water / coast / desert / day / night
  snow / ice / snglnt / hi_elev / antarctic / sh_ocean ...
```

每一个像素的"身份证"。后续光谱检验路径根据这组标志决定走哪条决策树分支。NDSI 雪检测（第 284-363 行）还包含 5 重假雪过滤（卷云/耀斑/冰云/水云/近红外亮度）。

#### 部件 4：光谱检验引擎 (`algorithm/tests/` 目录)

```
18 条检验路径:
  land_day.py      → land_day_standard / coast / desert / desert_coast
  ocean_day.py     → ocean_day
  land_nite.py     → land_nite
  ocean_nite.py    → ocean_nite
  polar_day.py     → polar_day_land / coast / desert / ocean / snow
  polar_nite.py    → polar_nite_land / ocean / snow
  snow_tests.py    → day_snow / nite_snow / antarctic_day
  restoral.py      → 8 个后处理修复
```

每条路径函数签名统一为：
```python
def xxx(pxldat, ..., thresholds, testbits, qa_bits) -> tuple[confdnc, nmtests, nbands]
```
- 读入 25 通道数据 + 阈值参数
- **原地修改** `testbits`（6 字节）和 `qa_bits`（10 字节）
- 返回 (置信度, 检验数, 使用波段数)

#### 部件 5：S 曲线置信度 (`confidence.py`)

```python
conf_test(val, locut, midpt, hicut, power) → [0.0, 1.0]
```

所有光谱检验的统一"打分函数"。4 个参数控制曲线形状：
- `locut` — 确定有云的下界
- `hicut` — 确定晴空的上界
- `midpt` — 50% 置信度转折点
- `power` — S 形弯曲程度

`encode_confidence()` 将连续置信度编码为 2-bit 离散掩膜：>0.99→确定晴空, >0.95→可能晴空, >0.66→可能多云, ≤0.66→有云。

#### 部件 6：位操作 (`bitops.py`)

```python
set_bit(testbits, bit_num)    # 置位
clear_bit(testbits, bit_num)  # 清零
check_bit(testbits, bit_num)  # 读取
fill_bit_pixel(...)           # 质量位组装
proc_path(...)                # 处理路径编码
convert_cloud_mask(...)       # 48-bit → 4-level 云掩膜
```

47 个命名位（0-47）编码了每一步检验的触发状态、通过的检验、最终掩膜值。10 字节 QA 位（80-bit）存储额外质量控制信息。

#### 部件 7：空间分析 (`spatial.py`)

```python
tview(vza, lat)                    # APOLLO BTD 11-12μm 阈值 (双线性查表)
get_regional_mean(data_3x3)        # 3x3 邻域均值
get_regional_std(data_3x3)         # 3x3 邻域标准差
check_reg_uniformity(data_3x3)     # 均匀性判定
```

云检测不只依赖单像素信息，还利用空间上下文。

#### 部件 8：云量统计 (`cloud_amount.py`)

```
5×5 像素盒 (1km×1km × 25 = 5km×5km) → 盒内 cloud_ratio → 0-100%
```

将 1km 云掩膜降分辨率到 5km 网格，统计每个 5×5 盒内的云像素占比。

#### 部件 9：输出写入器 (`output/writer.py`)

```python
write_cloud_mask()       → CLM 产品 (HDF5 Group: Cloud_Mask_1km)
write_cloud_amount()     → CLA 产品 (HDF5 Group: Cloud_Amount_5km)
write_combined_product() → 融合产品 (单文件包含两组 Group)
```

输出符合 CF-1.8 约定的 HDF5 格式，包含全局属性（机构/数据源/处理时间）。

### 数据流全景

```
┌──────────────┐   ┌──────────────┐   ┌──────────────┐
│ L1b HDF5     │   │ GEO HDF5     │   │ NWP GFS 二进制│
│ 25通道 DN值   │   │ 经纬度/角度    │   │ 地表温度/气压   │
└──────┬───────┘   └──────┬───────┘   └──────┬───────┘
       │                  │                  │
       ▼                  ▼                  ▼
  DN → 物理量          角度缩放            最近邻插值
  VIS: 反射率          sza/vza/glint       → 像素网格
  IR: Planck逆变换
       │                  │                  │
       └──────────────────┼──────────────────┘
                          │
              ┌───────────▼───────────┐
              │  run_cloud_mask_swath │
              │  逐像素循环 (4M pixels) │
              │                      │
              │  每个像素:             │
              │  ① 地表分类           │
              │  ② 分派光谱检验 (4组)   │
              │  ③ S曲线置信度         │
              │  ④ 后处理修复 (8种)    │
              │  ⑤ 位编码输出          │
              └───────────┬───────────┘
                          │
              ┌───────────▼───────────┐
              │   云量统计 (5km网格)    │
              │   5×5 盒 → 0-100%    │
              └───────────┬───────────┘
                          │
              ┌───────────▼───────────┐
              │   HDF5 输出文件        │
              │   Cloud_Mask_1km/     │
              │   Cloud_Amount_5km/   │
              └───────────────────────┘
```

### 核心设计决策

| 设计决策 | 原因 |
|---------|------|
| **Fortran 算法保持原样** | 确保与业务系统 MOD35 的一致性，所有阈值参数不变 |
| **Numba JIT 而不是纯 NumPy** | 逐像素循环有大量分支判断，无法向量化，JIT 是最佳折中 |
| **YAML 阈值替代 namelist** | 可读性、可版本管理、支持递归覆盖 |
| **双后端架构** | 开发/调试用纯 Python，生产用 C++/Fortran + OpenMP |
| **48+80 bit 位编码** | 兼容 MODIS MOD35 标准格式，每个检验都可追溯 |
| **原位修改 testbits/qa_bits** | 避免内存分配开销，与 Fortran 调用模式一致 |

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
│   ├── test_algorithms.py           # Algorithm component tests
│   └── test_pipeline_e2e.py         # E2E test with real FY-3D data
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

### Processing Pipeline

```
Input: FY-3D MERSI-II L1b (HDF5) + GEO (HDF5) + NWP (GRIB2)
  │
  ├─ 1. Read L1b data (25 channels: 19 VIS + 6 IR)
  ├─ 2. Read GEO data (lat, lon, SZA, VZA, land/sea flag)
  ├─ 3. Read NWP data (surface temp, pressure, precip water)
  │
  ▼
Per-pixel loop (2048 × 2000):
  │
  ├─ Step 1: DN → Physical units
  │    VIS: reflectance = (coef0 + coef1×DN + coef2×DN²) × 0.01 / cos(SZA) / esd²
  │    IR:  DN IS radiance (mW/m²/sr/cm⁻¹) for FY-3D
  │         BT = c2×ν / ln(c1×ν³ / (1e-5×R) + 1)   [Planck inversion]
  │         BT_corrected = (BT_raw - TCI) / TCS       [TCS/TCI correction]
  │
  ├─ Step 2: Classify surface type
  │    Land / Water / Coast / Desert / Snow / Ice / Polar (>60°)
  │    Day / Night (SZA > 85°)
  │
  ├─ Step 3: Dispatch spectral tests
  │    ├─ Daytime Land:  land_day_standard / land_day_coast / land_day_desert
  │    ├─ Daytime Ocean: ocean_day
  │    ├─ Nighttime:     land_nite / ocean_nite
  │    ├─ Polar:         polar_day_* / polar_nite_*
  │    └─ Snow/Ice:      day_snow / nite_snow / antarctic_day
  │
  │    Each test has 4 groups:
  │      Group 1: IR threshold + PFMFT/NFMFT + SST test
  │      Group 2: BTD tests (11-12μm, 11-4μm, 8-11μm)
  │      Group 3: Visible tests (0.64μm, 0.86μm, ratio)
  │      Group 4: NIR test (1.38μm thin cirrus)
  │
  │    Confidence = geometric mean of active group minimums
  │
  ├─ Step 4: Restoral tests (post-processing)
  │    1. chk_land_restoral     — restore clear if surface temp ≈ BT
  │    2. chk_coast_restoral    — restore clear for coastal pixels
  │    3. chk_sunglint_restoral — sun glint clear-sky restoral
  │    4. chk_shallow_water     — shallow water correction
  │    5. chk_spatial_var       — boost confidence for uniform scenes
  │    6. chk_thin_cirrus_ir    — flag thin cirrus (no confidence change)
  │    7. chk_shadow            — cloud shadow detection
  │    8. chk_cloud_adj         — non-cloud obstruction check (dust/smoke)
  │
  ├─ Step 5: Encode output
  │    confidence → 2-bit: (>0.99→3, >0.95→2, >0.66→1, ≤0.66→0)
  │    Assemble 48-bit testbits + 80-bit QA bits
  │
  ▼
Output: Cloud_Mask (48-bit), QA (80-bit), Cloud_Mask_Value (0-3), Confidence (0-1)
  │
  └─ Compute 5km Cloud Amount (5×5 pixel boxes → 0-100%)
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
# Run all unit tests
PYTHONPATH=src python -m pytest tests/test_algorithms.py tests/test_cloud_mask.py -v

# Run E2E test with real FY-3D data (~12 min)
PYTHONPATH=src python -m pytest tests/test_pipeline_e2e.py::TestPipelineE2E::test_full_orbit_with_nwp -v

# Run E2E test without NWP (uses dummy NWP)
PYTHONPATH=src python -m pytest tests/test_pipeline_e2e.py::TestPipelineE2E::test_full_orbit -v
```

### Test Types

- **Unit tests** (`test_algorithms.py`): confidence functions, bit operations, spatial analysis
- **Integration tests** (`test_cloud_mask.py`): pixel-level cloud mask processing
- **E2E tests** (`test_pipeline_e2e.py`): full orbit with real L1b/GEO/NWP data, output verification

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

- **Single orbit processing**: ~12 minutes (2048×2000 pixels, with NWP interpolation)
- **Memory usage**: ~4-8 GB for single orbit
- **Acceleration**: Numba JIT compilation for hot loops

## Migration from Fortran

This Python implementation is a port of the original Fortran system (retrieval_system_V3.1_cldmask):

| Component | Original Fortran | Python Implementation |
|-----------|------------------|----------------------|
| Lines of code | ~37,000 | ~5,000 |
| Configuration | Namelist (.nml) | YAML |
| Data I/O | HDF5 Fortran | h5py |
| NWP processing | wgrib2 shell scripts | xarray + cfgrib |
| Ancillary data | HDF4 Fortran | pyhdf |
| RTM (PFAAST/PLoD) | Included | Deferred (future) |

### What's Preserved

- All cloud detection algorithms (land_day, ocean_day, land_nite, ocean_nite, polar, snow)
- All threshold values (from mersi_ii3d_v8.yaml)
- Bit layout and encoding (48-bit testbits, 80-bit QA)
- Decision tree logic and surface classification rules
- Restoral tests (land, coast, sunglint, spatial var, thin cirrus, shadow)
- Planck function with FY-3D TCS/TCI correction
- APOLLO tview lookup table for thin cirrus BTD threshold

### What's Improved

- Modular, maintainable code structure
- Comprehensive unit tests + E2E tests with real data
- YAML configuration
- Python ecosystem integration
- Better error handling

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

### V3.1.0 (2026-05-20)

- Initial Python release
- Port from Fortran (retrieval_system_V3.1_cldmask)
- All cloud mask algorithms implemented
- E2E test with real FY-3D L1b/GEO/NWP data
- YAML configuration
- CLI interface
- HDF5 output
- Fixed IR Planck formula (correct c1=2hc², conversion 1e-5, FY-3D TCS/TCI)
- Fixed restoral tests (chk_cloud_adj, chk_sunglint_restoral matched to Fortran)

---

**Note**: This is a research codebase. For production use, please validate against reference datasets.
