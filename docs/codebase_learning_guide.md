# FY-3D MERSI-II Cloud Mask 代码学习指南

## 第一层：理解系统在做什么

**核心问题**：给定一颗卫星（FY-3D）扫描地球的一轨数据，判断每个像素是"云"还是"晴空"。

**输出**：每个像素一个 4 级分类：
- `0` = 确定有云 (CM_CLOUDY)
- `1` = 可能有云 (CM_PROB_CLOUDY)
- `2` = 可能晴空 (CM_PROB_CLEAR)
- `3` = 确定晴空 (CM_CLEAR)

**数据来源**：FY-3D MERSI-II 卫星，25 个光谱通道（19 个可见光 + 6 个红外），分辨率 1km。

---

## 第二层：理解数据流（pipeline.py 是入口）

从 `pipeline.py:89` 的 `process_orbit()` 可以看到完整的 6 步流程：

```
输入: L1b HDF (25通道辐射数据) + GEO HDF (经纬度/角度)
  ↓ Step 1: 读卫星数据 (_read_satellite_data)
  ↓ Step 2: 读 NWP 数据（气象预报场）(_read_nwp_data)
  ↓ Step 3: 读辅助数据（海温、雪冰掩膜等）(_read_ancillary_data)
  ↓ Step 4: 运行云检测算法 (run_cloud_mask_swath) ← 核心
  ↓ Step 5: 计算云量（5km 网格）(compute_cloud_amount)
  ↓ Step 6: 写 HDF5 输出 (write_combined_product)
```

**双后端架构**：
- **Python/Numba 后端**（默认）：纯 Python + `@njit` JIT 编译，每轨约 12 分钟
- **C++/Fortran 后端**（可选）：OpenMP 并行，每轨约 10-30 秒，需要编译 `ext/`

运行时自动选择后端，逻辑在 `algorithm/native_backend.py`。

---

## 第三层：理解算法核心（最重要）

### 3.1 单像素处理流程

**入口**：`algorithm/cloud_mask.py:59` 的 `run_cloud_mask_pixel()`

每个像素的处理分 5 步：

```
Step 1: classify_pixel_surface()
        → 判断表面类型：陆地/海洋/海岸/沙漠/雪/冰
        → 判断光照条件：白天/夜晚
        → 判断地理位置：极地/南极/格林兰/新西兰
        → 输出: PixelFlags 对象 + 修正后的 pxldat

Step 2: proc_path()
        → 设置处理路径标记位到 testbits

Step 3: _dispatch_test()
        → 根据"表面类型 × 光照"路由到 18 种光谱测试之一
        → 输出: 置信度 (confdnc)、测试数量 (nmtests)、使用波段数 (nbands)

Step 4: 8 种 restoral 测试（后处理修正）
        → chk_land_restoral: 陆地晴空恢复
        → chk_coast_restoral: 海岸晴空恢复
        → chk_sunglint_restoral: 太阳耀光恢复
        → chk_shallow_water: 浅水区恢复
        → chk_spatial_var: 空间变异性检查
        → chk_thin_cirrus_ir: 薄卷云 IR 检查
        → chk_shadow: 阴影检测
        → chk_cloud_adj: 云邻接检查

Step 5: encode_confidence()
        → S曲线置信度 → 2位编码 → 4级云掩膜
        → 组装 testbits + qa_bits
```

### 3.2 测试路由表

`_dispatch_test()` (cloud_mask.py:275) 根据 surface type × lighting 选择不同的测试函数：

| 表面类型 | 白天 | 夜晚 |
|---------|------|------|
| **陆地（标准）** | `land_day_standard` | `land_nite` |
| **陆地（海岸）** | `land_day_coast` | `land_nite` |
| **陆地（沙漠）** | `land_day_desert` | `land_nite` |
| **陆地（沙漠+海岸）** | `land_day_desert_coast` | `land_nite` |
| **海洋** | `ocean_day` | `ocean_nite` |
| **雪/冰** | `day_snow` | `nite_snow` |
| **极地陆地** | `polar_day_land` | `polar_nite_land` |
| **极地海岸** | `polar_day_coast` | `polar_nite_land` |
| **极地沙漠** | `polar_day_desert` | `polar_nite_land` |
| **极地海洋** | `polar_day_ocean` | `polar_nite_ocean` |
| **极地雪** | `polar_day_snow` | `polar_nite_snow` |
| **南极** | `antarctic_day` | (走极地路径) |

测试函数位于 `algorithm/tests/` 目录下，每个文件对应一类表面/光照组合。

---

## 第四层：理解 S 曲线置信度（confidence.py）

### 4.1 conf_test() 函数

这是核心数学函数，位于 `confidence.py:15`。

**输入**：一个光谱测试值 `val` + 4 个阈值参数 `(locut, midpt, hicut, power)`

**输出**：0.0 ~ 1.0 的置信度值

```
val <= locut → 0.0（确定有云）
val >= hicut → 1.0（确定晴空）
中间区域 → S 形曲线插值

当 hicut < locut 时曲线翻转（值越小越晴朗）
```

**S 曲线形状**：
```
置信度
1.0 |                    ___________
    |                  /
    |                /
0.5 |              /
    |            /
    |__________/
0.0 |________________________
    locut    midpt    hicut    测试值
```

### 4.2 置信度编码

`encode_confidence()` 将连续置信度映射到 2 位离散值：

| 置信度范围 | 编码 | 含义 |
|-----------|------|------|
| > 0.99 | (1, 1) | 确定晴空 |
| > 0.95 | (0, 1) | 可能晴空 |
| > 0.66 | (1, 0) | 可能有云 |
| <= 0.66 | (0, 0) | 确定有云 |

### 4.3 阈值来源

每个光谱测试的 4 个阈值参数定义在 `config/thresholds/mersi_ii3d_v8.yaml`（790+ 个参数），按表面类型和光照条件组织。

---

## 第五层：理解表面分类（surface_classifier.py）

`classify_pixel_surface()` 函数位于 `surface_classifier.py:90`，是算法的第一步。

### 5.1 分类维度

1. **光照条件**：`sza > 85°` → 夜晚，否则白天
2. **地理位置**：
   - 极地：`|lat| > 60°`
   - 南极：`lat < -60°`
   - 格林兰：特殊经纬度范围 + 生态类型
3. **表面类型**（基于 land/sea flag `lsf`）：
   - `lsf=1,4` → 陆地
   - `lsf=0` → 海洋
   - `lsf=2` → 海岸
   - `lsf=3` → 浅湖
4. **沙漠检测**：IGBP 生态类型 + 区域规则（非洲/欧亚/澳大利亚）
5. **雪/冰检测**：
   - 辅助数据雪掩膜 + NDSI 指数
   - NDSI = (R_0.55 - R_2.13) / (R_0.55 + R_2.13)
   - 多重假雪过滤：薄卷云、太阳耀光、冰云、水云

### 5.2 PixelFlags 数据类

包含 30+ 个布尔标记，记录像素的所有表面属性：
- `land`, `water`, `coast`, `desert`
- `day`, `night`, `polar`, `antarctic`
- `snow`, `ice`, `snglnt`（太阳耀光）
- `hi_elev`（高海拔）、`bad_value`（坏数据）
- 等等

---

## 第六层：理解光谱测试（algorithm/tests/）

每个测试函数执行一组光谱判据，输出置信度。

### 6.1 测试函数签名（以 land_day_standard 为例）

```python
def land_day_standard(
    pxldat,      # 25 通道像素数据
    bt_clr,      # 7 通道晴空亮温（来自 RTM/NWP）
    vza,         # 卫星天顶角
    is_cold_sfc, # 是否冷表面
    hi_elev,     # 是否高海拔
    thresholds,  # 阈值字典
    testbits,    # 6 字节 testbits（原位修改）
    qa_bits,     # 10 字节 QA bits（原位修改）
) -> (confdnc, nmtests, nbands)
```

### 6.2 典型测试项目

每个测试路径包含多个光谱测试，例如陆地白天路径可能包含：

1. **可见光反射率测试**：R_0.64 > 阈值 → 有云
2. **红外亮温测试**：BT_11um < 阈值 → 有云
3. **双通道差值测试**：BT_11um - BT_12um (分裂窗) → 有云信号
4. **1.38μm 卷云测试**：R_1.38 > 阈值 → 薄卷云
5. **NDSI 雪测试**：排除雪地假阳性
6. **GEMI 植被指数**：排除植被干扰
7. **VRAT 比值测试**：R_0.86 / R_0.64 → 区分云和地表

每个测试通过 `conf_test()` 计算置信度，取所有测试的最小值作为该路径的最终置信度。

---

## 第七层：理解位操作（bitops.py）

### 7.1 testbits 布局（6 字节 = 48 位）

```
Byte 0 (bits 0-7):
  bit 0: processed（已处理）
  bit 1: conf_lsb（置信度低位）
  bit 2: conf_msb（置信度高位）
  bit 3: day（白天标记）
  bit 4: no_sunglint（无太阳耀光）
  bit 5: no_snow_ice（无雪冰）
  bit 6: coast（海岸）
  bit 7: desert（沙漠）

Byte 1 (bits 8-15):
  bit 8: nco（非云目标）
  bit 9: thin_cirrus_solar（薄卷云-太阳通道）
  bit 10: shadow（阴影）
  bit 11: thin_cirrus_ir（薄卷云-IR）
  bit 12: cloud_adj（云邻接）
  bit 14: pfmft（正假薄卷云标记）
  bit 15: nfmft（负假薄卷云标记）

Byte 2 (bits 16-23):
  bit 16: nir_138（1.38μm 近红外）
  bit 18: btd_11_12（11-12μm 差值）
  bit 19: btd_11_4（11-4μm 差值）
  bit 20: ref_064（0.64μm 反射率）
  bit 21: gemi（GEMI 植被指数）

Byte 3 (bits 24-31):
  bit 24: temporal（时间标记）
  bit 26: land_restoral（陆地恢复）
  bit 28: suspended_dust（悬浮尘埃）
```

### 7.2 qa_bits 布局（10 字节 = 80 位）

QA bits 包含更详细的质量信息，用于后续产品验证和调试。

### 7.3 convert_cloud_mask()

将 testbits 中的置信度位解码为 0-3 的云掩膜值。

---

## 第八层：理解空间分析（spatial.py）

### 8.1 3x3 邻域操作

空间分析使用 3x3 像素窗口，主要用于：

1. **`get_regional_mean()`**：计算 3x3 区域均值
2. **`get_regional_std()`**：计算 3x3 区域标准差
3. **`get_regional_diff()`**：计算中心像素与区域均值的差值

### 8.2 应用场景

- **空间变异性检查** (`chk_spatial_var`)：海洋上的均匀云层 → 低标准差 → 可能是云
- **阴影检测** (`chk_shadow`)：可见光暗 + 红外冷 → 可能是云阴影
- **云邻接检查** (`chk_cloud_adj`)：周围像素是云 → 当前像素可能是云边缘

---

## 第九层：理解配置系统

### 9.1 配置层次

```
config/default.yaml              → 运行时路径、处理选项
config/sensors/fy3d_mersi_ii.yaml → 波段波长、传感器几何
config/thresholds/mersi_ii3d_v8.yaml → 790+ 算法阈值（按表面类型/光照分组）
```

### 9.2 FY3Config 数据类

`config.py` 中定义的 `FY3Config` 数据类，支持 YAML 加载和递归合并覆盖。

---

## 第十层：理解测试框架

### 10.1 测试层次

```
tests/test_algorithms.py  → 单元测试（合成数据，快速）
tests/test_cloud_mask.py  → 集成测试（合成数据，快速）
tests/test_pipeline_e2e.py → 端到端测试（需要真实 FY-3D 数据，~12 分钟）
```

### 10.2 运行测试

```bash
# 单元/集成测试（无需真实数据）
PYTHONPATH=src python -m pytest tests/test_algorithms.py tests/test_cloud_mask.py -v

# 端到端测试（需要真实数据）
PYTHONPATH=src python -m pytest tests/test_pipeline_e2e.py -v
```

---

## 建议的阅读顺序

1. **`pipeline.py`** — 先看 `process_orbit()` 理解全流程
2. **`constants.py`** — 理解 25 个通道、位位置、物理常数
3. **`algorithm/cloud_mask.py`** — 重点看 `run_cloud_mask_pixel()` 和 `_dispatch_test()`
4. **`algorithm/surface_classifier.py`** — 理解表面分类逻辑
5. **`algorithm/confidence.py`** — 理解 S 曲线
6. **`algorithm/tests/land_day.py`** — 选一个具体的测试路径看细节
7. **`algorithm/bitops.py`** — 理解输出编码
8. **`algorithm/spatial.py`** — 空间分析（3x3 邻域）
9. **`config/thresholds/mersi_ii3d_v8.yaml`** — 浏览阈值配置
10. **`ext/cloudmask_pybind.cpp`** — 了解 C++/Fortran 后端接口（可选）

---

## 关键设计决策

1. **双后端架构**：Python 便于开发调试，C++/Fortran 用于生产性能
2. **S 曲线置信度**：比硬阈值更平滑，减少边界效应
3. **表面类型路由**：不同地表需要完全不同的检测策略
4. **位操作编码**：紧凑存储所有测试结果，便于调试和追溯
5. **Fortran 移植保真**：所有常量和逻辑严格对应原始 Fortran 代码
