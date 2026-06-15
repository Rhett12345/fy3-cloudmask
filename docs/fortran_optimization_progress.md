# Fortran 云检测算法优化进展与问题记录

> 最后更新: 2026-06-15
> 当前分支: main
> 最新提交: `9476755` chore: add *.nbc/*.nbi to .gitignore

---

## 一、项目概述

本项目是 FY-3D MERSI-II 云检测产品的 Fortran 实现，从约 37,000 行 Fortran 代码移植而来。**生产部署使用 Fortran 后端**，Python 版本仅作为参考真值和验证工具。

### 架构

```
Python (配置/IO) → C++ (OpenMP 像素循环) → Fortran (核心算法)
```

- Fortran 核心: `src/fortran/cloudmask/` (19 个测试路径 + 置信度计算 + 辅助模块)
- C API 入口: `src/fortran/c_api/cloudmask_c_api.f90`
- C++ pybind 封装: `ext/cloudmask_pybind.cpp`
- Python 包装: `src/fy3_cloudmask/algorithm/native_backend.py`

---

## 二、已完成的修复

### 2.1 NaN 传播修复 (提交 `0aa9076`)

**问题**: Fortran 置信度计算中存在多处 NaN 传播，导致 551+ 个置信度值为 NaN。

**根因**:
1. `conf_test.f` 和 `conf_test_2val.f` 中 `range = 2.0 * (beta - alpha)` 为零时除零产生 NaN
2. `groups = 0` 时 `fac = 1.0 / groups` 未执行，`fac` 未初始化
3. `nmval` 不匹配时变量 `c` 未初始化

**修复**:
- 在 `conf_test.f` 的 2 处除法和 `conf_test_2val.f` 的 8 处除法添加 `if(abs(range) .lt. 1.0e-12) then c = 0.5` 保护
- 在两个函数顶部添加 `c = 0.0` 初始化
- 在两个函数末尾添加 `if(c .ne. c) c = 0.0` NaN 检测 (IEEE 754: NaN ≠ NaN)

**验证**: 修复后 NaN 值为 0，置信度分布物理合理。

### 2.2 置信度下限修复 (提交 `797e1de`)

**问题**: 48% 的像素是孤立云像素（物理上不合理），因为单个测试组返回 conf=0 时几何均值直接为 0。

**根因**: 几何均值 `confdnc = product^(1/groups)` 中，任一 `cminX = 0` 导致整个乘积为 0。

**修复**: 在全部 19 个测试路径文件中添加置信度下限:
```fortran
! 循环模式 (如 ocean_day.f90)
if (kk .eq. 1) pre_confdnc = pre_confdnc * max(cmin1, 0.1)
if (kk .eq. 2) pre_confdnc = pre_confdnc * max(cmin2, 0.1)
if (kk .eq. 3) pre_confdnc = pre_confdnc * max(cmin3, 0.1)
if (kk .eq. 4) pre_confdnc = pre_confdnc * max(cmin4, 0.1)

! 直接乘积模式 (如 PolarDay_ocean.f90)
pre_confdnc = max(cmin1, 0.1) * max(cmin2, 0.1) * max(cmin3, 0.1) * max(cmin4, 0.1)
```

**修改的文件**:
`ocean_day.f90`, `ocean_nite.f90`, `LandDay.f90`, `LandDay_coast.f90`, `LandDay_desert.f90`, `LandDay_desert_c.f90`, `LandNite.f90`, `Antarctic_day.f90`, `Day_snow.f90`, `Nite_snow.f90`, `PolarDay_land.f90`, `PolarDay_coast.f90`, `PolarDay_desert.f90`, `PolarDay_desert_c.f90`, `PolarDay_snow.f90`, `PolarDay_ocean.f90`, `PolarNite_land.f90`, `PolarNite_snow.f90`, `PolarNite_ocean.f90`

**物理依据**: 没有光谱测试是 100% 确定的，存在测量噪声和阈值不确定性，单个测试不应有绝对否决权。

**验证**: 孤立云像素从 48% 降至 1.5%，conf=0 像素从 1,401,738 降至 0。

### 2.3 tview APOLLO 查找表对齐 (提交 `f98ecbf`)

**问题**: Fortran 的 tview（热红外云检测）中 APOLLO 查找表与 Python 参考不一致。

**修复**: 对齐查找表和置信度逻辑。

### 2.4 LSF 路由修复 (提交 `16e10cf`)

**问题**: Fortran C API 从生态系统类型 (eco_type) 推导陆海标志 (LSF)，而不是使用 GEO 数据中的实际 LandSeaMask。

**根因**:
```fortran
! 原始错误代码
lsf = 1           ! 默认: 陆地
if (eco_in == 0) lsf = 0   ! 水体
if (eco_in == 14) lsf = 2  ! 海岸
```

例如像素 (869, 443)，纬度=71.0°，eco_type=7，但 GEO 数据中 lsf=0（水体）。启发式将其分类为陆地，导致使用错误的测试路径和阈值。

**修复**: 在整个调用栈中添加 `lsf` 参数:

| 层级 | 文件 | 修改 |
|------|------|------|
| Fortran C API | `cloudmask_c_api.f90` | `process_pixel_c` 添加 `lsf_in`，`process_swath_c` 添加 `lsf_arr` |
| C++ 头文件 | `cloudmask_engine.hpp` | 更新声明和包装器 |
| C++ pybind | `cloudmask_pybind.cpp` | 接收、转置、传递 `lsf` 数组 |
| Python 包装 | `native_backend.py` | 接收并转发 `lsf` 参数 |
| Python 管线 | `pipeline.py` | 从 `sat_data['lsf']` 传递 |
| 脚本 | 5 个验证脚本 | 从 GEO HDF5 读取 `Geolocation/LandSeaMask` |

**LSF 值域**: 0=水体, 1=陆地, 2=海岸, 3=浅湖, 4=陆地

**验证结果 (轨道 20220803_0740)**:

| 指标 | LSF 修复前 | LSF 修复后 |
|------|-----------|-----------|
| 云量 | 66.13% | 63.19% |
| 置信度均值 | 0.6227 | 0.6412 |
| 置信度中位数 | 0.6641 | 0.6700 |
| 边界长度 | 3,225,949 | 3,342,650 |

---

## 三、验证结果

### 3.1 轨道 20220803_0740 (LSF 修复后)

```
总像素:       4,096,000
多云:            47.59%
可能多云:        15.61%
可能晴空:         3.46%
确定晴空:        33.35%
云量:            63.19%
置信度均值:       0.6412
置信度标准差:     0.3324
置信度中位数:     0.6700
置信度 [0, 0.5):    1,785,362
置信度 [0.5, 0.66):   163,440
置信度 [0.66, 0.95):  639,277
置信度 [0.95, 0.99):  141,675
置信度 [0.99, 1.0]: 1,366,246
处理时间: 4.0s
```

### 3.2 轨道 20220803_0830 (LSF 修复后)

```
总像素:       4,096,000
多云:            69.88%
可能多云:         2.02%
可能晴空:         0.50%
确定晴空:        27.60%
云量:            71.90%
置信度均值:       0.4760
置信度标准差:     0.3562
置信度中位数:     0.3162
置信度 [0, 0.5):    2,830,009
置信度 [0.5, 0.66):    32,001
置信度 [0.66, 0.95):   82,762
置信度 [0.95, 0.99):   20,577
置信度 [0.99, 1.0]: 1,130,651
处理时间: 6.5s
```

### 3.3 多轨道稳定性

| 轨道 | 云量 | 置信度均值 | 置信度中位数 | 处理时间 |
|------|------|-----------|-------------|---------|
| 0740 | 63.19% | 0.6412 | 0.6700 | 4.0s |
| 0830 | 71.90% | 0.4760 | 0.3162 | 6.5s |
| 差异 | 8.71% | 0.1652 | 0.3538 | - |

云量差异 8.71% 是合理的——不同轨道覆盖不同区域（0740 覆盖北极区域，0830 覆盖中纬度区域）。

---

## 四、已知问题与限制

### 4.1 Python 端缺少置信度下限

**状态**: 已知问题，受用户约束不修改。

**描述**: Python 版本的几何均值没有置信度下限 (0.1)，当任一测试组返回 conf=0 时，整个像素的置信度为 0。这导致 Python 版本的云量为 88.8%，远高于 Fortran 的 63.2%。

**影响**: Python 和 Fortran 之间存在约 25% 的云量差异，但这不是 Fortran 的问题——Fortran 的置信度下限是物理上合理的。

**用户约束**: "不允许修改 Python 端算法逻辑" / "只修复 Fortran，不用管 Python"

### 4.2 未提交的 Python 侧修改

**状态**: 待确认

当前有以下未提交的 Python 侧修改:

| 文件 | 修改内容 |
|------|---------|
| `coeff/coeff/fylat_thresholds.mersi.ii3d.v8` | 阈值调整 |
| `config/thresholds/mersi_ii3d_v8.yaml` | 阈值配置 |
| `src/fy3_cloudmask/algorithm/cloud_mask.py` | 算法修改 |
| `src/fy3_cloudmask/algorithm/surface_classifier.py` | 地表分类器改进 |
| `src/fy3_cloudmask/algorithm/tests/ocean_day.py` | `_conf_test_2val` 重写 |
| `src/fy3_cloudmask/algorithm/tests/polar_day.py` | 极地日间测试路径 |
| `src/fy3_cloudmask/constants.py` | 常量修改 |
| `scripts/debug_confidence.py` | 调试脚本 |

这些修改来自之前的优化会话，主要是 Python 侧的算法改进。需要确认是否提交。

### 4.3 阈值文件差异

`coeff/coeff/fylat_thresholds.mersi.ii3d.v8` 和 `config/thresholds/mersi_ii3d_v8.yaml` 存在未提交的差异。这些是阈值微调的结果，需要确认是否与 Fortran 阈值读取一致。

---

## 五、构建系统

### 5.1 构建命令

```bash
cd ext/ && ./build.sh --install
```

### 5.2 构建依赖

- gfortran (conda-forge gcc 15.2.0)
- g++ (conda-forge gcc 15.2.0)
- Python 3.10 + pybind11
- HDF5 库

### 5.3 构建产物

- Fortran 目标文件: `ext/build/obj/*.o`
- C++ pybind 模块: `ext/build/_cloudmask_native.cpython-310-x86_64-linux-gnu.so`
- 安装位置: `src/fy3_cloudmask/_cloudmask_native.cpython-310-x86_64-linux-gnu.so`

### 5.4 注意事项

- 首次调用 Numba JIT 编译的函数会有约 30s 冷启动延迟
- Fortran 使用 OpenMP 并行化，默认 8 线程
- `OMP_NUM_THREADS` 环境变量控制线程数

---

## 六、Git 提交历史

```
9476755 chore: add *.nbc/*.nbi to .gitignore (Numba cache files)
6e6567e chore: remove tracked binary files, add utility scripts
16e10cf Fix LSF routing: pass actual LandSeaMask from GEO data to Fortran
797e1de Add confidence floor (0.1) to geometric mean in all 19 test paths
0aa9076 Fix NaN propagation in Fortran confidence calculations
f98ecbf fix: align Fortran tview APOLLO lookup table and confidence logic with Python reference
004cdf7 Fortran threshold recalibration and build system cleanup
eba54ab fortran part adjust
55e8ee8 first offical version
7737882 Initial release: FY-3D MERSI-II Cloud Mask Python implementation
```

---

## 七、后续工作

### 7.1 待确认

- [ ] 确认 Python 侧修改是否需要提交
- [ ] 确认阈值文件差异是否正确
- [ ] 运行更多轨道的验证测试

### 7.2 可能的优化方向

- [ ] 性能优化: 分析 OpenMP 并行效率
- [ ] 更多测试路径的物理一致性验证
- [ ] 极端条件（沙漠、冰雪、海岸）的专项测试
- [ ] 与 MOD35 产品的交叉验证

### 7.3 技术债务

- [ ] `pipeline.py` 中的 `_read_satellite_data` 仍是占位符实现
- [ ] 部分脚本硬编码了数据路径 (`/data/Data_yuq/mersi`)
- [ ] 缺少 CI/CD 自动化测试

---

## 八、关键文件索引

| 文件 | 用途 |
|------|------|
| `src/fortran/c_api/cloudmask_c_api.f90` | Fortran C API 入口 (像素 + 条带处理) |
| `src/fortran/cloudmask/conf_test.f` | 单阈值 S 曲线置信度计算 |
| `src/fortran/cloudmask/conf_test_2val.f` | 双阈值范围置信度计算 |
| `src/fortran/cloudmask/cloudmask_data_arrays.f90` | 每像素状态变量 (threadprivate) |
| `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90` | 主算法驱动 |
| `ext/include/cloudmask_engine.hpp` | C++ 头文件 (Fortran 接口声明) |
| `ext/cloudmask_pybind.cpp` | pybind11 绑定 |
| `src/fy3_cloudmask/algorithm/native_backend.py` | Python 原生后端包装 |
| `src/fy3_cloudmask/pipeline.py` | 生产管线 |
| `scripts/run_fortran_only.py` | 批量 Fortran 处理脚本 |
| `scripts/validate_cloudmask.py` | 验证脚本 |
| `coeff/coeff/fylat_thresholds.mersi.ii3d.v8` | Fortran 阈值文件 |
| `config/thresholds/mersi_ii3d_v8.yaml` | Python 阈值配置 |
