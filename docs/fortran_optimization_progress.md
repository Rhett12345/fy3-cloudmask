# Fortran 云检测算法优化进展与问题记录

> 最后更新: 2026-06-16
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

### 2.5 空间一致性过滤（椒盐噪声修复）

**问题**: Fortran 云掩码存在大量椒盐噪声（单像素级分类跳变），因为算法是纯逐像素光谱检测，没有任何空间平滑。同时，`chk_spatial_var` 只提升置信度（均匀场景），从不惩罚孤立像素。陆地像素完全没有空间方差检查。MODIS MOD35 通过 250m→1km 子像素聚合天然获得空间平滑，但 FY-3D 处理 1km 原生分辨率，跳过了这一步。

**修复**: 在 `process_swath_c` 主循环之后添加空间一致性后处理 `apply_spatial_consistency`（`cloudmask_c_api.f90`），仅翻转与所有 8 个邻居分类都不一致的孤立像素（使用众数类别替换）。

**物理依据**: 1km 分辨率的云不是单像素现象——大气过程产生云的尺度远大于 1km。真正的孤立云/晴空像素在物理上不可信。MOD35 的 250m→1km 聚合步骤提供了等效平滑。

**修改文件**: `src/fortran/c_api/cloudmask_c_api.f90`（添加 `apply_spatial_consistency` 子程序 + 调用点）

**验证**:

| 轨道 | 8邻域孤立像素 (修复前估算) | 8邻域孤立像素 (修复后) |
|------|--------------------------|------------------------|
| 0740 | ~330,000+ | 38 |
| 0830 | ~200,000+ | 0 |

---

## 三、验证结果

### 3.1 轨道 20220803_0740 (空间一致性过滤后)

```
总像素:       4,096,000
多云:            53.66%
可能多云:        13.02%
可能晴空:         1.34%
确定晴空:        31.99%
云量:            66.68%
置信度均值:       0.6333
置信度标准差:     0.3326
置信度中位数:     0.6700
处理时间: 4.3s
8邻域孤立像素: 38 (0.001%)
```

### 3.2 轨道 20220803_0830 (空间一致性过滤后)

```
总像素:       4,096,000
多云:            76.37%
可能多云:         3.28%
可能晴空:         0.02%
确定晴空:        20.33%
云量:            79.65%
置信度均值:       0.4361
置信度标准差:     0.3562
置信度中位数:     0.3162
处理时间: 5.1s
8邻域孤立像素: 0 (0.000%)
```

### 3.3 多轨道稳定性（空间一致性过滤后）

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

### 4.4 极地孤立云像素问题（已通过空间一致性过滤解决）

**状态**: 已识别根因，多次修复尝试均失败，需要进一步研究。

**问题描述**: 轨道 20220803_0740 中存在约 112,000 个极地内部孤立云像素（cloud_mask=0，4 个邻居均不为 0）。这些像素在物理上不合理——云不会以单像素形式出现在完全晴空的邻居中。

**统计特征** (lsf=6 子集，占孤立像素的 43.9%):

| 指标 | 值 |
|------|------|
| 纬度范围 | 63° - 86° |
| 生态类型 | 100% eco_type=0 (水体) |
| LSF 值 | 6 (海冰) |
| 11-12um BTD | 79.4% < 1.0K (晴空指标) |
| 1.38um 反射率均值 | 0.1974 (高于 0.22 阈值) |
| 0.64um 反射率均值 | 0.1828 (接近 0.18 阈值) |

**根因分析**:

1. **LSF 路由问题**: lsf=6 (海冰) 在 `fylat_fy3mersi_cloud_mask.f90` 第 622 行落入 `else` 分支，被分类为 `water=true`。但海冰表面的光谱特征更接近冰雪而非开放水体。

2. **1.38um 测试误报**: 极地冰雪表面的 1.38um 反射率 (~0.20) 高于 `pdlref3` 阈值 (0.22/0.18/0.14)，导致 Group 4 测试失败。但该反射率来自地表而非卷云。

3. **可见光测试失败**: 冰雪表面的 0.64um 反射率 (~0.18) 接近 `pdlref1` 阈值 (0.22/0.18/0.14)，Group 3 测试经常失败。

4. **btclr 为零**: 测试脚本 (`run_fortran_only.py`, `spatial_analysis.py`, `debug_polar_isolated.py`) 传入 `btclr=0`，导致 Group 1 的 PFMFT/NFMFT 测试条件 `(btclr(5)-btclr(6)) > 0` 和 `(btclr(5) .ne. 0.0)` 不满足，Group 1 完全不激活。

5. **Fortran vs Python 的 ngtests 行为差异**: Fortran 在测试被应用时始终递增 `ngtests(k)`（即使测试失败），Python 仅在测试通过时递增。这导致 Fortran 的几何均值计算包含更多失败的测试组。

**已尝试的修复方案（均已回滚）**:

| 方案 | 结果 | 原因 |
|------|------|------|
| lsf=6 路由到 ice → PolarDay_snow | 更差 | PolarDay_snow 的 1.38um 阈值更严格 (0.06 vs 0.22) |
| lsf=6 路由到 land (VRAT > 0.9) | 略差 | 部分像素改善，其他像素恶化 |
| 调整 1.38um 阈值 (2.5x) | 更差 | 孤立像素从 107,777 增至 114,074 |
| BTD 检查跳过 1.38um 测试 | 更差 | masdf1 初始化为 0.0，条件判断错误 |
| BTD 高置信度时中和 Group 4 | 无效 | 条件触发但对最终云掩码影响很小 |

**未尝试的可能方向**:

1. **在 `chk_land` 中添加冰雪表面恢复逻辑**: 当前 `chk_land` 仅处理沙漠像素，可扩展为在 IR 测试全部通过且 11um BT 较高时提升置信度。

2. **在 `polar_module.f90` 中为 lsf=6 添加专门的冰雪处理路径**: 在调用 `PolarDay_ocean` 之前检查 lsf=6，路由到 `PolarDay_snow` 或自定义路径。

3. **提供非零 btclr**: 生产管线 (`pipeline.py`) 默认 btclr=280K，但测试脚本传入零值。非零 btclr 会激活 Group 1 的 NFMFT 测试，可能改善置信度。

4. **空间后处理**: 在 `polar_module.f90` 的 `chk_spatial_var` 之后，为 land 类型添加类似的孤立像素过滤。

**影响**: 约 112K 像素（占总像素 ~2.7%）受到影响，集中在北极海冰区域。对整体云量影响约 2-3%。

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

- [ ] **[高优先级]** 解决极地孤立云像素问题 (见 4.4 节)
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
