# FY-3D 云检测代码 Bug 完整汇总

**当前版本**: v3.6.0
**GitHub 状态**: main 分支，未推送
**核心目标**: 消除椒盐噪声，提升与 MYD35 的一致性
**测试日期**: 2022-08-03（MERSI + NWP binary + MYD35 三者齐全）
**最后更新**: 2026-06-23（v3.6.0 数组布局修复完成）

---

## 版本历史

| 版本 | 主要变更 |
|------|----------|
| v3.3.3 | 基线版本 |
| v3.3.4~v3.3.7 | Goal 1: 恢复原始 Fortran 计算逻辑 |
| v3.4.0 | 移除 Python 算法后端，仅保留 Fortran |
| v3.4.1 | 禁用 pfmft/nfmft 测试 (btclr 缺失) |
| v3.4.2 | 禁用 APOLLO 11-12um 测试 (MERSI-II 不适用) |
| v3.4.3 | 修复 3 个关键 bug + 椒盐噪声缓解 |
| v3.4.4 | 修复 10 个 Fortran 计算 bug |
| v3.4.5 | 修复 relaz 和沙漠检测 |
| v3.5.0 | 重新启用 APOLLO + PFMFT 测试 |
| v3.5.1 | 根因修复: NWP 加载 + ocean PFMFT conf_test 恢复 + btclr sfctmp 估算 |
| v3.5.2 | 修正 IR 波数 (从 HDF5 Effect_Center_WaveLength 读取) + TBB 线性修正 |
| v3.5.3 | 清理死代码 planck_module.f90 |
| v3.5.4 | btclr 按表面类型区分 (Ocean BTD_clr=0.5K, Land BTD_clr=1.9K) |
| v3.6.0 | **修复 pybind 数组布局转换下标写反 (根因0)** + 两套定标系数反演 |

---

## 当前性能 (v3.6.0 vs MYD35, 2022-08-03, 0740 orbit, onboard)

```
  Cloud fraction: 95.2% (FY3D) — 阵列修复后整体云量略有下降
  孤立椒盐像素: 0.42% (v3.5.4 预计 >> 1%)
  连通域中位数: 3 像素 (云), 3 像素 (晴)
  Confidence >=0.66: 51.8%, >=0.95: 5.9%
  转置前向/往返验证: 全部 PASS
  Recal vs Onboard: CF 差异 < 1%
```

> **注意**: MYD35 对比尚未完成（数据量过大需优化匹配算法）。
> v3.5.4 的性能指标（Accuracy 49.98%, HSS -0.010）已不具参考价值，
> 因其基于错位数据。v3.6.0 需重新评估全量 MYD35 对比。

---

## 椒盐噪声根因权重总览

| 优先级 | 根因 | 严重度 | 状态 | 影响范围 |
|--------|------|--------|------|----------|
| 1 | pybind 数组布局转换下标写反 | **CRITICAL** | ✅ 已修复 (v3.6.0) | 全部输入场空间错位，椒盐噪声结构性根因 |
| 2 | PFMFT/APOLLO 阈值不适配 MERSI-II | HIGH | 未修复 | 系统性偏云 + 阈值附近随机跳变 |
| 3 | 后处理滤波过弱且 bit/QA 不同步 | MEDIUM-HIGH | 未修复 | 无法消除 2-3 像元斑块，bit 不一致 |
| 4 | snow_mask 全零 | HIGH | 未修复 | 夜间/高纬/冬季雪冰误判 |
| 5 | 正式 pipeline relaz 全零占位 | MEDIUM | 未修复 | 白天海洋几何相关路径 |
| 6 | btclr 经验估算（非 RTM） | MEDIUM | 未修复 | PFMFT/APOLLO 物理一致性 |
| 7 | 原始 Fortran apply_spatial_filter() 不生效 | LOW | 已确认 | native 后端不走该路径 |

---

## 根因 0 (NEW)：pybind 数组布局转换下标写反 [CRITICAL] -- ✅ 已修复 (v3.6.0)

**文件**: ext/cloudmask_pybind.cpp
**修复日期**: 2026-06-23
**影响版本**: 全部使用 native 后端的版本（v3.4.0 至今）

### 问题描述

Python 传入数组 shape 为 `(nElem=2048, nLine=2000)`，C-order (row-major) 连续内存。
Fortran 侧期望 `(nElem, nLine)` 列主序 (column-major) 内存。
C++ transpose helper 的源下标和目标下标**全部写反**。

### 2D 转换错误 (line 43)

**当前代码**:
```cpp
dst[i * nLine + j] = src[j * nElem + i];
```

**正确写法**:
```cpp
dst[j * nElem + i] = src[i * nLine + j];
```

**推演** (以 4x3 小矩阵为例):
- C-order `(4,3)` 内存: `[a00, a01, a02, a10, a11, a12, a20, a21, a22, a30, a31, a32]`
- 正确 Fortran 列主序: `[a00, a10, a20, a30, a01, a11, a21, a31, a02, a12, a22, a32]`
- 当前代码输出: `[a00, a01, a02, a10, a11, a12, a20, a21, a22, a30, a31, a32]` -- 就是原始 C-order，没有转置

对非方阵 `(2048, 2000)`，这不是简单未转置，而是带有周期性剪切的空间错位重排。

### 3D 转换错误 (line 66)

**当前代码**:
```cpp
const T* src_px = src + (j * nElem + i) * K;   // 错误
dst[(k * nLine + j) * nElem + i] = src_px[k];
```

**正确写法**:
```cpp
const T* src_px = src + (i * nLine + j) * K;   // 修正
dst[(k * nLine + j) * nElem + i] = src_px[k];
```

3D 目标下标 `(k * nLine + j) * nElem + i` 是正确的 Fortran 列主序，但源下标用反了。

### 输出 reshape (line 183) -- 正确

```cpp
return np.attr("reshape")(flat, py::make_tuple(nElem, nLine), py::arg("order") = "F");
```

输出 reshape `order='F'` 正确解读了 Fortran 列主序输出。但输入已经被错误搬运，输出的空间对应关系自然是错的。

### 影响链条

1. `ref_vis`, `tbb_ir`, `lat/lon`, `sza/vza`, `lsf`, `eco`, `snow_mask`, `btclr` 全部被同样错误地搬运
2. Fortran `process_swath_c` 按 `(ielem, iline)` 正常读取这些数组，但数据已经是错位的
3. 3x3 邻域 `indat_local` (line 418-436) 基于错位后的空间网格取邻域，"邻域不再是真实邻域"
4. 输出用 `order='F'` reshape 回 Python，但只能正确解释 Fortran 输出，不能修复输入已经错位
5. **这是椒盐噪声的结构性根因**: 云/晴判别基于错位的光谱+空间信息，输出自然表现为随机跳变

### 修复方案

```cpp
// 2D: Python C-order (nElem, nLine) -> Fortran column-major (nElem, nLine)
dst[j * nElem + i] = src[i * nLine + j];

// 3D: Python C-order (nElem, nLine, K) -> Fortran column-major (nElem, nLine, K)
const T* src_px = src + ((i * nLine + j) * K);
dst[(k * nLine + j) * nElem + i] = src_px[k];
```

### 验证测试

修复前必须先写最小测试，不跑云检测，只验证 C++ helper 语义：
- 2D: 构造 `arr[i, j] = 10000*i + j`，转换后 Fortran 侧 `A(i,j)` 应严格等于原始 `arr[i,j]`
- 3D: 构造 `arr[i,j,k] = 10000*i + 10*j + k`，同理
- 验证输出 reshape 回 Python 后坐标不变

---

## 根因 1：find_matched_files.py 配置错误导致 NWP 永远不加载 [CRITICAL] ✅ 已修复

**文件**: scripts/find_matched_files.py

两个配置错误：
```python
MERSI_ROOT = Path('/data/Data_yuq/mersi_test')   # 错误！应为 /data/Data_yuq/mersi
NWP_PATTERN = 'fnl_{date}_{hh}_00'               # 错误！应为 gfs0p25_41L_{date}_{hh}_00
```

**已修复**: v3.5.1

---

## 根因 2：btclr 全零导致 PFMFT Group 1 失效 [CRITICAL] ✅ 已修复

**文件**: scripts/run_fortran_only.py

**已修复**: v3.5.1 (sfctmp 估算), v3.5.4 (按表面类型区分: Ocean BTD_clr=0.5K, Land BTD_clr=1.9K)

---

## 根因 3：ocean 路径 PFMFT conf_test 被注释 [CRITICAL] ✅ 已修复

**已修复**: v3.5.1

---

## 根因 4：snow_mask 全零 [HIGH] 未修复

```python
# scripts/run_fortran_only.py:337
snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
```

雪/冰检测完全依赖 NDSI 可见光方法。夜间和高纬度冬季无法检测积雪。

Fortran 侧夜间路径 (`cloudmask_c_api.f90:658-660`):
```fortran
else
    ! Night: use ancillary map only
    ice  = map_ice
    snow = map_snow
end if
```

当 `snow_mask_val=0` 时 `map_snow=.false., map_ice=.false.`，积雪地表被当普通陆地处理。

**依赖**: NISE 雪冰掩码产品，或季节/纬度 fallback。

---

## 新发现的 Bug (v3.5.1~v3.5.4 诊断过程中发现)

### BUG-A：IR 波数硬编码错误 [CRITICAL] ✅ 已修复

**修复**: v3.5.2 -- 从 HDF5 `Effect_Center_WaveLength` 属性实时读取

### BUG-B：缺少 TBB 线性修正 [HIGH] ✅ 已修复

**修复**: v3.5.2

### BUG-C：Fortran planck_module.f90 死代码 [LOW] ✅ 已清理

**清理**: v3.5.3

### BUG-D：APOLLO tview 查找表 100% 回退 [HIGH] 未修复

tview 返回值恒 < 0.1，导致 APOLLO 全部回退到静态阈值 dfthrsh=3.0K。
输入参数（sec(VZA)、BT11）均在表范围内，怀疑 tview 表数据未正确加载。

**影响**: ocean_day/ocean_nite Group 2 的 11-12um APOLLO 测试失效

### BUG-E：PFMFT 阈值与 MERSI-II 实际分布不匹配 [HIGH] 未修复

**文件**: coeff/coeff/fylat_thresholds.mersi.ii3d.v8:265

```
pfmft_ocean : 1.9, 1.8, 1.7, 1
```

实际 PFMFT 信号: P50=-0.76K, P90=+1.31K -> 93% 海洋像素 PFMFT 触底 (c=0.1)

**需要**: 阈值重标定为 [0.2, 0.6]K 量级。**但必须先修复根因 0 数组布局问题，再基于正确数据重标定。**

### BUG-F (NEW)：正式 pipeline relaz 全零占位 [MEDIUM] 未修复

**文件**: src/fy3_cloudmask/pipeline.py:221

```python
relaz = np.zeros_like(lat, dtype=np.float32)  # Placeholder
```

验证脚本 `run_fortran_only.py:162` 已正确计算 `relaz = |saa - vaa|`，但正式 pipeline 没接上。
相对方位角影响几何相关测试、耀斑/可见光路径，白天海洋尤其敏感。

### BUG-G (NEW)：btclr 仍是经验估算，非 RTM [MEDIUM] 未修复

**文件**: scripts/run_fortran_only.py:223-234

```python
bt[ocean, 4] = sfctmp[ocean] - 5.6
bt[ocean, 5] = sfctmp[ocean] - 6.1
bt[land, 4] = sfctmp[land] - 9.9
bt[land, 5] = sfctmp[land] - 11.8
```

比全零好，但不是按廓线、视角、水汽、地表发射率计算的 clear-sky BT。PFMFT/APOLLO 依赖 btclr，估算误差直接进入判别函数。

pipeline.py 里更严重，直接 `np.full((*lat.shape, 7), 280.0, dtype=np.float32)`。

### BUG-H (NEW)：后处理不更新 bit/QA 数组 [MEDIUM] 未修复

**文件**: src/fortran/c_api/cloudmask_c_api.f90:477-482

后处理只改了 `out_cloud_mask` 和 `out_confidence`：
```fortran
call smooth_conf_reclassify(out_confidence, out_cloud_mask, nElem, nLine)
call apply_spatial_consistency(out_cloud_mask, nElem, nLine)
```

`out_cm_bitarray` / `out_qa_bitarray` 没有同步回写。HDF5 里的 `Cloud_Mask_Value` 和 `TestBits`/`QA_Bits` 不一致，影响验证和诊断。

### BUG-I (NEW)：原始 Fortran apply_spatial_filter() 在 native 后端不生效 [LOW] 已确认

**文件**: src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90:1325

原始模块的 `convert_cloud_mask()` 末尾调用了 `apply_spatial_filter()`，但 native 后端走的是 C API 的 `process_swath_c`，逐像元调用 `process_pixel_c`，后处理用的是 `smooth_conf_reclassify` + `apply_spatial_consistency`。原始模块的 `apply_spatial_filter()` 在 native 流程中从未被执行。

---

## 关于 APOLLO 的纠正说明

goals.md 旧版本认为 ocean APOLLO 测试会产生垃圾阈值，需要禁用（FIX-D）。**这是错误的。**

APOLLO 有内置回退机制：
```fortran
if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
  dfthrsh = do11_12hi(1)    ! 回退到静态阈值 3.0K
else
  dfthrsh = diftemp          ! 用 APOLLO 自适应阈值
end if
```

当 btclr=0 时，tview 返回 diftemp < 0.1，自动回退到静态阈值。但 v3.5.4 修复波数+TBB 后 btclr 正常，仍需排查 tview 表数据加载。

---

## 其他已知 Bug

### 已修复 (v3.4.3~v3.5.4)

| BUG | 描述 | 修复版本 |
|-----|------|----------|
| BUG-1 | IR_WAVENUMBERS 硬编码错误 (909.458->933.358) | v3.5.2 |
| BUG-2 | relaz 从 saa/vaa 计算，不再传零 | v3.4.5 |
| BUG-4 | set_confdnc 添加 set_bit(testbits,0) | v3.4.4 |
| BUG-5 | 沙漠判断改为完整原始逻辑 | v3.4.5 |
| BUG-6 | fill_bit_pixel 独立输出参数 | v3.4.4 |
| BUG-7 | smooth_conf_reclassify 排除未处理像素 | v3.4.4 |
| BUG-8 | PolarNite_snow/ocean_nite 统一用 pxldat(20) | v3.4.4 |
| BUG-9 | ocean_nite groups==0 分支修复 | v3.4.4 |
| BUG-12 | ocean_nite 取消 c4/c6 注释 | v3.4.4 |
| BUG-13 | LandNite 取消 c6 注释 | v3.4.4 |
| BUG-14 | LandNite i4=0 死代码清理 | v3.4.4 |
| BUG-15 | LandDay_desert cirrus_vis 覆盖修复 | v3.4.4 |
| -- | find_matched_files.py MERSI_ROOT/NWP_PATTERN 错误 | v3.5.1 |
| -- | btclr 全零 (sfctmp 估算 -> 表面类型区分) | v3.5.1, v3.5.4 |
| -- | ocean_day/ocean_nite PFMFT conf_test 注释 | v3.5.1 |
| -- | 缺少 TBB_Trans 线性修正 | v3.5.2 |
| -- | Fortran planck_module 死代码 | v3.5.3 |

### 未修复

| BUG | 描述 | 严重度 | 状态 |
|-----|------|--------|------|
| ~~根因0~~ | ~~pybind 数组布局转换下标写反~~ | ~~CRITICAL~~ | ✅ v3.6.0 |
| 根因4 | snow_mask 全零 | HIGH | 缺 NISE 数据 |
| BUG-D | APOLLO tview 表 100% 回退 | HIGH | 待查 tview 表加载 |
| BUG-E | PFMFT 阈值不匹配 (1.7~1.9K vs 实际 -0.8~1.3K) | HIGH | 需先修根因0再重标定 |
| BUG-F | pipeline relaz 全零占位 | MEDIUM | 需接入 GEO 计算逻辑 |
| BUG-G | btclr 经验估算（非 RTM） | MEDIUM | 需 RTM 或更精细估算 |
| BUG-H | 后处理不更新 bit/QA 数组 | MEDIUM | 需同步回写 |
| BUG-16 | nfmft 阈值范围过宽（-23~-22K，实际值 -1~1K） | MEDIUM | 需重新标定 |
| BUG-17 | ocean_nite 11-4um 阈值过窄 (2.25K) | MEDIUM | 需重新标定 |
| BUG-18 | GEMI 公式除零风险 (LandDay s1->1.0) | MEDIUM | 代码修改 |
| BUG-I | 原始 Fortran apply_spatial_filter() 不生效 | LOW | native 后端不走该路径 |

---

## 修复计划

### v3.6.0 -- 数组布局修复 + 验证 (当前最高优先级)

**目标**: 修复 pybind 数组布局转换错误，验证椒盐噪声是否消除

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-0a | 修复 2D transpose 下标: `dst[j*nElem+i] = src[i*nLine+j]` | ext/cloudmask_pybind.cpp:43 |
| FIX-0b | 修复 3D transpose 源下标: `src + ((i*nLine+j)*K)` | ext/cloudmask_pybind.cpp:66 |
| FIX-0c | 编写最小验证测试 (2D/3D 坐标一致性) | tests/test_array_layout.py |
| FIX-0d | 重建 native 后端，跑同一轨对比 | -- |

**验收标准**:
- 构造 `arr[i,j] = 10000*i + j` 经 C++ 转换后 Fortran 侧读取一致
- Cloud_Mask_Value 图像从随机颗粒变为有连续云团
- 相邻像元云/晴跳变率显著下降
- lsf/eco/sza/vza/bt11/bt12 输出抽样空间连续

### v3.6.1 -- 阈值重新标定

> 必须在 v3.6.0 修复数组布局之后执行。当前分布数据可能受错位影响。

**目标**: 将 MERSI-II 阈值从 MODIS 标定值调整到 MERSI-II 实际观测值

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-E1 | PFMFT ocean 阈值 [1.7,1.9] -> [0.2,0.6]K | pfmft_nfmft_thr.inc |
| FIX-E2 | PFMFT land 阈值重标定 | pfmft_nfmft_thr.inc |
| FIX-D | APOLLO tview 表加载修复 | tview.f |
| FIX-6 | nfmft 阈值重标定 | pfmft_nfmft_thr.inc |
| FIX-7 | ocean_nite 11-4um 阈值调宽 | ocean_nite_thr.inc (或 coeff/) |

**前置条件**: v3.6.0 完成后，基于正确空间数据重新统计 PFMFT/NFMFT/APOLLO 分布。

### v3.6.2 -- 输入完整性 + 后处理

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-4 | snow_mask 传入真实 NISE 数据（或季节/纬度 fallback） | scripts/run_fortran_only.py, pipeline.py |
| FIX-F | pipeline relaz 接入 GEO 计算逻辑 | src/fy3_cloudmask/pipeline.py |
| FIX-G | btclr 改用 RTM 或更精细估算 | scripts/run_fortran_only.py, pipeline.py |
| FIX-H | 后处理同步更新 cm_bitarray/qa_bitarray | src/fortran/c_api/cloudmask_c_api.f90 |
| FIX-5 | GEMI 除零保护 | -- |

### v3.6.3 -- 后处理增强（可选）

当前 `apply_spatial_consistency` 的 `max_count >= 7` 过于保守。
修复根因 0 后椒盐应大幅减少，再评估是否需要增强后处理：
- 连通域面积阈值（去除小斑块）
- 3x3/5x5 majority filter（保留云边界）
- 如修改 cloud_mask，同步更新 cm_bitarray 前 3 个关键 bit

---

## 验证方案

### 测试数据

**日期**: 2022-08-03
- MERSI: `/data/Data_yuq/mersi/20220803/`（3 轨道：0740, 0830, 0920）
- NWP: `/data/nwp/20220803/ORG/gfs0p25_41L_20220803_{00..24}_00`（9 时次，1.1GB/个）
- MYD35: `/data/Data_yuq/aqua_modis/MYD35_L2/20220803/`

### 快速验证

```bash
PYTHONPATH=src python scripts/run_fortran_only.py --date 20220803 --output /data/Data_yuq/fy3_cloud
```

### v3.6.0 验证清单 (2026-06-23)

- [x] 数组布局单元测试通过 (2D/3D 坐标一致性) -- `verify_transpose(10,8)` 全部 PASS
- [x] Cloud_Mask_Value 空间图像无椒盐 -- 孤立盐/胡椒像素仅 0.42%
- [x] 相邻像元跳变率 < 5%（当前预计 > 30%）-- 水平 8.17%, 垂直 8.29% (含真实云边界)
- [x] confidence 分布不再大量堆积在阈值边界 -- 分布平滑: >=0.66: 51.8%, >=0.95: 5.9%
- [x] lsf/eco/sza/vza 空间连续性目视检查通过 -- 转置前向验证通过
- [x] 两套定标系数 (onboard + recal) 均成功运行 3 轨

### 完整验证

- 7 个测试日期（20220803~20250302）全轨道验证
- 与 MYD35 对比：Agreement、HSS、混淆矩阵
- BTD(11-12) 逐通道统计
- 空间分布图：检查椒盐是否消除
