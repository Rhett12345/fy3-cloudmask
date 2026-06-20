# FY-3D 云检测代码 Bug 完整汇总

**当前版本**: v3.4.3
**GitHub 状态**: main 分支，已推送

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

---

## v3.4.3 已修复的问题

### ✅ FIX-1 | convert_cloud_mask_value 未处理像素返回值

**状态**: 已修复 (v3.4.3)
**原 BUG-4 的部分修复**: `btest(tb(1), 0)` 读取 bit0 判断是否已处理，但 `set_confdnc.f` 从未设置 bit0。
**当前修复**: `b0==0` 时返回 `cm_val=5` (fill value) 而非 `cm_val=0` (cloudy)。

**剩余问题**: `set_confdnc.f` 仍然不设置 bit0，依赖下游 `convert_cloud_mask_value` 的特判。如果未来有其他代码读取 bit0，仍会出错。

---

### ✅ FIX-2 | ocean_nite groups==0 分支

**状态**: 部分修复 (v3.4.3)
**当前代码**:
```fortran
if (groups .gt. 0) then
    fac = 1.0 / groups
else
    confdnc = 1.0
end if
confdnc = pre_confdnc**fac   ! ← 这行仍然会执行
confdnc = max(confdnc, 0.1)
```

**剩余问题**: 当 groups==0 时，`confdnc=1.0` 被设置后，下一行 `confdnc = pre_confdnc**fac` 又会覆盖它。如果 fac 未初始化或为 0.0，结果碰巧正确 (x^0=1.0)，但逻辑不严谨。应改为：
```fortran
if (groups .gt. 0) then
    fac = 1.0 / groups
    confdnc = pre_confdnc**fac
else
    confdnc = 1.0
end if
confdnc = max(confdnc, 0.1)
```

---

### ✅ FIX-3 | 空间一致性滤波

**状态**: 已修复 (v3.4.3)
- 从 0/8 (仅孤立像素) 改为 7/8 多数投票
- 添加了 3x3 中值滤波 (`smooth_conf_reclassify`)
- 添加了 cmin floor = 0.1 和 confdnc floor = 0.1

---

## 未修复的 Bug

### BUG-1 | IR 波数表错误 | 致命 | 椒盐直接根因

**状态**: ❌ 未修复
**文件**: `scripts/run_fortran_only.py`，第 39 行

**当前值**:
```python
IR_WAVENUMBERS = np.array([2643.4359, 2471.654, 1382.621, 1168.182, 933.364, 836.941])
#                                                                      ^^^^^^^ 10.7μm，错误
```

**正确值** (应从 HDF 文件头 CenterWavenum 属性读取):
```python
IR_WAVENUMBERS = np.array([2643.4359, 2471.654, 1382.621, 1168.182, 909.458, 836.941])
#                          3.8μm      4.05μm    7.3μm     8.5μm    11μm     12μm
```

**影响**: 11μm 通道 BT 计算偏低约 10-15K，所有依赖 BT11 的测试全部接收到错误输入。
**修复版本**: v3.4.4

---

### BUG-2 | relaz（相对方位角）全部传零 | 致命

**状态**: ❌ 未修复
**文件**: `scripts/run_fortran_only.py`，第 295 行

**当前值**:
```python
relaz=np.ascontiguousarray(np.zeros_like(geo['sza']).astype(np.float32)),
```

**修复**:
```python
# 在 read_geo_data 里计算 relaz
relaz = np.abs(geo['saa'] - geo['vaa'])  # 或 (saa - vaa) mod 360
```

**影响**: 耀斑检测结果不可靠。
**修复版本**: v3.4.4

---

### BUG-3 | snow_mask 全部传零 | 严重

**状态**: ❌ 未修复
**文件**: `scripts/run_fortran_only.py`，第 301 行

**当前值**:
```python
snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
```

**修复**: 传入真实 NISE 雪冰掩码数据。
**影响**: 极地和高纬度夜间像素积雪/海冰判断失效。
**修复版本**: v3.4.5 (需要额外数据源)

---

### BUG-4 | set_bit 编号与 btest 读取不一致 | 严重

**状态**: ⚠️ 部分修复 (v3.4.3 通过特判绕过)
**文件**: `src/fortran/cloudmask/set_confdnc.f`

**当前状态**: `set_confdnc.f` 仍未设置 bit0，但 `convert_cloud_mask_value` 已通过 `b0==0 → cm_val=5` 特判绕过。

**建议**: 在 `set_confdnc.f` 中添加 `set_bit(testbits, 0)` 使逻辑自洽：
```fortran
if(confdnc .gt. 0.99) then
    call set_bit(testbits, 0)   ! 已处理标志
    call set_bit(testbits, 1)
    call set_bit(testbits, 2)
else if(confdnc .gt. 0.95) then
    call set_bit(testbits, 0)
    call set_bit(testbits, 2)
else if(confdnc .gt. 0.66) then
    call set_bit(testbits, 0)
    call set_bit(testbits, 1)
else
    call set_bit(testbits, 0)   ! 有云也是已处理
end if
```

**修复版本**: v3.4.4

---

### BUG-5 | compute_pixel_flags 沙漠判断逻辑错误 | 严重

**状态**: ❌ 未修复
**文件**: `src/fortran/c_api/cloudmask_c_api.f90`

**当前值**:
```fortran
desert = .false.
if (eco_int >= 7 .and. eco_int <= 10) desert = .true.
if (eco_int == 16) desert = .true.
```

**修复**: 迁移原始 `get_pxldat` 的完整沙漠判断逻辑。
**影响**: 大量像素走错处理路径。
**修复版本**: v3.4.5

---

### BUG-6 | fill_bit_pixel 输入输出传了同一数组 | 严重

**状态**: ❌ 未修复
**文件**: `src/fortran/c_api/cloudmask_c_api.f90`

**当前值**:
```fortran
call fill_bit_pixel(..., testbits, qa_bits, testbits, qa_bits)
```

**修复**: 使用独立的输出变量。
**影响**: 位操作结果随机。
**修复版本**: v3.4.4

---

### BUG-7 | smooth_conf_reclassify 未排除未处理像素 | 中等

**状态**: ❌ 未修复
**文件**: `src/fortran/c_api/cloudmask_c_api.f90`

**当前值**:
```fortran
if (conf_val >= 0.0) then   ! 未处理像素 confidence=0.0 也被纳入滤波
    n = n + 1
    window(n) = conf_val
end if
```

**修复**:
```fortran
if (conf_val >= 0.0 .and. cloud_mask(i+di, j+dj) /= 5) then
    n = n + 1
    window(n) = conf_val
end if
```

**影响**: 边界像素的 confidence 被向下拖拽。
**修复版本**: v3.4.4

---

### BUG-8 | PolarNite_snow 和 ocean_nite 用错了 3.8μm 通道 | 严重

**状态**: ❌ 未修复
**文件**: `src/fortran/cloudmask/PolarNite_snow.f90` 和 `src/fortran/cloudmask/ocean_nite.f90`

**当前值**:
```fortran
masir4 = pxldat(21)   ! 4.05μm (错误)
```

**修复**: 改为 `masir4 = pxldat(20)` (3.8μm)
**影响**: 11-4μm BTD 和 4-12μm 测试使用错误通道。
**修复版本**: v3.4.4

---

### BUG-9 | ocean_nite groups==0 分支逻辑不严谨 | 中等

**状态**: ⚠️ 部分修复 (v3.4.3)
**文件**: `src/fortran/cloudmask/ocean_nite.f90`

**当前代码**:
```fortran
if (groups .gt. 0) then
    fac = 1.0 / groups
else
    confdnc = 1.0
end if
confdnc = pre_confdnc**fac   ! ← 这行仍然执行
```

**修复**: 将 `confdnc = pre_confdnc**fac` 移入 if 块内。
**修复版本**: v3.4.4

---

### BUG-10 | APOLLO 11-12μm 测试被 .false. 硬禁用 | 待定

**状态**: ⚠️ 故意禁用 (v3.4.2)
**原因**: APOLLO 查找表为 MODIS 设计 (2-10K BTD 范围)，MERSI-II 的 11-12μm BTD 范围更小 (0-4K)，导致 98.8% 像素触发测试。

**后续**: 需要为 MERSI-II 重新标定 APOLLO 查找表，或改用 MERSI-II 专用阈值。
**修复版本**: v3.5.0 (需要重新标定)

---

### BUG-11 | pfmft 和 nfmft 全部注释 | 待定

**状态**: ⚠️ 故意禁用 (v3.4.1)
**原因**: btclr (晴空亮温) 缺失，传入全零。

**后续**:
1. 提供 btclr 数据后解注释
2. 重新标定 `nfmft_land` 阈值 (-23~-22K 远离典型值)

**修复版本**: v3.5.0 (需要 NWP RTM 数据)

---

### BUG-12 | ocean_nite 三光谱和 11-4μm 置信度注释掉了 | 中等

**状态**: ❌ 未修复
**文件**: `src/fortran/cloudmask/ocean_nite.f90`

**当前值**:
```fortran
! cmin2 = min(cmin2, c4)       ! 三光谱，被注释 by wuxiao
! cmin2 = min(cmin2, c6)       ! 11-4μm，被注释 by wuxiao
```

**修复**: 取消注释。
**修复版本**: v3.4.4

---

### BUG-13 | LandNite 7.3-11μm 置信度注释掉了 | 中等

**状态**: ❌ 未修复
**文件**: `src/fortran/cloudmask/LandNite.f90`

**当前值**:
```fortran
! cmin2 = min(cmin2, c6)   ! 被注释
```

**修复**: 取消注释。
**修复版本**: v3.4.4

---

### BUG-14 | LandNite 12-3.7μm 测试是死代码 | 低

**状态**: ❌ 未修复
**文件**: `src/fortran/cloudmask/LandNite.f90`

**当前值**:
```fortran
i4 = 0
if (i4 == 1) then   ! 永远为假
```

**修复**: 决定是否启用，否则删除。
**修复版本**: v3.4.4

---

### BUG-15 | LandDay_desert 薄卷云标志被立即覆盖 | 低

**状态**: ❌ 未修复
**文件**: `src/fortran/cloudmask/LandDay_desert.f90`

**当前值**:
```fortran
cirrus_vis = .true.
cirrus_vis = .false.   ! 立即覆盖
```

**修复**: 删除第二行。
**修复版本**: v3.4.4

---

## 修复计划

### v3.4.4 (修复直接导致椒盐的 Bug) ✅ 已完成

| BUG | 描述 | 文件 | 状态 |
|-----|------|------|------|
| BUG-1 | IR_WAVENUMBERS 波数值 | scripts/run_fortran_only.py | ✅ |
| BUG-4 | set_confdnc 加 set_bit(testbits,0) | src/fortran/cloudmask/set_confdnc.f | ✅ |
| BUG-6 | fill_bit_pixel 独立输出参数 | src/fortran/c_api/cloudmask_c_api.f90 | ✅ |
| BUG-7 | smooth_conf_reclassify 排除未处理像素 | src/fortran/c_api/cloudmask_c_api.f90 | ✅ |
| BUG-8 | PolarNite_snow/ocean_nite 统一用 pxldat(20) | PolarNite_snow.f90, ocean_nite.f90 | ✅ |
| BUG-9 | ocean_nite groups==0 分支 | ocean_nite.f90 | ✅ |
| BUG-12 | ocean_nite 取消 c4/c6 注释 | ocean_nite.f90 | ✅ |
| BUG-13 | LandNite 取消 c6 注释 | LandNite.f90 | ✅ |
| BUG-14 | LandNite i4=0 死代码清理 | LandNite.f90 | ✅ |
| BUG-15 | LandDay_desert cirrus_vis 覆盖 | LandDay_desert.f90 | ✅ |

**v3.4.4 MYD35 验证结果 (2020-03-08)**:

| Orbit | 精度 | FY3D云量 | MYD35云量 | POD_cld | POD_clr | 区域 |
|-------|------|----------|-----------|---------|---------|------|
| 1345 | 27.6% | 97.9% | 26.6% | 97.7% | 2.1% | 南极 (-79~-54) |
| 1435 | 44.9% | 90.6% | 43.6% | 90.6% | 9.5% | 北极 (56~82) |
| 1525 | 70.2% | 84.4% | 79.1% | 84.5% | 16.0% | 南极 (-85~-58) |

**v3.4.4 结论**:
- BUG-1 修复 (11μm波数纠正) 后 BT11 变暖约10-15K，但极地云量仍然过高
- Orbit 1525 与 MYD35 一致性尚可 (70%)，但 1345/1435 极地严重过判云
- 主要残留问题：pfmft/nfmft 禁用、APOLLO 禁用、snow_mask=0、沙漠判断不完整
- 极地对晴空判识能力极弱 (POD_clr < 16%)，需要后续版本解决

### v3.4.5 (改善精度，需要额外数据)

| BUG | 描述 | 文件 |
|-----|------|------|
| BUG-2 | relaz 传真实值 | scripts/run_fortran_only.py |
| BUG-3 | snow_mask 传真实值 | scripts/run_fortran_only.py |
| BUG-5 | 沙漠判断用离散列表 | src/fortran/c_api/cloudmask_c_api.f90 |

### v3.5.0 (改善检测灵敏度，需重新标定)

| BUG | 描述 | 文件 |
|-----|------|------|
| BUG-10 | APOLLO 测试重新标定 | 所有 .f90 测试路径 |
| BUG-11 | pfmft/nfmft 解注释 + 重新标定阈值 | 所有 .f90 测试路径 |

---

## 当前工作目录状态

```
On branch main
Your branch is up to date with 'origin/main'.

Changes not staged for commit:
  modified:   .claude/settings.local.json
  modified:   ext/build/fortran_modules/cloudmask_c_api_mod.mod
  modified:   plans/goals.md
  deleted:    plans/original_vs_new_diff.md

Untracked files:
  scripts/diag_pfmft_nfmft.py
  scripts/diag_run_fortran.py
  scripts/process_20200308_btclr.py
  scripts/process_20200308_f90.py
  scripts/test_btclr_fix.py
```
