# FY3 Cloud Mask 迁移代码审查报告

**审查日期**: 2026-05-28
**审查范围**: `fy3_cloudmask/src/fy3_cloudmask/` (Python 迁移版本)
**对照基准**: 原始 Fortran 版本 (retrieval_system_V3.1_cldmask)
**审查目标**: 识别迁移过程中的 bug、逻辑偏差和性能/阈值设计问题

---

## 一、严重 Bug（会导致结果错误）

### BUG-001: `fill_bit_pixel` 缺少 `desert` 参数

**文件**: `src/fy3_cloudmask/algorithm/bitops.py:94`

**问题描述**:
Fortran 版本签名:
```fortran
subroutine fill_bit_pixel(nmtests, nbands, bad_value, bad_geo,
                          snglnt, desert, testbits, qa_bits, cm, qa)
```
Python 版本签名:
```python
def fill_bit_pixel(nmtests, nbands, bad_geo, snglnt, testbits, qa_bits):
```
`desert` 参数被完全省略。Fortran 中 `desert` 用于影响沙漠路径下的质量标记逻辑。

**影响**: 沙漠像素的质量标记可能不正确。

**修复建议**: 添加 `desert: bool` 参数，并在函数内部使用。

---

### BUG-002: `cloud_mask.py` 中 surface type bits 重复设置

**文件**: `src/fy3_cloudmask/algorithm/cloud_mask.py:252-258`

**问题描述**:
`proc_path()` 已经根据 `land/coast/desert` 设置了 bit 6 和 bit 7。但在 `run_cloud_mask_pixel()` 的 Step 5 中又重复设置了这些 bit:

```python
# proc_path 已经设置了这些 bit
proc_path(result.testbits, ...)

# 这里又重复设置
if flags.land:
    set_bit(result.testbits, BIT_COAST)   # bit 6
    set_bit(result.testbits, BIT_DESERT)  # bit 7
elif flags.coast:
    set_bit(result.testbits, BIT_COAST)
elif flags.desert:
    set_bit(result.bits, BIT_DESERT)
```

虽然 OR 操作使结果相同，但说明对 Fortran 流程理解有误。Fortran 只在 `proc_path` 中设置一次。

**影响**: 代码冗余，可能掩盖更深层的逻辑错误。

**修复建议**: 删除 `cloud_mask.py:252-258` 中重复的 surface type bit 设置代码。

---

### BUG-003: 缺少 `set_quality_A` 调用

**文件**: `src/fy3_cloudmask/algorithm/cloud_mask.py`

**问题描述**:
Fortran 在 `fill_bit_pixel` 之前调用了 `set_quality_A(nmtests, nbands, lsf, qa_bits)`，它设置以下 QA bit:
- QA bits 48-51: 测试数量编码（>4 → 50+51, >2 → 51, >0 → 50）
- QA bit 64: 生态系统文件标记（始终设置）
- QA bits 70-71: 陆地/海洋掩码标记（lsf==-1 时设置）

Python 版本完全没有实现这个函数。

**影响**: QA 字节 6-9 中的元数据信息丢失，影响产品质量评估和后续处理。

**修复建议**: 在 `bitops.py` 中实现 `set_quality_A` 函数，并在 `cloud_mask.py` 中调用。

---

### BUG-004: `land_day_coast` 错误地委托给 `land_day_standard`

**文件**: `src/fy3_cloudmask/algorithm/tests/land_day.py:274`

**问题描述**:
函数注释说"使用沿海特定阈值"，但实际直接调用 `land_day_standard`：
```python
def land_day_coast(...):
    # Same logic as standard land but with coastal thresholds
    return land_day_standard(pxldat, bt_clr, vza, is_cold_sfc, hi_elev,
                            thresholds, testbits, qa_bits)
```

配置文件中 `land_day_coast` 有独立的阈值：
```yaml
land_day_coast:
  ref064: [0.22, 0.18, 0.14, 1.0]  # 与 land_day 不同
  btd_11_12: 3.0
```

但 `land_day_standard` 使用的是 `land_day` 的阈值 `[0.24, 0.20, 0.16]`。

**影响**: 沿海地区使用了错误的 0.64um 反射率阈值，可能导致云检测偏差。

**修复建议**: 实现独立的 `land_day_coast` 函数，使用 `thresholds.get('land_day_coast', {})` 获取阈值。

---

### BUG-005: `land_day_desert_coast` 错误地委托给 `land_day_desert`

**文件**: `src/fy3_cloudmask/algorithm/tests/land_day.py:408`

**问题描述**:
与 BUG-004 类似，`land_day_desert_coast` 直接调用 `land_day_desert`，使用了 `land_day_desert` 的阈值而非 `land_day_desert_coast` 的。

关键阈值差异：
| 阈值 | desert | desert_coast |
|------|--------|--------------|
| btd_11_4_hi | [-3.0, -5.0, -7.0] | [2.0, 0.0, -2.0] |
| btd_11_4_lo | [-25.0, -23.0, -21.0] | [-20.0, -18.0, -16.0] |
| ref086 | [0.42, 0.39, 0.36] | [0.34, 0.30, 0.26] |

**影响**: 沿海沙漠区域的云检测阈值完全错误。

**修复建议**: 实现独立的 `land_day_desert_coast` 函数，使用 `thresholds.get('land_day_desert_coast', {})`。

---

### BUG-006: 处理顺序与 Fortran 不一致

**文件**: `src/fy3_cloudmask/algorithm/cloud_mask.py:163-231`

**问题描述**:
Fortran 的处理顺序：
```
dispatch test → shadows → noncld_obs_chk → thin_ci_chk_ir → proc_path
→ set_unused_bits → set_confdnc → set_quality_A → fill_bit_pixel
```

Python 的处理顺序：
```
dispatch test → chk_land_restoral → chk_coast_restoral → chk_sunglint_restoral
→ chk_shallow_water → chk_spatial_var → chk_thin_cirrus_ir → chk_shadow
→ chk_cloud_adj → encode_confidence → fill_bit_pixel → set_unused_bits
→ convert_cloud_mask
```

关键差异：
1. Fortran 的 `shadows` 只在 `(.not.water .and. .not.coast .and. day .and. .not.polar .and. confdnc>=0.66)` 时执行
2. Python 的 `chk_shadow` 只要有 3x3 数据就执行，缺少前置条件
3. Fortran 的 `noncld_obs_chk` 只在 `(land .and. day .and. .not. snow)` 时执行
4. Python 的 `chk_cloud_adj` 是空操作

**影响**: 在不应该执行某些测试的条件下错误执行。

**修复建议**: 严格按照 Fortran 的调用顺序和前置条件重构。

---

## 二、中等 Bug（功能不完整）

### BUG-007: `chk_sunglint_restoral` 是空操作

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:224-231`

**问题描述**:
```python
def chk_sunglint_restoral(confdnc, pxldat, refang, snglnt, thresholds,
                          testbits, qa_bits):
    if not snglnt:
        return confdnc
    return confdnc  # 直接返回，什么都不做
```

注释说"Fortran 版本不 cap confidence"，但实际上 Fortran 的 `chk_sunglint.f90` 会执行清除天空恢复测试，可以将 confidence 提升到 0.96。

**影响**: 太阳耀斑区域的清除天空像素无法被恢复为晴天。

**修复建议**: 实现 Fortran `chk_sunglint.f90` 中的清除天空恢复逻辑。

---

### BUG-008: `chk_cloud_adj` 是空操作

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:352-357`

**问题描述**:
```python
def chk_cloud_adj(confdnc, cm_array, row, col, n_rows, n_cols,
                  thresholds, testbits, qa_bits):
    return confdnc  # 直接返回
```

注释说"Fortran noncld_obs_chk.f90 不修改 confidence"，但 Fortran 实际上会：
1. 检查 11-12um BTD 来检测沙尘/烟雾
2. 清除 test bit 28（suspended dust）如果检测到非云遮挡物
3. 设置 smoke 标志

**影响**: 非云遮挡物（沙尘/烟雾）检测功能完全缺失。

**修复建议**: 实现 Fortran `noncld_obs_chk.f90` 中的 IR BTD 测试逻辑。

---

### BUG-009: `chk_spatial_var` 逻辑与 Fortran 不同

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:271-319`

**问题描述**:
Python 版本的逻辑是"如果空间均匀（std <= threshold），提升 confidence"。但 Fortran 的 `chk_spatial_var` 是一个独立的空间变异性测试，使用不同的阈值和逻辑。

Python 使用 bit 25，而 Fortran 使用的 bit 编号可能不同。

**影响**: 空间变异性检查的行为与原始算法不一致。

**修复建议**: 对照 Fortran `chk_spatial_var.f` 重新实现。

---

### BUG-010: Shadow 检测缺少前置条件

**文件**: `src/fy3_cloudmask/algorithm/cloud_mask.py:220`

**问题描述**:
Fortran 只在以下条件满足时检测阴影：
```fortran
if(.not.water .and. .not.coast .and. day .and.
   .not.polar .and. confdnc.ge.0.66) then
    call shadows(pxldat, shadow, visusd, qa_bits)
end if
```

Python 版本：
```python
if indat_3x3_vis is not None and indat_3x3_11um is not None:
    confdnc = chk_shadow(confdnc, pxldat, indat_3x3_vis, indat_3x3_11um, ...)
```

缺少的条件：非水域、非沿海、白天、非极区、confdnc >= 0.66。

**影响**: 在水域、沿海、极区、夜间等不应该检测阴影的地方错误执行，可能导致误判。

**修复建议**: 添加完整的前置条件检查。

---

### BUG-011: `land_nite` 缺少 `ice`/`snow` 参数和门控

**文件**: `src/fy3_cloudmask/algorithm/tests/land_nite.py`

**问题描述**:
Fortran 的 `LandNite` 接收 `ice` 和 `snow` 参数，并在某些测试中检查这些条件。Python 版本完全没有这些参数。

Fortran 签名：
```fortran
subroutine LandNite(pxldat, plat, vza, ice, snow, coast, tbadj,
                    desert, hi_elev, sh_lake, sfctmp, eco_type,
                    nmtests, testbits, qa_bits, confdnc, ptwp,
                    btclr, is_cold_sfc)
```

**影响**: 冰雪覆盖的夜间陆地像素可能使用错误的测试路径。

**修复建议**: 添加 `ice` 和 `snow` 参数，并在相应测试中使用。

---

### BUG-012: `chk_thin_cirrus_ir` 前置条件不完整

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:360`

**问题描述**:
Fortran 调用 `thin_ci_chk_ir` 的条件是：
```fortran
if ( (.not. snow) .and. (.not. ice) ) then
    call thin_ci_chk_ir(pxldat, vza, cirrus_ir, qa_bits, testbits)
end if
```

Python 版本没有检查 snow/ice 条件。

**影响**: 在冰雪区域可能错误地执行薄卷云检测。

**修复建议**: 添加 snow/ice 前置条件检查。

---

## 三、性能/阈值设计问题

### PERF-001: `ocean_nite` 的 `btd_11_4` 阈值方向需确认

**文件**: 配置 `config/thresholds/mersi_ii3d_v8.yaml` → `ocean_nite.btd_11_4`

**问题描述**:
阈值配置为 `[1.25, 1.0, -1.0, 1.0]`，Python 代码中 `if mas11_4 <= no11_4lo[1]`（即 <= 1.0）时设置 testbit。这表示 11-4um BTD 越小（更负）越可能是云。

需要对照 Fortran 原始阈值确认方向是否正确。

---

### PERF-002: `chk_shallow_water` 硬编码 max_confidence

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:266`

**问题描述**:
```python
max_conf = thr.get('max_confidence', 0.95)
confdnc = min(confdnc, max_conf)
```

浅水区域永远无法达到"confident clear"（需要 > 0.99）。这是设计意图还是过度保守？

---

### PERF-003: `conf_test` 忽略 `nmval` 参数

**文件**: `src/fy3_cloudmask/algorithm/confidence.py:15`

**问题描述**:
Fortran 的 `conf_test` 有 `nmval` 参数（1 或 2 个阈值），用于区分单阈值和双阈值测试。Python 版本完全忽略此参数，只处理单阈值情况。

如果存在双阈值测试（如 `conf_test_2val`），行为会不同。

**影响**: 潜在的双阈值测试行为偏差。

---

### PERF-004: `land_day.py` 中 `visusd` 门控不完整

**文件**: `src/fy3_cloudmask/algorithm/tests/land_day.py:161-232`

**问题描述**:
Fortran 的 Group 2（11-4um BTD）和 Group 3（可见光）测试都在 `if (visusd)` 块内。Python 版本中 11-4um BTD 测试没有 `visusd` 门控，而 Fortran 中它是受 `visusd` 保护的。

---

## 四、代码质量问题

### CODE-001: `compute_group_confidence` 中的无效过滤

**文件**: `src/fy3_cloudmask/algorithm/confidence.py:146`

```python
active = [c for c in group_confidences if c < 1.0 or True]  # all groups count
```

`or True` 使过滤条件永远为真，这是死代码。

---

### CODE-002: `@njit` 装饰器与 dict 操作不兼容

**文件**: 多个文件

多个使用 `@njit(cache=True)` 装饰的函数内部调用了 Python dict 的 `.get()` 方法。Numba 的 nopython 模式不支持 Python dict 操作，这些函数实际上不会被 JIT 编译，会回退到 object mode 或直接报错。

---

### CODE-003: `chk_sunglint_restoral` 中 `refang` 参数未使用

**文件**: `src/fy3_cloudmask/algorithm/tests/restoral.py:193`

函数接收 `refang` 参数但从未使用。

---

## 五、对照验证建议

修复上述 bug 后，建议进行以下验证：

1. **Bit-by-bit 对比**: 使用相同的输入数据，对比 Fortran 和 Python 输出的 `testbits` 和 `qa_bits` 数组
2. **阈值一致性**: 逐个检查 YAML 配置文件中的阈值是否与 Fortran 的 `.inc` 文件一致
3. **处理路径覆盖**: 确保所有 surface type × lighting 组合都有测试覆盖
4. **边界条件**: 测试极区、沿海、沙漠、冰雪等特殊场景

---

## 六、优先级排序

| 优先级 | Bug ID | 描述 | 影响范围 |
|--------|--------|------|----------|
| P0 | BUG-004 | land_day_coast 使用错误阈值 | 沿海白天 |
| P0 | BUG-005 | land_day_desert_coast 使用错误阈值 | 沿海沙漠白天 |
| P0 | BUG-006 | 处理顺序不一致 | 全局 |
| P0 | BUG-010 | Shadow 缺少前置条件 | 水域/极区/夜间 |
| P1 | BUG-001 | fill_bit_pixel 缺少 desert 参数 | 沙漠 |
| P1 | BUG-003 | 缺少 set_quality_A | QA 元数据 |
| P1 | BUG-007 | sunglint restoral 空操作 | 太阳耀斑区 |
| P1 | BUG-008 | cloud adj 空操作 | 沙尘/烟雾检测 |
| P2 | BUG-002 | surface type bits 重复设置 | 代码质量 |
| P2 | BUG-009 | spatial_var 逻辑不同 | 空间变异性 |
| P2 | BUG-011 | land_nite 缺少 ice/snow | 夜间冰雪 |
| P2 | BUG-012 | thin_cirrus_ir 缺少条件 | 冰雪区卷云 |
