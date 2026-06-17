# 原始代码 vs 当前 Fortran 后端 — 计算差异分析

**对照基准**: `retrieval_system_V3.1_cldmask/src/cloudmask/`
**审查对象**: `fy3_cloudmask/src/fortran/cloudmask/`
**审查范围**: 纯计算逻辑差异，不含注释、格式、架构差异

---

## 前提

阈值文件（`.inc`）**完全一致**，无差异。差异全部集中在算法逻辑层面。

---

## 差异 1：置信度下限（新增逻辑）

**涉及文件**: 全部 19 个测试路径文件

原始代码直接相乘：
```fortran
pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4
```

新版加了 0.1 下限：
```fortran
pre_confdnc = max(cmin1, 0.1) * max(cmin2, 0.1) * max(cmin3, 0.1) * max(cmin4, 0.1)
```

**效果**: 单组置信度为 0 时不再拖垮整像素。原来几何平均可以算出 0，现在最低 `0.1^(1/n)`。这是修复 48% 孤立云像素的核心改动。

---

## 差异 2：groups==0 的防御性分支（新增逻辑）

**涉及文件**: 全部测试路径文件

原始代码：
```fortran
if (groups .gt. 0) fac = 1.0 / groups
confdnc = pre_confdnc**fac
```

新版：
```fortran
if (groups .gt. 0) then
    fac = 1.0 / groups
    confdnc = pre_confdnc**fac
else
    confdnc = 1.0    ! 或 fac = 0.0，视文件而定
end if
```

**效果**: 无测试组时置信度默认 1.0（晴空），而非对 pre_confdnc 做未定义的幂运算。

---

## 差异 3：LandDay_desert / LandDay_desert_c 测试组重构（结构性变化）

**涉及文件**: `LandDay_desert.f90`, `LandDay_desert_c.f90`

原始代码有 4 组：
- Group 1: pfmft + nfmft (11-12um BTD)
- Group 2: 11-12um thin cirrus + 11-4um fog/low cloud
- Group 3: 0.86um 反射率
- Group 4: 1.38um 卷云

新版只有 3 组：
- Group 1: 11-12um + 11-4um（无 pfmft/nfmft）
- Group 2: 0.86um + 1.38um（合并为一组）
- Group 3: （无）

同时删除了 `include 'pfmft_nfmft_thr.inc'`，移除了 pfmft/nfmft 相关变量（`tv11_12`, `cosvza`, `schi` 等）。

**效果**: 沙漠路径丢失 pfmft/nfmft 测试，1.38um 从独立组变为子测试，分组权重改变。

---

## 差异 4：nmtests/ngtests 计数位置修正（逻辑修正）

**涉及文件**: `LandDay.f90`, `ocean_day.f90`, `LandDay_coast.f90`

原始代码在条件判断**前**无条件递增：
```fortran
nmtests = nmtests + 1
call set_qa_bit(qa_bits,20)
if (masv66.le.dlref1(2)) then
    call set_bit(testbits,20)
    ngtests(3) = ngtests(3) + 1
end if
```

新版在条件**通过后**才递增：
```fortran
call set_qa_bit(qa_bits,20)
if (masv66.le.dlref1(2)) then
    nmtests = nmtests + 1
    ngtests(3) = ngtests(3) + 1
end if
```

**效果**: `nmtests` 更准确，影响 `fill_bit_pixel` 的质量等级判定（`nmtests < 3` → 质量 4，`< 7` → 质量 6，`>= 7` → 质量 7）。

---

## 差异 5：chk_land / chk_coast restoral 置信度上限降低

**涉及文件**: `chk_land.f90`, `chk_coast.f90`

```fortran
! 原始: confdnc = 1.0
! 新版: confdnc = 0.97
```

**效果**: 恢复后的像素从"确定晴空"(>0.99) 降为"可能晴空"(0.95-0.99)，在 `set_confdnc` 中编码结果不同。

---

## 差异 6：chk_land restoral 空间均匀性门控（新增逻辑）

**涉及文件**: `land_module.f90`（白天和夜晚路径都加了）

新增 BT11 3x3 标准差检查，std < 1.5K 才允许 restoral：

```fortran
bt11_sum = 0.0
bt11_n = 0
do dj = 1, 3
  do di = 1, 3
    bt11_val = indat_11um(di, dj)
    if (abs(bt11_val - bad_data) .gt. 0.1) then
      bt11_sum = bt11_sum + bt11_val
      bt11_n = bt11_n + 1
    end if
  end do
end do
if (bt11_n .ge. 5) then
  bt11_mean = bt11_sum / real(bt11_n)
  bt11_std = sqrt(bt11_std / real(bt11_n))
  if (bt11_std .lt. 1.5) then
    call chk_land(...)
  end if
end if
```

**效果**: 云边界区域（BT11 变化大）的晴空恢复被跳过。

---

## 差异 7：check_reg_uniformity 容差放宽

**涉及文件**: `fylat_fy3mersi_cloud_mask.f90`

原始代码：任一不一致就 `uniform = .false.`

新版：允许最多 2 个不一致邻域（>= 7/9 匹配）：
```fortran
if (nmismatch .gt. 2) uniform = .false.
```

**效果**: 更多像素被标记为 uniform，使 `chk_spatial_var` 和 restoral 在更多场景下生效。

---

## 差异 8：空间一致性后处理过滤器（新增子程序）

**涉及文件**: `fylat_fy3mersi_cloud_mask.f90`

原始代码：无后处理。

新版在主循环后添加了 `apply_spatial_filter()`，对孤立像素（周围 8 邻域无同类）做多数投票修正：

```fortran
subroutine apply_spatial_filter()
    ! 遍历所有像素（跳过边缘）
    ! 对每个像素，统计 3x3 邻域 cloud mask 值
    ! 如果 same_count == 0（周围无同类），修正为多数类别
end subroutine
```

**效果**: 消除与 8 邻域完全不一致的孤立像素。

---

## 差异 9：PolarDay nfmft btclr 零值保护（新增条件）

**涉及文件**: `PolarDay_land.f90`, `PolarDay_coast.f90`

新增条件：
```fortran
(btclr(5) .ne. 0.0 .or. btclr(6) .ne. 0.0)
```

**效果**: 无晴空亮温参考数据时跳过 nfmft 测试，避免错误判定。

---

## 影响程度汇总

| # | 差异 | 影响范围 | 影响程度 |
|---|------|----------|----------|
| 1 | 置信度下限 0.1 | 全部路径 | **高** — 防止零置信度 |
| 2 | groups==0 防御 | 全部路径 | 低 — 边界情况 |
| 3 | 沙漠路径重构 | LandDay_desert, desert_c | **高** — 丢失 pfmft/nfmft |
| 4 | nmtests 计数修正 | LandDay, ocean_day, coast | 中 — 影响质量标记 |
| 5 | restoral 上限 0.97 | chk_land, chk_coast | 中 — 影响编码等级 |
| 6 | chk_land 空间门控 | land_module | 中 — 减少边界 restoral |
| 7 | uniformity 容差 | check_reg_uniformity | 低 — 扩大 uniform 范围 |
| 8 | 空间一致性滤波器 | 主循环后处理 | **高** — 新增去噪 |
| 9 | btclr 零值保护 | PolarDay_land, coast | 低 — 防御性检查 |
