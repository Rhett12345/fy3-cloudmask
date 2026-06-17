# Fortran 云掩膜椒盐噪声根因分析

## 背景

Python 版本产出的云掩膜没有椒盐噪声，但 Fortran 版本存在明显的椒盐噪声。本文档基于源码验证，逐项分析 6 个根因。

---

## 原因 1：逐像素独立处理，无空间上下文

**文件**: `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90:124-126`

```fortran
line_loop_1: do iline= 1, sat%nLine
element_loop_1: do ielem= 1, sat%nElem
```

每个像素的云判识完全独立。虽然提取了 3×3 窗口 `indat(1:3,1:3,1:25)`，但仅用于：
- `check_reg_uniformity`：判断邻域地表类型一致性（不参与云判识）
- `spatial_var`：仅在海洋白天 + uniform=true + 非冰面时触发

**`spatial_var` 的局限** (`src/fortran/cloudmask/spatial_var.f`):
- 仅对海洋像素生效
- 仅检查 11μm 单一波段
- 仅在 `uniform=true` 时才触发（uniform 条件极其严格）
- 只能提升置信度（从 0.67 提升到 0.96），不能降低

**结论**: 结构性根因，但不单独导致噪声。需要与原因 2、3、4 结合才会产生噪声。

---

## 原因 2：S 曲线阈值边界效应

**文件**: `src/fortran/cloudmask/conf_test.f`, `src/fortran/cloudmask/set_confdnc.f`

S 曲线将光谱测试值映射到 [0,1] 置信度。在阈值中点（50% confidence）附近，光谱值的微小波动导致置信度大幅跳变。

```fortran
! set_confdnc.f:47-54
if(confdnc .gt. 0.99) then       ! 晴空 + 高质量
  call set_bit(testbits,1)      
  call set_bit(testbits,2)
else if(confdnc .gt. 0.95) then  ! 晴空
  call set_bit(testbits,2) 
else if(confdnc .gt. 0.66) then  ! 可能晴空
  call set_bit(testbits,1) 
end if 
! confdnc <= 0.66 → 有云
```

只有 3 个切分点（0.66、0.95、0.99），将连续值强制二值化。在地表光谱不均匀区域（混合植被/裸土），相邻像素 confdnc 在 0.93~0.97 之间波动 → 判为不同类别。

**结论**: 算法设计层面的根因。MODIS MOD35 原始算法依赖 250m→1km 子像素聚合来平滑，Fortran 缺失这一层。

---

## 原因 3：chk_land 恢复测试独立翻转像素

**文件**: `src/fortran/cloudmask/chk_land.f90:79-82`, `src/fortran/cloudmask/land_module.f90:79-82`

```fortran
if(.not. (snow .or. ice)) then
  if(confdnc .le. 0.95) then
     call chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,testbits)
  end if
end if
```

**问题**:
1. 如果 IR 通道通过晴空测试 + 11μm BT 足够高 → 直接改 confdnc 为 0.96（晴空），不考虑邻域
2. md1（3.8-3.959μm）检查被注释掉（line 164-165），恢复条件更宽松
3. 被周围云包围的像素如果恰巧通过 IR 测试，会被"恢复"为晴空 → 形成孤立晴空点（椒盐噪声中的"盐"）

**结论**: 椒盐噪声中"盐"（孤立晴空点）的主要来源。

---

## 原因 4：PFMFT/NFMFT 测试被注释掉的组贡献 — ⚠️ 最严重

这是一个**跨 8 个文件**的系统性问题：

| 文件 | PFMFT cmin1 | NFMFT cmin1 | 影响 |
|------|------------|------------|------|
| `LandDay.f90` (line 208, 227) | ✅ 活跃 | ✅ 活跃 | 白天陆地正常 |
| `LandDay_coast.f90` (line 182, 201) | ✅ 活跃 | ✅ 活跃 | 白天海岸正常 |
| **`LandNite.f90`** (line 194, 212) | ❌ 注释 | ❌ 注释 | 夜间陆地 G1 只有表面温度测试 |
| **`PolarNite_land.f90`** (line 184, 202) | ❌ 注释 | ❌ 注释 | 极夜陆地同上 |
| **`Nite_snow.f90`** (line 168, 186) | ❌ 注释 | ❌ 注释 | 夜间雪面 G1 严重削弱 |
| **`PolarNite_snow.f90`** (line 175, 195) | ❌ 注释 | ❌ 注释 | 极夜雪面同上 |
| **`PolarDay_snow.f90`** (line 163, 184) | ❌ 注释 | ❌ 注释 | 极昼雪面 G1 严重削弱 |
| **`Antarctic_day.f90`** (line 164, 182) | ❌ 注释 | ❌ 注释 | 南极白天同上 |

**LandNite.f90 典型示例** (line 194-195, 212-213):
```fortran
! PFMFT test — 测试代码完整执行，c1 被计算，但：
  !      cmin1 = min(cmin1,c1)       ← 被注释掉了
  !      ngtests(1) = ngtests(1) + 1 ← 被注释掉了

! NFMFT test — 同理：
    !    cmin1 = min(cmin1,c2)       ← 被注释掉了
    !    ngtests(1) = ngtests(1) + 1 ← 被注释掉了
```

**后果**：夜间路径的 Group 1（高厚云测试组）只有表面温度测试（c9）贡献置信度 → 置信度系统性偏高 → 夜间噪声比白天更严重。

**额外问题**：LandDay.f90 line 329，11-4μm BTD 测试的 cmin2 也被注释掉：
```fortran
!         cmin2 = min(cmin2,c4)       ← 被注释掉了
!         ngtests(2) = ngtests(2) + 1 ← 被注释掉了
```

**结论**: 置信度计算的结构性缺陷，影响所有非白天陆地路径。

---

## 原因 5：空间均匀性测试条件太严格

**文件**: `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90:1041-1207`

`check_reg_uniformity` 的 3×3 窗口要求 **9 个像素全部满足**：
1. eco_type 完全一致
2. land/water/coast 分类完全一致
3. snow 状态完全一致

额外限制：
- 边缘像素 → 直接 `uniform = .false.`
- 雪/冰中心像素 → 直接 `uniform = .false.`
- 海岸线混合 → 直接 `uniform = .false.` + 强制 coast=true

**后果**：地表过渡带、海岸线、雪线等区域的 uniform 几乎永远为 false，chk_spatial_var（唯一利用空间信息的机制）在这些区域完全不起作用。

**结论**: 空间上下文利用不足，最需要约束的区域反而最缺乏约束。

---

## 原因 6：无后处理滤波

**Fortran 独立主程序** (`fylat_fy3mersi_cloud_mask.f90`):
- 主循环结束后直接输出，没有任何后处理
- 无多数票滤波、无形态学操作、无连通域分析

**Fortran C API** (`cloudmask_c_api.f90:477-481`):
- `commit e2b2edc` (2026-06-16) 新增 `apply_spatial_consistency`
- 仅翻转 0/8 匹配的孤立像素（不是真正的 median filter）
- 独立 Fortran 程序路径不受益

**Python**:
- 完全没有空间一致性后处理
- `scripts/spatial_analysis.py` 中的 `numpy_salt_pepper` 仅用于诊断

**结论**: 独立 Fortran 程序路径缺少兜底机制。

---

## 噪声产生机制总结

```
逐像素独立处理（原因1）
    +
S曲线硬阈值二值化（原因2）
    +
chk_land 独立恢复（原因3）
    +
PFMFT/NFMFT组贡献缺失（原因4）→ 置信度偏高，夜间尤甚
    +
uniform条件过严，空间均匀性测试不生效（原因5）
    ↓
椒盐噪声
    ↓
无后处理滤波（原因6，独立Fortran程序路径）→ 噪声直接输出
```

## 修复优先级

| 优先级 | 原因 | 修复方案 | 影响范围 |
|--------|------|---------|---------|
| P0 | 原因 4 | 恢复 8 个文件中被注释的 cmin/ngtests 行 | 夜间 + 雪面路径 |
| P1 | 原因 6 | 在独立 Fortran 程序中也加入 `apply_spatial_consistency` | 所有路径 |
| P2 | 原因 3 | chk_land 增加邻域一致性检查 | 白天陆地 |
| P3 | 原因 5 | 放宽 uniform 条件（允许 1-2 个不一致邻域） | 海洋 |
| P4 | 原因 1/2 | 结构性改进（需更大范围重构） | 所有路径 |
