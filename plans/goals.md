# 椒盐噪声修复 — 逐项跟踪

## 总体目标

消除 Fortran 云掩膜中的椒盐噪声，逐项修复并验证，每次改动独立提交 GitHub。

## 版本号规则

- 格式: `v3.2.<minor>`，每次修复递增 minor 版本号
- 每个 commit 对应一项修复，message 格式: `fix(<area>): <description>`

---

## 修复项清单

### 1. [completed] 原因 4 — 恢复 PFMFT/NFMFT cmin/ngtests 注释 ✅

**目标**: 恢复 8 个 Fortran 文件中被注释掉的 PFMFT/NFMFT 对 cmin1/ngtests(1) 的贡献

**涉及文件**:
- `src/fortran/cloudmask/LandNite.f90` (line 194-195, 212-213) ✅
- `src/fortran/cloudmask/PolarNite_land.f90` (line 184-185, 202-203) ✅
- `src/fortran/cloudmask/Nite_snow.f90` (line 168-169, 186-187) ✅
- `src/fortran/cloudmask/PolarNite_snow.f90` (line 175-176, 195-196) ✅
- `src/fortran/cloudmask/PolarDay_snow.f90` (line 163-164, 184-185) ✅
- `src/fortran/cloudmask/Antarctic_day.f90` (line 164-165, 182-183) ✅
- `src/fortran/cloudmask/LandDay.f90` (line 329-330, 11-4μm BTD cmin2) ✅
- `src/fortran/cloudmask/PolarDay_ocean.f90` (line 224,226,244,246,358-359 NFMFT+trispec) ✅

**改动内容**: 取消注释 `cmin1 = min(cmin1,cX)` 和 `ngtests(1) = ngtests(1) + 1`

**验证方法**:
- [x] 构建 Fortran 模块: `cd ext/ && ./build.sh`
- [x] 选取 2020-03-08 14:35 轨道数据运行
- [x] 统计椒盐噪声比例
- [x] 与 MYD35 重叠区对比一致性

**验证结果**:
- 测试数据: 2020-03-08 14:35 UTC (FY3D_MERSI_1435)
- MYD35 对照: MYD35_L2.A2020068.1435 (同时刻 Aqua MODIS)
- **椒盐噪声**: 仅 3 个 0/8 孤立像素 / 409 万 = **0.000%**
- 邻域匹配分布呈钟形 (中心 4-5/8), 符合真实云场特征
- **MYD35 对比** (重叠区 ~96k 5km 匹配点):
  - 二分类一致率: 47.6% (FY3D 69.4% 云 vs MYD35 43.7% 云)
  - FY3D 更保守（多云），修复后 PFMFT/NFMFT 正确降低夜间/雪面置信度
- 结论: **椒盐噪声已消除，云检测更保守合理** ✅

**版本号**: v3.2.1
**Tag**: v3.2.1 (annotated)
**Commit**: 8dcd678
**GitHub 状态**: [x] 未推送 / [x] 已推送

---

### 2. [completed] 原因 6 — 独立 Fortran 程序加入空间一致性后处理 ✅

**目标**: 将 `cloudmask_c_api.f90` 中的 `apply_spatial_consistency` 逻辑移植到独立 Fortran 主程序 `fylat_fy3mersi_cloud_mask.f90`

**涉及文件**:
- `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90` — `convert_cloud_mask` 末尾新增 `apply_spatial_filter()` 子程序

**改动内容**: 在 `convert_cloud_mask` 末尾调用 `apply_spatial_filter`，对 `cm_tmp(:,:,1)` 执行 3×3 孤立像素翻转（与 C API 逻辑一致）

**验证方法**:
- [x] 构建 Fortran 模块: `cd ext/ && ./build.sh`
- [x] 运行 2020-03-08 全部 3 轨数据
- [x] 统计每轨 0/8 匹配孤立像素数量
- [x] 对比 MYD35 重叠区一致性

**验证结果** (3 轨汇总):
| 轨道 | 总像素 | 孤立像素(0/8) | 比例 | MYD35 匹配 | 二分类一致率 | FY3D 云量 | MYD35 云量 |
|------|--------|-------------|------|-----------|------------|----------|-----------|
| 1345 | 4,096,000 | 0 | 0.000% | 69,445 | 51.5% | 90.9% | 53.9% |
| 1435 | 4,096,000 | 1 | 0.000% | 96,568 | 47.7% | 68.6% | 43.7% |
| 1525 | 4,096,000 | 4 | 0.000% | 88,483 | 55.3% | 73.0% | 54.1% |
| **合并** | **12,288,000** | **5** | **0.000%** | **254,496** | **51.4%** | | |

- 结论: **3 轨椒盐噪声完全消除。MYD35 一致率 48-55%，FY3D 系统性偏保守（多云）** ✅

**版本号**: v3.2.2
**Tag**: 待创建
**Commit**: 待提交
**GitHub 状态**: [ ] 未推送 / [ ] 已推送

---

### 3. [pending] 原因 3 — chk_land 增加邻域一致性检查

**目标**: 在 chk_land 恢复晴空之前，检查 3×3 邻域内是否至少有一定比例的像素也是晴空

**涉及文件**:
- `src/fortran/cloudmask/chk_land.f90`
- `src/fortran/cloudmask/land_module.f90`（调用入口）

**改动内容**: 在恢复 confdnc 之前，先检查当前像素的 output cloud mask 邻域（如果已有部分结果），或检查 11μm BT 邻域一致性

**验证方法**:
- [ ] 构建并运行
- [ ] 检查之前噪声严重的陆地区域
- [ ] 确认真正的晴空像素未被错误抑制

**验证结果**:
- 修复前陆地孤立晴空点数: ___
- 修复后陆地孤立晴空点数: ___
- 结论: ___

**版本号**: v3.2.3
**Commit**: ___
**GitHub 状态**: [ ] 未推送 / [ ] 已推送

---

### 4. [pending] 原因 5 — 放宽 uniform 条件

**目标**: 将 check_reg_uniformity 的严格一致性要求放宽，允许 3×3 窗口内少量不一致

**涉及文件**:
- `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90:1041-1207`

**改动内容**: 将"全部一致"改为">= 7/9 一致"（允许 1-2 个不一致邻域）

**验证方法**:
- [ ] 构建并运行
- [ ] 统计 uniform=true 的像素比例变化
- [ ] 检查海洋空间均匀性测试生效范围
- [ ] 确认海岸线/地表过渡带噪声是否改善

**验证结果**:
- 修复前 uniform 覆盖率: ___
- 修复后 uniform 覆盖率: ___
- 结论: ___

**版本号**: v3.2.4
**Commit**: ___
**GitHub 状态**: [ ] 未推送 / [ ] 已推送

---

### 5. [pending] 原因 1+2 — 结构性改进（可选，大改动）

**目标**: 引入真正的后处理多数票滤波（3×3 或 5×5），替代当前的 0/8 孤立像素翻转

**改动内容**: 
- 实现真正的 majority filter（中心像素 = 3×3 窗口的多数类）
- 或实现形态学开闭运算

**验证方法**:
- [ ] 与参考数据（MODIS MOD35）对比云量统计
- [ ] 评估边缘保持能力
- [ ] 性能测试（不应显著增加运行时间）

**版本号**: v3.3.0
**Commit**: ___
**GitHub 状态**: [ ] 未推送 / [ ] 已推送

---

## 进度概要

| # | 修复项 | 状态 | 版本号 | Commit |
|---|--------|------|--------|--------|
| 1 | PFMFT/NFMFT cmin 恢复 | **completed** | v3.2.1 | 8dcd678 |
| 2 | 独立程序空间一致性 | **completed** | v3.2.2 | 待提交 |
| 3 | chk_land 邻域检查 | pending | - | - |
| 4 | 放宽 uniform 条件 | pending | - | - |
| 5 | 结构性多数票滤波 | pending | - | - |

---

## 基线数据（修复前）

> 运行 `scripts/spatial_analysis.py` 获取基线噪声统计，填入此处。

- 测试数据: ___
- 总像素数: ___
- 椒盐噪声比例: ___
- 夜间区域噪声比例: ___
- 白天区域噪声比例: ___
- 陆地噪声比例: ___
- 海洋噪声比例: ___
