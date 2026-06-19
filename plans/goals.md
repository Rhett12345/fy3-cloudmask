# 云检测算法修复目标

**创建时间**: 2026-06-17
**当前版本**: v3.3.3 (commit `47ae0dc`, branch `main`, tag v3.3.3)
**基准代码**: `retrieval_system_V3.1_cldmask/src/cloudmask/`

---

## 版本规划总览

每个子任务独立分支、独立版本号。流程：分支 → 实现 → 验证 → 合并 main → 打 tag。

| 版本 | Task | 内容 | 优先级 |
|------|------|------|--------|
| v3.3.4 | 1.1 | 恢复 LandDay_desert 测试组（pfmft/nfmft） | CRITICAL |
| v3.3.5 | 1.2 | 恢复 chk_land/chk_coast restoral 置信度上限 | HIGH |
| v3.3.6 | 1.3 | 恢复 nmtests/ngtests 计数位置 | MEDIUM |
| v3.3.7 | 1.4 | 移除 chk_land restoral 空间均匀性门控 | MEDIUM |
| v3.3.8 | 1.5 | 恢复 check_reg_uniformity 原始容差 | LOW |
| v3.4.0 | 2.2 | 添加空间一致性后处理过滤器 | CONDITIONAL |
| v3.4.1 | 2.3 | 放松海洋空间变异性测试条件 | CONDITIONAL |

> Task 1.6（置信度下限）和 Task 1.7（PolarDay 保护）不产生新版本，详见下方说明。

---

## Goal 1: 恢复与原始代码的计算一致性

> 优先级最高。新 Fortran 后端相对于原始代码存在 9 处计算差异，需逐一修复。

### Task 1.1 → v3.3.4: 恢复 LandDay_desert / LandDay_desert_c 测试组结构 [CRITICAL]

- **问题**: 原始代码有 4 组测试（含 pfmft/nfmft），新版重构为 3 组，丢失了 pfmft/nfmft 测试
- **影响**: 沙漠区域云检测能力下降，大量本应检测到的云被漏判
- **文件**: `src/fortran/cloudmask/LandDay_desert.f90`, `LandDay_desert_c.f90`
- **操作**:
  1. 恢复 `include 'pfmft_nfmft_thr.inc'`
  2. 恢复 Group 1: pfmft + nfmft (11-12um BTD)
  3. 将 0.86um 和 1.38um 恢复为独立组（Group 3 和 Group 4）
  4. 恢复相关变量声明（`tv11_12`, `cosvza`, `schi` 等）
- **验证**: 对比原始代码的 4 组结构，逐行确认一致
- **分支**: `fix/restore-desert-test-groups`
- **tag**: `v3.3.4`

### Task 1.2 → v3.3.5: 恢复 chk_land / chk_coast restoral 置信度上限 [HIGH]

- **问题**: 原始代码 `confdnc = 1.0`，新版改为 `0.97`
- **影响**: 恢复后的像素从"确定晴空"(>0.99) 降为"可能晴空"(0.95-0.99)，编码结果不同
- **文件**: `src/fortran/cloudmask/chk_land.f90`, `chk_coast.f90`
- **操作**: 将 `confdnc = 0.97` 改回 `confdnc = 1.0`
- **分支**: `fix/restore-restoral-confidence`
- **tag**: `v3.3.5`

### Task 1.3 → v3.3.6: 恢复 nmtests/ngtests 计数位置 [MEDIUM]

- **问题**: 原始代码在条件判断**前**无条件递增 `nmtests`，新版改为条件通过后递增
- **影响**: `nmtests` 值偏小，影响 `fill_bit_pixel` 质量等级判定（`<3` → 质量 4，`<7` → 质量 6，`>=7` → 质量 7）
- **文件**: `src/fortran/cloudmask/LandDay.f90`, `ocean_day.f90`, `LandDay_coast.f90`
- **操作**: 将 `nmtests = nmtests + 1` 移回 `if` 块之前
- **分支**: `fix/restore-nmtests-count`
- **tag**: `v3.3.6`

### Task 1.4 → v3.3.7: 移除 chk_land restoral 空间均匀性门控 [MEDIUM]

- **问题**: 新增了 BT11 3x3 标准差 < 1.5K 的门控，原始代码无此限制
- **影响**: 云边界区域的晴空恢复被跳过，加剧边界噪声
- **文件**: `src/fortran/cloudmask/land_module.f90`
- **操作**: 移除 `bt11_std .lt. 1.5` 条件判断及其相关计算代码
- **分支**: `fix/remove-spatial-gate`
- **tag**: `v3.3.7`

### Task 1.5 → v3.3.8: 恢复 check_reg_uniformity 原始容差 [LOW]

- **问题**: 原始代码任一不一致即 `uniform = .false.`，新版允许最多 2 个不一致
- **影响**: 更多像素被标记为 uniform，扩大了 chk_spatial_var 和 restoral 的执行范围
- **文件**: `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90` (check_reg_uniformity 子程序)
- **操作**: 将 `nmismatch .gt. 2` 改回 `nmismatch .gt. 0`
- **分支**: `fix/restore-uniformity-tolerance`
- **tag**: `v3.3.8`

### Task 1.6: 置信度下限 0.1 [暂不处理]

- **问题**: 新版对所有路径添加了 `max(cminX, 0.1)` 下限，原始代码无此限制
- **决策**: **暂不移除**。此改动解决了原始代码的缺陷（零置信度导致 48% 孤立像素）。待 Goal 2 空间滤波器生效后，再评估是否移除
- **不产生新版本**

### Task 1.7: PolarDay btclr 零值保护 [保留]

- **问题**: 新版添加了 `(btclr(5) .ne. 0.0 .or. btclr(6) .ne. 0.0)` 条件，原始代码无此保护
- **决策**: **保留**。合理的防御性编程，不改变正常路径行为
- **不产生新版本**

---

## Goal 2: 消除椒盐噪声（向 MYD35 靠齐）

> 在 Goal 1 全部完成后执行。先诊断 Goal 1 修复后椒盐噪声是否仍存在。

### 诊断: 评估 Goal 1 修复效果

- **操作**:
  1. 用 2020-03-08 数据运行 v3.3.8 Fortran 版本
  2. 运行 `scripts/validate_spatial.py` 统计孤立像素比例
  3. 对比 v3.3.3 的指标（fully_isolated: 16.5%, near_isolated: 44.5%）
- **判断标准**:
  - `fully_isolated` < 5% → 椒盐噪声已解决，跳过 Task 2.2/2.3
  - `fully_isolated` 5-15% → 只需 Task 2.2
  - `fully_isolated` > 15% → 需要 Task 2.2 + 2.3

### Task 2.2 → v3.4.0: 添加空间一致性后处理过滤器 [CONDITIONAL]

- **前提**: 诊断结果 `fully_isolated` >= 5%
- **文件**: `src/fortran/cloudmask/fylat_fy3mersi_cloud_mask.f90`
- **操作**: 在主循环结束后添加 MYD35 风格的 3x3 多数投票后处理
- **逻辑**:
  - 遍历所有像素（跳过边缘 1 像素）
  - 统计 3x3 邻域 cloud mask 各等级计数
  - 如果当前像素与周围 >= 7/8 邻域不一致，修正为多数类别
  - 保留原始 testbits 供追溯
- **验证**: 统计修复前后 `salt_pepper.fully_isolated` 变化
- **分支**: `fix/spatial-consistency-filter`
- **tag**: `v3.4.0`

### Task 2.3 → v3.4.1: 放松海洋空间变异性测试条件 [CONDITIONAL]

- **前提**: 诊断结果 `fully_isolated` >= 15%（或 Task 2.2 后海洋区域仍有噪声）
- **文件**: `src/fortran/cloudmask/water_module.f90`
- **操作**: 在 `water_day` 和 `water_nite` 中移除 `uniform` 门控
- **分支**: `fix/relax-ocean-spatial-test`
- **tag**: `v3.4.1`

---

## 分支管理

```
main (v3.3.3)
 ├── fix/restore-desert-test-groups      → merge → tag v3.3.4
 ├── fix/restore-restoral-confidence     → merge → tag v3.3.5
 ├── fix/restore-nmtests-count           → merge → tag v3.3.6
 ├── fix/remove-spatial-gate             → merge → tag v3.3.7
 ├── fix/restore-uniformity-tolerance    → merge → tag v3.3.8
 ├── fix/spatial-consistency-filter      → merge → tag v3.4.0 (conditional)
 └── fix/relax-ocean-spatial-test        → merge → tag v3.4.1 (conditional)
```

### 提交规范

```
fix(scope): description

- 详细说明修改内容
- 引用原始代码对应位置
- 引用 plans/original_vs_new_diff.md 差异编号
```

### 合并流程

1. 从 `main` 拉分支
2. 实现修改 → 本地编译验证
3. 用 2020-03-08 数据运行，确认无回归
4. PR 合并到 `main`
5. 打 tag（如 `v3.3.4`）
6. 继续下一个 Task

### 回退策略

如果某个版本引入回归：
- `git revert <commit>` 或 `git reset --hard <tag>` 回到上一个版本
- 单独回退该 Task，不影响其他已合并的修复

---

## 进度跟踪

| 版本 | Task | 状态 | 分支 | 验证结果 |
|------|------|------|------|----------|
| v3.3.4 | 1.1 恢复沙漠测试组 | ✅ 已完成 | main | 编译通过，恢复 pfmft/nfmft + 4组结构，保留 0.1 下限和 groups==0 防御 |
| v3.3.5 | 1.2 恢复 restoral 上限 | ✅ 已完成 | main | 编译通过，chk_land/chk_coast confdnc 0.97→1.0 |
| v3.3.6 | 1.3 恢复 nmtests 计数 | ✅ 已完成 | main | 编译通过, LandDay/ocean_day/LandDay_coast nmtests 5+5+4 处恢复 |
| v3.3.7 | 1.4 移除空间门控 | ✅ 已完成 | main | land_module day/night 移除 BT11 std < 1.5K gate |
| v3.3.7 | 1.5 恢复 uniformity 容差 | ✅ 已完成 | main | nmismatch > 2 → nmismatch > 0 |
| v3.3.7 | 1.6 移除置信度下限 | ✅ 已完成 | main | 19 个文件 max(cminX,0.1) → cminX |
| — | 1.7 PolarDay 保护 (保留) | ✅ | - | - |
| v3.4.0 | 2.2 空间一致性滤波 | ⬜ 待定 | - | 诊断后决定 |
| v3.4.1 | 2.3 海洋空间测试 | ⬜ 待定 | - | 诊断后决定 |
