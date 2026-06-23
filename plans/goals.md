# FY-3D 云检测代码 Bug 完整汇总

**当前版本**: v3.5.4
**GitHub 状态**: main 分支，已推送
**核心目标**: 消除椒盐噪声，提升与 MYD35 的一致性
**测试日期**: 2022-08-03（MERSI + NWP binary + MYD35 三者齐全）

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

---

## 当前性能 (v3.5.4 vs MYD35, 2022-08-03)

```
  Accuracy: 49.98%  |  HSS: -0.010  |  FY3_CF: 96%  |  MYD35_CF: 50%
  Ocean BTD median: +0.66K (MODIS reference: +0.5~+3.0K) ✓
```

BTD 已修正到物理合理范围，但云检测精度仍 ~50%。**阈值重标定是下一阶段的核心任务。**

---

## 椒盐噪声根因诊断 (v3.5.0 确定性分析)

### 现象

Fig 2 MYD35 验证统计：Agreement 46.5%，HSS 0.085
- truth:Clr -> pred:Cld = 41%（晴空误判为云）
- truth:Cld -> pred:Clr = 38%（云误判为晴空）
- 双向对称误判 -> 算法置信度在 0.5 附近随机跳动

### 根因 1：find_matched_files.py 配置错误导致 NWP 永远不加载 [CRITICAL] ✅ 已修复

**文件**: scripts/find_matched_files.py

两个配置错误：
```python
MERSI_ROOT = Path('/data/Data_yuq/mersi_test')   # 错误！应为 /data/Data_yuq/mersi
NWP_PATTERN = 'fnl_{date}_{hh}_00'               # 错误！应为 gfs0p25_41L_{date}_{hh}_00
```

实际数据位置：
- MERSI: `/data/Data_yuq/mersi/{date}/FY3D_MERSI_GBAL_L1_{date}_{time}_1000M_MS.HDF`
- NWP: `/data/nwp/{date}/ORG/gfs0p25_41L_{date}_{hh}_00`

**已修复**: v3.5.1

### 根因 2：btclr 全零导致 PFMFT Group 1 失效 [CRITICAL] ✅ 已修复

**文件**: scripts/run_fortran_only.py

btclr 全部传零 → btclr(5)-btclr(6)=0 > 0.5 永远为假 → LandDay/LandNite 的 Group 1 永远不触发。

**已修复**: v3.5.1 (sfctmp 估算), v3.5.4 (按表面类型区分: Ocean BTD_clr=0.5K, Land BTD_clr=1.9K)

### 根因 3：ocean 路径 PFMFT conf_test 被注释 [CRITICAL] ✅ 已修复

**ocean_day.f90** / **ocean_nite.f90** PFMFT conf_test 调用和 cmin1/ngtests 更新被注释。

**已修复**: v3.5.1

### 根因 4：snow_mask 全零 [HIGH] 未修复

```python
snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
```

雪/冰检测完全依赖 NDSI 可见光方法。夜间和高纬度冬季无法检测积雪。依赖: NISE 雪冰掩码产品。

---

## 新发现的 Bug (v3.5.1~v3.5.4 诊断过程中发现)

### BUG-A：IR 波数硬编码错误 [CRITICAL] ✅ 已修复

**文件**: scripts/run_fortran_only.py, src/fy3_cloudmask/constants.py

| 通道 | 文档 λ(μm) | 正确 ν(cm⁻¹) | 原代码 ν(cm⁻¹) | 偏差 |
|------|-----------|-------------|---------------|------|
| B24  | 10.714    | 933.358     | 909.458       | -23.9 |
| B25  | 11.948    | 836.960     | 836.941       | -0.02 |

**修复**: v3.5.2 — 从 HDF5 `Effect_Center_WaveLength` 属性实时读取

### BUG-B：缺少 TBB 线性修正 [HIGH] ✅ 已修复

**文件**: scripts/run_fortran_only.py

公式: Tbb = A * Te + B，系数从 HDF5 `TBB_Trans_Coefficient_A/B` 读取。

老代码使用 TCS/TCI 非线性修正（已废弃）。两个修正叠加使 Ocean BTD 从 +0.11K → +0.66K。

**修复**: v3.5.2

### BUG-C：Fortran planck_module.f90 死代码 [LOW] ✅ 已清理

planck_module 在云检测流程中从未被调用（BT 转换全部在 Python 端完成）。

**清理**: v3.5.3

### BUG-D：APOLLO tview 查找表 100% 回退 [HIGH] 未修复

tview 返回值恒 < 0.1，导致 APOLLO 全部回退到静态阈值 dfthrsh=3.0K。
输入参数（sec(VZA)、BT11）均在表范围内，怀疑 tview 表数据未正确加载。

**影响**: ocean_day/ocean_nite Group 2 的 11-12um APOLLO 测试失效（全部用静态阈值 3.0K）

### BUG-E：PFMFT 阈值与 MERSI-II 实际分布不匹配 [HIGH] 未修复

PFMFT ocean 阈值: locut=1.7K, hicut=1.9K
实际 PFMFT 信号: P50=-0.76K, P90=+1.31K
→ 93% 海洋像素 PFMFT 触底 (c=0.1)

**需要**: 阈值重标定为 [0.2, 0.6]K 量级

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
| BUG-1 | IR_WAVENUMBERS 硬编码错误 (909.458→933.358) | v3.5.2 |
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
| — | find_matched_files.py MERSI_ROOT/NWP_PATTERN 错误 | v3.5.1 |
| — | btclr 全零 (sfctmp 估算 → 表面类型区分) | v3.5.1, v3.5.4 |
| — | ocean_day/ocean_nite PFMFT conf_test 注释 | v3.5.1 |
| — | 缺少 TBB_Trans 线性修正 | v3.5.2 |
| — | Fortran planck_module 死代码 | v3.5.3 |

### 未修复

| BUG | 描述 | 严重度 | 状态 |
|-----|------|--------|------|
| 根因4 | snow_mask 全零 | HIGH | 缺 NISE 数据 |
| BUG-D | APOLLO tview 表 100% 回退 | HIGH | 待查 tview 表加载 |
| BUG-E | PFMFT 阈值不匹配 (1.7~1.9K vs 实际 -0.8~1.3K) | HIGH | 需重新标定 |
| BUG-16 | nfmft 阈值范围过宽（-23~-22K，实际值 -1~1K） | MEDIUM | 需重新标定 |
| BUG-17 | ocean_nite 11-4um 阈值过窄 (2.25K) | MEDIUM | 需重新标定 |
| BUG-18 | GEMI 公式除零风险 (LandDay s1→1.0) | MEDIUM | 代码修改 |

---

## 修复计划

### v3.6.0 (阈值重新标定 — 当前最高优先级)

**目标**: 将 MERSI-II 阈值从 MODIS 标定值调整到 MERSI-II 实际观测值

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-E1 | PFMFT ocean 阈值 [1.7,1.9] → [0.2,0.6]K | pfmft_nfmft_thr.inc |
| FIX-E2 | PFMFT land 阈值重标定 | pfmft_nfmft_thr.inc |
| FIX-D | APOLLO tview 表加载修复 | tview.f |
| FIX-6 | nfmft 阈值重标定 | pfmft_nfmft_thr.inc |
| FIX-7 | ocean_nite 11-4um 阈值调宽 | ocean_nite_thr.inc (或 coeff/) |

### v3.6.1 (代码安全)

| 修复项 | 描述 |
|--------|------|
| FIX-4 | snow_mask 传入真实 NISE 数据 |
| FIX-5 | GEMI 除零保护 |

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

### 完整验证

- 7 个测试日期（20220803~20250302）全轨道验证
- 与 MYD35 对比：Agreement、HSS、混淆矩阵
- BTD(11-12) 逐通道统计
- 空间分布图：检查椒盐是否消除
