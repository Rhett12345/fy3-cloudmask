# FY-3D 云检测代码 Bug 完整汇总

**当前版本**: v3.5.0
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

---

## 椒盐噪声根因诊断 (v3.5.0 确定性分析)

### 现象

Fig 2 MYD35 验证统计：Agreement 46.5%，HSS 0.085
- truth:Clr -> pred:Cld = 41%（晴空误判为云）
- truth:Cld -> pred:Clr = 38%（云误判为晴空）
- 双向对称误判 -> 算法置信度在 0.5 附近随机跳动

### 根因 1：find_matched_files.py 配置错误导致 NWP 永远不加载 [CRITICAL]

**文件**: scripts/find_matched_files.py

两个配置错误：
```python
MERSI_ROOT = Path('/data/Data_yuq/mersi_test')   # 错误！应为 /data/Data_yuq/mersi
NWP_PATTERN = 'fnl_{date}_{hh}_00'               # 错误！应为 gfs0p25_41L_{date}_{hh}_00
```

实际数据位置：
- MERSI: `/data/Data_yuq/mersi/{date}/FY3D_MERSI_GBAL_L1_{date}_{time}_1000M_MS.HDF`
- NWP: `/data/nwp/{date}/ORG/gfs0p25_41L_{date}_{hh}_00`（20220803 等日期有预处理好的二进制格式）

连锁效应：
1. `find_matched_triplets()` 对所有日期返回空列表
2. `run_fortran_only.py` 中 NWP 数据从未加载
3. btclr 全部传零 → 触发后续所有级联失败

**这是所有问题的源头。修复此文件后，btclr 从 NWP 真实数据读取，不需要 sfctmp 估算。**

### 根因 2：btclr 全零导致 PFMFT Group 1 失效 [CRITICAL]

**文件**: scripts/run_fortran_only.py
```python
btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
```

**LandDay.f90** 的 PFMFT 触发条件：
```fortran
(btclr(5)-btclr(6)) > pfmft_btd_min(1)   ! pfmft_btd_min = 0.5
```

btclr 全零 → 0.0 > 0.5 永远为假 → LandDay 的 Group 1 永远没有测试执行 → ngtests(1) = 0，Group 1 不参与几何平均。

**LandNite.f90** 同理，btclr(5)-btclr(6) > 0.5 永远为假。

### 根因 3：ocean 路径 PFMFT conf_test 被注释 [CRITICAL]

**ocean_day.f90** 第240-243行：
```fortran
!     call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4),   &
!                    pfmft_ocean(2),1,c2)
!   cmin1 = min(cmin1,c2)          ! annotation by minmin (to close this threshold)
!   ngtests(1) = ngtests(1) + 1
```

ocean 的 PFMFT 触发条件用的是 `(masir11-masir12) < pfmft_btd_min(1)`，不依赖 btclr，所以能触发。
但 conf_test 调用和 cmin1/ngtests 更新全部被注释 → testbit 设了，nmtests 加了，但置信度贡献为零。

**ocean_nite.f90** 同理。

### 根因 4：snow_mask 全零导致雪面误判为云 [HIGH]

**文件**: scripts/run_fortran_only.py
```python
snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
```

雪/冰检测完全依赖 NDSI 可见光方法。夜间和高纬度冬季无法检测积雪，雪面高反射率被可见光测试判为云。

---

## 关于 APOLLO 的纠正说明

goals.md 旧版本认为 ocean APOLLO 测试会产生垃圾阈值，需要禁用（FIX-D）。**这是错误的。**

APOLLO 有内置回退机制：
```fortran
! ocean_day.f90 第430行
if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
  dfthrsh = do11_12hi(1)    ! 回退到静态阈值 3.0K
else
  dfthrsh = diftemp          ! 用 APOLLO 自适应阈值
end if
```

当 btclr=0 时，`tview` 查找表返回 `diftemp < 0.1`，自动回退到静态阈值 3.0K。APOLLO 不会产生垃圾值。禁用 APOLLO 反而会让 ocean 路径丢失自适应阈值能力。

LandDay/LandNite 已用 `.false.` 门控禁用 APOLLO，这也是不必要的（但影响较小，因为陆地有静态阈值兜底）。

---

## 其他已知 Bug

### 已修复 (v3.4.3~v3.5.0)

| BUG | 描述 | 修复版本 |
|-----|------|----------|
| BUG-1 | IR_WAVENUMBERS[4] 从 933.364 改为 909.458 | v3.4.4 |
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

### 未修复

| BUG | 描述 | 严重度 | 状态 |
|-----|------|--------|------|
| 根因1 | find_matched_files.py 配置错误 | CRITICAL | ✅ 已修复 |
| 根因2 | btclr 全零（根因1的下游后果） | CRITICAL | ✅ 已修复 (sfctmp估算) |
| 根因3 | ocean PFMFT conf_test 注释 | CRITICAL | ✅ 已修复 |
| 根因4 | snow_mask 全零 | HIGH | 缺 NISE 数据 |
| BUG-16 | nfmft 阈值范围过宽（-23~-22K，实际值 -1~1K） | MEDIUM | 需重新标定 |
| BUG-17 | ocean_nite 11-4um 阈值过窄 (2.25K) | MEDIUM | 需重新标定 |
| BUG-18 | GEMI 公式除零风险 (s1→1.0) | MEDIUM | 待修 |

---

## 修复计划

### v3.5.1 (消除椒盐噪声 - 当前最高优先级)

**目标**: 修复根因，将 Agreement 从 46% 提升到 65%+

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-1 | 修正 find_matched_files.py 的 MERSI_ROOT 和 NWP_PATTERN | scripts/find_matched_files.py |
| FIX-2 | ocean_day PFMFT 取消 conf_test 注释 | ocean_day.f90 |
| FIX-3 | ocean_nite PFMFT 取消 conf_test 注释 | ocean_nite.f90 |

#### FIX-1: 修正文件匹配配置

```python
# scripts/find_matched_files.py
MERSI_ROOT = Path('/data/Data_yuq/mersi')    # 修正路径
NWP_PATTERN = 'gfs0p25_41L_{date}_{hh}_00'   # 修正文件名模式
```

修复后：
- `find_matched_triplets()` 能正确找到 L1B+GEO+NWP 三元组
- NWP 数据正常加载，btclr 从 NWP 真实数据读取
- PFMFT 触发条件 btclr(5)-btclr(6) > 0.5 正常满足
- APOLLO 查找表拿到合理参考温度
- Group 1 重新参与置信度计算

#### FIX-2/3: ocean PFMFT 置信度恢复

```fortran
! ocean_day.f90 和 ocean_nite.f90，取消注释：
call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4), &
               pfmft_ocean(2),1,c2)
cmin1 = min(cmin1,c2)
ngtests(1) = ngtests(1) + 1
```

**预期效果**:
- Group 1 在所有路径都产生置信度贡献
- Confidence 分布从 0.5 附近集中变为更明确的 0/1 分布
- Agreement 预期从 46% 提升到 65%+

### v3.5.2 (补充数据源)

| 修复项 | 描述 | 依赖 |
|--------|------|------|
| FIX-4 | snow_mask 传入真实 NISE 数据 | NISE 雪冰掩码产品 |
| FIX-5 | GEMI 除零保护 | 代码修改 |

### v3.6.0 (阈值重新标定)

| 修复项 | 描述 | 依赖 |
|--------|------|------|
| FIX-6 | nfmft 阈值重新标定 (-23~-22K -> 合理范围) | MERSI-II 统计分析 |
| FIX-7 | ocean_nite 11-4um 阈值调宽 (2.25K -> 5K+) | 验证数据 |

---

## 验证方案

### 测试数据

**日期**: 2022-08-03
- MERSI: `/data/Data_yuq/mersi/20220803/`（3 轨道：0740, 0830, 0920）
- NWP: `/data/nwp/20220803/ORG/gfs0p25_41L_20220803_{00..24}_00`（9 时次，1.1GB/个）
- MYD35: `/data/Data_yuq/aqua_modis/MYD35_L2/20220803/`（9 时次覆盖 MERSI 时刻）

### 快速验证 (v3.5.1 后)

1. 运行 `python scripts/run_fortran_only.py --date 20220803`
2. 检查日志中 NWP 是否正确加载（不应跳过）
3. 检查 btclr(5) 是否为合理温度值（250-310K，非零）
4. 检查 PFMFT 触发率是否 > 0%
5. 检查 Confidence 分布是否不再集中在 0.5 附近
6. 与 MYD35 对比：Agreement 是否提升到 65%+

### 完整验证

- 7 个测试日期全轨道验证
- 与 MYD35 对比：Agreement、HSS、混淆矩阵
- 空间分布图：检查椒盐是否消除
