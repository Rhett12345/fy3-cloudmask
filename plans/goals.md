# FY-3D 云检测代码 Bug 完整汇总

**当前版本**: v3.5.0
**GitHub 状态**: main 分支，已推送
**核心目标**: 消除椒盐噪声，提升与 MYD35 的一致性

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

### 根因 1：PFMFT 触发条件导致 Group 1 完全失效 [CRITICAL]

**LandDay.f90** 的 PFMFT 触发条件：
```fortran
(btclr(5)-btclr(6)) > pfmft_btd_min(1)   ! pfmft_btd_min = 0.5
```

btclr 全部传零（见根因 2），所以 0.0 > 0.5 永远为假。
-> LandDay 的 Group 1（热红外裂窗测试）永远没有测试执行
-> ngtests(1) = 0，cmin1 = 1.0，Group 1 不参与几何平均

**ocean_day.f90** 的 PFMFT 更糟糕：
- 触发条件用的是 (masir11-masir12) < pfmft_btd_min（不用 btclr），所以能触发
- 但 conf_test 调用和 cmin1/ngtests 更新全部被注释掉
- 结果：testbit 设了，nmtests 加了，但置信度贡献为零

**ocean_nite.f90** 的 PFMFT：同 ocean_day，conf_test 调用被注释，cmin1/ngtests 未更新

**LandNite.f90** 的 PFMFT：触发条件用 btclr(5)-btclr(6) > 0.5，同 LandDay 永远不触发

**结论：所有路径的 Group 1（热红外裂窗）实际不产生置信度贡献。**

### 根因 2：btclr 全零导致 APOLLO 和 PFMFT 同时崩溃 [CRITICAL]

**文件**: scripts/run_fortran_only.py
```python
btclr=np.ascontiguousarray(np.zeros((n_elem, n_line, 7), dtype=np.float32)),
```

连锁效应：
1. PFMFT 触发条件 btclr(5)-btclr(6) > 0.5 永远为假 -> Group 1 失效
2. APOLLO 查找表用 btclr(5) 作参考温度 -> 0K 查找返回极端值
3. conf_test 对极端阈值的处理：abs(range) < 1e-12 时返回 c = 0.5
4. 大量像素的 Group 2 confidence 落在 0.5 附近
5. 几何平均后 confidence 在 0.5-0.7 区间随机分布
6. set_confdnc 在 cloudy/prob_cloudy 边界随机跳动 -> 椒盐

### 根因 3：snow_mask 全零导致雪面误判为云 [HIGH]

**文件**: scripts/run_fortran_only.py
```python
snow_mask=np.ascontiguousarray(np.zeros((n_elem, n_line), dtype=np.int8)),
```

雪/冰检测完全依赖 NDSI 可见光方法。夜间和高纬度冬季无法检测积雪，雪面高反射率被可见光测试判为云。

### 根因 4：ocean 路径 APOLLO 测试仍启用但阈值不适用 [HIGH]

APOLLO 查找表为 MODIS 设计（BTD 范围 2-10K），MERSI-II 的 11-12um BTD 范围 0-4K。
- 在 270K/天底时 APOLLO 返回 4.0K 阈值
- 静态阈值是 3.0K
- 3-4K 范围的像素：APOLLO 判晴空，静态判有云 -> 边界噪声

LandDay/LandNite 已用 .false. 禁用 APOLLO，但 ocean_day/ocean_nite 仍然启用。

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
| 根因1 | PFMFT Group1 全路径失效 | CRITICAL | 待修 |
| 根因2 | btclr 全零 | CRITICAL | 待修 |
| 根因3 | snow_mask 全零 | HIGH | 缺数据 |
| 根因4 | ocean APOLLO 阈值不适用 | HIGH | 待修 |
| BUG-16 | ocean_day PFMFT conf_test 注释 | CRITICAL | 待修 |
| BUG-17 | ocean_nite PFMFT conf_test 注释 | CRITICAL | 待修 |
| BUG-18 | GEMI 公式除零风险 (s1->1.0) | MEDIUM | 待修 |
| BUG-19 | ocean_nite 11-4um 阈值过窄 (2.25K) | MEDIUM | 待标定 |

---

## 修复计划

### v3.5.1 (消除椒盐噪声 - 当前最高优先级)

**目标**: 修复 4 个根因，将 Agreement 从 46% 提升到 65%+

| 修复项 | 描述 | 文件 |
|--------|------|------|
| FIX-A | btclr 用 sfctmp 估算填充 | scripts/run_fortran_only.py |
| FIX-B | ocean_day PFMFT 取消 conf_test 注释 | ocean_day.f90 |
| FIX-C | ocean_nite PFMFT 取消 conf_test 注释 | ocean_nite.f90 |
| FIX-D | ocean 路径禁用 APOLLO (加 .false. 门控) | ocean_day.f90, ocean_nite.f90 |

#### FIX-A: btclr 估算方案

```python
# scripts/run_fortran_only.py
sfctmp_arr = nwp_interp['tsfc'].astype(np.float32)
btclr_arr = np.zeros((n_elem, n_line, 7), dtype=np.float32)
btclr_arr[:, :, 4] = sfctmp_arr          # btclr(5): 11um 晴空BT ~ 地表温度
btclr_arr[:, :, 5] = sfctmp_arr - 1.0    # btclr(6): 12um 晴空BT ~ sfctmp - 1K
btclr_arr[:, :, 0] = sfctmp_arr - 25.0   # btclr(1): 3.8um 晴空BT
```

效果：
- btclr(5) - btclr(6) = 1.0 > 0.5 -> PFMFT 正常触发
- APOLLO 查找表拿到合理参考温度 -> 阈值不再返回垃圾值
- Group 1 重新参与置信度计算

#### FIX-B/C: ocean PFMFT 置信度恢复

```fortran
! ocean_day.f90 和 ocean_nite.f90，取消注释：
call conf_test(tv11_12,pfmft_ocean(1),pfmft_ocean(3),pfmft_ocean(4), &
               pfmft_ocean(2),1,c2)
cmin1 = min(cmin1,c2)
ngtests(1) = ngtests(1) + 1
```

#### FIX-D: ocean APOLLO 禁用

```fortran
! ocean_day.f90 和 ocean_nite.f90，在 APOLLO 调用前加门控：
if (.false.) then   ! APOLLO disabled for MERSI-II
  call tview(1,schi,r24,diftemp)
  ...
end if
dfthrsh = do11_12hi(1)   ! 直接用静态阈值
```

**预期效果**:
- Group 1 重新产生置信度贡献（pfmft 测试生效）
- Group 2 的 11-12um 测试不再受 APOLLO 垃圾阈值干扰
- Confidence 分布从 0.5 附近集中变为更明确的 0/1 分布
- Agreement 预期从 46% 提升到 65%+

### v3.5.2 (补充数据源)

| 修复项 | 描述 | 依赖 |
|--------|------|------|
| FIX-E | snow_mask 传入真实 NISE 数据 | NISE 雪冰掩码产品 |
| FIX-F | btclr 改用 RTM 计算值替代 sfctmp 估算 | NWP RTM 模块 |

### v3.6.0 (阈值重新标定)

| 修复项 | 描述 | 依赖 |
|--------|------|------|
| FIX-G | nfmft 阈值重新标定 (-23~-22K -> 合理范围) | MERSI-II 统计分析 |
| FIX-H | ocean_nite 11-4um 阈值调宽 (2.25K -> 5K+) | 验证数据 |
| FIX-I | GEMI 除零保护 | 代码修改 |
| FIX-J | APOLLO 表为 MERSI-II 重新标定 | 统计分析 |

---

## 验证方案

### 快速验证 (v3.5.1 后)

1. 在 process_pixel_c 中添加诊断输出：
```fortran
if (ielem_in == 1000 .and. iline_in == 1000) then
    write(0,*) 'DEBUG: confdnc=', confdnc, ' mask=', out_mask
    write(0,*) 'DEBUG: btclr(5)=', btclr_in(5), ' btclr(6)=', btclr_in(6)
end if
```

2. 运行 2020-03-08 数据，检查：
   - btclr(5) 不再为 0（应为 250-310K）
   - PFMFT 触发率 > 0%（之前为 0%）
   - Confidence 分布不再集中在 0.5 附近
   - Agreement 提升到 65%+

### 完整验证

- 7 个测试日期（20220803~20250302）全轨道验证
- 与 MYD35 对比：Agreement、HSS、混淆矩阵
- 空间分布图：检查椒盐是否消除
