# 项目 5 执行笔记

## 阶段 0 进展

### ✅ 环境
- `~/Final Project/p5-vcc-null-audit`（临时位置，等连了专用文件夹再 `mv`）
- uv venv, Python 3.12.16；`cell-eval2 0.22.0` + numba 0.68.1（CPU 路径）
- 官方 tutorial 子集：`vendor/cell-eval2/docs/data/H1-VCC-2025-training.h5ad`
  600 细胞 × 1000 基因，原始计数，6 组各 100 细胞（non-targeting + MED12/SRC/STAT1/TET1/TMSB4X）

### ✅ 硬门槛：手算 expr_mae 与官方吻合到浮点精度
最大绝对偏差 **2.8e-17**，相对偏差 4.8e-16。脚本：`scripts/00_hard_gate.py`

**正确的伪bulk 约定**（读 `metrics/delta.py` 与 `metrics/discrimination.py` 的文档得到）：
```
pseudobulk_g = log1p(bulk_target_sum * group_gene_count_sum_g / group_total_count)
bulk_target_sum = 50000
```
即：组内按基因把**原始计数求和**，除以组内**总计数**，乘 5e4，再 log1p。

**走过的两条弯路**（记下来免得再踩）：
1. 逐细胞 CPM 归一化再对 log1p 取均值 → 偏高 4–9%
2. 逐细胞归一化再求均值再 log1p → 偏差降到 ~2%，仍不对
真正的做法是**计数加权**的组级伪bulk，不是逐细胞归一化后的平均。

**另外两点**：
- `filter_gene_min_cpm_cell=5.0` 只作用于 DE（`distances.py:107` 明确说它 never shrinks the panel），expr_mae 不做基因过滤
- `resolve_target_sum` 对对照子集返回 22376（对照中位文库大小），但 v2 counts/counts 运行下
  每个 delta 类指标用的是 `bulk_lognorm`，常数是 `bulk_target_sum=5e4`，不是这个值

### B0 基线的首次结果（5MB 子集，profile=vcc）
B0_boot（每个被扰动细胞 = 随机抽一个对照细胞，即「什么都没变」但保留真实方差）：

| 指标 | 结果 |
|---|---|
| `de_wilcoxon_overlap` | 全部 0.000（符合预期：预测不出任何差异表达基因）|
| `expr_mae` | 0.0388 – 0.0580 |
| `pds_l1` | 0.00 / 1.00 / 0.50 / 0.50 / 0.25（5 个扰动，方差极大）|

`pds_l1` 在「什么都没预测」的情况下居然出现 1.00 和 0.50 —— 子集只有 5 个扰动，判别
基本靠运气。这条不是发现，是样本量太小的提醒：**PDS 必须在全量 300 个靶点上看。**

---

## 读源码的重要发现（改变项目的问题表述）

### 1. 2026 的指标已经针对尺度/深度做了加固
`PROFILES['vcc2026']` = `pds_cosine`, `expr_mse_unbiased`, `expr_mse_unbiased_capped`,
`expr_distance_unbiased`, `expr_real_mass_ratio`, `expr_mse_unbiased_capped_norm`,
`de_wilcoxon_lfc_nmae`, `de_wilcoxon_direction_fidelity_yield_raw`,
`de_wilcoxon_direction_reach_raw`, `de_wilcoxon_sig_jaccard`

查 `MetricSpec.scoring.scored`：
- **计分**：`expr_mse_unbiased_capped_norm`（= `expr_mse_unbiased_capped` / `expr_distance_unbiased`，
  `agg='ratio_of_sums'`，clamp 到 [0,1]）、`pds_cosine`、`de_wilcoxon_sig_jaccard` 等
- **不计分（诊断用）**：`expr_mse_unbiased`、`expr_real_mass_ratio`

含义：MSE 被**除以真实数据自身的距离尺度**再截断到 [0,1] —— 这是针对「预测幅度刷分」的
尺度不变化修正。`expr_mse_unbiased` 还带 `needs_moments=True`，用矩校正减掉有限样本噪声。

**这正是我们「经验零分布」的解析版本。** 项目的问题因此应当重新表述为：

> 官方用**解析/矩校正**来扣除噪声底；我们用**经验零分布**（N1 惰性靶点 + N2 对照伪扰动）
> 来扣除。两者给出的「无效应」定义一致吗？2026 的加固真的让计分指标对测序深度和
> 预测尺度不变了吗？

这比「他们有漏洞」既更准确也更有价值。

### 2. 官方自带噪声注入器，但不在计分路径上
`cell_eval2.noise.noise_block(X, *, kind, level, seed, shard_idx)`，两种 kind：
- `gaussian`：非零元乘 lognormal；源码注释自己说 `exp(N(0,level))` 均值 >1，
  会**抬高期望文库大小**（level=0.5 时约 13%），是有意为之的尺度扰动
- `downsample`：二项稀释，`keep = 1 - level`

源码注释明确写：`noise=` 不进入任何计分路径（无 CLI 开关，默认 None，调用方只有
`tools/scale` 和 `tools/rowstore`）。所以这**不是**可利用的提交漏洞，而是现成的实验工具 ——
我们的 X2（稀释）、X4（加噪）可以直接用官方实现，可信度更高。

### 3. 一个已记录但故意不修的数值缺陷（诚实记录，不要当成评分漏洞）
`downsample` 分支的整数检查用了 numpy 默认**相对**容差的 `np.allclose`，源码注释实测
`50000.5` 能通过检查、随后被 `astype(np.int64)` 静默截断为 `50000` —— 正是那条报错信息
本该防止的事。注释说明已于 2026-08-19 裁定保留，理由就是它不进计分路径，并要求
不要重新作为新 bug 提交。

**我们不把它写成「发现的漏洞」。** 它是内部工具的已知边界，且官方已公开记录。

### 4. 提交端另有尺度闸门
`norm.check_scale_limit` —— 「Reject submissions whose per-cell totals exceed
max_counts_per_cell」，默认 `max_counts_per_cell=1e6`。所以单纯放大计数会被提交端拒掉，
X3（幅度缩放）必须在这个闸门以内做。

### 5. 整条竞赛标度可以本地复现
CLI 提供 `--preset vcc2026`、`--anchor`（重复锚点 = 分数 1 那一端，默认 5 次拆分平均）、
`--ceiling`（自拆分 + Spearman-Brown 数据天花板）、`baseline` 子命令（通用响应基线 =
分数 0 那一端）。**不需要提交到排行榜就能得到官方标度下的分数。**

### 6. 目标基因泄漏：官方自己在警告
跑子集时 cell-eval2 直接打印警告：5 个标签里只有 1 个解析到了实测基因，另外 4 个
**未做排除**地参与计分，各自保留了自己的转录本 —— 一个仅凭标签就能预测的坐标，可能
抬高分数。官方实测数字：未解析的靶点对一个「如实再现敲低」的模型读出 `pds_cosine` 1.0000，
而已解析的读 0.7885。默认 `DiscriminationParams(exclude_target_gene=True, exclusion_scope='panel')`。

**这是一个可直接做成实验的条目**：靶基因自身那一列对 PDS 的贡献有多大。官方给了
1.0000 vs 0.7885 这组对照，但只在他们的 panel 上测过一次。

---

## 阻塞项：全量数据下载
`storage.googleapis.com`、`zenodo.org`、`huggingface.co` 在**云端容器和 Mac mini VM 两边**
都被出站允许清单拒绝（只有 PyPI 和 GitHub 可达）。`gsutil`/`gcloud` 虽然装了，同样被挡。
→ 需要用户自行下载，见 README 的说明。

---

## 阶段 0 结果（全量训练集：221,273 细胞 × 18,080 基因，150 扰动 + 38,176 对照，48 批次）

工程约束：device_bash 的 VM 只有 **3 GB 内存、4 核**（不是 Mac mini 的 16 GB）。
后台进程会在回复结束后被回收 → 所有长任务都写成可分段续跑。

### 事实 1 —— 预注册预测作废
150 个靶点全部在面板内且全部表达：最低 CPM 2.84（SSBP1），中位 105，<CPM 5 仅 3 个。
惰性靶点 = 0 → **N1 无法构建**。另：Spearman(‖δ‖, 靶基因表达量) = +0.231，
「低表达 → 无效应」很弱；KAT2A 对照 CPM 3.3，却 ‖δ‖=12.5（2120 细胞）。

### 事实 3 —— 预注册预测错误（猜 40–60%）
| 量 | 原始 | 去噪后 |
|---|---|---|
| 平均方向解释的平方和/信号能量 | 12.4% | 13.3% |
| 最佳 rank-1 | 13.7% | — |
| 前 5 主方向 | 47.3% | — |
| 特异性中位 | 0.978 | 0.974 |

原因一：训练集**完全不含**核糖体/剪接体/蛋白酶体/Pol II/Integrator 靶点（我预注册的
「强但通用」类别在本 split 中不存在；Integrator 预测随之作废）。
原因二：细胞数混杂（见下）。

### 事实 2 —— N2 零分布
- 池 6000 对照细胞，参考 bulk = 其余 32,176 个对照，严格不相交
- 14 个对数网格点 × B=40；拟合 **null² = a/n + b**，a=2375.24，b=0.4495
- **留出验证**：不用 n=100、1000 拟合 → 预测误差 +0.29%、−0.41%
- **池偏差**：理论 b = a/6000 + a/32176 = 0.4697，实测 ‖池bulk−参考‖² = 0.4635，
  拟合 0.4495，三者吻合。真实扰动的噪声无池项 → 校正 null_true² = a/n + b − a/6000。
  不校正时 n=4760 处零分布**虚高 31%**
- **去混杂成功**：Spearman(‖δ‖, n) = −0.545 → Spearman(m_p, n) = **+0.074 (p=0.37)**

### 两个被检验的「底」假说
m_p 最小 1.149，无一接近 1，说明零分布漏了某个普遍方差源。

**批次 —— 否定。** 每个扰动覆盖中位 48/48 批次，单批次占比 ≤9%；批次偏移中位 0.096，
仅为零分布的 6%；批次匹配前后 m_p 几乎一字不差。实验设计非常平衡。

**靶基因自身的敲低 —— 证实。**
- 150/150 靶点被压低，**0 个上调**；剩余表达中位 **5.9%**（0.7%–89.3%），log 降中位 −1.50
- 靶基因一列占 ‖δ‖² 中位 10.1%；>25% 有 44 个，>50% 有 10 个，最大 89.2%
- 去掉后：最小 m 1.149 → **1.041**；m<1.1 有 4 个，<1.2 有 9 个；>2 从 111 降到 94
- Spearman(剩余表达, m_excl) = +0.13 (n.s.)：**压得多狠不预测效应多大**

### 值得跟进的条目（尚未下结论）
1. **「敲低成功但转录组不响应」**：LAD1 (m_excl 1.041)、RAB3B 1.042、ANTXR1 1.061、
   DZIP3 1.093 —— 均 ~95% 敲低。MKI67 也在其中（1.28），与「MKI67 对增殖非必需」的文献一致，
   可作合理性检查。注意：无转录组响应 ≠ 冗余，也可能只是该基因在 hESC 中无功能需求
2. **「敲低很弱但效应很强」**：HIRA 剩余 89.3% 却 m_excl 4.4；ARID1A 剩余 50% 却 4.9；
   STAG2 剩余 18% 却 7.5。候选解释：剂量敏感、或 CRISPRi 波及邻近/双向启动子（脱靶）。
   需查基因组邻居后再说
3. **新对抗预测 X7（审计用）**：「对照 + 只把靶基因压低 1.5 个 log 单位」—— 不含任何
   生物学，却能捕获中位 10% 的 δ 能量。需逐个核对官方各指标是否排除靶基因
   （PDS 默认 `exclude_target_gene=True, exclusion_scope='panel'`，其余指标待查）

### 信号最强的扰动（去混杂后 m_p）
SMARCA4 14.1、TADA1 12.5、DHX36 12.1、KAT2A 11.5、SMARCA5 10.9、DNMT1 10.6、
METTL3 10.4、METTL14 10.2、PRDM14 9.8、KDM1A 9.0、SOX2 8.5 —— 染色质/表观与多能性因子。
注意 TADA1、KAT2A 原始 ‖δ‖ 不在前列，靠大细胞数才显出来；PRDM14、KDM1A 原始排名靠前
主要因细胞少。**排序在校正前后显著改变。**
