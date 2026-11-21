"""阶段0 事实1（惰性靶点计数）与事实3（共同响应方向解释的方差）。"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
import sys
S = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"results/stage0"
z = np.load(S/"group_sums.npz", allow_pickle=True)
gene_sums, cats, genes = z["gene_sums"], list(z["cats"]), list(z["genes"])
group_cells = z["group_cells"]
CTRL = "non-targeting"
ci = cats.index(CTRL)
TS = 50000.0

gidx = {g: i for i, g in enumerate(genes)}
perts = [c for c in cats if c != CTRL]
print(f"组数 {len(cats)}（扰动 {len(perts)} + 对照）；对照细胞 {group_cells[ci]:,}")
print(f"扰动细胞数 中位 {np.median([group_cells[cats.index(p)] for p in perts]):.0f} "
      f"范围 {min(group_cells[cats.index(p)] for p in perts)}–{max(group_cells[cats.index(p)] for p in perts)}")

# ---------- 事实 1：靶基因在对照中的表达量 ----------
Pc = gene_sums[ci]; cpm_ctrl = 1e6*Pc/Pc.sum()
rows = []
for p in perts:
    i = gidx.get(p)
    rows.append({"target": p, "in_panel": i is not None,
                 "ctrl_cpm": float(cpm_ctrl[i]) if i is not None else np.nan,
                 "n_cells": int(group_cells[cats.index(p)])})
f1 = pd.DataFrame(rows).sort_values("ctrl_cpm", na_position="first")
f1.to_csv(S/"fact1_target_expression.csv", index=False)

n_miss = int((~f1.in_panel).sum())
print(f"\n=== 事实1：靶基因在 H1 对照中的表达 ===")
print(f"不在基因面板里（无法评估，也无法被 CRISPRi 的读出捕捉）: {n_miss}")
ex = f1[f1.in_panel]
for thr, name in [(0.5,"CPM<0.5 近乎不表达"),(1,"CPM<1"),(5,"CPM<5 低表达"),(10,"CPM<10")]:
    print(f"  {name:22s}: {int((ex.ctrl_cpm<thr).sum()):3d} / {len(ex)}")
print(f"  表达量中位数 CPM = {ex.ctrl_cpm.median():.1f}")
print("\n最低表达的 12 个靶点:")
print(ex.head(12).to_string(index=False))

# ---------- 事实 3：共同响应方向 ----------
def bulk(P):
    t = P.sum()
    return np.log1p(TS*P/t) if t > 0 else np.zeros_like(P)

bc = bulk(Pc)
D = np.vstack([bulk(gene_sums[cats.index(p)]) - bc for p in perts])   # (150, 18080)
norms = np.linalg.norm(D, axis=1)
tot = float((D**2).sum())

mean_d = D.mean(0); v_mean = mean_d/np.linalg.norm(mean_d)
proj_mean = D @ v_mean
ve_mean = float((proj_mean**2).sum())/tot

# 最佳 rank-1（第一主方向，不中心化）
U, sv, Vt = np.linalg.svd(D, full_matrices=False)
ve_pc1 = float(sv[0]**2)/tot
v_pc1 = Vt[0]
cos_mean_pc1 = float(abs(v_mean @ v_pc1))

print(f"\n=== 事实3：共同响应方向 ===")
print(f"扰动数 {D.shape[0]}，基因数 {D.shape[1]}")
print(f"||δ|| 中位 {np.median(norms):.3f}  范围 {norms.min():.3f}–{norms.max():.3f}")
print(f"平均方向 v̂ 解释的平方和比例      = {ve_mean:6.1%}   ← 预注册猜测 40–60%")
print(f"最佳 rank-1（第一主方向）        = {ve_pc1:6.1%}   （上界）")
print(f"cos(平均方向, 第一主方向)        = {cos_mean_pc1:.4f}")
print(f"前5个主方向累计                  = {float((sv[:5]**2).sum())/tot:6.1%}")

# 每个扰动的 (m_p 用 ||δ|| 代替, s_p)
s_p = np.sqrt(np.maximum(norms**2 - proj_mean**2, 0))/norms
out = pd.DataFrame({"target": perts, "delta_norm": norms,
                    "proj_on_vmean": proj_mean, "specificity_s": s_p,
                    "n_cells": [group_cells[cats.index(p)] for p in perts]})
out = out.merge(f1[["target","ctrl_cpm","in_panel"]], on="target", how="left")
out.sort_values("delta_norm", ascending=False).to_csv(S/"fact3_magnitude_specificity.csv", index=False)
print(f"\n特异性 s_p: 中位 {np.median(s_p):.3f} 范围 {s_p.min():.3f}–{s_p.max():.3f}")
print("\n幅度最大的 10 个扰动:")
print(out.sort_values("delta_norm", ascending=False).head(10).to_string(index=False))
print("\n幅度最小的 10 个扰动:")
print(out.sort_values("delta_norm").head(10).to_string(index=False))
