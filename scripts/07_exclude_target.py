"""假说：m_p 的底（最小 1.15）来自靶基因自身的敲低那一列。去掉它重算。"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]; S = ROOT/"results/stage0"
TS, A = 50000.0, 2375.24
z = np.load(S/"group_sums.npz", allow_pickle=True)
cats = list(z["cats"]); genes = list(z["genes"]); gs = z["gene_sums"]; gc = z["group_cells"]
ci = cats.index("non-targeting"); gidx = {g:i for i,g in enumerate(genes)}
bm = pd.read_csv(S/"fact2_batch_matched.csv").set_index("target")
def bulk(P): return np.log1p(TS*P/P.sum())
bc = bulk(gs[ci])

rows = []
for p in [c for c in cats if c != "non-targeting"]:
    d = bulk(gs[cats.index(p)]) - bc; t = gidx[p]
    full = float(np.linalg.norm(d)); dt = float(d[t])
    excl = float(np.sqrt(max(full**2 - dt**2, 0)))
    null = bm.loc[p, "null"]
    # 剩余表达比例：用计数尺度 (expm1) 比，避免 log1p 在低计数处失真
    rem = float(np.expm1(bulk(gs[cats.index(p)])[t]) / np.expm1(bc[t]))
    rows.append(dict(target=p, n_cells=int(gc[cats.index(p)]), ctrl_cpm=1e6*gs[ci][t]/gs[ci].sum(),
                     d_target=dt, remaining=rem, target_share=dt**2/full**2,
                     m_full=full/null, m_excl=excl/null))
d = pd.DataFrame(rows); d.to_csv(S/"fact2_exclude_target.csv", index=False)

print("=== 靶基因自身的敲低 ===")
print(f"  d_target (log): 中位 {d.d_target.median():+.2f}  范围 {d.d_target.min():+.2f} – {d.d_target.max():+.2f}")
print(f"  剩余表达比例:    中位 {d.remaining.median():.1%}  范围 {d.remaining.min():.1%} – {d.remaining.max():.1%}")
print(f"  靶基因上调(>0)的扰动数: {int((d.d_target>0).sum())}")
print(f"  靶基因占 ||δ||² 的比例: 中位 {d.target_share.median():.1%}  最大 {d.target_share.max():.1%}")
for t in (0.25, 0.5, 0.75):
    print(f"    占比 > {t:.0%}: {int((d.target_share>t).sum())} 个扰动")

print("\n=== 去掉靶基因后的 m_p ===")
for col in ("m_full", "m_excl"):
    print(f"  {col}: 中位 {d[col].median():.2f}  最小 {d[col].min():.3f}  "
          f"<1.1 {int((d[col]<1.1).sum())}  <1.2 {int((d[col]<1.2).sum())}  >2 {int((d[col]>2).sum())}  >3 {int((d[col]>3).sum())}")
r, pv = spearmanr(d.m_excl, d.n_cells); print(f"  Spearman(m_excl, n_cells) = {r:+.3f} (p={pv:.2f})")
r, pv = spearmanr(d.remaining, d.m_excl); print(f"  Spearman(剩余表达, m_excl)  = {r:+.3f} (p={pv:.1e})  ← 压得越狠效应越大？")
r, pv = spearmanr(d.ctrl_cpm, d.remaining); print(f"  Spearman(靶基因表达量, 剩余) = {r:+.3f} (p={pv:.1e})")

print("\n去掉靶基因后最弱的 15 个（全转录组几乎无响应，但敲低本身成功）:")
print(d.sort_values("m_excl")[["target","n_cells","ctrl_cpm","remaining","target_share","m_full","m_excl"]].head(15).round(3).to_string(index=False))
print("\n敲低最差的 8 个（剩余比例最高）:")
print(d.sort_values("remaining", ascending=False)[["target","n_cells","ctrl_cpm","remaining","m_excl"]].head(8).round(3).to_string(index=False))
