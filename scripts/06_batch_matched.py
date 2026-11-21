"""批次匹配：零分布最小 m_p=1.15、无一接近 1，怀疑批次组成不匹配。

对每个扰动构造「同批次组成」的对照 bulk（按该扰动在各批次的细胞比例重加权各批次对照计数），
比较 δ_raw（对全体对照）与 δ_matched（对批次匹配对照）。
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np, h5py, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]; S = ROOT/"results/stage0"
H5 = ROOT/"data/raw/adata_Training.h5ad"; TS, ROWS = 50000.0, 4000
A_FIT = 2375.24   # 来自 05：null² ≈ a/n + a/n_ref
t0 = time.time()
z = np.load(S/"group_sums.npz", allow_pickle=True)
cats = list(z["cats"]); codes = z["codes"]; bcodes = z["bcodes"].astype(int)
gene_sums = z["gene_sums"]; bcats = list(z["bcats"])
ci = cats.index("non-targeting"); nb = len(bcats); G = gene_sums.shape[1]

out = S/"ctrl_batch_sums.npy"
if out.exists():
    cbs = np.load(out)
else:
    cbs = np.zeros((nb, G))
    with h5py.File(H5, "r") as f:
        ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]; N = len(codes)
        for s in range(0, N, ROWS):
            e = min(s+ROWS, N); sel = np.flatnonzero(codes[s:e] == ci)
            if not sel.size: continue
            lo, hi = int(ip[s]), int(ip[e]); d = Xd[lo:hi].astype(np.float64); idx = Xi[lo:hi]
            off = ip[s:e+1]-lo
            for r in sel:
                a, b = int(off[r]), int(off[r+1]); np.add.at(cbs[bcodes[s+r]], idx[a:b], d[a:b])
    np.save(out, cbs)
print(f"对照按批次累计完成 {time.time()-t0:.1f}s", flush=True)

ctrl_n_b = np.bincount(bcodes[codes == ci], minlength=nb)
def bulk(P): return np.log1p(TS*P/P.sum())
bc_all = bulk(gene_sums[ci])

perts = [c for c in cats if c != "non-targeting"]
rows = []
for p in perts:
    gi = cats.index(p); m = codes == gi
    cnt = np.bincount(bcodes[m], minlength=nb).astype(float); n = cnt.sum()
    w = cnt/n
    ok = ctrl_n_b > 0
    if (w[~ok] > 0).any(): print(f"  !! {p} 有批次没有对照细胞", flush=True)
    # 每批次对照的「每细胞平均计数」按扰动的批次比例加权 -> 批次匹配对照 bulk
    per_cell = np.zeros_like(cbs); per_cell[ok] = cbs[ok]/ctrl_n_b[ok, None]
    Pm = (w[:, None]*per_cell).sum(0)
    bm = bulk(Pm)
    n_ref_eff = 1.0/np.sum(w[ok]**2/ctrl_n_b[ok])
    bp = bulk(gene_sums[gi])
    d_raw = np.linalg.norm(bp-bc_all); d_mat = np.linalg.norm(bp-bm)
    off_b = np.linalg.norm(bm-bc_all)
    null = np.sqrt(A_FIT/n + A_FIT/n_ref_eff)
    rows.append(dict(target=p, n_cells=int(n), n_batches=int((cnt>0).sum()),
                     top_batch_frac=float(w.max()), n_ref_eff=float(n_ref_eff),
                     batch_offset=float(off_b), delta_raw=float(d_raw), delta_matched=float(d_mat),
                     null=float(null), m_raw=float(d_raw/null), m_matched=float(d_mat/null)))
d = pd.DataFrame(rows)
d.to_csv(S/"fact2_batch_matched.csv", index=False)

print("\n=== 批次结构 ===")
print(f"  每个扰动覆盖的批次数: 中位 {d.n_batches.median():.0f} 范围 {d.n_batches.min()}–{d.n_batches.max()}（共 {nb}）")
print(f"  最大单批次占比: 中位 {d.top_batch_frac.median():.2f} 范围 {d.top_batch_frac.min():.2f}–{d.top_batch_frac.max():.2f}")
print(f"  批次偏移 ||bm−bc||: 中位 {d.batch_offset.median():.3f} 范围 {d.batch_offset.min():.3f}–{d.batch_offset.max():.3f}")
print(f"  批次偏移/零分布: 中位 {(d.batch_offset/d.null).median():.2f}")

print("\n=== m_p：全体对照参考 vs 批次匹配参考 ===")
for col in ("m_raw","m_matched"):
    print(f"  {col:9s}: 中位 {d[col].median():.2f}  最小 {d[col].min():.3f}  "
          f">1.5 {int((d[col]>1.5).sum())}  >2 {int((d[col]>2).sum())}  >3 {int((d[col]>3).sum())}")
r,p = spearmanr(d.m_matched, d.n_cells); print(f"  Spearman(m_matched, n_cells) = {r:+.3f} (p={p:.2f})")
print("\n批次匹配后最弱的 12 个:")
print(d.sort_values("m_matched")[["target","n_cells","n_batches","batch_offset","delta_raw","delta_matched","null","m_raw","m_matched"]].head(12).round(3).to_string(index=False))
print("\n批次匹配后最强的 8 个:")
print(d.sort_values("m_matched",ascending=False)[["target","n_cells","delta_matched","null","m_matched"]].head(8).round(2).to_string(index=False))
