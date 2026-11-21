"""机制检验：噪声预测为何能拿 FIDELITY（未做机会校正的方向一致性）。

单细胞层面：每细胞归一化到 1e4、log1p，各组取均值，与真实对照比。
- q：真实扰动中 |Δ| 最大的 TOP 个基因（不含靶基因）里下调所占比例
- 对每个预测：其 TOP 基因里下调的比例，以及在真实 TOP 基因上的方向一致率
"""
import sys
from pathlib import Path
import numpy as np, anndata as ad, pandas as pd
from scipy import sparse
O = Path(sys.argv[1]); REF = Path(sys.argv[2]); TOP = 50; CTRL = "non-targeting"

def group_meanlog(a):
    X = a.X.tocsr().astype(np.float64)
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
    L = sparse.diags(1e4/tot) @ X; L.data = np.log1p(L.data)
    lab = a.obs["target_gene"].astype(str).to_numpy()
    return {g: np.asarray(L[lab == g].mean(0)).ravel() for g in np.unique(lab)}

ref = ad.read_h5ad(REF); gidx = {g: i for i, g in enumerate(ref.var_names)}
R = group_meanlog(ref); c = R[CTRL]; perts = [p for p in R if p != CTRL]; del ref
def top_idx(d, p):
    d = d.copy(); d[gidx[p]] = 0; return np.argsort(-np.abs(d))[:TOP]
real_top = {p: top_idx(R[p]-c, p) for p in perts}
q = np.mean([np.mean((R[p]-c)[real_top[p]] < 0) for p in perts])
print(f"真实扰动 TOP{TOP} 基因中下调比例 q = {q:.3f}")
rows = []
for f in sorted((O/"preds").glob("*.h5ad")):
    P = group_meanlog(ad.read_h5ad(f))
    down = np.mean([np.mean((P[p]-c)[top_idx(P[p]-c, p)] < 0) for p in perts])
    agree = np.mean([np.mean(np.sign((P[p]-c)[real_top[p]]) == np.sign((R[p]-c)[real_top[p]])) for p in perts])
    mag = np.median([np.abs(P[p]-c).mean() for p in perts])
    rows.append(dict(variant=f.stem, pred_top_down=round(down,3), agree_on_real_top=round(agree,3), median_abs_shift=round(mag,4)))
df = pd.DataFrame(rows); df.to_csv(O/"sign_check.csv", index=False)
print(df.to_string(index=False))
