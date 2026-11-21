"""阶段3：构建能在 3 GB 内存里评测的参考集。

保留全部扰动与全部基因；每个扰动抽 K 个细胞（不足则全取），对照抽 C 个。
用法: python 09_build_reference.py K C [seed] [NPERT]  （NPERT：只取细胞数≥K 的随机 NPERT 个扰动）
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np, h5py, pandas as pd, anndata as ad
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
H5 = ROOT/"data/raw/adata_Training.h5ad"
K, C = int(sys.argv[1]), int(sys.argv[2]); SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 0
NPERT = int(sys.argv[4]) if len(sys.argv) > 4 else 0
OUT = ROOT/(f"results/s3/real_K{K}_C{C}" + (f"_P{NPERT}" if NPERT else "") + (f"_s{SEED}" if SEED else "") + ".h5ad"); OUT.parent.mkdir(parents=True, exist_ok=True)
t0 = time.time(); rng = np.random.default_rng(SEED)

with h5py.File(H5, "r") as f:
    dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    tcat = dec(f["obs/target_gene/categories"][:]); tcode = f["obs/target_gene/codes"][:]
    bcat = dec(f["obs/batch/categories"][:]); bcode = f["obs/batch/codes"][:]
    gcat = dec(f["obs/guide_id/categories"][:]); gcode = f["obs/guide_id/codes"][:]
    cells = dec(f["obs/_index"][:]); genes = dec(f["var/_index"][:])
    gene_id = dec(f["var/gene_id"][:])
    ci = tcat.index("non-targeting")
    counts = np.bincount(tcode, minlength=len(tcat))
    keep_g = set(range(len(tcat)))
    if NPERT:
        elig = [g for g in range(len(tcat)) if g != ci and counts[g] >= K]
        keep_g = set(rng.choice(elig, size=NPERT, replace=False).tolist()) | {ci}
    sel = []
    for g in sorted(keep_g):
        rows = np.flatnonzero(tcode == g); k = C if g == ci else K
        sel.append(rows if rows.size <= k else rng.choice(rows, size=k, replace=False))
    sel = np.sort(np.concatenate(sel))
    ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]
    blocks = []
    # 按连续块读，块内用 scipy 做行选择
    STEP = 2000
    for s in range(0, len(tcode), STEP):
        e = min(s+STEP, len(tcode)); m = sel[(sel >= s) & (sel < e)] - s
        if not m.size: continue
        lo, hi = int(ip[s]), int(ip[e])
        blk = sparse.csr_matrix((Xd[lo:hi], Xi[lo:hi], (ip[s:e+1]-lo).astype(np.int64)),
                                shape=(e-s, len(genes)))
        blocks.append(blk[m]); del blk
X = sparse.vstack(blocks, format="csr").astype(np.float32)
obs = pd.DataFrame({"target_gene": pd.Categorical([tcat[i] for i in tcode[sel]]),
                    "batch": pd.Categorical([bcat[i] for i in bcode[sel]]),
                    "guide_id": pd.Categorical([gcat[i] for i in gcode[sel]])},
                   index=[cells[i] for i in sel])
var = pd.DataFrame({"gene_id": gene_id}, index=genes)
a = ad.AnnData(X=X, obs=obs, var=var)
a.write_h5ad(OUT, compression=None)
vc = obs.target_gene.value_counts()
print(f"{a.shape} nnz={X.nnz:,} ≈{X.nnz*8/1e6:.0f}MB 写入 {OUT.name}  {time.time()-t0:.1f}s")
print(f"对照 {vc['non-targeting']}；扰动细胞数 最少 {vc.drop('non-targeting').min()} 中位 {vc.drop('non-targeting').median():.0f}")
