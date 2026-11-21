"""事实2 / N2：对照伪扰动经验零分布（按细胞数配平），可分段续跑。

- 池 6000 个对照细胞（与参考 bulk 严格不相交），缓存到 results/stage0/pool.npz
- 零分布在对数间隔的网格 n 上计算，每个 n 的结果单独存盘，重跑自动跳过
- 用法: python 04_null_n2.py [预算秒数]
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np, h5py
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
S = ROOT/"results/stage0"; H5 = ROOT/"data/raw/adata_Training.h5ad"
CACHE = S/"null_grid"; CACHE.mkdir(exist_ok=True)
TS, B, SEED, POOL, ROWS = 50000.0, 40, 0, 6000, 2000
BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 150
t0 = time.time()

z = np.load(S/"group_sums.npz", allow_pickle=True)
cats = list(z["cats"]); codes = z["codes"]; gene_sums = z["gene_sums"]
ci = cats.index("non-targeting"); n_genes = gene_sums.shape[1]

pool_f = S/"pool.npz"
if pool_f.exists():
    Xc = sparse.load_npz(pool_f)
else:
    rng0 = np.random.default_rng(SEED)
    pool_rows = np.sort(rng0.choice(np.flatnonzero(codes == ci), size=POOL, replace=False))
    in_pool = np.zeros(codes.size, bool); in_pool[pool_rows] = True
    blocks = []
    with h5py.File(H5, "r") as f:
        ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]
        N = int(f["X"].attrs["shape"][0])
        for s in range(0, N, ROWS):
            e = min(s+ROWS, N); sel = np.flatnonzero(in_pool[s:e])
            if not sel.size: continue
            lo, hi = int(ip[s]), int(ip[e])
            blk = sparse.csr_matrix((Xd[lo:hi], Xi[lo:hi], (ip[s:e+1]-lo).astype(np.int64)),
                                    shape=(e-s, n_genes))
            blocks.append(blk[sel]); del blk
    Xc = sparse.vstack(blocks, format="csr"); del blocks
    sparse.save_npz(pool_f, Xc, compressed=False)
print(f"池 {Xc.shape} nnz {Xc.nnz:,}  {time.time()-t0:.1f}s", flush=True)

pool_sums = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
ref = gene_sums[ci] - pool_sums
def bulk(P):
    t = P.sum(); return np.log1p(TS*P/t)
bc = bulk(ref)
np.save(S/"ref_bulk.npy", bc)

grid = sorted(set(np.unique(np.round(np.geomspace(33, 4760, 12)).astype(int)).tolist()
                  + [100, 1000]))              # 100、1000 两个点之后用于验证插值
for n in grid:
    out = CACHE/f"n{n}.npy"
    if out.exists(): continue
    if time.time()-t0 > BUDGET:
        print("预算用完，下次续跑", flush=True); break
    rng = np.random.default_rng(SEED + n)
    norms = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, POOL, size=n)
        if n <= 600:
            u, c = np.unique(idx, return_counts=True)
            P = np.asarray(Xc[u].T @ c.astype(np.float64)).ravel()
        else:
            w = np.bincount(idx, minlength=POOL).astype(np.float64)
            P = Xc.T @ w
        norms[b] = np.linalg.norm(bulk(P) - bc)
    np.save(out, norms)
    print(f"  n={n:>4}  零中位 {np.median(norms):.3f}  p95 {np.quantile(norms,.95):.3f}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
done = sorted(int(p.stem[1:]) for p in CACHE.glob("n*.npy"))
print(f"已完成 {len(done)}/{len(grid)}: {done}", flush=True)
