"""阶段3：生成不含（或只含指定）生物学信息的对抗预测，写成与参考集同结构的 h5ad。

所有预测：
- 与参考集同 obs 标签、同每组细胞数；对照组直接复制参考集的对照
- 扰动组的「底细胞」从训练集中**不在参考集里**的对照细胞有放回抽取
- 输出为整数计数（乘法缩放后随机舍入，保持期望）

用法: python 11_adversarial.py <参考集.h5ad> <输出目录> [变体...]
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np, h5py, pandas as pd, anndata as ad
from scipy import sparse
from cell_eval2.noise import noise_block

ROOT = Path(__file__).resolve().parents[1]
H5 = ROOT/"data/raw/adata_Training.h5ad"
REF, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)
ALL = ["X1","X2_keep75","X2_keep50","X2_keep25","X4_gauss05","X3_a1","X3_a2","X7","X7plus","X5","X6"]
WANT = sys.argv[3:] or ALL
POOL, TS, KD_REMAIN, SEED = 4000, 50000.0, 0.059, 0
t0 = time.time(); rng = np.random.default_rng(SEED)
CTRL = "non-targeting"

ref = ad.read_h5ad(REF)
lab = ref.obs["target_gene"].astype(str).to_numpy()
genes = ref.var_names.to_list(); gidx = {g: i for i, g in enumerate(genes)}; G = len(genes)
is_c = lab == CTRL
perts = sorted(set(lab) - {CTRL})
Xr = ref.X.tocsr()
print(f"参考集 {ref.shape}，扰动 {len(perts)}", flush=True)

# ---- 预测用对照池：训练集中不在参考集里的对照
ref_ids = set(ref.obs_names)
with h5py.File(H5, "r") as f:
    dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]
    tcat = dec(f["obs/target_gene/categories"][:]); tcode = f["obs/target_gene/codes"][:]
    ids = dec(f["obs/_index"][:]); ci = tcat.index(CTRL)
    cand = np.array([i for i in np.flatnonzero(tcode == ci) if ids[i] not in ref_ids])
    pick = np.sort(rng.choice(cand, size=POOL, replace=False))
    ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]
    rows_d, rows_i, ptr = [], [], [0]
    for r in pick:
        a, b = int(ip[r]), int(ip[r+1])
        rows_d.append(Xd[a:b]); rows_i.append(Xi[a:b]); ptr.append(ptr[-1] + b - a)
Pool = sparse.csr_matrix((np.concatenate(rows_d), np.concatenate(rows_i), np.array(ptr)),
                         shape=(POOL, G)).astype(np.float32)
assert not (set(ids[i] for i in pick) & ref_ids)
print(f"预测对照池 {Pool.shape}（与参考集不相交）{time.time()-t0:.0f}s", flush=True)

# ---- 参考集伪bulk 与 δ（X3/X6 用）
def bulk(P): return np.log1p(TS * P / P.sum())
P_c = np.asarray(Xr[is_c].sum(0)).ravel().astype(np.float64); bc = bulk(P_c)
delta = {}
for p in perts:
    m = lab == p; delta[p] = bulk(np.asarray(Xr[m].sum(0)).ravel().astype(np.float64)) - bc
mu = np.zeros(G)                                   # 平均响应，不含各自靶基因（同官方 baseline）
cnt = np.zeros(G)
for p in perts:
    d = delta[p].copy(); keep = np.ones(G, bool); keep[gidx[p]] = False
    mu[keep] += d[keep]; cnt[keep] += 1
mu /= np.maximum(cnt, 1)

def stoch_round(M):
    M = M.tocsr(); fl = np.floor(M.data)
    M.data = (fl + (rng.random(M.data.size) < (M.data - fl))).astype(np.float32)
    M.eliminate_zeros(); return M

def scale_cols(M, factor):          # 每个基因乘一个因子
    return M @ sparse.diags(factor.astype(np.float32))

def base_cells(n):                  # n 个底细胞（有放回）
    return Pool[rng.integers(0, POOL, size=n)]

def knock_target(M, p):            # 只把靶基因按二项稀释压到 KD_REMAIN
    M = M.tocsc(copy=True); j = gidx[p]; s, e = M.indptr[j], M.indptr[j+1]
    M.data[s:e] = rng.binomial(M.data[s:e].astype(np.int64), KD_REMAIN).astype(np.float32)
    M = M.tocsr(); M.eliminate_zeros(); return M

def build(variant):
    blocks, labels = [Xr[is_c]], [CTRL]*int(is_c.sum())
    shuf = None
    if variant == "X5":
        shuf = dict(zip(perts, rng.permutation(perts)))   # 组 p 拿组 shuf[p] 的真实细胞
    for k, p in enumerate(perts):
        n = int((lab == p).sum())
        if variant == "X5":
            src = np.flatnonzero(lab == shuf[p]); M = Xr[rng.choice(src, size=n, replace=True)]
        else:
            M = base_cells(n)
            if variant.startswith("X2_keep"):
                M = noise_block(M, kind="downsample", level=1 - int(variant[-2:])/100, seed=SEED, shard_idx=k)
            elif variant == "X4_gauss05":
                M = stoch_round(noise_block(M, kind="gaussian", level=0.5, seed=SEED, shard_idx=k))
            elif variant in ("X3_a1", "X3_a2", "X7plus"):
                a = 2.0 if variant == "X3_a2" else 1.0
                M = stoch_round(scale_cols(M, np.exp(a*mu)))
            elif variant == "X6":
                M = stoch_round(scale_cols(M, np.exp(rng.permutation(delta[p]))))
            if variant in ("X7", "X7plus"):
                M = knock_target(M, p)
        blocks.append(M.tocsr().astype(np.float32)); labels += [p]*n
    X = sparse.vstack(blocks, format="csr")
    obs = pd.DataFrame({"target_gene": pd.Categorical(labels)}, index=[f"pred_{i}" for i in range(len(labels))])
    return ad.AnnData(X=X, obs=obs, var=ref.var.copy())

for v in WANT:
    t = time.time(); a = build(v)
    assert np.allclose(a.X.data, np.rint(a.X.data)), v
    a.write_h5ad(OUT/f"{v}.h5ad")
    print(f"{v:11s} {a.shape} nnz={a.X.nnz:,} 每细胞总计数中位={np.median(np.asarray(a.X.sum(1)).ravel()):,.0f}  {time.time()-t:.0f}s", flush=True)
print(f"全部完成 {time.time()-t0:.0f}s", flush=True)
