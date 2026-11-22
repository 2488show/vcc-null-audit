"""ST 训练输入（v2）：全部 200 个扰动都参与训练，但**排除 VAL 参考集里的细胞**。

ST 用 one-hot 编码扰动，无法预测训练时未见的扰动（基因嵌入文件需从 HuggingFace 取，不可达）。
因此改为「留出细胞」设计：模型见过这些扰动，但没见过被当作真值的那批细胞。
留出 5 个 train split 的扰动作为 Lightning 的 val/test，不影响 50 个 validation 扰动。
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np, h5py, pandas as pd, anndata as ad
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/"results/s3/realVAL_K50_C1000.h5ad"
OUT = ROOT/"data/st2/h1"; OUT.mkdir(parents=True, exist_ok=True)
K, C, STEP, CTRL = 150, 3000, 2000, "non-targeting"
SRC = [("train", ROOT/"data/raw/adata_Training.h5ad"),
       ("val",   ROOT/"data/raw/validation/adata_Validation.h5ad")]
rng = np.random.default_rng(7); t0 = time.time()
dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]

ref_ids = set(ad.read_h5ad(REF, backed="r").obs_names)
print(f"VAL 参考集细胞 {len(ref_ids):,}（将从训练数据中排除）", flush=True)

blocks, obs_rows, genes_ref, val_perts, train_perts = [], [], None, [], []
for which, f5 in SRC:
    with h5py.File(f5, "r") as f:
        tcat = dec(f["obs/target_gene/categories"][:]); tcode = f["obs/target_gene/codes"][:]
        bcat = dec(f["obs/batch/categories"][:]);       bcode = f["obs/batch/codes"][:]
        ids = dec(f["obs/_index"][:]); genes = dec(f["var/_index"][:])
        if genes_ref is None: genes_ref = genes
        assert genes == genes_ref
        ci = tcat.index(CTRL)
        keep = np.array([i for i in range(len(ids)) if ids[i] not in ref_ids])
        kmask = np.zeros(len(ids), bool); kmask[keep] = True
        if which == "val": val_perts = [c for c in tcat if c != CTRL]
        else: train_perts = [c for c in tcat if c != CTRL]
        sel = []
        for g in range(len(tcat)):
            if g == ci and which == "val": continue
            rows = np.flatnonzero((tcode == g) & kmask); k = C if g == ci else K
            if rows.size == 0: continue
            sel.append(rows if rows.size <= k else rng.choice(rows, size=k, replace=False))
        sel = np.sort(np.concatenate(sel))
        assert not (set(ids[i] for i in sel) & ref_ids), "与参考集有重叠"
        ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]
        for s in range(0, len(tcode), STEP):
            e = min(s+STEP, len(tcode)); m = sel[(sel >= s) & (sel < e)] - s
            if not m.size: continue
            lo, hi = int(ip[s]), int(ip[e])
            blk = sparse.csr_matrix((Xd[lo:hi], Xi[lo:hi], (ip[s:e+1]-lo).astype(np.int64)),
                                    shape=(e-s, len(genes)))
            blocks.append(blk[m].astype(np.float32)); del blk
        obs_rows += [(ids[i], tcat[tcode[i]], bcat[bcode[i]], which) for i in sel]
        print(f"{which}: {sel.size} 细胞  {time.time()-t0:.0f}s", flush=True)

X = sparse.vstack(blocks, format="csr"); del blocks
obs = pd.DataFrame(obs_rows, columns=["cell","target_gene","batch","split_src"]).set_index("cell")
obs["cell_type"] = pd.Categorical(["H1"]*len(obs))
# cell-load 要求按 (cell_type, perturbation) 连续分组
obs = obs.sort_values("target_gene", kind="stable")
X = X[[list(pd.DataFrame(obs_rows, columns=["cell","t","b","s"]).set_index("cell").index).index(c) for c in obs.index]] \
    if False else X[np.argsort(pd.Categorical(pd.DataFrame(obs_rows)[1]).codes, kind="stable")]
for c in ("target_gene","batch","split_src"): obs[c] = pd.Categorical(obs[c])
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes_ref))
a.write_h5ad(OUT/"H1.h5ad", compression="gzip")
print(f"\n{a.shape} nnz={X.nnz:,}  {time.time()-t0:.0f}s")
print(obs.groupby('split_src', observed=True).size().to_string())

holdout = sorted(train_perts)[:5]
q = ", ".join(f'"{p}"' for p in holdout)
(ROOT/"data/st2/vcc.toml").write_text(f'''[datasets]
h1 = "{OUT}/"

[training]
h1 = "train"

[zeroshot]

[fewshot]

[fewshot."h1.H1"]
val = [{q}]
test = [{q}]
''')
(ROOT/"data/st2/val_perts.txt").write_text("\n".join(sorted(val_perts)))
print(f"Lightning 的 val/test 留出 5 个 train 扰动: {holdout}")
print(f"待预测的 50 个 validation 扰动已全部参与训练（但细胞与参考集不重叠）")
