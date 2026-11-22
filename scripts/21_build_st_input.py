"""为 ST 训练合并一个 h5ad：train 的 150 个扰动 + validation 的 50 个扰动 + 对照。

validation 的 50 个扰动在 cell-load 的 [fewshot] 里被留出作 test，模型训练时见不到它们。
用法: python 21_build_st_input.py [每扰动细胞数] [对照细胞数]
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np, h5py, pandas as pd, anndata as ad
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
K = int(sys.argv[1]) if len(sys.argv) > 1 else 150
C = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
OUT = ROOT/"data/st/h1"; OUT.mkdir(parents=True, exist_ok=True)
SRC = [("train", ROOT/"data/raw/adata_Training.h5ad"),
       ("val",   ROOT/"data/raw/validation/adata_Validation.h5ad")]
CTRL, STEP = "non-targeting", 2000
rng = np.random.default_rng(0); t0 = time.time()
dec = lambda a: [x.decode() if isinstance(x, bytes) else str(x) for x in a]

blocks, obs_rows, genes_ref, val_perts = [], [], None, []
for which, f5 in SRC:
    with h5py.File(f5, "r") as f:
        tcat = dec(f["obs/target_gene/categories"][:]); tcode = f["obs/target_gene/codes"][:]
        bcat = dec(f["obs/batch/categories"][:]);       bcode = f["obs/batch/codes"][:]
        gcat = dec(f["obs/guide_id/categories"][:]);    gcode = f["obs/guide_id/codes"][:]
        ids = dec(f["obs/_index"][:]); genes = dec(f["var/_index"][:])
        if genes_ref is None: genes_ref = genes
        assert genes == genes_ref, "两个文件的基因面板不一致"
        ci = tcat.index(CTRL)
        if which == "val": val_perts = [c for c in tcat if c != CTRL]
        sel = []
        for g in range(len(tcat)):
            if g == ci and which == "val": continue      # 对照只从 train 取一次
            rows = np.flatnonzero(tcode == g); k = C if g == ci else K
            sel.append(rows if rows.size <= k else rng.choice(rows, size=k, replace=False))
        sel = np.sort(np.concatenate(sel))
        ip = f["X/indptr"][:]; Xd, Xi = f["X/data"], f["X/indices"]
        for s in range(0, len(tcode), STEP):
            e = min(s+STEP, len(tcode)); m = sel[(sel >= s) & (sel < e)] - s
            if not m.size: continue
            lo, hi = int(ip[s]), int(ip[e])
            blk = sparse.csr_matrix((Xd[lo:hi], Xi[lo:hi], (ip[s:e+1]-lo).astype(np.int64)),
                                    shape=(e-s, len(genes)))
            blocks.append(blk[m].astype(np.float32)); del blk
        obs_rows += [(ids[i], tcat[tcode[i]], bcat[bcode[i]], gcat[gcode[i]], which) for i in sel]
        print(f"{which}: 取 {sel.size} 细胞  {time.time()-t0:.0f}s", flush=True)

X = sparse.vstack(blocks, format="csr"); del blocks
obs = pd.DataFrame(obs_rows, columns=["cell", "target_gene", "batch", "guide_id", "split_src"]).set_index("cell")
obs["cell_type"] = pd.Categorical(["H1"]*len(obs))
for c in ("target_gene", "batch", "guide_id", "split_src"): obs[c] = pd.Categorical(obs[c])
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes_ref))
a.write_h5ad(OUT/"H1.h5ad", compression="gzip")
print(f"\n{a.shape} nnz={X.nnz:,} -> {OUT/'H1.h5ad'}  {time.time()-t0:.0f}s")
print(obs.groupby('split_src', observed=True).size().to_string())
print(f"留出（fewshot test）的 validation 扰动: {len(val_perts)} 个")
(ROOT/"data/st/val_perts.txt").write_text("\n".join(sorted(val_perts)))

toml = ROOT/"data/st/vcc.toml"
quoted = ", ".join(f'"{p}"' for p in sorted(val_perts))
toml.write_text(f'''[datasets]
h1 = "{OUT}/"

[training]
h1 = "train"

[zeroshot]

[fewshot]

[fewshot."h1.H1"]
val = [{quoted}]
test = [{quoted}]
''')
print(f"TOML -> {toml}")
