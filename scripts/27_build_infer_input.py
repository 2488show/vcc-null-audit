"""构建 ST 推理输入：VAL 参考集的对照细胞，按训练时的归一化约定处理。

训练数据的约定（已核对）：每细胞 normalize_total 到 53418.5，再 log1p；X_hvg 是其 2000 个高变基因子集。
同时写出 TSV，指定要为哪 50 个 validation 扰动各生成多少细胞。
"""
from __future__ import annotations
from pathlib import Path
import numpy as np, anndata as ad, pandas as pd, h5py
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/"results/s3/realVAL_K50_C1000.h5ad"
TRAIN = ROOT/"data/st2/h1_hvg/H1.h5ad"
OUT = ROOT/"data/st2/infer"; OUT.mkdir(parents=True, exist_ok=True)
TS, CTRL, NCELL = 53418.5, "non-targeting", 50

ref = ad.read_h5ad(REF)
with h5py.File(TRAIN, "r") as f:
    hv = f["var/highly_variable"][:]
    tr_genes = [x.decode() if isinstance(x, bytes) else str(x) for x in f["var/gene_name"][:]]
assert list(ref.var_names) == tr_genes, "基因面板与训练数据不一致"
print(f"高变基因 {hv.sum()}")

m = (ref.obs["target_gene"].astype(str) == CTRL).to_numpy()
X = ref.X[m].tocsr().astype(np.float64)
tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
L = sparse.diags(TS/tot) @ X
L.data = np.log1p(L.data)
L = L.tocsr().astype(np.float32)

obs = pd.DataFrame({
    "target_gene": pd.Categorical([CTRL]*int(m.sum())),
    "batch": pd.Categorical(ref.obs["batch"][m].astype(str).to_numpy()),
    "cell_type": pd.Categorical(["H1"]*int(m.sum())),
}, index=ref.obs_names[m])
a = ad.AnnData(X=L, obs=obs, var=pd.DataFrame(index=ref.var_names))
a.obsm["X_hvg"] = np.asarray(L[:, hv].todense(), dtype=np.float32)
p = OUT/"infer_input.h5ad"; a.write_h5ad(p)

# cell-load 按 h5py 读 var/gene_name
with h5py.File(p, "r+") as f:
    if "gene_name" not in f["var"]:
        f["var"].create_dataset("gene_name", data=np.array(list(ref.var_names), dtype=object),
                                dtype=h5py.string_dtype(encoding="utf-8"))

val_perts = sorted((ROOT/"data/st2/val_perts.txt").read_text().split())
tsv = OUT/"perts.tsv"
pd.DataFrame({"perturbation": val_perts, "num_cells": [NCELL]*len(val_perts)}).to_csv(tsv, sep="\t", index=False)
print(f"{a.shape} 对照细胞，obsm X_hvg {a.obsm['X_hvg'].shape} -> {p.name}")
print(f"TSV: {len(val_perts)} 个扰动 × {NCELL} 细胞 -> {tsv.name}")
