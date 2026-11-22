"""把 ST 的预测整理成与 VAL 参考集同结构的整数计数 h5ad，供 cell-eval2 打分。"""
from __future__ import annotations
from pathlib import Path
import numpy as np, anndata as ad, pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT/"results/s3/realVAL_K50_C1000.h5ad"
P = ROOT/"data/st2/infer/st_pred.h5ad"
OUT = ROOT/"results/s3/VAL/preds"; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)

ref = ad.read_h5ad(REF); pred = ad.read_h5ad(P)
X = pred.X.toarray() if sparse.issparse(pred.X) else np.asarray(pred.X)
print(f"预测 {pred.shape}；X 范围 {X.min():.3f}–{X.max():.1f}；是否整数 {np.allclose(X, np.rint(X))}")
print(f"每细胞总和 中位 {np.median(X.sum(1)):,.0f}")
lab = pred.obs["target_gene"].astype(str).to_numpy()
print("组数:", len(set(lab)), "| 每组细胞数 中位", int(np.median(pd.Series(lab).value_counts())))

# ST 的输出在 log1p 空间（训练输入即 log1p，每细胞归一化到 53418.5）：
# 先 expm1 还原为计数尺度，再随机舍入成整数，保持期望。
X = np.expm1(np.clip(X, 0, None))
fl = np.floor(X); X = fl + (rng.random(X.shape) < (X - fl))
X = X.astype(np.float32)
a = ad.AnnData(X=sparse.csr_matrix(X),
               obs=pd.DataFrame({"target_gene": pd.Categorical(lab)},
                                index=[f"st_{i}" for i in range(len(lab))]),
               var=pd.DataFrame(index=ref.var_names))
assert np.allclose(a.X.data, np.rint(a.X.data))
assert list(a.var_names) == list(ref.var_names)
missing = set(ref.obs["target_gene"].astype(str)) - set(lab)
print("参考集中缺少预测的扰动:", missing or "无")
a.write_h5ad(OUT/"STATE.h5ad")
print(f"写出 {a.shape} -> {OUT/'STATE.h5ad'}；每细胞总计数 中位 {np.median(np.asarray(a.X.sum(1))):,.0f}")
