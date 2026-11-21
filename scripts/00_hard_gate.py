"""硬门槛：在官方 5MB 子集上独立手算 expr_mae，与 cell-eval2 输出对账。

对不上就不往下走。参见项目规划阶段 0。
"""
from __future__ import annotations
import subprocess, sys, json
from pathlib import Path
import numpy as np
import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUB = ROOT / "vendor/cell-eval2/docs/data/H1-VCC-2025-training.h5ad"
OUT = ROOT / "results/gate"
OUT.mkdir(parents=True, exist_ok=True)
PERT, CTRL = "target_gene", "non-targeting"
rng = np.random.default_rng(0)

real = ad.read_h5ad(SUB)
X = real.X.toarray().astype(np.float32)
is_ctrl = (real.obs[PERT] == CTRL).to_numpy()
print(f"real {real.shape}, controls {is_ctrl.sum()}")

# --- B0_boot：每个被扰动的细胞 = 随机抽一个对照细胞（「什么都没变」，但保留真实方差）
ctrl_idx = np.flatnonzero(is_ctrl)
Xp = X.copy()
pert_rows = np.flatnonzero(~is_ctrl)
Xp[pert_rows] = X[rng.choice(ctrl_idx, size=pert_rows.size, replace=True)]
pred = ad.AnnData(X=Xp, obs=real.obs.copy(), var=real.var.copy())

real_p = ROOT / "results/gate/real.h5ad"
pred_p = ROOT / "results/gate/pred_b0boot.h5ad"
real.write_h5ad(real_p); pred.write_h5ad(pred_p)

# --- 跑官方评测
cmd = [str(ROOT/".venv/bin/cell-eval2"), "run", "-ap", str(pred_p), "-ar", str(real_p),
       "-o", str(OUT/"ce2_vcc"), "--profile", "vcc", "--input-type", "counts",
       "--pert-col", PERT, "--control", CTRL]
print("\n$", " ".join(cmd[1:]))
r = subprocess.run(cmd, capture_output=True, text=True)
print(r.stdout[-2500:]); print(r.stderr[-2500:], file=sys.stderr)
if r.returncode != 0:
    print(f"!! cell-eval2 退出码 {r.returncode}"); sys.exit(1)

# --- 手算 expr_mae：用官方归一化，自己做伪bulk 与 MAE
from cell_eval2 import norm
tsum = norm.resolve_target_sum(real[is_ctrl].copy(), input_type="counts", target_sum=None)
print(f"\n官方解析出的 target_sum = {tsum}")

def lognorm(A, target_sum):
    tot = A.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    return np.log1p(A / tot * target_sum)

Lr, Lp = lognorm(X, tsum), lognorm(Xp, tsum)
perts = [p for p in real.obs[PERT].unique() if p != CTRL]
rows = []
for p in perts:
    m = (real.obs[PERT] == p).to_numpy()
    rows.append({"perturbation": p,
                 "expr_mae_hand": float(np.abs(Lp[m].mean(0) - Lr[m].mean(0)).mean())})
hand = pd.DataFrame(rows).sort_values("perturbation")
hand.to_csv(OUT/"hand_expr_mae.csv", index=False)

# --- 对账
res = pd.read_csv(OUT/"ce2_vcc/results.csv")
print("\nresults.csv 列:", list(res.columns))
print(res.head(20).to_string())
print("\n手算:"); print(hand.to_string(index=False))
print(f"\n手算均值 = {hand['expr_mae_hand'].mean():.6f}")
