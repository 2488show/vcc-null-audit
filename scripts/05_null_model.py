"""拟合零分布 null² = a/n + b，核对池大小偏差，留出验证插值，计算 m_p 与去噪后的事实3。"""
from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]; S = ROOT/"results/stage0"
POOL, NREF, TS = 6000, 38176-6000, 50000.0
grid = {int(p.stem[1:]): np.load(p) for p in (S/"null_grid").glob("n*.npy")}
ns = np.array(sorted(grid)); med = np.array([np.median(grid[n]) for n in ns])

def fit(nn, mm):
    A = np.c_[1/nn, np.ones_like(nn, dtype=float)]
    (a, b), *_ = np.linalg.lstsq(A, mm**2, rcond=None); return a, b

# --- 1. 留出验证：不用 n=100 和 1000 拟合，去预测它们
keep = ~np.isin(ns, [100, 1000])
a_h, b_h = fit(ns[keep].astype(float), med[keep])
print("=== 留出验证插值 ===")
for n in (100, 1000):
    pred = np.sqrt(a_h/n + b_h); obs = np.median(grid[n])
    print(f"  n={n:>4}: 预测 {pred:.4f}  实测 {obs:.4f}  相对误差 {pred/obs-1:+.2%}")

# --- 2. 全网格拟合
a, b = fit(ns.astype(float), med)
resid = med - np.sqrt(a/ns + b)
print(f"\n=== 拟合 null² = a/n + b ===\n  a = {a:.2f}   b = {b:.4f}   最大残差 {np.abs(resid).max():.4f}")

# --- 3. 池大小偏差核对：理论 b = a/POOL + a/NREF
b_theory = a/POOL + a/NREF
Xc = sparse.load_npz(S/"pool.npz")
bc = np.load(S/"ref_bulk.npy")
P = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
pool_vs_ref = float(np.linalg.norm(np.log1p(TS*P/P.sum()) - bc))
print(f"\n=== 池大小偏差核对 ===")
print(f"  拟合的 b                = {b:.4f}")
print(f"  理论 a/{POOL} + a/{NREF} = {b_theory:.4f}")
print(f"  实测 ||池bulk − 参考||² = {pool_vs_ref**2:.4f}")
# 校正：真实扰动的噪声只有 a/n + a/NREF，扣掉池项 a/POOL
def null_true(n): return np.sqrt(np.maximum(a/n + b - a/POOL, 1e-12))
def null_boot(n): return np.sqrt(a/n + b)
for n in (33, 1000, 4760):
    print(f"  n={n:>4}: 自助法 {null_boot(n):.3f} → 校正后 {null_true(n):.3f} "
          f"（虚高 {null_boot(n)/null_true(n)-1:+.1%}）")

# --- 4. m_p 与去噪的幅度/特异性
d = pd.read_csv(S/"fact3_magnitude_specificity.csv")
d["null"] = null_true(d.n_cells.astype(float))
d["m_p"] = d.delta_norm/d.null
d["signal_energy"] = np.maximum(d.delta_norm**2 - d.null**2, 0)
d["signal_norm"] = np.sqrt(d.signal_energy)
spec_e = np.maximum(d.signal_energy - d.proj_on_vmean**2, 0)
d["s_denoised"] = np.where(d.signal_energy > 0, np.sqrt(spec_e)/np.sqrt(d.signal_energy.clip(lower=1e-12)), np.nan)
d.sort_values("m_p", ascending=False).to_csv(S/"fact2_mp.csv", index=False)

print("\n=== 去混杂检查：校正后还跟细胞数相关吗 ===")
for col in ("delta_norm", "m_p", "signal_norm", "specificity_s", "s_denoised"):
    r, p = spearmanr(d[col], d.n_cells, nan_policy="omit")
    print(f"  Spearman({col:13s}, n_cells) = {r:+.3f}  (p={p:.1e})")

print("\n=== 信噪比 m_p 的分布 ===")
for t in (1.0, 1.1, 1.5, 2, 3, 5):
    print(f"  m_p > {t:<4}: {int((d.m_p>t).sum()):3d}/150")
print(f"  中位 m_p = {d.m_p.median():.2f}")

tot_sig = d.signal_energy.sum()
print("\n=== 去噪后的事实3 ===")
print(f"  共同方向占信号能量比例 = {float((d.proj_on_vmean**2).sum()/tot_sig):.1%}"
      f"   （未去噪是 12.4%）")
print(f"  去噪特异性 s 中位 = {d.s_denoised.median():.3f}  （未去噪 0.978）")
q = pd.qcut(d.n_cells, 5, labels=[f"Q{i+1}" for i in range(5)])
print("\n按细胞数五分位（校正后）:")
print(d.groupby(q, observed=True).agg(n=("n_cells","median"), raw=("delta_norm","median"),
      m_p=("m_p","median"), signal=("signal_norm","median"), s_dn=("s_denoised","median")).round(3).to_string())

print("\n=== m_p 前 15 ===")
print(d.sort_values("m_p", ascending=False)[["target","n_cells","delta_norm","null","m_p","signal_norm","s_denoised"]].head(15).round(3).to_string(index=False))
print("\n=== m_p 后 10 ===")
print(d.sort_values("m_p")[["target","n_cells","delta_norm","null","m_p"]].head(10).round(3).to_string(index=False))
