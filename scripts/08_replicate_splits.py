"""用训练集拟合的零分布（三个 split 共用同一批对照细胞）在 validation/test 上复现去混杂。"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[1]
A, B, POOL = 2375.24, 0.4495, 6000
def null(n): return np.sqrt(np.maximum(A/n + B - A/POOL, 1e-12))
rows = []
for sp, d in [("train","results/stage0"),("validation","results/stage0_validation"),("test","results/stage0_test")]:
    f = pd.read_csv(ROOT/d/"fact3_magnitude_specificity.csv")
    m = f.delta_norm/null(f.n_cells.astype(float))
    r_raw,_ = spearmanr(f.delta_norm, f.n_cells); r_m,p_m = spearmanr(m, f.n_cells)
    rows.append(dict(split=sp, n_perts=len(f), raw_vs_n=round(r_raw,3), mp_vs_n=round(r_m,3), p=round(p_m,3),
                     mp_median=round(float(m.median()),2), mp_min=round(float(m.min()),3),
                     gt2=int((m>2).sum()), gt3=int((m>3).sum())))
    f.assign(m_p=m).to_csv(ROOT/d/"mp.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
