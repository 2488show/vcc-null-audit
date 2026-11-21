"""检查 ||δ|| 与细胞数的混杂，并核对预注册的基因类别是否在 150 个靶点里。"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

ROOT = Path(__file__).resolve().parents[1]
import sys
S = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT/"results/stage0"
d = pd.read_csv(S/"fact3_magnitude_specificity.csv")

print("=== ||δ|| 与细胞数的关系（噪声混杂检查）===")
r_s, p_s = spearmanr(d.delta_norm, d.n_cells)
r_p, _   = pearsonr(d.delta_norm, 1/np.sqrt(d.n_cells))
print(f"Spearman(||δ||, n_cells)      = {r_s:+.3f}  (p={p_s:.2e})")
print(f"Pearson (||δ||, 1/sqrt(n))    = {r_p:+.3f}   ← 若接近 +1，||δ|| 主要在测噪声")
r_s2,_ = spearmanr(d.specificity_s, d.n_cells)
print(f"Spearman(s_p, n_cells)        = {r_s2:+.3f}")
r_s3,_ = spearmanr(d.delta_norm, d.ctrl_cpm)
print(f"Spearman(||δ||, 靶基因表达量) = {r_s3:+.3f}")

# 按细胞数分层看 ||δ||
d["n_bin"] = pd.qcut(d.n_cells, 5, labels=[f"Q{i+1}" for i in range(5)])
print("\n按细胞数五分位看幅度与特异性:")
print(d.groupby("n_bin", observed=True).agg(
    n_cells_med=("n_cells","median"), delta_norm_med=("delta_norm","median"),
    s_med=("specificity_s","median"), k=("target","size")).to_string())

# ---- 预注册基因类别核对 ----
targets = set(d.target)
groups = {
 "核心多能性TF": ["POU5F1","SOX2","NANOG","PRDM14","FOXH1","VENTX"],
 "核糖体": [t for t in targets if t.startswith(("RPL","RPS","MRPL","MRPS"))],
 "剪接体": ["SF3B1","SNRPD1","SF3A1","PRPF8","SNRNP200","U2AF1"],
 "RNA Pol II": [t for t in targets if t.startswith("POLR")],
 "蛋白酶体": [t for t in targets if t.startswith(("PSMA","PSMB","PSMC","PSMD"))],
 "Integrator": [t for t in targets if t.startswith("INTS")],
 "染色质/表观": ["KDM1A","SMARCA4","DNMT1","METTL3","METTL14","KAT2A","TADA1","HIRA","SUPT4H1"],
 "谱系TF(预期无效应)": ["MYOD1","NEUROD1","FOXA2","FOXA1","GATA4","PAX6"],
}
print("\n=== 预注册基因类别是否在 150 个训练靶点里 ===")
for name, gl in groups.items():
    present = sorted(set(gl) & targets)
    print(f"\n{name}: 命中 {len(present)}/{len(set(gl))}")
    if present:
        sub = d[d.target.isin(present)].sort_values("delta_norm", ascending=False)
        print(sub[["target","delta_norm","specificity_s","n_cells","ctrl_cpm"]].to_string(index=False))
    else:
        print("  （无）")
