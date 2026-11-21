"""把官方标度拆开看：每个指标的原始值——0 端（通用响应基线）、1 端（重复锚点）、各对抗预测。"""
import sys
from pathlib import Path
import pandas as pd, polars as pl
O = Path(sys.argv[1])
M = {"pds_cosine":"PDS","expr_mse_unbiased_capped_norm":"MSE","de_wilcoxon_lfc_nmae":"NMAE",
     "de_wilcoxon_direction_fidelity_yield_raw":"FIDELITY","de_wilcoxon_direction_reach_raw":"REACH",
     "de_wilcoxon_sig_jaccard":"JACCARD","expr_real_mass_ratio":"mass_ratio(诊断)"}
def mean_row(csv):
    d = pd.read_csv(csv); d = d[d.statistic == "mean"]
    return {v: float(d[k].iloc[0]) for k, v in M.items() if k in d}
rows = {"0端 baseline": mean_row(O/"bundle/baseline_agg.csv")}
anc = pl.read_parquet(O/"bundle/anchor_agg.parquet")
print("anchor 列:", anc.columns); print(anc.head(14))
for r in sorted(O.glob("runs/*/agg_results.csv")):
    rows[r.parent.name] = mean_row(r)
df = pd.DataFrame(rows).T
df.to_csv(O/"raw_values.csv"); pd.set_option("display.width", 220)
print(df.round(4).to_string())
