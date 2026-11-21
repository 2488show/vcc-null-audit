"""汇总各对抗预测的 vcc2026 官方标度分数（from_replicate = 竞赛分数）。"""
import sys
from pathlib import Path
import pandas as pd
O = Path(sys.argv[1]); rows = []
for s in sorted(O.glob("runs/*/score.csv")):
    d = pd.read_csv(s); d["variant"] = s.parent.name; rows.append(d[["variant","metric","from_baseline","from_replicate"]])
df = pd.concat(rows); df.to_csv(O/"scores_long.csv", index=False)
w = df.pivot(index="variant", columns="metric", values="from_replicate")
short = {"pds_cosine":"PDS","expr_mse_unbiased_capped_norm":"MSE","de_wilcoxon_lfc_nmae":"NMAE",
         "de_wilcoxon_direction_fidelity_yield_raw":"FIDELITY","de_wilcoxon_direction_reach_raw":"REACH",
         "de_wilcoxon_sig_jaccard":"JACCARD","avg_score":"AVG"}
w = w.rename(columns=short)[[c for c in short.values() if c in w.rename(columns=short).columns]]
w.to_csv(O/"scores_wide.csv")
pd.set_option("display.width", 200)
print(w.round(3).to_string())
