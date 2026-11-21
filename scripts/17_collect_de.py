"""(a) 机会校正版方向指标 vs 竞赛计分的 _raw 版，以及显著基因数（de 档位）。"""
import sys
from pathlib import Path
import pandas as pd
O = Path(sys.argv[1])
COLS = {"de_wilcoxon_nsig_counts_real":"n_sig_real","de_wilcoxon_nsig_counts_pred":"n_sig_pred",
        "de_wilcoxon_direction_fidelity_raw":"fidelity_raw","de_wilcoxon_direction_fidelity":"fidelity_corrected",
        "de_wilcoxon_direction_fidelity_yield_raw":"FY_raw(竞赛计分)","de_wilcoxon_direction_fidelity_yield":"FY_corrected",
        "de_wilcoxon_direction_coverage":"coverage","de_wilcoxon_direction_reach_raw":"reach_raw",
        "de_wilcoxon_direction_reach":"reach_corrected","de_wilcoxon_sig_jaccard":"jaccard",
        "de_wilcoxon_precision":"precision","de_wilcoxon_sig_recall":"recall"}
rows = {}
for f in sorted(O.glob("runs_de/*/agg_results.csv")):
    d = pd.read_csv(f)
    if f.parent.name == sorted(p.parent.name for p in O.glob("runs_de/*/agg_results.csv"))[0]:
        print("statistic 行:", d.statistic.tolist()); print("缺失列:", [c for c in COLS if c not in d.columns])
    for stat in ("mean", "median"):
        r = d[d.statistic == stat]
        rows[(f.parent.name, stat)] = {v: float(r[k].iloc[0]) for k, v in COLS.items() if k in d.columns}
df = pd.DataFrame(rows).T
df.to_csv(O/"de_profile_values.csv"); pd.set_option("display.width", 250)
print(df.xs("mean", level=1).round(3).to_string())
print("\n中位数（显著基因数）:")
print(df.xs("median", level=1)[["n_sig_real","n_sig_pred"]].round(0).to_string())
