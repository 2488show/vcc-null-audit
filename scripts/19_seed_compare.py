"""跨种子汇总：每个变体在各指标上的三次独立复现值、均值、极差。

种子 0 = results/s3/K50（第一轮），种子 1/2 = results/s3/K50_s{1,2}。
每个种子独立重抽参考集、重建标度（baseline dispersed 发射与 anchor 拆分都换种子）、重生成预测。
"""
from pathlib import Path
import pandas as pd, numpy as np
ROOT = Path(__file__).resolve().parents[1]
SHORT = {"pds_cosine":"PDS","expr_mse_unbiased_capped_norm":"MSE","de_wilcoxon_lfc_nmae":"NMAE",
         "de_wilcoxon_direction_fidelity_yield_raw":"FIDELITY","de_wilcoxon_direction_reach_raw":"REACH",
         "de_wilcoxon_sig_jaccard":"JACCARD","avg_score":"AVG"}
dirs = {0: ROOT/"results/s3/K50", 1: ROOT/"results/s3/K50_s1", 2: ROOT/"results/s3/K50_s2"}
rows = []
for seed, d in dirs.items():
    for f in sorted(d.glob("runs/*/score.csv")):
        s = pd.read_csv(f)
        for _, r in s.iterrows():
            if r.metric in SHORT:
                rows.append(dict(seed=seed, variant=f.parent.name, metric=SHORT[r.metric],
                                 score=r.from_replicate))
if not rows: print("还没有结果"); raise SystemExit
d = pd.DataFrame(rows)
d.to_csv(ROOT/"results/s3/seed_scores_long.csv", index=False)
common = sorted(set.intersection(*[set(d[d.seed==s].variant) for s in d.seed.unique()]))
d = d[d.variant.isin(common)]
pd.set_option("display.width", 240)
print(f"种子: {sorted(d.seed.unique())}；变体: {common}\n")
for m in ["FIDELITY","AVG","JACCARD","MSE","NMAE","REACH","PDS"]:
    sub = d[d.metric==m]
    if sub.empty: continue
    w = sub.pivot(index="variant", columns="seed", values="score")
    w["均值"] = w.mean(1); w["极差"] = w.iloc[:, :len(dirs)].max(1) - w.iloc[:, :len(dirs)].min(1)
    print(f"=== {m} ===");  print(w.round(3).to_string()); print()
summ = d.groupby(["metric","variant"]).score.agg(["mean","std","min","max"]).round(3)
summ.to_csv(ROOT/"results/s3/seed_summary.csv")
print("按指标的跨变体极差（衡量指标本身的种子敏感性）:")
print(d.groupby("metric").score.agg(lambda x: round(x.max()-x.min(),3)).to_string())
