#!/bin/zsh
# 串行流水线：等 bundle 完成 → 逐个预测 run → score → 汇总。可重跑，已完成的跳过。
# 用法: 12_run_preds.sh <参考集.h5ad> <K目录>
cd "${0:A:h}/.."
R=$1; O=$2; CE=.venv-mac/bin/cell-eval2
C=(--preset vcc2026 --pert-col target_gene --control non-targeting --input-type counts)
until grep -q BUNDLE_DONE $O/progress 2>/dev/null; do sleep 20; done
echo "$(date +%T) bundle 就绪" >> $O/runs_progress
for f in $O/preds/*.h5ad; do
  v=${f:t:r}; d=$O/runs/$v; mkdir -p $d
  if [ ! -f $d/agg_results.csv ]; then
    /usr/bin/time -l $CE run -ap $f -ar $R -o $d $C > $d.log 2>&1
    echo "$(date +%T) run $v exit=$?" >> $O/runs_progress
  fi
  if [ -f $d/agg_results.csv ] && [ ! -f $d/score.csv ]; then
    $CE score --user-agg $d/agg_results.csv --real-bundle $O/bundle -o $d/score.csv > $d/score.log 2>&1
    echo "$(date +%T) score $v exit=$?" >> $O/runs_progress
  fi
done
.venv-mac/bin/python scripts/13_collect.py $O >> $O/runs_progress 2>&1
echo "$(date +%T) ALL_DONE" >> $O/runs_progress
