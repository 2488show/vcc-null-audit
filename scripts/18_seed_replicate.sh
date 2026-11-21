#!/bin/zsh
# 多种子复现：每个种子独立重抽参考集、重建官方标度、重生成预测并评测。
# 全部可重跑，已完成的步骤跳过。用法: 18_seed_replicate.sh <种子> [种子...]
cd "${0:A:h}/.."
PY=.venv-mac/bin/python
VARIANTS=(X1 X2_keep25 X3_a1 X4_gauss05 X5 X7plus)
for S in "$@"; do
  R=results/s3/real_K50_C1000_s$S.h5ad; O=results/s3/K50_s$S; L=$O/seed_progress
  mkdir -p $O
  [ -f $R ] || $PY scripts/09_build_reference.py 50 1000 $S >> $L 2>&1
  echo "$(date +%T) seed $S 参考集就绪" >> $L
  grep -q BUNDLE_DONE $O/progress 2>/dev/null || scripts/10_scale.sh $R $O $S
  echo "$(date +%T) seed $S 标度 exit=$?" >> $L
  [ -d $O/preds ] || P5_SEED=$S $PY scripts/11_adversarial.py $R $O/preds $VARIANTS >> $L 2>&1
  scripts/12_run_preds.sh $R $O
  echo "$(date +%T) seed $S ALL_DONE" >> $L
done
echo "$(date +%T) SEEDS_DONE" >> results/s3/seeds_progress
