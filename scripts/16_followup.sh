#!/bin/zsh
# 后续验证链（可重跑）：
#  (a) 用 de 档位重跑若干对抗预测，取得机会校正版 fidelity_yield 与显著基因数
#  (b) K 敏感性：40 个扰动 × 200 细胞 + 2000 对照，重建标度并评测 6 个预测
cd "${0:A:h}/.."
CE=.venv-mac/bin/cell-eval2; PY=.venv-mac/bin/python
C=(--preset vcc2026 --pert-col target_gene --control non-targeting --input-type counts)
O=results/s3/K50; R=results/s3/real_K50_C1000.h5ad; LOG=$O/followup_progress
mkdir -p $O/runs_de
for v in X3_a1 X1 X2_keep25 X4_gauss05 X5 X6 X3_a2; do
  d=$O/runs_de/$v
  [ -f $d/agg_results.csv ] || $CE run -ap $O/preds/$v.h5ad -ar $R -o $d $C --profile de > $d.log 2>&1
  echo "$(date +%T) de $v exit=$?" >> $LOG
done
echo "$(date +%T) DE_DONE" >> $LOG
R2=results/s3/real_K200_C2000_P40.h5ad; O2=results/s3/K200P40
[ -f $R2 ] || $PY scripts/09_build_reference.py 200 2000 0 40 >> $LOG 2>&1
grep -q BUNDLE_DONE $O2/progress 2>/dev/null || scripts/10_scale.sh $R2 $O2
echo "$(date +%T) scale K200P40 exit=$?" >> $LOG
[ -d $O2/preds ] || $PY scripts/11_adversarial.py $R2 $O2/preds X1 X2_keep25 X3_a1 X4_gauss05 X5 X7plus >> $LOG 2>&1
scripts/12_run_preds.sh $R2 $O2
echo "$(date +%T) FOLLOWUP_DONE" >> $LOG
