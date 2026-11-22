#!/bin/zsh
# 真实模型对照的参考集：validation split 的 50 个扰动（模型训练时未见），
# 建标度，并在同一参考上重做对抗预测，保证与模型分数可比。可重跑。
cd "${0:A:h}/.."
PY=.venv-mac/bin/python
R=results/s3/realVAL_K50_C1000.h5ad; O=results/s3/VAL; L=$O/val_progress
mkdir -p $O
[ -f $R ] || P5_SRC=$PWD/data/raw/validation/adata_Validation.h5ad P5_TAG=VAL \
  $PY scripts/09_build_reference.py 50 1000 0 >> $L 2>&1
echo "$(date +%T) 参考集就绪" >> $L
grep -q BUNDLE_DONE $O/progress 2>/dev/null || scripts/10_scale.sh $R $O 0
echo "$(date +%T) 标度 exit=$?" >> $L
[ -d $O/preds ] || P5_SRC=$PWD/data/raw/validation/adata_Validation.h5ad \
  $PY scripts/11_adversarial.py $R $O/preds X1 X2_keep25 X3_a1 X4_gauss05 X5 >> $L 2>&1
scripts/12_run_preds.sh $R $O
echo "$(date +%T) VAL_DONE" >> $L
