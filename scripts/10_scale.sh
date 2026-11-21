#!/bin/zsh
# 构建官方 vcc2026 标度：0 端 baseline + 1 端重复锚点 bundle（macOS 原生环境）
set -e
cd "${0:A:h}/.."
R=$1; O=$2; mkdir -p $O; rm -f $O/progress
CE=.venv-mac/bin/cell-eval2
C=(--preset vcc2026 --pert-col target_gene --control non-targeting --input-type counts)
/usr/bin/time -l $CE baseline -ar $R -o $O/baseline --save-pred $O/baseline_pred.h5ad $C > $O/baseline.log 2>&1
echo BASELINE_DONE >> $O/progress
/usr/bin/time -l $CE prep-real-bundle --real $R --baseline $O/baseline_pred.h5ad -o $O/bundle --force $C > $O/bundle.log 2>&1
echo BUNDLE_DONE >> $O/progress
