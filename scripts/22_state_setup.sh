#!/bin/zsh
# 装 arc-state（独立 venv，避免和 cell-eval2 的依赖冲突）并做 HVG 预处理
cd "${0:A:h}/.."
export UV_CACHE_DIR="$PWD/.uv-cache-mac"
L=logs/state_setup.log
uv venv --python 3.12 .venv-state >> $L 2>&1
uv pip install --python .venv-state/bin/python arc-state >> $L 2>&1
echo "=== 安装完成 exit=$? ===" >> $L
.venv-state/bin/state --help >> $L 2>&1
echo "=== preprocess_train ===" >> $L
.venv-state/bin/state tx preprocess_train \
  --adata data/st/h1/H1.h5ad \
  --output data/st/h1_hvg/H1.h5ad \
  --num_hvgs 2000 >> $L 2>&1
echo "=== preprocess exit=$? ===" >> $L
