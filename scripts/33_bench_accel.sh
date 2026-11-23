#!/bin/zsh
# 用法: 33_bench_accel.sh <mps|cpu> <steps>  —— 短跑测吞吐，输出到 runs/bench_<accel>
cd "${0:A:h}/.."
A=$1; N=${2:-80}; rm -rf runs/bench_$A
START=$(date +%s)
STATE_ACCEL=$A PYTORCH_ENABLE_MPS_FALLBACK=1 .venv-state/bin/python scripts/32_state_accel.py tx train \
  data.kwargs.toml_config_path=$PWD/data/st2/vcc.toml data.kwargs.embed_key=X_hvg \
  data.kwargs.pert_col=target_gene data.kwargs.cell_type_key=cell_type data.kwargs.batch_col=batch \
  data.kwargs.control_pert=non-targeting data.kwargs.output_space=all \
  training.max_steps=$N training.batch_size=8 model=state \
  output_dir=$PWD/runs/bench_$A name=vcc_h1 > logs/bench_$A.log 2>&1
echo "BENCH $A steps=$N exit=$? seconds=$(( $(date +%s) - START ))" >> logs/bench_$A.log
