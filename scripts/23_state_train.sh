#!/bin/zsh
# 在本机 MPS 上训练 ST。可随时停，last.ckpt 保留当前进度。
cd "${0:A:h}/.."
STEPS=${1:-40000}
PYTORCH_ENABLE_MPS_FALLBACK=1 .venv-state/bin/state tx train \
  data.kwargs.toml_config_path=$PWD/data/st/vcc.toml \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.output_space=all \
  training.max_steps=$STEPS \
  training.batch_size=8 \
  model=state \
  output_dir=$PWD/runs/state \
  name=vcc_h1 > logs/state_train.log 2>&1
echo "TRAIN_EXIT=$?" >> logs/state_train.log
