#!/bin/zsh
# ST 预处理 + 在 MPS 上训练。可随时停，last.ckpt 保留进度。
cd "${0:A:h}/.."
S=.venv-state/bin/state; L=logs/state2.log; mkdir -p data/st2/h1_hvg
[ -f data/st2/h1_hvg/H1.h5ad ] || $S tx preprocess_train --adata data/st2/h1/H1.h5ad \
  --output data/st2/h1_hvg/H1.h5ad --num_hvgs 2000 >> $L 2>&1
echo "PREPROCESS_EXIT=$?" >> $L
PYTORCH_ENABLE_MPS_FALLBACK=1 $S tx train \
  data.kwargs.toml_config_path=$PWD/data/st2/vcc.toml \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=batch \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.output_space=all \
  training.max_steps=${1:-40000} \
  training.batch_size=8 \
  model=state \
  output_dir=$PWD/runs/state2 \
  name=vcc_h1 >> $L 2>&1
echo "TRAIN_EXIT=$?" >> $L
