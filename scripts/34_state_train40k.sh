#!/bin/zsh
# 阶段 5：把 State ST 训练到官方的 40,000 步（阶段 4 只有 2,000 步）。
# 在本机 MPS 上跑：免费 Colab T4 只有 2 核 CPU，数据加载成为瓶颈（实测约 0.12 步/秒，GPU 利用率 0%），
# 比本机还慢。除 max_steps、accelerator 和输出目录外，参数与 scripts/25_state_train2.sh 相同。
# num_workers 从默认 12 降到 4：12 个进程在 16 GB 机器上把系统压进 swap（实测 9.5 GB），反而更慢。
# 可随时停（screen -S p5train -X quit），再次运行会从 runs/state40k 的 last.ckpt 续训。
cd "${0:A:h}/.."
L=logs/state40k.log
echo "=== start $(date '+%F %T') ===" >> $L
STATE_ACCEL=mps PYTORCH_ENABLE_MPS_FALLBACK=1 .venv-state/bin/python scripts/32_state_accel.py tx train \
  data.kwargs.toml_config_path=$PWD/data/st2/vcc.toml \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=batch \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.output_space=all \
  data.kwargs.num_workers=${NW:-4} \
  training.max_steps=40000 \
  training.batch_size=8 \
  model=state \
  output_dir=$PWD/runs/state40k \
  name=vcc_h1 >> $L 2>&1
echo "TRAIN_EXIT=$? $(date '+%F %T')" >> $L
