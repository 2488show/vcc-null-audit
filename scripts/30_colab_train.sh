#!/usr/bin/env bash
# 在 Colab T4 上把 State ST 训练到 40,000 步（阶段 4 只跑了 2,000 步）。
# 除 max_steps 与设备外，参数与 scripts/25_state_train2.sh 完全相同。
# num_workers 从默认 12 降到 2：Colab 标准机型只有 2 核 12 GB，默认值会被 OOM 杀掉（只影响数据加载，不影响训练）。
# 可重入：output_dir 里已有 last.ckpt 时 state 会从它续训。
set -uo pipefail
R=/content/p5; cd $R; mkdir -p logs sync runs
if [ ! -x .venv/bin/state ]; then
  uv venv --python 3.12 .venv
  uv pip install --python .venv/bin/python "arc-state==0.11.1" "cell-load==0.10.4" "anndata==0.13.4"
fi
export MPLBACKEND=Agg
cat > vcc.toml <<'T'
[datasets]
h1 = "/content/p5/data/h1_hvg/"

[training]
h1 = "train"

[zeroshot]

[fewshot]

[fewshot."h1.H1"]
val = ["ACAT2", "ACVR1B", "AKT2", "ANTXR1", "ARID1A"]
test = ["ACAT2", "ACVR1B", "AKT2", "ANTXR1", "ARID1A"]
T
.venv/bin/python -c "import torch,importlib.metadata as m;print('torch',torch.__version__,torch.cuda.is_available(),'arc-state',m.version('arc-state'),'cell-load',m.version('cell-load'),'anndata',m.version('anndata'))"

# 每当 last.ckpt 更新，就做一份稳定副本并写出标记，供本机拉取
( while true; do
    C=runs/state40k/vcc_h1/checkpoints/last.ckpt
    if [ -f "$C" ] && [ "$C" -nt sync/last.ckpt ] 2>/dev/null || { [ -f "$C" ] && [ ! -f sync/last.ckpt ]; }; then
      sleep 20; cp "$C" sync/last.tmp && mv sync/last.tmp sync/last.ckpt
      echo "$(date -u +%FT%TZ) $(stat -c %s sync/last.ckpt) $(sha256sum sync/last.ckpt | cut -c1-16)" > sync/last.marker
    fi
    sleep 60
  done ) &

.venv/bin/state tx train \
  data.kwargs.toml_config_path=$R/vcc.toml \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.batch_col=batch \
  data.kwargs.control_pert=non-targeting \
  data.kwargs.output_space=all \
  data.kwargs.num_workers=${NW:-2} \
  training.max_steps=${1:-40000} \
  training.batch_size=8 \
  model=state \
  output_dir=$R/runs/state40k \
  name=vcc_h1
echo "TRAIN_EXIT=$? $(date -u +%FT%TZ)"
sleep 150   # 让最后一个 checkpoint 完成同步副本
echo "ALL_DONE"
