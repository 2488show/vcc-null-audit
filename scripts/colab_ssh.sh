#!/bin/zsh
# 通过一条复用的主连接访问 Colab 会话（colab ssh 同一时间只允许一个连接；
# 用 ControlMaster 复用后，命令、rsync 可以并行走同一条连接）。
# 用法: scripts/colab_ssh.sh <会话名> [ssh 参数/命令...]
S=$1; shift
exec ssh -o "ProxyCommand=/Users/vmtakx/.local/bin/colab ssh --proxy-mode -s $S -i $HOME/.ssh/id_ed25519.pub" \
  -o ControlMaster=auto -o "ControlPath=/tmp/cm-$S" -o ControlPersist=12h \
  -o ServerAliveInterval=30 -o ServerAliveCountMax=6 \
  -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o LogLevel=ERROR \
  root@colab-$S "$@"
