#!/bin/zsh
# 守护循环：训练进程若被中断（休眠、应用重启、被杀），自动从 last.ckpt 续训，直到出现 final.ckpt。
# 每 5 分钟往 logs/state40k_rate.log 记一行：时间、已到步数、swap 用量。
cd "${0:A:h}/.."
( while true; do
    echo "$(date '+%F %T') step=$(tail -n 1 runs/state40k/vcc_h1/version_0/metrics.csv 2>/dev/null | cut -d, -f4) swap=$(sysctl -n vm.swapusage | awk '{print $6}')" >> logs/state40k_rate.log
    sleep 300
  done ) &
MON=$!
for i in {1..40}; do
  [ -f runs/state40k/vcc_h1/checkpoints/final.ckpt ] && break
  echo "=== attempt $i $(date '+%F %T') ===" >> logs/state40k_loop.log
  zsh scripts/34_state_train40k.sh
  sleep 20
done
kill $MON
echo "LOOP_END $(date '+%F %T') final=$(ls runs/state40k/vcc_h1/checkpoints/final.ckpt 2>/dev/null)" >> logs/state40k_loop.log
