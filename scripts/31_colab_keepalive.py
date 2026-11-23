# 在 Colab kernel 里运行（colab exec -f ... --timeout 0），防止运行时因 kernel 空闲被回收。
import os, time
t0 = time.time(); L = "/content/p5/logs/train.log"
while time.time() - t0 < 14 * 3600:
    done = os.path.exists(L) and "ALL_DONE" in open(L, errors="replace").read()[-2000:]
    print(f"alive {int((time.time()-t0)/60)} min done={done}", flush=True)
    if done: break
    time.sleep(60)
