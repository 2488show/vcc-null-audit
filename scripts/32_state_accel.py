"""state CLI 包装器：允许用环境变量 STATE_ACCEL 指定 Lightning 的 accelerator。

arc-state 0.11.1 的训练入口写死了「有 CUDA 用 gpu，否则 cpu」，Apple Silicon 的 MPS 用不上。
这里只改 Trainer 的 accelerator 参数，并沿用 28_state_cli.py 对本机自产 checkpoint 的放行
（续训时要加载 last.ckpt）。模型、数据、损失一概不动。
"""
import os, sys, numpy as np, torch
try:
    from numpy._core.multiarray import scalar
except ImportError:
    from numpy.core.multiarray import scalar
torch.serialization.add_safe_globals([scalar, np.dtype, np.ndarray])
_orig_load = torch.load
def _load(*a, **k):
    k["weights_only"] = False   # 仅用于本机自产 checkpoint
    return _orig_load(*a, **k)
torch.load = _load

ACCEL = os.environ.get("STATE_ACCEL")
if ACCEL:
    import lightning.pytorch as pl
    _init = pl.Trainer.__init__
    def _patched(self, *a, **k):
        if k.get("accelerator") == "cpu":
            k["accelerator"] = ACCEL
        print(f"[32_state_accel] Trainer accelerator = {k.get('accelerator')}", flush=True)
        return _init(self, *a, **k)
    pl.Trainer.__init__ = _patched

if __name__ == "__main__":   # macOS 用 spawn 起数据加载子进程，会重新导入本文件
    from state.__main__ import main
    sys.exit(main())
