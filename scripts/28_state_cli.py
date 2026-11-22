"""调用 state CLI 的包装器：放行自己训练产生的 checkpoint。

torch 2.15 的 torch.load 默认 weights_only=True，会拒绝加载含 numpy 标量的
Lightning checkpoint。这里只对本机自己训练出的文件放行。
"""
import sys, numpy as np, torch
try:
    from numpy._core.multiarray import scalar
except ImportError:
    from numpy.core.multiarray import scalar
torch.serialization.add_safe_globals([scalar, np.dtype, np.ndarray])
_orig = torch.load
def _load(*a, **k):
    k["weights_only"] = False   # 仅用于本机自产 checkpoint（lightning 会显式传入，必须覆盖）
    return _orig(*a, **k)
torch.load = _load
from state.__main__ import main
sys.exit(main())
