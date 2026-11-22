"""给预处理后的 h5ad 补一个 var/gene_name 普通字符串数据集。

anndata 0.13 把 var/_index 写成 nullable-string-array 分组，而 cell-load 用 h5py
按数据集读取（var/gene_name → var/gene_name_index → var/gene_name/categories → var/_index），
四条路径全部失败。这里只新增数据集，不改 column-order，所以 anndata 读取不受影响。
"""
import sys
from pathlib import Path
import h5py, numpy as np

p = Path(sys.argv[1])
with h5py.File(p, "r+") as f:
    v = f["var"]
    if "gene_name" in v:
        print("已存在 var/gene_name，跳过"); raise SystemExit
    node = v["_index"]
    raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
    names = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in raw], dtype=object)
    v.create_dataset("gene_name", data=names, dtype=h5py.string_dtype(encoding="utf-8"))
    print(f"写入 var/gene_name：{len(names)} 条，前 3 个 {list(names[:3])}")
