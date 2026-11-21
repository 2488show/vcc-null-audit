"""流式扫过全量 CSR，累计每组的按基因计数总和、组总计数、每细胞文库大小。

只需一次 I/O。输出 results/stage0/group_sums.npz
"""
from __future__ import annotations
import time
from pathlib import Path
import numpy as np, h5py

ROOT = Path(__file__).resolve().parents[1]
H5 = ROOT/"data/raw/adata_Training.h5ad"
OUT = ROOT/"results/stage0"; OUT.mkdir(parents=True, exist_ok=True)
ROWS = 4000   # 每块细胞数

t0 = time.time()
with h5py.File(H5, "r") as f:
    cats = [c.decode() if isinstance(c, bytes) else str(c)
            for c in f["obs/target_gene/categories"][:]]
    codes = f["obs/target_gene/codes"][:].astype(np.int32)
    genes = [g.decode() if isinstance(g, bytes) else str(g)
             for g in f["var/_index"][:]]
    bcats = [c.decode() if isinstance(c, bytes) else str(c)
             for c in f["obs/batch/categories"][:]]
    bcodes = f["obs/batch/codes"][:].astype(np.int16)

    n_cells, n_genes = f["X"].attrs["shape"]
    n_groups = len(cats)
    indptr = f["X/indptr"][:]
    Xd, Xi = f["X/data"], f["X/indices"]

    gene_sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    group_cells = np.bincount(codes, minlength=n_groups).astype(np.int64)
    cell_totals = np.zeros(n_cells, dtype=np.float64)
    cell_nnz = np.zeros(n_cells, dtype=np.int32)

    print(f"cells={n_cells} genes={n_genes} groups={n_groups} nnz={indptr[-1]:,}", flush=True)
    for start in range(0, n_cells, ROWS):
        stop = min(start+ROWS, n_cells)
        lo, hi = int(indptr[start]), int(indptr[stop])
        d = Xd[lo:hi].astype(np.float64)
        idx = Xi[lo:hi]
        off = indptr[start:stop+1] - lo
        for r in range(stop-start):
            a, b = int(off[r]), int(off[r+1])
            if b > a:
                g = codes[start+r]
                np.add.at(gene_sums[g], idx[a:b], d[a:b])
                cell_totals[start+r] = d[a:b].sum()
                cell_nnz[start+r] = b-a
        if (start//ROWS) % 10 == 0:
            el = time.time()-t0
            pct = stop/n_cells
            print(f"  {stop:>7}/{n_cells} ({pct:5.1%}) {el:6.1f}s "
                  f"eta {el/max(pct,1e-9)*(1-pct):6.1f}s", flush=True)

np.savez_compressed(
    OUT/"group_sums.npz",
    gene_sums=gene_sums, group_cells=group_cells,
    cell_totals=cell_totals, cell_nnz=cell_nnz,
    codes=codes, bcodes=bcodes,
    cats=np.array(cats, dtype=object), genes=np.array(genes, dtype=object),
    bcats=np.array(bcats, dtype=object),
)
print(f"done in {time.time()-t0:.1f}s -> {OUT/'group_sums.npz'}", flush=True)
