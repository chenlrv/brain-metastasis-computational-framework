"""Write a memory-light cluster prior (clusters.csv) for decontX.

decontX only runs its own UMAP + cell-type estimation when no `z` is supplied, and
that internal step blows up (std::bad_alloc) on the large slices. Supplying `z`
skips it entirely.

A scanpy Leiden prior was not used: its neighbour-graph step OOMs at this cell count.
TruncatedSVD + MiniBatchKMeans never builds an N x N graph and stays well inside
RAM, and the same partition is applied to all six slices. decontX only needs a
coarse population partition, not a biologically tuned clustering.

Row order/count matches the exported counts.mtx exactly (same kept-cell mask).

Run: conda run -n thesis_research python analysis/03_data_quality/contamination_decontx/write_cluster_labels.py
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_decontx as m  # noqa: E402

N_COMPS = 30
N_CLUSTERS = 25
SEED = 0


def cluster_slice(name, path):
    work = m.workdir_for(name)
    cfile = os.path.join(work, "clusters.csv")
    if not os.path.exists(os.path.join(work, "counts.mtx")):
        print(f"{name}: no counts.mtx, skipping")
        return
    adata, _, _, _ = m._load(path)
    keep = np.load(os.path.join(work, "kept_cells.npy"))
    adata = adata[keep]
    X = sp.csr_matrix(adata.X)

    # library-size normalise -> log1p (sparse throughout)
    totals = np.asarray(X.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    X = sp.diags(np.median(totals) / totals) @ X
    X = sp.csr_matrix(X)
    X.data = np.log1p(X.data)

    emb = TruncatedSVD(n_components=N_COMPS, random_state=SEED).fit_transform(X)
    z = MiniBatchKMeans(n_clusters=N_CLUSTERS, random_state=SEED,
                        batch_size=10000, n_init=3).fit_predict(emb)

    pd.DataFrame({"cluster": z.astype(int)}).to_csv(cfile, index=False)
    sizes = np.bincount(z, minlength=N_CLUSTERS)
    print(f"{name}: {X.shape[0]:,} cells -> {len(np.unique(z))} clusters "
          f"(min {sizes[sizes > 0].min():,}, max {sizes.max():,})  -> {cfile}")


if __name__ == "__main__":
    for name, path in m.SLICES.items():
        cluster_slice(name, path)
    print("done.")
