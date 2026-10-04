"""DecontX partition sensitivity -- does the probe-anomaly conclusion depend on the
population groups (z) supplied to DecontX?

DecontX is run on slice 1 with four partitions and, for each, the drop in the
percentage of cells positive for tdTomato, GFP, Lyve1 and Meg3 is reported:

  k=10, k=50   k-means on the same 30 truncated-SVD components as the thesis run
               (write_cluster_labels.py, seed 0), with k changed
  k=25         the thesis partition (score_genes/decontx_clusters/slice_1_clusters.csv)
  leiden5      the earlier five-cluster Leiden partition
               (score_genes/decontx_clusters/sensitivity/slice_1_leiden5_clusters.csv)

Positivity is at least one count, before and after correction, over the cells
in the exported count matrix; corrected counts are rounded as in assembly.

Requires the DecontX export of slice 1 (run_decontx.py export), i.e.
resources/cache/decontx/slice_1_work/{counts.mtx, kept_cells.npy}, and Rscript.
DecontX runs are cached under resources/cache/decontx/sensitivity/ and skipped
when their output already exists.

Saves: thesis_plots/decontx_partition_sensitivity.csv

Run (repository root on PYTHONPATH):
  conda run -n thesis_research python analysis/03_data_quality/contamination_decontx/partition_sensitivity.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD

ROOT = PROJECT_ROOT_STR
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_decontx as m  # noqa: E402
import write_cluster_labels as w  # noqa: E402

SLICE = "slice_1"
WORK = m.workdir_for(SLICE)
OUT_PARENT = os.path.join(m.OUT_DIR, "sensitivity")
CLUSTERS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "decontx_clusters")
OUT_CSV = os.path.join(ROOT, "thesis_plots", "decontx_partition_sensitivity.csv")
PROBES = ["tdTomato", "GFP", "Lyve1", "Meg3"]
K_EXTRA = (10, 50)


def kmeans_labels(adata, ks):
    """Same preprocessing and seed as write_cluster_labels.py, for each k."""
    X = sp.csr_matrix(adata.X)
    totals = np.asarray(X.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    X = sp.csr_matrix(sp.diags(np.median(totals) / totals) @ X)
    X.data = np.log1p(X.data)
    emb = TruncatedSVD(n_components=w.N_COMPS, random_state=w.SEED).fit_transform(X)
    return {k: MiniBatchKMeans(n_clusters=k, random_state=w.SEED, batch_size=10000,
                               n_init=3).fit_predict(emb).astype(int) for k in ks}


def prepare(label, z):
    """One parent dir per partition, so decontx_model.R processes only that run."""
    parent = os.path.join(OUT_PARENT, label)
    work = os.path.join(parent, SLICE + "_work")
    os.makedirs(work, exist_ok=True)
    pd.DataFrame({"cluster": np.asarray(z, dtype=int)}).to_csv(
        os.path.join(work, "clusters.csv"), index=False)
    if not os.path.exists(os.path.join(work, "counts.mtx")):
        shutil.copyfile(os.path.join(WORK, "counts.mtx"), os.path.join(work, "counts.mtx"))
    return parent, work


def main():
    if not os.path.exists(os.path.join(WORK, "counts.mtx")):
        raise FileNotFoundError(f"{WORK}/counts.mtx missing -- run "
                                "analysis/03_data_quality/contamination_decontx/run_decontx.py export first")
    adata, _, _, _ = m._load(m.SLICES[SLICE])
    adata = adata[np.load(os.path.join(WORK, "kept_cells.npy"))]
    var = list(adata.var_names)

    partitions = {f"k={k}": z for k, z in kmeans_labels(adata, K_EXTRA).items()}
    partitions["k=25"] = pd.read_csv(os.path.join(CLUSTERS, f"{SLICE}_clusters.csv"))["cluster"]
    partitions["leiden5"] = pd.read_csv(
        os.path.join(CLUSTERS, "sensitivity", f"{SLICE}_leiden5_clusters.csv"))["cluster"]

    rscript = m.find_rscript()
    idx = [var.index(p) for p in PROBES]
    raw = mmread(os.path.join(WORK, "counts.mtx")).tocsr()[idx]       # genes x cells
    raw_pos = 100 * (raw > 0).mean(axis=1).A1

    rows = []
    for label in ("k=10", "k=25", "k=50", "leiden5"):
        z = np.asarray(partitions[label])
        if len(z) != adata.n_obs:
            raise ValueError(f"{label}: {len(z)} labels != {adata.n_obs} cells")
        parent, work = prepare(label.replace("=", ""), z)
        if not os.path.exists(os.path.join(work, "decontx_counts.mtx")):
            print(f"{label}: running DecontX ...")
            subprocess.run([rscript, m.R_SCRIPT, parent], check=True)
        dec = mmread(os.path.join(work, "decontx_counts.mtx")).tocsr()[idx]
        dec.data = np.rint(dec.data)
        pos = 100 * (dec > 0).mean(axis=1).A1
        cont = pd.read_csv(os.path.join(work, "decontx_contamination.csv"))["contamination"]
        row = {"partition": label, "n_clusters": len(np.unique(z)),
               "median_contamination": float(np.median(cont))}
        for p, r, c in zip(PROBES, raw_pos, pos):
            row[f"{p}_raw_pct"] = r
            row[f"{p}_drop_pp"] = r - c
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    show = ["partition", "n_clusters", "median_contamination"] + [f"{p}_drop_pp" for p in PROBES]
    print(f"{SLICE}: drop in % positive cells (percentage points)")
    print(out[show].round(3).to_string(index=False))
    print("\nSaved:", OUT_CSV)


if __name__ == "__main__":
    main()
