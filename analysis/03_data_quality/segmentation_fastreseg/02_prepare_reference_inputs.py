"""FastReseg inputs for the five slice-1 FOVs (451, 512, 514, 515, 523).

FastReseg (Wu et al. 2025; NanoString R package) needs
  * a cells x genes count matrix for the whole dataset, with a cluster label per
    cell, from which it estimates one reference expression profile per cluster;
  * the transcript table of the FOVs to re-segment, with coordinates, z-plane,
    cell ID and gene.
Counts and clusters are those of all QC-passed slice-1 cells; the clusters are
the vendor's own InSituType cell typing shipped with the export, so no
clustering choice is made here. Transcripts are the raw vendor transcripts of
the five FOVs (automation/outputs/segmentation/tx_by_fov/), restricted to the 958
panel genes (negative and system-control probes removed).

Writes only to automation/outputs/segmentation_fastreseg/inputs/ and refuses to run
if that folder already holds files.

Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/02_prepare_reference_inputs.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

ROOT = pathlib.Path(PROJECT_ROOT_STR)
COUNTS_H5AD = ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad"
TX_DIR = ROOT / "automation/outputs/segmentation/tx_by_fov"
OUT = ROOT / "automation/outputs/segmentation_fastreseg/inputs"
FOVS = [451, 512, 514, 515, 523]
CLUST_COL = "RNA_Basic.run_Cell.Typing.InSituType.1_1_clusters"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if any(OUT.iterdir()):
        raise FileExistsError(f"{OUT} is not empty; refusing to overwrite")

    a = ad.read_h5ad(COUNTS_H5AD)
    cell_ids = a.obs["cell"].astype(str).to_numpy()          # e.g. c_1_514_7
    if len(set(cell_ids)) != len(cell_ids):
        raise ValueError("cell IDs are not unique")
    genes = list(a.var_names)
    X = sp.csr_matrix(a.X).astype(np.int32)
    scipy.io.mmwrite(str(OUT / "counts.mtx"), X)               # cells x genes
    pd.Series(cell_ids).to_csv(OUT / "cells.csv", index=False, header=["cell"])
    pd.Series(genes).to_csv(OUT / "genes.csv", index=False, header=["gene"])
    pd.DataFrame({"cell": cell_ids, "clust": a.obs[CLUST_COL].astype(str).to_numpy()}) \
        .to_csv(OUT / "clust.csv", index=False)
    print(f"counts: {X.shape[0]:,} cells x {X.shape[1]} genes; "
          f"{a.obs[CLUST_COL].nunique()} vendor clusters")

    gene_set = set(genes)
    info = []
    for fov in FOVS:
        tx = pd.read_csv(TX_DIR / f"tx_fov{fov}.csv",
                         usecols=["fov", "cell_ID", "cell", "x_global_px", "y_global_px", "z", "target"])
        n0 = len(tx)
        tx = tx[tx["target"].isin(gene_set)].copy()
        tx["transcript_id"] = [f"t_{fov}_{i}" for i in range(len(tx))]
        path = OUT / f"tx_fov{fov}.csv"
        tx[["transcript_id", "cell", "cell_ID", "x_global_px", "y_global_px", "z", "target"]] \
            .to_csv(path, index=False)
        info.append({"file_path": str(path).replace("\\", "/"), "slide": 1, "fov": fov,
                     "stage_X": 0.0, "stage_Y": 0.0})
        print(f"FOV {fov}: {len(tx):,} of {n0:,} transcripts on panel genes; "
              f"{tx.loc[tx.cell_ID != 0, 'cell'].nunique():,} segmented cells")
    pd.DataFrame(info).to_csv(OUT / "transDF_fileInfo.csv", index=False)
    print(f"written to {OUT}")


if __name__ == "__main__":
    main()
