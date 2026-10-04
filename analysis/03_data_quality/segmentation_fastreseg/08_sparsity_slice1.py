"""Per-cell sparsity of the vendor segmentation, slice 1, QC-passed non-tumor cells.

Reports (1) the share of intracellular transcripts that are the only copy of
their gene in their cell, and (2) the share of Lyve1-positive cells carrying
exactly one Lyve1 transcript. These support the statement that per-cell
expression profiles are too sparse for cell-level transcript reassignment.
Inputs: the slice-1 transcript files of 05_prepare_slice1_inputs.py (vendor cell
assignment). Writes only agents/outputs/segmentation_fastreseg_slice1/eval_qcpassed/sparsity.json
(must not exist).

Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/08_sparsity_slice1.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import json
import pathlib

import anndata as ad
import pandas as pd

ROOT = pathlib.Path(PROJECT_ROOT_STR)
B = ROOT / "agents/outputs/segmentation_fastreseg_slice1"
OUT = B / "eval_qcpassed" / "sparsity.json"


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} exists; refusing to overwrite")
    a = ad.read_h5ad(ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad", backed="r")
    tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    qc = set(a.obs["cell"].astype(str))
    a.file.close()
    info = pd.read_csv(B / "inputs" / "transDF_fileInfo.csv")
    total = single = lyve_pos = lyve_single = 0
    for fov in info.fov:
        o = pd.read_csv(B / "inputs" / f"tx_fov{fov}.csv", usecols=["cell", "cell_ID", "target"])
        o = o[(o.cell_ID != 0) & ~o.cell.isin(tumor) & o.cell.isin(qc)]
        c = o.groupby(["cell", "target"]).size()
        total += int(c.sum())
        single += int((c == 1).sum())
        ly = c[c.index.get_level_values("target") == "Lyve1"]
        lyve_pos += len(ly)
        lyve_single += int((ly == 1).sum())
    res = {"n_transcripts": total, "pct_single_copy_transcripts": 100 * single / total,
           "n_lyve1_positive_cells": lyve_pos, "pct_lyve1_cells_single_transcript": 100 * lyve_single / lyve_pos}
    json.dump(res, open(OUT, "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
