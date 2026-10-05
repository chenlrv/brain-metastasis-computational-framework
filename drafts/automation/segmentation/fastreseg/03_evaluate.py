"""Thesis read-outs before and after FastReseg, five slice-1 FOVs, non-tumor cells.

"Before" is the vendor segmentation: every intracellular transcript of the
FastReseg input files (01_prep_inputs.py), assigned to its vendor cell. "After"
is FastReseg's result: the per-FOV updated transcript tables written by
02_run_fastreseg.R, which hold each remaining intracellular transcript with its
updated_cellID; transcripts FastReseg trimmed to extracellular space are absent
from these tables. Any difference is therefore due to FastReseg alone. Tumor cells are
excluded using the final Stage-3 calls (with_tumor_prediction_final/). The
read-outs are those of agents/segmentation/seg_metrics.py, as used for the
earlier reassignment test: myeloid marker purity, the GFP-tdTomato and
GFP-Cx3cr1 correlations, and Lyve1 breadth.

Writes only to agents/outputs/segmentation_fastreseg/eval/ (must be empty).

Run: conda run -n thesis_research python agents/segmentation/fastreseg/03_evaluate.py
"""
import json
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd

ROOT = pathlib.Path(r"D:/thesis-research")
sys.path.insert(0, str(ROOT / "agents/segmentation"))
from seg_metrics import compute_all  # noqa: E402

BASE = ROOT / "agents/outputs/segmentation_fastreseg"
RES = BASE / "fastreseg_out"
INP = BASE / "inputs"
OUT = BASE / "eval_v2"   # eval/ holds a first run whose "before" omitted the trimmed transcripts
TUMOR_H5AD = ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad"
FOVS = [451, 512, 514, 515, 523]   # FastReseg numbers its output files 1..5 in this order


def matrix(cell, gene, cells, genes):
    """Cells x genes count matrix from one row per transcript."""
    ci = pd.Index(cells).get_indexer(cell)
    gi = pd.Index(genes).get_indexer(gene)
    ok = (ci >= 0) & (gi >= 0)
    m = np.zeros((len(cells), len(genes)), dtype=np.int32)
    np.add.at(m, (ci[ok], gi[ok]), 1)
    return m


def summary(r):
    p = r["purity"]
    return {"n_cells": r["n_cells"], "purity": p["purity_own_frac"], "n_purity_cells": p["n_valid"],
            "gfp_tdtomato_r": r.get("F1_full_r"), "gfp_cx3cr1_r": r.get("F2_full_r"),
            "lyve1_pct_pos": r.get("F3_lyve_pct_pos"), "lyve1_pct_lacking_bam": r.get("F3_lyve_pct_lacking_bam"),
            "median_total_counts": r["median_total_counts"]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if any(OUT.iterdir()):
        raise FileExistsError(f"{OUT} is not empty; refusing to overwrite")

    t = ad.read_h5ad(TUMOR_H5AD, backed="r")
    tumor_cells = set(t.obs.loc[t.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    genes = list(t.var_names)
    t.file.close()

    rows, pooled = [], {"before": [], "after": []}
    for k, fov in enumerate(FOVS, start=1):
        tx = pd.read_csv(RES / f"{k}_updated_transDF.csv",
                         usecols=["UMI_transID", "target", "UMI_cellID", "updated_cellID"])
        if not tx["UMI_transID"].str.startswith(f"t_{fov}_").all():
            raise ValueError(f"file {k} is not FOV {fov}")
        orig = pd.read_csv(INP / f"tx_fov{fov}.csv", usecols=["transcript_id", "cell", "cell_ID", "target"])
        orig = orig[orig["cell_ID"] != 0]                        # intracellular transcripts only
        before_cells = sorted(set(orig["cell"]))
        after_cells = sorted(set(tx["updated_cellID"].dropna()))
        mb = matrix(orig["cell"], orig["target"], before_cells, genes)
        ma = matrix(tx["updated_cellID"].dropna(), tx.loc[tx["updated_cellID"].notna(), "target"],
                    after_cells, genes)
        n_trimmed = len(orig) - int(tx["updated_cellID"].notna().sum())
        # non-tumor cells with at least one transcript, as in the earlier test;
        # a cell created by FastReseg inherits the tumor status of its parent cell
        parent = lambda c: "_".join(c.split("_")[:4])
        kb = np.array([parent(c) not in tumor_cells for c in before_cells]) & (mb.sum(1) > 0)
        ka = np.array([parent(c) not in tumor_cells for c in after_cells]) & (ma.sum(1) > 0)
        rb, ra = summary(compute_all(mb[kb], genes)), summary(compute_all(ma[ka], genes))
        n_moved = int((tx["UMI_cellID"] != tx["updated_cellID"]).sum())
        rows.append({"fov": fov, "n_intracellular": len(orig), "n_trimmed": n_trimmed,
                     "n_changed_cell": n_moved, "pct_changed_cell": 100 * n_moved / len(orig),
                     **{f"before_{key}": v for key, v in rb.items()},
                     **{f"after_{key}": v for key, v in ra.items()}})
        pooled["before"].append(mb[kb]); pooled["after"].append(ma[ka])
        print(f"FOV {fov}: {len(orig):,} intracellular transcripts, {n_trimmed:,} trimmed, "
              f"{n_moved:,} moved to another cell; non-tumor cells {kb.sum():,} -> {ka.sum():,}; "
              f"purity {rb['purity']:.3f} -> {ra['purity']:.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_fov_before_after.csv", index=False)
    pooled_res = {s: summary(compute_all(np.vstack(m), genes)) for s, m in pooled.items()}
    json.dump(pooled_res, open(OUT / "pooled_before_after.json", "w"), indent=2)

    print("\nmean over FOVs:  purity {:.3f} -> {:.3f}".format(df.before_purity.mean(), df.after_purity.mean()))
    for s in ("before", "after"):
        p = pooled_res[s]
        print(f"pooled {s:6s}: cells {p['n_cells']:,}, purity {p['purity']:.3f} (n={p['n_purity_cells']:,}), "
              f"GFP-tdTomato r {p['gfp_tdtomato_r']:.3f}, GFP-Cx3cr1 r {p['gfp_cx3cr1_r']:.3f}, "
              f"Lyve1+ {p['lyve1_pct_pos']:.1f}% ({p['lyve1_pct_lacking_bam']:.0f}% lacking Mrc1/Cd163)")
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
