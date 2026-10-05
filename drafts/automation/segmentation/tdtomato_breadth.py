"""tdTomato breadth before and after each segmentation-correction test.

The tdTomato anomaly (Figure 17) is its breadth: it is detected in as many cells
as the pan-myeloid transcripts Cx3cr1 and Csf1r, or more. If that breadth came
from spillover out of genuinely tdTomato-positive cells, correcting
misassignment should shrink it relative to Cx3cr1 and Csf1r. For the five
slice-1 FOVs, non-tumor cells (final calls), this reports the percentage of
cells with >= 1 count of each gene:
  * FastReseg -- before = vendor assignment of the FastReseg input transcripts,
    after = FastReseg's updated_cellID (agents/outputs/segmentation_fastreseg/);
  * permissive reassignment -- reassign_v3.py's configuration, recomputed with
    its own functions (the procedure is deterministic).
Writes only to agents/outputs/segmentation_tdtomato_breadth/ (must not exist).

Run: conda run -n thesis_research python agents/segmentation/tdtomato_breadth.py
"""
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import reassign_v3  # noqa: E402,F401  (sets the v3 configuration on reassign_v2)
import reassign_v2 as rv  # noqa: E402

OUT = rv.ROOT / "agents/outputs/segmentation_tdtomato_breadth"
FR = rv.ROOT / "agents/outputs/segmentation_fastreseg"
GENES = ["tdTomato", "Cx3cr1", "Csf1r"]


def pct_pos(m, genes, cells, tumor):
    keep = np.array([c.rsplit("_g", 1)[0] not in tumor for c in cells]) & (m.sum(1) > 0)
    return {g: 100 * float((m[keep, genes.index(g)] > 0).mean()) for g in GENES}, int(keep.sum())


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} already exists; refusing to overwrite")
    OUT.mkdir(parents=True)
    a = ad.read_h5ad(rv.CACHE, backed="r")
    genes = list(a.var_names)
    tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    a.file.close()
    gi = {g: i for i, g in enumerate(genes)}

    rows = []
    agg = {k: [] for k in ("fr_before", "fr_after", "pr_before", "pr_after")}
    for k, fov in enumerate(rv.FOVS, start=1):
        # FastReseg
        orig = pd.read_csv(FR / "inputs" / f"tx_fov{fov}.csv", usecols=["cell", "cell_ID", "target"])
        orig = orig[orig.cell_ID != 0]
        upd = pd.read_csv(FR / "fastreseg_out" / f"{k}_updated_transDF.csv",
                          usecols=["UMI_transID", "target", "updated_cellID"]).dropna()
        assert upd.UMI_transID.str.startswith(f"t_{fov}_").all()
        for tag, cell, tgt in (("fr_before", orig.cell, orig.target), ("fr_after", upd.updated_cellID, upd.target)):
            cells = sorted(set(cell))
            m = np.zeros((len(cells), len(genes)))
            ci = pd.Index(cells).get_indexer(cell); gidx = pd.Index(genes).get_indexer(tgt)
            ok = gidx >= 0
            np.add.at(m, (ci[ok], gidx[ok]), 1)
            p, n = pct_pos(m, genes, cells, tumor)
            agg[tag].append((p, n))
            rows.append({"method": "FastReseg", "stage": tag.split("_")[1], "fov": fov, "n_cells": n, **p})
        # permissive reassignment (v3 configuration)
        d = rv.load_fov(fov, gi)
        after = rv.reassign(d, len(genes), **rv.CANON)
        for tag, assign in (("pr_before", d["home"]), ("pr_after", after)):
            m = rv.counts(assign, d["gene"], len(d["cells"]), len(genes))
            p, n = pct_pos(m, genes, list(d["cells"]), tumor)
            agg[tag].append((p, n))
            rows.append({"method": "permissive", "stage": tag.split("_")[1], "fov": fov, "n_cells": n, **p})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_fov.csv", index=False)
    pooled = []
    for tag, items in agg.items():
        n = sum(x[1] for x in items)
        p = {g: sum(x[0][g] * x[1] for x in items) / n for g in GENES}   # cell-weighted
        pooled.append({"method": "FastReseg" if tag.startswith("fr") else "permissive",
                       "stage": tag.split("_")[1], "n_cells": n, **p})
    pooled = pd.DataFrame(pooled)
    pooled.to_csv(OUT / "pooled.csv", index=False)
    print(pooled.round(2).to_string(index=False))
    print("\nper FOV:\n", df.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
