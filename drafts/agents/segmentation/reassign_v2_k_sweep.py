"""Sweep of the number of candidate cells (K) for the transcript reassignment test.

Runs reassign_v2.py's procedure unchanged except for K, the maximum number of
nearest cells (within the 14.4 um radius) considered as candidates for each
transcript: K = 5, 6 (the configuration used), 10, 15, 20. All other parameters
stay at the canonical values. Reported on FOV 514, as the other sweeps, and
pooled over the five FOVs.

Writes only to agents/outputs/segmentation_reassign_v2_ksweep/ (must not exist).

Run: conda run -n thesis_research python agents/segmentation/reassign_v2_k_sweep.py
"""
import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import reassign_v2 as rv  # noqa: E402
from seg_metrics import compute_all  # noqa: E402

OUT = rv.ROOT / "agents/outputs/segmentation_reassign_v2_ksweep"
KS = [5, 6, 10, 15, 20]


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} already exists; refusing to overwrite")
    OUT.mkdir(parents=True)
    a = ad.read_h5ad(rv.CACHE, backed="r")
    genes = list(a.var_names)
    tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    a.file.close()
    gi = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)
    data = {fov: rv.load_fov(fov, gi) for fov in rv.FOVS}

    rows = []
    for k in KS:
        rv.K_NEIGH = k                       # reassign() reads the module-level K_NEIGH
        before, after, moved, n_tx = [], [], 0, 0
        for fov, d in data.items():
            assign = rv.reassign(d, n_genes, **rv.CANON)
            mb, rb = rv.readouts(d["home"], d, n_genes, genes, tumor)
            ma, ra = rv.readouts(assign, d, n_genes, genes, tumor)
            before.append(mb); after.append(ma)
            intra = d["home"] >= 0
            moved += int((intra & (assign != d["home"])).sum()); n_tx += len(assign)
            if fov == rv.SWEEP_FOV:
                p514 = (rb["purity"]["purity_own_frac"], ra["purity"]["purity_own_frac"])
        pb = compute_all(np.vstack(before), genes)["purity"]["purity_own_frac"]
        pa = compute_all(np.vstack(after), genes)["purity"]["purity_own_frac"]
        rows.append({"k": k, "fov514_before": p514[0], "fov514_after": p514[1],
                     "fov514_delta": p514[1] - p514[0], "pooled_before": pb, "pooled_after": pa,
                     "pooled_delta": pa - pb, "pct_moved_between_cells": 100 * moved / n_tx})
        r = rows[-1]
        print(f"K={k:2d}: FOV 514 purity {r['fov514_before']:.3f} -> {r['fov514_after']:.3f} "
              f"({r['fov514_delta']:+.3f}); pooled {pb:.3f} -> {pa:.3f} ({pa - pb:+.3f}); "
              f"{r['pct_moved_between_cells']:.2f}% of transcripts moved between cells", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "k_sweep.csv", index=False)
    print(f"written to {OUT}")


if __name__ == "__main__":
    main()
