"""Transcript reassignment test with up to FIVE candidate cells per transcript.

Runs reassign_v2.py unchanged except K_NEIGH = 5 (the thesis text uses five
candidates) and writes to a new folder, agents/outputs/segmentation_reassign_v2_k5/
(must not exist). Also reports how often a transcript's own cell is among its
candidates under this setting. The six-candidate run in segmentation_reassign_v2/
is left untouched.

Run: conda run -n thesis_research python agents/segmentation/reassign_v2_k5.py
"""
import json
import pathlib
import sys

import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import reassign_v2 as rv  # noqa: E402

rv.K_NEIGH = 5
rv.OUT = rv.ROOT / "agents/outputs/segmentation_reassign_v2_k5"


def own_cell_coverage():
    a = ad.read_h5ad(rv.CACHE, backed="r")
    gi = {g: i for i, g in enumerate(a.var_names)}
    a.file.close()
    n_intra = n_missing = 0
    for fov in rv.FOVS:
        d = rv.load_fov(fov, gi)
        _, idx = d["tree"].query(d["xy"], k=rv.K_NEIGH,
                                 distance_upper_bound=rv.CANON["radius_um"] / rv.PX_UM)
        intra = d["home"] >= 0
        own = (idx == d["home"][:, None]).any(1)
        n_intra += int(intra.sum())
        n_missing += int((intra & ~own).sum())
    return {"k": rv.K_NEIGH, "n_intracellular": n_intra, "own_cell_not_candidate": n_missing,
            "pct_own_cell_candidate": 100 * (1 - n_missing / n_intra)}


if __name__ == "__main__":
    rv.main()                                   # refuses if the output folder exists
    cov = own_cell_coverage()
    json.dump(cov, open(rv.OUT / "own_cell_coverage.json", "w"), indent=2)
    print(f"own cell among candidates: {cov['pct_own_cell_candidate']:.2f}% of "
          f"{cov['n_intracellular']:,} intracellular transcripts")
