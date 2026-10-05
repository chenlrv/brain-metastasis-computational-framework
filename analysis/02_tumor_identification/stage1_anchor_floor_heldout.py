"""Anchor floor on held-out data, at the adopted candidate thresholds (0.2 / 0.08).

Companion to stage1_threshold_heldout.py, whose anchor sweep ran at the candidate
setting selected from its grid. Here the anchor floor (0.30-0.50) is varied at the
thesis candidate setting; at each floor the classifier is retrained and false calls
are counted on control-slice look-alikes withheld from training (out-of-fold and
leave-one-control-slice-out), with the same functions as stage1_threshold_heldout.py.
The 0.40 row must reproduce the published numbers; otherwise the script stops.

Writes only to thesis_plots/stage1_sensitivity_heldout_v2_anchor/ (must not exist).
Run: conda run -n thesis_research python analysis/02_tumor_identification/stage1_anchor_floor_heldout.py
"""
import os

for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")

import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import stage1_threshold_heldout as s1  # noqa: E402

OUT_DIR = s1.BASE_DIR / "thesis_plots" / "stage1_sensitivity_heldout_v2_anchor"


def main():
    if OUT_DIR.exists():
        raise FileExistsError(f"{OUT_DIR} already exists; refusing to overwrite")
    REF = s1.reference_superset()
    CAND = s1.candidate_superset(REF[3])
    floor, margin, _ = s1.ADOPTED
    rows = []
    for a in s1.ANCHOR_FLOORS:
        r = s1.evaluate(floor, margin, a, REF, CAND)
        r["anchor_purity_pct"] = s1.anchor_purity(a)
        rows.append(r)
        print(f"  anchor floor {a:.2f} done", flush=True)
    df = pd.DataFrame(rows)
    base = df[df["anchor_floor"] == s1.ADOPTED[2]].iloc[0]
    expected = dict(n_lookalikes=630, n_anchors=863, oof_fp=13, loso_fp_L34=12, loso_fp_L321=18, calls_all4=20873)
    bad = {k: (base[k], v) for k, v in expected.items() if base[k] != v}
    if bad:
        raise SystemExit(f"published numbers NOT reproduced (got, expected): {bad}")
    print("sanity check passed: published numbers reproduced exactly", flush=True)
    OUT_DIR.mkdir()
    df.to_csv(OUT_DIR / "sweep_anchor_floor_adopted_candidates.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(df[["anchor_floor", "n_anchors", "anchor_purity_pct", "oof_fp", "loso_fp",
                  "loso_fp_L321", "loso_fp_L34", "calls_heldout", "calls_all4"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
