"""Held-out specificity of the Stage-2 classifiers in sham-injected control tissue.

The control-slice tumor candidates scored in Stage 3 are the same 630 cells that
form the negative training class (both are selected by the candidate filter of
identify_tumor_cells.py), so "zero tumor calls in the control slices" measures
fit to the training negatives, not specificity. This script estimates the
specificity on control cells the classifier did not see, in two ways:

  A  Out-of-fold: 5-fold stratified CV on the joint reference pool (same folds
     and seed as fig08_classifier_comparison.py). Every look-alike is scored by a
     model that did not train on it.
  B  Leave-one-control-slice-out: the anchors are kept, the look-alikes of one
     control slice are used for training, and the look-alikes of the other
     control slice are scored. Run in both directions. Both control slices come
     from the same animal (mouse 1), so this tests generalization across slides
     (L321 vs L34), not across animals.

Models are fitted and applied with fit_classifiers / score_candidates from
fig09_14_spatial_refinement.py, so the decision rules (P(tumor) > 0.5; the
LogReg + KNN hybrid calls tumor only if both signals exceed 0.5) match the
candidate-retention figures exactly. Within each fold the PCA used by the KNN
hybrid is fitted on the training cells only.

Saves: thesis_plots/control_specificity/oof_false_positives.csv
       thesis_plots/control_specificity/loso_false_positives.csv

Run: conda run -n thesis_research python analysis/02_tumor_identification/table03_heldout_specificity.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import os

# Single-threaded BLAS/OpenMP: the multi-threaded runtimes deadlock on this machine.
for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")

import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fig09_14_spatial_refinement import fit_classifiers, score_candidates  # noqa: E402
from thesis_research.pipeline.cell_type_annotation.tumor_cells.identify_tumor_cells import (  # noqa: E402
    _get_healthy_ref_ids,
    _get_tumor_ref_ids,
)

BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
SLIDE_CACHE = BASE_DIR / "resources" / "cache"
OUT_DIR = BASE_DIR / "thesis_plots" / "control_specificity"

RANDOM_STATE = 42
N_SPLITS = 5
CONTROL_SLICE = {"L321": 3, "L34": 4}
MODEL_ORDER = ["LogReg", "LogReg + PCA", "LogReg + KNN (PCA)", "Random Forest", "XGBoost"]


def build_reference_pool():
    """X, y (1 = tumor anchor) and, for each cell, its slide.

    Same construction as figure_3_model_comparison.build_reference_pool, but the
    slide of each cell is kept so the look-alikes can be split by control slice.
    """
    slide_ids = ["L321", "L34"]
    healthy = set().union(*(_get_healthy_ref_ids(s) for s in slide_ids))
    tumor = set().union(*(_get_tumor_ref_ids(s) for s in slide_ids))
    adatas = {}
    for sid in slide_ids:
        adata = ad.read_h5ad(SLIDE_CACHE / f"sample_{sid}_adata.h5ad")
        sub = adata[np.asarray(adata.obs_names.isin(healthy | tumor))].copy()
        del adata
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        adatas[sid] = sub
    joint = ad.concat(adatas, join="inner", label="slide_id")
    X = joint.X.toarray() if hasattr(joint.X, "toarray") else np.asarray(joint.X)
    y = np.asarray(joint.obs_names.isin(tumor)).astype(int)
    slide = joint.obs["slide_id"].astype(str).to_numpy()
    print(f"reference pool: {len(y):,} cells "
          f"({int(y.sum()):,} anchors, {int((y == 0).sum()):,} look-alikes)")
    return X, y, slide


def out_of_fold(X, y, slide):
    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    is_tumor = {m: np.zeros(len(y), dtype=bool) for m in MODEL_ORDER}
    for k, (tr, va) in enumerate(cv.split(X, y), start=1):
        print(f"  fold {k}/{N_SPLITS}", flush=True)
        calls = score_candidates(X[va], fit_classifiers(X[tr], y[tr]))
        for m in MODEL_ORDER:
            is_tumor[m][va] = calls[m]

    rows = []
    neg = y == 0
    for m in MODEL_ORDER:
        row = {"model": m}
        for sid, s in CONTROL_SLICE.items():
            mask = neg & (slide == sid)
            row[f"slice{s}_fp"] = int(is_tumor[m][mask].sum())
            row[f"slice{s}_n"] = int(mask.sum())
        row["total_fp"] = int(is_tumor[m][neg].sum())
        row["total_n"] = int(neg.sum())
        row["fp_rate_pct"] = round(100 * row["total_fp"] / row["total_n"], 2)
        row["anchors_retained_pct"] = round(100 * is_tumor[m][y == 1].mean(), 2)
        rows.append(row)
    return pd.DataFrame(rows)


def leave_one_slice_out(X, y, slide):
    rows = []
    for train_sid, test_sid in (("L321", "L34"), ("L34", "L321")):
        train = (y == 1) | ((y == 0) & (slide == train_sid))
        test = (y == 0) & (slide == test_sid)
        print(f"  train look-alikes: slice {CONTROL_SLICE[train_sid]}; "
              f"test: slice {CONTROL_SLICE[test_sid]}", flush=True)
        calls = score_candidates(X[test], fit_classifiers(X[train], y[train]))
        for m in MODEL_ORDER:
            fp, n = int(calls[m].sum()), int(test.sum())
            rows.append({
                "model": m,
                "train_control_slice": CONTROL_SLICE[train_sid],
                "test_control_slice": CONTROL_SLICE[test_sid],
                "n_train_lookalikes": int(((y == 0) & (slide == train_sid)).sum()),
                "n_test": n,
                "fp": fp,
                "fp_rate_pct": round(100 * fp / n, 2),
            })
    return pd.DataFrame(rows)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    X, y, slide = build_reference_pool()

    print("\nA. Out-of-fold (5-fold stratified CV)")
    oof = out_of_fold(X, y, slide)
    oof.to_csv(OUT_DIR / "oof_false_positives.csv", index=False)
    print(oof.to_string(index=False))

    print("\nB. Leave-one-control-slice-out")
    loso = leave_one_slice_out(X, y, slide)
    loso.to_csv(OUT_DIR / "loso_false_positives.csv", index=False)
    print(loso.to_string(index=False))
    print(f"\nSaved to {OUT_DIR}")


if __name__ == "__main__":
    main()
