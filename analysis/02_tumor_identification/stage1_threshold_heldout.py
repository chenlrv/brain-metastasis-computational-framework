"""Stage-1 thresholds justified on held-out data, with the classifier retrained per setting.

stage1_threshold_sensitivity.py kept the Stage-3 classifier fixed while the
candidate thresholds were varied, so at the adopted values every control-slice
candidate was a training cell and "zero control calls" was in-sample. Here every
setting is run as the pipeline would run it: the look-alikes (control-slice
candidates) are redefined, XGBoost is retrained, and specificity is measured on
look-alikes the classifier did not see.

Settings
  Grid  candidate score_tumor floor {0.10..0.30} x delta_score margin {0..0.15}
        (anchors fixed at score_tumor >= 0.4, margin > 0.08).
  C     anchor floor {0.30..0.50} at the candidate setting selected from the grid.

Per setting
  * held-out false calls in tumor-free tissue: out-of-fold (5-fold, seed 42)
    predictions for the look-alikes; leave-one-control-slice-out as a check;
  * tumor yield on slices 2 and 6, which never contribute training cells;
  * estimated false-discovery proportion (FDP) of the yield: the per-cell false
    call rate of each control slice (3 for slide L321, 4 for L34) times the
    non-called cells of the tumor slice on the same slide, over the tumor calls.

Decision rules (fixed before the run)
  Grid: among settings with estimated FDP <= 5%, the one with the most estimated
        true tumor calls (calls minus estimated false calls).
  C:    among anchor floors with anchor purity >= 99.5%, the lowest out-of-fold
        false calls; floors within one Poisson standard error of it are tied, and
        the lowest (largest anchor set) of those is chosen.

Sanity check: at 0.2 / 0.08 / 0.4 the published numbers must be reproduced
(630 look-alikes, 863 anchors, 13 out-of-fold and 12 + 18 held-out false
positives, 20,873 tumor calls); otherwise the script stops.

Writes only to thesis_plots/stage1_sensitivity_heldout_v2/ (must not exist; the
earlier stage1_sensitivity_heldout/ belongs to stage1_threshold_sensitivity_heldout.py).
Run: conda run -n thesis_research python analysis/02_tumor_identification/stage1_threshold_heldout.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import os

for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")

import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from scipy.stats import chi2
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from xgboost import XGBClassifier

BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
SLIDE_CACHE = BASE_DIR / "resources" / "cache"
OUT_DIR = BASE_DIR / "thesis_plots" / "stage1_sensitivity_heldout_v2"

RANDOM_STATE, N_SPLITS, PROB_THRESH = 42, 5, 0.5
SLICES = {1: "L321", 2: "L321", 3: "L321", 4: "L34", 5: "L34", 6: "L34"}
CONTROL = {"L321": 3, "L34": 4}          # look-alike source per slide
ANCHOR_SLICE = {"L321": 1, "L34": 5}     # anchor source per slide
HELDOUT_TUMOR = {"L321": 2, "L34": 6}    # tumor slices with no training cells
ANCHOR_DELTA = 0.08

FLOORS = [0.10, 0.15, 0.20, 0.25, 0.30]
MARGINS = [0.00, 0.04, 0.06, 0.08, 0.10, 0.15]
ANCHOR_FLOORS = [0.30, 0.35, 0.40, 0.45, 0.50]
ADOPTED = (0.20, 0.08, 0.40)
FDP_CAP, PURITY_MIN = 0.05, 99.5


def annot(s):
    slide = SLICES[s]
    df = pd.read_csv(BASE_DIR / "outputs" / "cell_annotation" / slide / "05" / str(s)
                     / f"slice_{s}_final_cell_annotations_refined_tabula_brain_tumor_avinoam.csv")
    n_cells = len(df)
    df = df[df["predicted_cell_type"] == "Tumor"].copy()
    df["next_best"] = df[["score_brain_struct", "score_brain_immune"]].max(axis=1)
    df["delta_score"] = (df["score_tumor"] - df["next_best"]).abs()
    df["cell_barcode"] = df["cell_barcode"].astype(str)
    return df, n_cells


ANNOT, N_CELLS = {}, {}
for _s in SLICES:
    ANNOT[_s], N_CELLS[_s] = annot(_s)


def candidates(s, floor, margin):
    d = ANNOT[s]
    return d[(d["score_tumor"] > floor) & (d["delta_score"] > margin) & (d["score_tumor"] > d["next_best"])]


def anchors(s, floor):
    d = ANNOT[s]
    return d[(d["score_tumor"] >= floor) & (d["delta_score"] > ANCHOR_DELTA) & (d["score_tumor"] > d["next_best"])]


def dense(x):
    return (x.toarray() if issparse(x) else np.asarray(x)).astype(np.float32)


def reference_superset():
    """Normalized expression for every cell any setting could train on, in pipeline order."""
    look = set().union(*(set(candidates(CONTROL[sl], min(FLOORS), min(MARGINS))["cell_barcode"]) for sl in CONTROL))
    anch = set().union(*(set(anchors(ANCHOR_SLICE[sl], min(ANCHOR_FLOORS))["cell_barcode"]) for sl in ANCHOR_SLICE))
    adatas = {}
    for slide in ["L321", "L34"]:
        a = ad.read_h5ad(SLIDE_CACHE / f"sample_{slide}_adata.h5ad")
        sub = a[np.asarray(a.obs_names.isin(look | anch))].copy()
        del a
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        adatas[slide] = sub
    joint = ad.concat(adatas, join="inner", label="slide_id")
    names = pd.Index(joint.obs_names.astype(str))
    print(f"reference superset: {len(names):,} cells", flush=True)
    return dense(joint.X), names, joint.obs["slide_id"].astype(str).to_numpy(), list(joint.var_names)


def candidate_superset(var_names):
    """Normalized expression for every tumor-slice candidate any setting could select."""
    out = {}
    for s in (1, 2, 5, 6):
        ids = set(candidates(s, min(FLOORS), min(MARGINS))["cell_barcode"])
        a = ad.read_h5ad(SLIDE_CACHE / f"slice_{s}_adata.h5ad")
        sub = a[np.asarray(a.obs_names.isin(ids))].copy()
        del a
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sub = sub[:, var_names].copy()
        out[s] = (pd.Index(sub.obs_names.astype(str)), dense(sub.X))
        print(f"  slice {s}: {len(ids):,} candidate cells at the loosest setting", flush=True)
    return out


def xgb():
    return XGBClassifier(random_state=RANDOM_STATE, n_jobs=4, verbosity=0, eval_metric="logloss")


def poisson_ci(k):
    lo = 0.0 if k == 0 else chi2.ppf(0.025, 2 * k) / 2
    return lo, chi2.ppf(0.975, 2 * (k + 1)) / 2


def evaluate(floor, margin, a_floor, REF, CAND):
    X_all, names, slide_of, _ = REF
    look_ids = {sl: set(candidates(CONTROL[sl], floor, margin)["cell_barcode"]) for sl in CONTROL}
    anch_ids = set().union(*(set(anchors(ANCHOR_SLICE[sl], a_floor)["cell_barcode"]) for sl in ANCHOR_SLICE))
    is_look = names.isin(look_ids["L321"] | look_ids["L34"])
    is_anch = names.isin(anch_ids) & ~is_look
    mask = np.asarray(is_look | is_anch)
    X, y, sl = X_all[mask], np.asarray(is_anch[mask]).astype(int), slide_of[mask]
    neg = y == 0

    cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    p = cross_val_predict(xgb(), X, y, cv=cv, method="predict_proba")[:, 1]
    oof_fp = {s: int(((p > PROB_THRESH) & neg & (sl == s)).sum()) for s in CONTROL}
    loso_fp = {}
    for train_s, test_s in (("L321", "L34"), ("L34", "L321")):
        tr = (y == 1) | (neg & (sl == train_s))
        te = neg & (sl == test_s)
        m = xgb().fit(X[tr], y[tr])
        loso_fp[test_s] = int((m.predict_proba(X[te])[:, 1] > PROB_THRESH).sum())

    full = xgb().fit(X, y)
    calls = {}
    for s, (cnames, cX) in CAND.items():
        sel = np.asarray(cnames.isin(set(candidates(s, floor, margin)["cell_barcode"])))
        calls[s] = int((full.predict_proba(cX[sel])[:, 1] > PROB_THRESH).sum()) if sel.any() else 0

    row = dict(cand_floor=floor, margin=margin, anchor_floor=a_floor,
               n_lookalikes=int(neg.sum()), n_anchors=int((y == 1).sum()),
               oof_fp=sum(oof_fp.values()), oof_fp_L321=oof_fp["L321"], oof_fp_L34=oof_fp["L34"],
               loso_fp=sum(loso_fp.values()), loso_fp_L321=loso_fp["L321"], loso_fp_L34=loso_fp["L34"],
               **{f"calls_slice{s}": c for s, c in calls.items()})
    row["calls_all4"] = sum(calls.values())
    row["calls_heldout"] = sum(calls[HELDOUT_TUMOR[s]] for s in HELDOUT_TUMOR)
    for tag, fp in (("oof", oof_fp), ("loso", loso_fp)):
        e_fp = sum(fp[s] / N_CELLS[CONTROL[s]] * (N_CELLS[HELDOUT_TUMOR[s]] - calls[HELDOUT_TUMOR[s]])
                   for s in CONTROL)
        k = sum(fp.values())
        scale = e_fp / k if k else sum((N_CELLS[HELDOUT_TUMOR[s]] - calls[HELDOUT_TUMOR[s]]) / N_CELLS[CONTROL[s]]
                                      for s in CONTROL) / 2
        lo, hi = poisson_ci(k)
        row[f"{tag}_est_false_calls"] = e_fp
        row[f"{tag}_fdp"] = e_fp / max(row["calls_heldout"], 1)
        row[f"{tag}_fdp_hi95"] = hi * scale / max(row["calls_heldout"], 1)
        row[f"{tag}_est_true_calls"] = row["calls_heldout"] - e_fp
    return row


def select_grid(df, tag):
    ok = df[df[f"{tag}_fdp"] <= FDP_CAP]
    return ok.loc[ok[f"{tag}_est_true_calls"].idxmax()]


def anchor_purity(a_floor):
    fp = sum(len(anchors(CONTROL[sl], a_floor)) for sl in CONTROL)
    tot = sum(len(anchors(s, a_floor)) for s in (1, 2, 5, 6))
    return 100.0 * tot / max(tot + fp, 1)


def main():
    if OUT_DIR.exists():
        raise FileExistsError(f"{OUT_DIR} already exists; refusing to overwrite")
    REF = reference_superset()
    CAND = candidate_superset(REF[3])

    base = evaluate(*ADOPTED, REF, CAND)
    expected = dict(n_lookalikes=630, n_anchors=863, oof_fp=13, loso_fp_L34=12, loso_fp_L321=18,
                    calls_all4=20873, calls_slice1=4267, calls_slice2=6272, calls_slice5=5039, calls_slice6=5295)
    bad = {k: (base[k], v) for k, v in expected.items() if base[k] != v}
    if bad:
        raise SystemExit(f"published numbers NOT reproduced (got, expected): {bad}")
    print("sanity check passed: published numbers reproduced exactly", flush=True)
    OUT_DIR.mkdir()

    grid = []
    for f in FLOORS:
        for m in MARGINS:
            grid.append(base if (f, m) == ADOPTED[:2] else evaluate(f, m, ADOPTED[2], REF, CAND))
            print(f"  grid floor {f:.2f} margin {m:.2f} done", flush=True)
    grid = pd.DataFrame(grid)
    grid.to_csv(OUT_DIR / "grid_candidate_floor_margin.csv", index=False)
    pick = select_grid(grid, "oof")
    pick_loso = select_grid(grid, "loso")
    f_star, m_star = float(pick["cand_floor"]), float(pick["margin"])

    rows = []
    for a in ANCHOR_FLOORS:
        r = base if (f_star, m_star, a) == ADOPTED else evaluate(f_star, m_star, a, REF, CAND)
        r = dict(r)
        r["anchor_purity_pct"] = anchor_purity(a)
        rows.append(r)
        print(f"  anchor floor {a:.2f} done", flush=True)
    anc = pd.DataFrame(rows)
    anc.to_csv(OUT_DIR / "sweep_anchor_floor.csv", index=False)
    elig = anc[anc["anchor_purity_pct"] >= PURITY_MIN]
    best = elig["oof_fp"].min()
    tied = elig[elig["oof_fp"] <= best + np.sqrt(max(best, 1))]
    a_star = float(tied["anchor_floor"].min())

    with pd.option_context("display.width", 250, "display.max_columns", 60):
        cols = ["cand_floor", "margin", "n_lookalikes", "oof_fp", "loso_fp", "calls_heldout",
                "oof_est_false_calls", "oof_fdp", "oof_fdp_hi95", "oof_est_true_calls", "loso_fdp", "calls_all4"]
        print("\n=== Grid (anchor floor 0.40) ===")
        print(grid[cols].round(4).to_string(index=False))
        print("\n=== Anchor floor at selected candidate setting ===")
        print(anc[["anchor_floor", "n_anchors", "anchor_purity_pct", "oof_fp", "loso_fp", "calls_heldout",
                   "oof_fdp", "calls_all4"]].round(4).to_string(index=False))
    summary = (f"Selected (out-of-fold rule): candidate floor {f_star}, margin {m_star}, anchor floor {a_star}\n"
               f"Leave-one-slice-out rule would select: floor {pick_loso['cand_floor']}, margin {pick_loso['margin']}\n"
               f"Adopted in thesis: {ADOPTED}\n")
    (OUT_DIR / "selection.txt").write_text(summary, encoding="utf-8")
    print("\n" + summary)


if __name__ == "__main__":
    main()
