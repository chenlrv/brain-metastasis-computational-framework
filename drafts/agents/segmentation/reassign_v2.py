"""Transcript reassignment test, version 2: every number of the thesis section from one run.

Same procedure as 04_reseg_reassign.py / 05_prior_sweep.py / 07_contamination_summary.py:
each transcript is scored against its vendor cell and the (up to) six nearest
cells whose centroid lies within the candidate radius,
    score = log p(gene | candidate profile) - distance / decay_length
            + home_advantage * [candidate is the vendor cell],
with profiles = log of pseudocount-smoothed per-cell gene proportions,
re-estimated from the current assignment in each of two rounds, and centroids
taken from the vendor assignment. Extracellular transcripts have no vendor cell
and join their best-scoring candidate if one lies within the radius.

Differences from the earlier scripts, all corrections:
  * genes are the 958 panel genes only; the earlier filter (prefixes negprb/
    falsecode/blank) did not match this panel's Negative*/SystemControl*
    probes, which therefore entered the profiles;
  * tumor cells are excluded with the final Stage-3 calls
    (with_tumor_prediction_final/);
  * the fraction of transcripts reassigned is measured against the vendor
    assignment (the earlier scripts compared round 2 with round 1);
  * all distances are set in micrometers and converted with the vendor pixel
    size (0.12028 um/px); the earlier radius of 60 px is 7.2 um, not 14.4 um.

Outputs (only to agents/outputs/segmentation_reassign_v2/, which must not exist):
  main_per_fov.csv      canonical configuration, per FOV, before/after read-outs
  main_pooled.json      canonical configuration, pooled over the five FOVs
  foreign_fraction.json per-cell fraction of vendor-assigned transcripts moved away
  sweep_fov514.csv      one-parameter-at-a-time sweeps on FOV 514 (Figure 20b)

Run: conda run -n thesis_research python agents/segmentation/reassign_v2.py
"""
import json
import pathlib
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = pathlib.Path(r"D:/thesis-research")
sys.path.insert(0, str(ROOT / "agents/segmentation"))
from seg_metrics import compute_all  # noqa: E402

TXDIR = ROOT / "agents/outputs/segmentation/tx_by_fov"
CACHE = ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad"
OUT = ROOT / "agents/outputs/segmentation_reassign_v2"
FOVS = [451, 512, 514, 515, 523]
SWEEP_FOV = 514

PX_UM = 0.12028        # vendor pixel size (median sqrt(Area.um2 / Area))
K_NEIGH = 6
PSEUDO = 1e-4
N_ROUNDS = 2
CANON = {"radius_um": 14.4, "decay_um": 3.6, "home_adv": 0.5}   # home_adv in nats: exp(0.5) = 1.65-fold
SWEEPS = {
    "radius_um": [7.2, 14.4, 24.1],
    "home_adv": [np.inf, 3.0, 1.5, 0.5, 0.0],            # inf = vendor cell always kept
    "decay_um": [1.8, 3.6, 7.2, 14.4, np.inf],           # inf = no distance penalty
}


def load_fov(fov, gene_index):
    df = pd.read_csv(TXDIR / f"tx_fov{fov}.csv",
                     usecols=["cell_ID", "cell", "x_global_px", "y_global_px", "target"])
    df["gene_idx"] = df["target"].map(gene_index).fillna(-1).astype(int)
    df = df[df["gene_idx"] >= 0].reset_index(drop=True)            # panel genes only
    intra = df["cell_ID"].to_numpy() != 0
    cell_ids = np.array(sorted(df.loc[intra, "cell"].unique()))
    home = np.full(len(df), -1)
    home[intra] = pd.Index(cell_ids).get_indexer(df.loc[intra, "cell"])
    cent = (df[intra].groupby("cell")[["x_global_px", "y_global_px"]].mean()
            .loc[cell_ids].to_numpy())
    return {"cells": cell_ids, "home": home, "gene": df["gene_idx"].to_numpy(),
            "xy": df[["x_global_px", "y_global_px"]].to_numpy(), "tree": cKDTree(cent)}


def counts(assign, gene, n_cells, n_genes):
    m = np.zeros((n_cells, n_genes))
    ok = assign >= 0
    np.add.at(m, (assign[ok], gene[ok]), 1.0)
    return m


def reassign(d, n_genes, radius_um, decay_um, home_adv):
    n_cells, home, gene = len(d["cells"]), d["home"], d["gene"]
    dist, idx = d["tree"].query(d["xy"], k=K_NEIGH, distance_upper_bound=radius_um / PX_UM)
    valid = idx < n_cells
    dist_um = dist * PX_UM
    cur = home.copy()
    for _ in range(N_ROUNDS):
        m = counts(cur, gene, n_cells, n_genes)
        tot = m.sum(1, keepdims=True); tot[tot == 0] = 1
        logp = np.log((m + PSEUDO) / (tot + PSEUDO * n_genes))
        best_s = np.full(len(gene), -np.inf); best_c = cur.copy()
        for k in range(K_NEIGH):
            vm = valid[:, k]
            cc = idx[vm, k]
            s = logp[cc, gene[vm]]
            if np.isfinite(decay_um):
                s = s - dist_um[vm, k] / decay_um
            is_home = (cc == home[vm]) & (home[vm] >= 0)
            s = s + np.where(is_home, home_adv, 0.0)
            rows = np.where(vm)[0]
            better = s > best_s[vm]
            best_s[rows[better]] = s[better]
            best_c[rows[better]] = cc[better]
        cur = best_c
    return cur


def readouts(assign, d, n_genes, genes, tumor):
    m = counts(assign, d["gene"], len(d["cells"]), n_genes)
    keep = np.array([c not in tumor for c in d["cells"]]) & (m.sum(1) > 0)
    return m[keep], compute_all(m[keep], genes)


def brief(r):
    p = r["purity"]
    return {"n_cells": r["n_cells"], "median_counts": r["median_total_counts"],
            "purity": p["purity_own_frac"], "n_purity_cells": p["n_valid"],
            "gfp_tdtomato_r": r["F1_full_r"], "gfp_cx3cr1_r": r["F2_full_r"],
            "lyve1_pct_pos": r["F3_lyve_pct_pos"], "lyve1_pct_lacking_bam": r["F3_lyve_pct_lacking_bam"]}


def move_stats(home, after, gene_n):
    intra = home >= 0
    return {"pct_changed": 100 * float((after != home).mean()),
            "pct_between_cells": 100 * float((intra & (after != home)).mean()),
            "pct_extracellular_captured": 100 * float((~intra & (after >= 0)).mean())}


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} already exists; refusing to overwrite")
    OUT.mkdir(parents=True)

    a = ad.read_h5ad(CACHE, backed="r")
    genes = list(a.var_names)                                         # 958 panel genes
    tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    a.file.close()
    gene_index = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)

    rows, pooled, foreign = [], {"before": [], "after": []}, []
    for fov in FOVS:
        t0 = time.time()
        d = load_fov(fov, gene_index)
        after = reassign(d, n_genes, **CANON)
        mb, rb = readouts(d["home"], d, n_genes, genes, tumor)
        ma, ra = readouts(after, d, n_genes, genes, tumor)
        pooled["before"].append(mb); pooled["after"].append(ma)
        # per non-tumor cell: share of its vendor-assigned transcripts moved elsewhere
        intra = d["home"] >= 0
        n_home = np.bincount(d["home"][intra], minlength=len(d["cells"]))
        n_away = np.bincount(d["home"][intra & (after != d["home"])], minlength=len(d["cells"]))
        nt = np.array([c not in tumor for c in d["cells"]]) & (n_home > 0)
        foreign.append(n_away[nt] / n_home[nt])
        rows.append({"fov": fov, "n_transcripts": len(d["gene"]), **move_stats(d["home"], after, n_genes),
                     **{f"before_{k}": v for k, v in brief(rb).items()},
                     **{f"after_{k}": v for k, v in brief(ra).items()}})
        r = rows[-1]
        print(f"FOV {fov}: {r['pct_changed']:.2f}% of transcripts changed cell "
              f"({r['pct_between_cells']:.2f}% between cells, {r['pct_extracellular_captured']:.2f}% "
              f"extracellular captured); purity {r['before_purity']:.3f} -> {r['after_purity']:.3f}; "
              f"median counts {r['before_median_counts']:.0f} -> {r['after_median_counts']:.0f} "
              f"({time.time() - t0:.0f}s)", flush=True)
    per_fov = pd.DataFrame(rows)
    per_fov.to_csv(OUT / "main_per_fov.csv", index=False)

    pooled_res = {s: brief(compute_all(np.vstack(m), genes)) for s, m in pooled.items()}
    pooled_res["config"] = {**CANON, "home_adv_fold": float(np.exp(CANON["home_adv"])),
                            "k_neighbours": K_NEIGH, "rounds": N_ROUNDS, "px_um": PX_UM}
    json.dump(pooled_res, open(OUT / "main_pooled.json", "w"), indent=2)
    fr = np.concatenate(foreign)
    json.dump({"n_cells": int(len(fr)), "median": float(np.median(fr)), "mean": float(fr.mean()),
               "p90": float(np.percentile(fr, 90))}, open(OUT / "foreign_fraction.json", "w"), indent=2)

    print("\nPOOLED (non-tumor cells, five FOVs)")
    for s in ("before", "after"):
        p = pooled_res[s]
        print(f"  {s:6s}: purity {p['purity']:.3f} over {p['n_purity_cells']:,} of {p['n_cells']:,} cells; "
              f"GFP-tdTomato r {p['gfp_tdtomato_r']:.3f}; GFP-Cx3cr1 r {p['gfp_cx3cr1_r']:.3f}; "
              f"Lyve1+ {p['lyve1_pct_pos']:.1f}% ({p['lyve1_pct_lacking_bam']:.0f}% lacking Mrc1/Cd163)")
    print(f"  mean of per-FOV purity: {per_fov.before_purity.mean():.3f} -> {per_fov.after_purity.mean():.3f}; "
          f"rose in {(per_fov.after_purity > per_fov.before_purity).sum()} of 5 FOVs")
    print(f"  per-cell share of transcripts moved away: median {np.median(fr) * 100:.1f}% ({len(fr):,} cells)")

    # one-parameter-at-a-time sweeps on FOV 514, others held at the canonical values
    d = load_fov(SWEEP_FOV, gene_index)
    _, rb = readouts(d["home"], d, n_genes, genes, tumor)
    base = rb["purity"]["purity_own_frac"]
    sweep = []
    for param, values in SWEEPS.items():
        for v in values:
            cfg = {**CANON, param: v}
            after = reassign(d, n_genes, **cfg)
            _, ra = readouts(after, d, n_genes, genes, tumor)
            p = ra["purity"]["purity_own_frac"]
            sweep.append({"param": param, "value": v, **cfg, "purity_before": base, "purity_after": p,
                          "delta_purity": p - base, **move_stats(d["home"], after, n_genes)})
            print(f"  sweep {param}={v}: purity {base:.3f} -> {p:.3f} (delta {p - base:+.3f}), "
                  f"{sweep[-1]['pct_changed']:.2f}% changed", flush=True)
    pd.DataFrame(sweep).to_csv(OUT / "sweep_fov514.csv", index=False)
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    main()
