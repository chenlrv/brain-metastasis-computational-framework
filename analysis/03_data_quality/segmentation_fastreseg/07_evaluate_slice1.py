"""Read-outs before and after FastReseg for ALL slice-1 FOVs, QC-passed non-tumor cells.

Same as slice1_03_evaluate.py, but restricted to the cells that passed the
low-transcript QC (the 120,689 QC-passed non-tumor cells used elsewhere in the
thesis), so the read-outs are comparable with Figures 17 and 18. A cell created
by FastReseg is kept if its parent cell passed QC.

"Before" = vendor segmentation: the intracellular transcripts of the input files
(05_prepare_slice1_inputs.py). "After" = FastReseg's updated tables
(06_run_fastreseg_slice1.R; transcripts it trimmed are absent). Tumor cells are
excluded with the final Stage-3 calls; a cell created by FastReseg inherits its
parent's status. A cell is positive for a gene if >= 1 transcript is assigned.

Reports, pooled over the slice:
  * Lyve1 positivity and the share of Lyve1-positive cells lacking Mrc1 and Cd163;
  * tdTomato positivity;
  * transcripts trimmed / moved, and the genes of the trimmed transcripts;
  * FastReseg's own cell typing of the non-tumor cells, and how many
    tdTomato-positive cells fall in the group(s) where a tdTomato transcript
    scores below -2 (needs figures/score_matrix.csv from the five-FOV run, whose
    reference is identical).
Writes only to agents/outputs/segmentation_fastreseg_slice1/eval/ (must not exist).

Run: conda run -n thesis_research python agents/segmentation/fastreseg/slice1_03_evaluate.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import json
import pathlib

import anndata as ad
import pandas as pd

ROOT = pathlib.Path(PROJECT_ROOT_STR)
B = ROOT / "agents/outputs/segmentation_fastreseg_slice1"
SCORES = ROOT / "agents/outputs/segmentation_fastreseg/figures/score_matrix.csv"
OUT = B / "eval_qcpassed"
GENES = ["Lyve1", "Mrc1", "Cd163", "tdTomato"]
QC_PASSED = set()


def positives(cell, target, tumor):
    """Per non-tumor cell: which of GENES it carries (>= 1 transcript)."""
    df = pd.DataFrame({"cell": cell.to_numpy(), "target": target.to_numpy()})
    parent = df.cell.str.rsplit("_g", n=1).str[0]
    df = df[~parent.isin(tumor) & parent.isin(QC_PASSED)]
    cells = df.cell.unique()
    pos = df[df.target.isin(GENES)].drop_duplicates().assign(v=True) \
        .pivot(index="cell", columns="target", values="v")
    pos = pos.reindex(index=cells, columns=GENES).fillna(False).astype(bool)
    return pos


def summarize(pos):
    lyve = pos.Lyve1
    return {"n_cells": len(pos), "lyve1_pct": 100 * lyve.mean(),
            "lyve1_lacking_bam_pct": 100 * (lyve & ~pos.Mrc1 & ~pos.Cd163).sum() / max(1, lyve.sum()),
            "tdtomato_pct": 100 * pos.tdTomato.mean()}


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} already exists; refusing to overwrite")
    OUT.mkdir(parents=True)
    a = ad.read_h5ad(ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad", backed="r")
    tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
    QC_PASSED.update(a.obs["cell"].astype(str))   # the cache holds only QC-passed cells
    a.file.close()
    info = pd.read_csv(B / "inputs" / "transDF_fileInfo.csv")
    low_groups = set(pd.read_csv(SCORES, index_col=0).loc["tdTomato"].pipe(lambda s: s[s < -2].index))

    before, after, typing, trimmed_genes, per_fov = [], [], [], [], []
    n_intra = n_trim = n_moved = 0
    for k, fov in enumerate(info.fov, start=1):
        orig = pd.read_csv(B / "inputs" / f"tx_fov{fov}.csv", usecols=["transcript_id", "cell", "cell_ID", "target"])
        orig = orig[orig.cell_ID != 0]
        upd = pd.read_csv(B / "fastreseg_out" / f"{k}_updated_transDF.csv",
                          usecols=["UMI_transID", "target", "UMI_cellID", "updated_cellID", "tLLR_maxCellType"])
        if not upd.UMI_transID.str.startswith(f"t_{fov}_").all():
            raise ValueError(f"output file {k} is not FOV {fov}")
        kept = upd.dropna(subset=["updated_cellID"])
        trim = orig[~orig.transcript_id.isin(set(kept.UMI_transID))]
        moved = int((kept.UMI_cellID != kept.updated_cellID).sum())
        n_intra += len(orig); n_trim += len(trim); n_moved += moved
        trimmed_genes.append(trim.target)
        before.append(positives(orig.cell, orig.target, tumor))
        after.append(positives(kept.updated_cellID, kept.target, tumor))
        typing.append(upd.drop_duplicates("UMI_cellID")[["UMI_cellID", "tLLR_maxCellType"]])
        per_fov.append({"fov": fov, "n_intracellular": len(orig), "n_trimmed": len(trim), "n_moved": moved})
        if k % 25 == 0:
            print(f"  {k}/{len(info)} FOVs", flush=True)

    pb, pa = pd.concat(before), pd.concat(after)
    typ = pd.concat(typing).rename(columns={"UMI_cellID": "cell", "tLLR_maxCellType": "group"})
    typ = typ[~typ.cell.isin(tumor) & typ.cell.isin(QC_PASSED)].set_index("cell")["group"]
    tdt_cells = pb.index[pb.tdTomato]
    tdt_in_low = int(typ.reindex(tdt_cells).isin(low_groups).sum())
    res = {"n_fovs": len(info), "n_intracellular": n_intra, "n_trimmed": n_trim, "n_moved": n_moved,
           "pct_trimmed": 100 * n_trim / n_intra,
           "max_pct_trimmed_per_fov": max(100 * r["n_trimmed"] / r["n_intracellular"] for r in per_fov),
           "before": summarize(pb), "after": summarize(pa),
           "trimmed_top_genes": pd.concat(trimmed_genes).value_counts().head(10).to_dict(),
           "trimmed_anomaly_genes": {g: int((pd.concat(trimmed_genes) == g).sum()) for g in GENES},
           "tdtomato_low_groups": sorted(low_groups),
           "group_share_nontumor_pct": (typ.value_counts(normalize=True) * 100).round(2).to_dict(),
           "tdtomato_pos_cells": int(len(tdt_cells)), "tdtomato_pos_in_low_groups": tdt_in_low}
    pd.DataFrame(per_fov).to_csv(OUT / "per_fov.csv", index=False)
    typ.to_csv(OUT / "cell_groups_nontumor.csv")
    json.dump(res, open(OUT / "summary.json", "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
