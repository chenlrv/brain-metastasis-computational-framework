"""Supplementary figure: where the 12 FastReseg cell groups (a-l) lie in slice 1.

One panel per group, ordered as in Figure 20 (by share of the QC-passed non-tumor
cells); the group's cells are drawn in color and all other cells in gray. Group
assignments are FastReseg's own (eval_qcpassed/cell_groups_nontumor.csv from
07_evaluate_slice1.py), the same as the top bar of Figure 20; tumor cells
(final Stage-3 calls) are not assigned a group and are drawn in gray. Coordinates
are the global tissue coordinates (CenterX/Y_global_px).

Writes only automation/outputs/segmentation_fastreseg_slice1/figures/figS1_cell_groups_slice1.png
(must not exist).
Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/figS1_cell_groups.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import pathlib

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = pathlib.Path(PROJECT_ROOT_STR)
B = ROOT / "automation/outputs/segmentation_fastreseg_slice1"
OUT = B / "figures" / "figS1_cell_groups_slice1.png"
PX_UM = 0.12028
GREY, COLOR = "#DADADA", "#1F5BFF"


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} exists; refusing to overwrite")
    a = ad.read_h5ad(ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad", backed="r")
    obs = a.obs[["cell", "CenterX_global_px", "CenterY_global_px"]].copy()
    a.file.close()
    obs["cell"] = obs["cell"].astype(str)
    groups = pd.read_csv(B / "eval_qcpassed" / "cell_groups_nontumor.csv")
    obs = obs.merge(groups, on="cell", how="left")
    if obs["group"].notna().sum() != len(groups):
        raise ValueError("not every grouped cell was found in slice 1")
    x = obs["CenterX_global_px"].to_numpy() * PX_UM / 1000.0
    y = obs["CenterY_global_px"].to_numpy() * PX_UM / 1000.0
    share = groups["group"].value_counts(normalize=True) * 100
    order = list(share.sort_values(ascending=False).index)

    fig, axes = plt.subplots(3, 4, figsize=(13, 8.2))
    for k, (g, ax) in enumerate(zip(order, axes.ravel())):
        sel = (obs["group"] == g).to_numpy()
        ax.scatter(x[~sel], y[~sel], s=0.05, c=GREY, linewidths=0, rasterized=True)
        ax.scatter(x[sel], y[sel], s=0.25, c=COLOR, linewidths=0, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"group {g}  ({share[g]:.1f}% of cells, n = {int(sel.sum()):,})", fontsize=9.5, pad=3)
        if k == 0:
            x0, x1 = ax.get_xlim(); y0, _ = ax.get_ylim()
            ax.plot([x1 - 1.4, x1 - 0.4], [y0 + 0.25, y0 + 0.25], "-", color="#333", lw=2)
            ax.text(x1 - 0.9, y0 + 0.38, "1 mm", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"saved {OUT}")
    print("order and shares:", {g: round(share[g], 1) for g in order})


if __name__ == "__main__":
    main()
