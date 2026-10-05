"""Two candidate figures on what FastReseg could recognize as misplaced (choose one).

Both use FastReseg's own transcript scores (score_matrix.csv, exported by
04_export_scores.R: ln of a gene's abundance in a cell type relative to its
highest abundance in any type) and FastReseg's own cell typing of the cells in
the five tested FOVs (tLLR_maxCellType in the updated transcript tables).
A transcript is treated as poorly fitting, i.e. recognizable as misplaced, only
if its score in its cell's type is below -2 (svmClass_score_cutoff default).

  A  heatmap: scores of the anomaly genes (Lyve1, tdTomato), the BAM markers
     (Mrc1, Cd163) and the three most-trimmed genes (Ttr, Apod, Ptgds) across the
     12 cell types; cells below -2 outlined; top bar = share of non-tumor cells
     of each type in the five FOVs.
  B  coverage: for every panel gene, the percentage of non-tumor cells in the
     five FOVs whose cell type gives that gene a score below -2, genes ranked.

Writes only new files to agents/outputs/segmentation_fastreseg/figures/.

Run: conda run -n thesis_research python agents/segmentation/fastreseg/05_figures.py
"""
import pathlib

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

ROOT = pathlib.Path("D:/thesis-research")
B = ROOT / "agents/outputs/segmentation_fastreseg"
OUT = B / "figures"
CUT = -2.0
GENES_A = ["Lyve1", "tdTomato", "Mrc1", "Cd163", "Ttr", "Apod", "Ptgds"]
LABELS_B = {"Lyve1": "#D62728", "tdTomato": "#D62728", "Mrc1": "#E07B00", "Cd163": "#E07B00",
            "Ttr": "#0072B2", "Apod": "#0072B2", "Ptgds": "#0072B2"}

for f in ("fig_A_reference_heatmap.png", "fig_B_gene_coverage.png", "cell_types_5fov.csv", "gene_coverage.csv"):
    if (OUT / f).exists():
        raise FileExistsError(f"{OUT / f} exists; refusing to overwrite")

score = pd.read_csv(OUT / "score_matrix.csv", index_col=0)

# FastReseg's cell typing of the non-tumor cells of the five FOVs
a = ad.read_h5ad(ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad", backed="r")
tumor = set(a.obs.loc[a.obs["pred_tumor_XGBoost"].to_numpy(dtype=bool), "cell"].astype(str))
a.file.close()
typed = []
for k in range(1, 6):
    t = pd.read_csv(B / "fastreseg_out" / f"{k}_updated_transDF.csv", usecols=["UMI_cellID", "tLLR_maxCellType"])
    typed.append(t.drop_duplicates("UMI_cellID"))
typed = pd.concat(typed).rename(columns={"UMI_cellID": "cell", "tLLR_maxCellType": "type"})
typed["tumor"] = typed.cell.isin(tumor)
typed.to_csv(OUT / "cell_types_5fov.csv", index=False)
share = typed.loc[~typed.tumor, "type"].value_counts(normalize=True).reindex(score.columns).fillna(0)

# per-gene coverage: share of non-tumor cells whose type scores the gene below the cutoff
coverage = (score < CUT).astype(float).mul(share, axis=1).sum(axis=1) * 100
pd.DataFrame({"min_score": score.min(axis=1), "fold_max_min": np.exp(-score.min(axis=1)),
              "n_types_below_cutoff": (score < CUT).sum(axis=1), "pct_cells_recognizable": coverage}) \
    .sort_values("pct_cells_recognizable").to_csv(OUT / "gene_coverage.csv")

# ---- A: heatmap -------------------------------------------------------------------
cols = list(share.sort_values(ascending=False).index)
M = score.loc[GENES_A, cols].to_numpy()
fig = plt.figure(figsize=(8.6, 4.4), dpi=200)
gs = fig.add_gridspec(2, 2, height_ratios=[0.9, 4], width_ratios=[30, 1], hspace=0.06, wspace=0.04)
axbar = fig.add_subplot(gs[0, 0])
axbar.bar(range(len(cols)), share[cols] * 100, color="#9A9A9A", width=0.8)
axbar.set_xlim(-0.5, len(cols) - 0.5)
axbar.set_xticks([])
axbar.set_ylabel("% of cells", fontsize=8)
axbar.tick_params(labelsize=7)
axbar.spines[["top", "right"]].set_visible(False)
ax = fig.add_subplot(gs[1, 0])
im = ax.imshow(M, cmap="viridis", vmin=-6.5, vmax=0, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if M[i, j] < CUT:
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="white", lw=1.8))
        ax.text(j, i, f"{M[i, j]:.1f}", ha="center", va="center", fontsize=6.5,
                color="white" if M[i, j] < -3.5 else "black")
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, fontsize=8)
ax.set_yticks(range(len(GENES_A)))
ax.set_yticklabels(GENES_A, fontsize=9)
ax.set_xlabel("cell type in the FastReseg reference (ordered by share of non-tumor cells in the five FOVs)",
              fontsize=8)
ax.axhline(3.5, color="white", lw=2.5)
cax = fig.add_subplot(gs[1, 1])
cb = fig.colorbar(im, cax=cax)
cb.set_label("transcript score (ln of abundance\nrelative to the type of highest abundance)", fontsize=7)
cb.ax.axhline(CUT, color="white", lw=1.5)
cb.ax.tick_params(labelsize=7)
fig.savefig(OUT / "fig_A_reference_heatmap.png", bbox_inches="tight", dpi=200)
plt.close(fig)

# ---- B: coverage per gene -----------------------------------------------------------
cov = coverage.sort_values()
fig, ax = plt.subplots(figsize=(7.6, 3.8), dpi=200)
x = np.arange(len(cov))
ax.plot(x, cov.to_numpy(), color="#333", lw=1.4)
n_zero = int((cov == 0).sum())
ax.axvspan(-0.5, n_zero - 0.5, color="#DDDDDD", zorder=0)
ax.text(n_zero / 2, 55, f"{n_zero} of {len(cov)} genes ({100 * n_zero / len(cov):.0f}%):\n"
        "never recognizable as misplaced", ha="center", fontsize=8, color="#444")
for g, color in LABELS_B.items():
    xi = int(np.where(cov.index == g)[0][0]); yi = cov[g]
    ax.plot(xi, yi, "o", color=color, ms=5, zorder=4)
    ax.annotate(f"{g} ({yi:.0f}%)", xy=(xi, yi), xytext=(xi - 40 if xi > 900 else xi + 15, yi + (8 if yi < 50 else -10)),
                fontsize=7.5, color=color, arrowprops=dict(arrowstyle="-", color=color, lw=0.6))
ax.set_xlim(-5, len(cov) + 5)
ax.set_ylim(-3, 103)
ax.set_xlabel("panel genes, ranked", fontsize=9)
ax.set_ylabel("non-tumor cells in which a transcript of\nthe gene could be recognized as misplaced (%)", fontsize=8.5)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "fig_B_gene_coverage.png", bbox_inches="tight", dpi=200)
plt.close(fig)

print("type shares (non-tumor, 5 FOVs):", (share * 100).round(1).to_dict())
print(f"genes never recognizable: {n_zero} of {len(cov)}")
print(cov.loc[list(LABELS_B)].round(1).to_dict())
print("written to", OUT)
