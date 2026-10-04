"""Refined heatmap figure (A2): what FastReseg's reference lets it recognize as misplaced.

Rows: the anomaly genes (Lyve1, tdTomato), the canonical BAM markers (Mrc1, Cd163)
and, for comparison, Ttr, the gene FastReseg trimmed most often. Columns: the 12
cell types of the reference, ordered by their share of the non-tumor cells of the
five FOVs (FastReseg's own cell typing; top bar). Colour: FastReseg's transcript
score, ln of the gene's abundance in a type relative to its highest abundance in
any type; outlined cells score below -2, where a transcript would be recognized as
misplaced. Inputs from 04_export_reference_scores.R and 05_figures.py (score_matrix.csv,
cell_types_5fov.csv). Writes one new file; refuses to overwrite.

Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/fig20_heatmap_base.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import pathlib

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

FIG = pathlib.Path(PROJECT_ROOT_STR + "/agents/outputs/segmentation_fastreseg/figures")
OUT = FIG / "fig_A7_reference_heatmap.png"   # A2-A6: earlier renders (A6: borders partly hidden)
EDGE = "#1F5BFF"   # border marking scores below the cutoff
CUT = -2.0
GENES = ["Lyve1", "tdTomato", "Mrc1", "Cd163", "Ttr"]
GROUPS = [("anomalies", 0, 1), ("BAM\nmarkers", 2, 3), ("comparison", 4, 4)]

if OUT.exists():
    raise FileExistsError(f"{OUT} exists; refusing to overwrite")
score = pd.read_csv(FIG / "score_matrix.csv", index_col=0)
types = pd.read_csv(FIG / "cell_types_5fov.csv")
share = types.loc[~types.tumor, "type"].value_counts(normalize=True).reindex(score.columns).fillna(0) * 100
cols = list(share.sort_values(ascending=False).index)
M = score.loc[GENES, cols].to_numpy()
M = np.where(np.abs(M) < 0.05, 0.0, M)              # avoid "-0.0" labels

fig = plt.figure(figsize=(8.8, 4.3), dpi=200)
gs = fig.add_gridspec(2, 2, height_ratios=[1, 4.2], width_ratios=[30, 1], hspace=0.08, wspace=0.04)
axb = fig.add_subplot(gs[0, 0])
axb.bar(range(len(cols)), share[cols], color="#9A9A9A", width=0.8)
for j, v in enumerate(share[cols]):
    if v >= 1:
        axb.text(j, v + 1, f"{v:.0f}", ha="center", fontsize=6.5, color="#555")
axb.set_xlim(-0.5, len(cols) - 0.5)
axb.set_ylim(0, share.max() * 1.25)
axb.set_xticks([])
axb.set_ylabel("% of cells", fontsize=7.5)
axb.tick_params(labelsize=7)
axb.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[1, 0])
im = ax.imshow(M, cmap="viridis", vmin=-6.5, vmax=0, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if M[i, j] < CUT:
            ax.add_patch(Rectangle((j - 0.42, i - 0.42), 0.84, 0.84, fill=False, ec=EDGE, lw=2.4, zorder=5))
        # two decimals near the cutoff, so -1.97 is not shown as an unoutlined "-2.0"
        txt = f"{M[i, j]:.2f}" if abs(M[i, j] - CUT) < 0.05 else f"{M[i, j]:.1f}"
        ax.text(j, i, txt, ha="center", va="center", fontsize=7,
                color="white" if M[i, j] < -3.5 else "black")
for y in (1.5, 3.5):
    ax.axhline(y, color="white", lw=3)
ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols, fontsize=8.5)
ax.set_yticks(range(len(GENES)))
ax.set_yticklabels(GENES, fontsize=9.5)
ax.set_xlabel("cell group (unsupervised InSituType cell typing from the vendor export)", fontsize=8.5)
for label, i0, i1 in GROUPS:
    ax.annotate("", xy=(-2.35, i0 - 0.4), xytext=(-2.35, i1 + 0.4), xycoords="data",
                arrowprops=dict(arrowstyle="-", lw=1, color="#555"), annotation_clip=False)
    ax.text(-2.5, (i0 + i1) / 2, label, ha="right", va="center", fontsize=7.5, color="#333")

cax = fig.add_subplot(gs[1, 1])
cb = fig.colorbar(im, cax=cax)
cb.set_label("transcript score\n(ln abundance relative to the\ngroup of highest abundance)", fontsize=7)
cb.ax.axhline(CUT, color=EDGE, lw=2)
cb.ax.tick_params(labelsize=7)
fig.text(0.5, -0.04, "blue border: score below −2 (≥ 7.4-fold below the maximum), where a transcript "
         "would be recognized as misplaced", ha="center", fontsize=7.5, color="#444")
fig.savefig(OUT, bbox_inches="tight", dpi=200)
plt.close(fig)
print("Saved:", OUT)
