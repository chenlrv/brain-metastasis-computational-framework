"""The three sets derived from the SingleR output, one figure per set.

  fig_set_candidates.png   tumor candidates, all six slices
                           (score_tumor > 0.2, margin > 0.08, best-class)
  fig_set_anchors.png      tumor anchors, slices 1 and 5 -- the same filters at
                           score_tumor >= 0.4, shown against the candidate pool
                           they are drawn from (positive training class)
  fig_set_lookalikes.png   healthy look-alikes, the two sham-injected control
                           slices, where every selected cell is a false positive
                           by construction (negative training class)

The set definitions are the ones used by the pipeline
(identify_tumor_cells.py: _get_tumor_candidates_ids, _get_tumor_ref_ids,
_get_healthy_ref_ids); they are re-applied here per slice so every section can
be drawn on one common scale.

Saves the three PNGs (200 dpi) and fig_singler_sets.csv with the counts.

Run: conda run -n thesis_research python analysis/02_tumor_identification/fig04_07_singler_sets.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = PROJECT_ROOT_STR
ANNOT = ROOT + "/outputs/cell_annotation/{slide}/05/{s}/slice_{s}_final_cell_annotations_refined_tabula_brain_tumor_avinoam.csv"
WTP = ROOT + "/resources/cache/with_tumor_prediction/slice_{}_adata.h5ad"
OUT = ROOT + "/thesis_plots/"

SLICES = {1: "L321", 2: "L321", 3: "L321", 4: "L34", 5: "L34", 6: "L34"}
MOUSE = {1: 2, 2: 3, 3: 1, 4: 1, 5: 3, 6: 2}     # Table 1
CONTROL = {3, 4}
FRAME_TUMOR, FRAME_CONTROL = "#C77400", "#2166AC"
ANCHOR_SLICES = [1, 5]          # one tumor-bearing slice per slide
CAND_FLOOR, ANCHOR_FLOOR, DELTA = 0.2, 0.4, 0.08
PX_UM = 0.12028

GREY, RAW, CAND, CAND_BG, ANCHOR, LOOKALIKE = (
    "#DCDEE0", "#D6604D", "#D6604D", "#F4C3B8", "#A50F15", "#2166AC")


def singler_sets(slice_id, slide_id):
    """Raw-label, candidate and anchor barcodes for one slice.

    raw    : every cell whose SingleR label is "Tumor", with no further filter
    cand   : the candidate filters (score_tumor > 0.2, margin > 0.08, best-class)
    anchor : the same filters at the stricter floor score_tumor >= 0.4
    """
    df = pd.read_csv(ANNOT.format(slide=slide_id, s=slice_id))
    t = df[df["predicted_cell_type"] == "Tumor"].copy()
    raw = set(t["cell_barcode"].astype(str))
    t["next_best"] = t[["score_brain_struct", "score_brain_immune"]].max(axis=1)
    t["delta"] = (t["score_tumor"] - t["next_best"]).abs()
    passes = t[(t["delta"] > DELTA) & (t["score_tumor"] > t["next_best"])]
    cand = set(passes.loc[passes["score_tumor"] > CAND_FLOOR, "cell_barcode"].astype(str))
    anchor = set(passes.loc[passes["score_tumor"] >= ANCHOR_FLOOR, "cell_barcode"].astype(str))
    return raw, cand, anchor


def load_all():
    data, rows = {}, []
    for s, slide in SLICES.items():
        raw, cand, anchor = singler_sets(s, slide)
        if s not in ANCHOR_SLICES:
            anchor = set()                      # anchors come from slices 1 and 5 only
        obs = ad.read_h5ad(WTP.format(s), backed="r").obs
        bc = obs.index.astype(str).to_numpy()
        x = obs["CenterX_global_px"].to_numpy() * PX_UM / 1000.0
        y = obs["CenterY_global_px"].to_numpy() * PX_UM / 1000.0
        data[s] = dict(x=x, y=y,
                       raw=np.isin(bc, list(raw)),
                       cand=np.isin(bc, list(cand)),
                       anchor=np.isin(bc, list(anchor)))
        rows.append(dict(slice=s, slide=slide,
                         type="control" if s in CONTROL else "tumor-bearing",
                         n_cells=len(bc), n_raw_tumor_label=len(raw),
                         n_candidates=len(cand), n_anchors=len(anchor),
                         n_lookalikes=len(cand) if s in CONTROL else 0))
        print(f"slice {s}: {len(raw):,} raw 'Tumor' labels, {len(cand):,} candidates, "
              f"{len(anchor):,} anchors of {len(bc):,} cells")
    return data, pd.DataFrame(rows)


def panel_span(data):
    """One shared data span, so every panel of every figure is at the same scale."""
    pad = 0.35
    span_x = max(d["x"].max() - d["x"].min() for d in data.values()) + 2 * pad
    span_y = max(d["y"].max() - d["y"].min() for d in data.values()) + 2 * pad
    return span_x, span_y


def draw(data, slices, layers, legend, out_name, span, mode, ncols=2, scalebar=True):
    """layers: list of (mask_key or None, color, size) drawn back to front."""
    span_x, span_y = span
    nrows = int(np.ceil(len(slices) / ncols))
    fig_w = 9.6 if ncols == 2 else 5.2
    panel_w = fig_w / ncols - 0.25
    fig_h = nrows * (panel_w * span_y / span_x + 0.52) + 0.6   # + two-line title strip
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), dpi=200,
                             squeeze=False)

    for ax, s in zip(axes.ravel(), slices):
        d = data[s]
        ax.scatter(d["x"], d["y"], s=0.14, c=GREY, linewidths=0, rasterized=True)
        for key, color, size in layers:
            m = d[key]
            if m.any():
                ax.scatter(d["x"][m], d["y"][m], s=size, c=color, linewidths=0,
                           rasterized=True)
        cx = (d["x"].min() + d["x"].max()) / 2
        cy = (d["y"].min() + d["y"].max()) / 2
        ax.set_xlim(cx - span_x / 2, cx + span_x / 2)
        ax.set_ylim(cy - span_y / 2, cy + span_y / 2)
        ax.set_aspect("equal")
        ax.set_anchor("N")
        ax.set_xticks([]); ax.set_yticks([])
        frame = FRAME_CONTROL if s in CONTROL else FRAME_TUMOR
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color(frame)
            sp.set_linewidth(1.6)
        ax.set_title(title_for(s, d, mode), fontsize=9, pad=4, color=frame,
                     linespacing=1.3)

    for ax in axes.ravel()[len(slices):]:
        ax.set_visible(False)

    if scalebar:
        ax0 = axes[0, 0]
        x1 = ax0.get_xlim()[1]; y0 = ax0.get_ylim()[0]
        ax0.plot([x1 - 1.5, x1 - 0.5], [y0 + 0.35, y0 + 0.35], "-", color="#333", lw=2.2)
        ax0.text(x1 - 1.0, y0 + 0.48, "1 mm", ha="center", va="bottom", fontsize=8)

    fig.legend(handles=legend, fontsize=9, frameon=False, ncol=min(len(legend), 3),
               loc="lower center", bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(h_pad=2.0, w_pad=2.0, rect=(0, 0.05, 1, 1))
    fig.savefig(OUT + out_name, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("Saved:", OUT + out_name)


def title_for(s, d, mode):
    kind = "sham-injected control" if s in CONTROL else "tumor-bearing"
    who = f"Slice {s} ({SLICES[s]}, mouse {MOUSE[s]}, {kind})"
    n_cand = int(d["cand"].sum())
    n_anchor = int(d["anchor"].sum())
    if mode == "raw":
        n_raw = int(d["raw"].sum())
        n_all = len(d["raw"])
        return (f"{who}\n{n_raw:,} labeled Tumor of {n_all:,} cells "
                f"({100 * n_raw / n_all:.1f}%)")
    n_all = len(d["cand"])
    if mode == "lookalikes":
        what = (f"{n_cand:,} look-alikes of {n_all:,} cells "
                f"({100 * n_cand / n_all:.1f}%)")
    elif mode == "anchors":
        what = (f"{n_anchor:,} anchors of {n_cand:,} candidates "
                f"({100 * n_anchor / n_cand:.1f}%)")
    else:
        what = (f"{n_cand:,} candidates of {n_all:,} cells "
                f"({100 * n_cand / n_all:.1f}%)")
    return f"{who}\n{what}"


def frame_handle(color, label):
    """Legend entry for the panel border, drawn as an empty box in its color."""
    return Patch(facecolor="none", edgecolor=color, linewidth=1.6, label=label)


def main():
    data, counts = load_all()
    counts.to_csv(OUT + "fig_singler_sets.csv", index=False)
    span = panel_span(data)

    # 0. the raw SingleR labels, before any filtering
    draw(data, list(SLICES), [("raw", RAW, 0.5)],
         [Patch(color=RAW, label="cell labeled Tumor by SingleR"),
          Patch(color=GREY, label="all other cells"),
          frame_handle(FRAME_TUMOR, "tumor-bearing section"),
          frame_handle(FRAME_CONTROL, "sham-injected control section")],
         "fig_set_raw_labels.png", span, mode="raw")

    # 1. candidates, all six slices
    draw(data, list(SLICES), [("cand", CAND, 0.6)],
         [Patch(color=CAND, label="tumor candidate"),
          Patch(color=GREY, label="all other cells"),
          frame_handle(FRAME_TUMOR, "tumor-bearing section"),
          frame_handle(FRAME_CONTROL, "sham-injected control section")],
         "fig_set_candidates.png", span, mode="candidates")

    # 2. anchors, against the candidate pool they are drawn from
    draw(data, ANCHOR_SLICES,
         [("cand", CAND_BG, 0.5), ("anchor", ANCHOR, 1.1)],
         [Patch(color=ANCHOR, label="tumor anchor (positive training class)"),
          Patch(color=CAND_BG, label="other tumor candidates"),
          Patch(color=GREY, label="all other cells"),
          frame_handle(FRAME_TUMOR, "tumor-bearing section")],
         "fig_set_anchors.png", span, mode="anchors")

    # 3. healthy look-alikes, the two control slices
    draw(data, sorted(CONTROL), [("cand", LOOKALIKE, 1.1)],
         [Patch(color=LOOKALIKE, label="healthy look-alike (negative training class)"),
          Patch(color=GREY, label="all other cells"),
          frame_handle(FRAME_CONTROL, "sham-injected control section")],
         "fig_set_lookalikes.png", span, mode="lookalikes")

    print()
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()

