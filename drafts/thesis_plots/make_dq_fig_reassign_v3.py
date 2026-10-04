"""Figure 20 (version 3) -- permissive transcript reassignment raises marker purity but leaves
the Lyve1 anomaly unchanged.

Values are read from the saved outputs of agents/segmentation/reassign_v3.py
(candidate radius 14.4 um, up to five candidates, decay length 6.3 um = median
cell radius, no current-cell advantage):
(a) myeloid marker purity before and after reassignment per FOV, with the mean;
(b) percentage of Lyve1-positive cells lacking both Mrc1 and Cd163, before and
    after, per FOV -- the probe anomaly the test was meant to explain;
(c) change in purity on FOV 514 under one-parameter-at-a-time sweeps.
Writes a new file only; earlier figure versions are untouched.

Run: conda run -n thesis_research python thesis_plots/make_dq_fig_reassign_v3.py
"""
import json
import pathlib

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES = pathlib.Path("D:/thesis-research/agents/outputs/segmentation_reassign_v3")
OUT = RES / "dq_fig_reassign_v3b.png"   # v3.png: first render, label layout since improved
GREY, BLUE, RED = "#9A9A9A", "#0072B2", "#D62728"

if OUT.exists():
    raise FileExistsError(f"{OUT} exists; refusing to overwrite")
per_fov = pd.read_csv(RES / "main_per_fov.csv")
sweep = pd.read_csv(RES / "sweep_fov514.csv")
canon = json.load(open(RES / "main_pooled.json"))["config"]


def spread(values, min_gap, n_iter=200):
    """Label positions: the values, pushed apart symmetrically until no two are closer
    than min_gap, so each label stays as near its point as possible."""
    pos = np.array(values, dtype=float)
    for _ in range(n_iter):
        order = np.argsort(pos)
        moved = False
        for a, b in zip(order[:-1], order[1:]):
            gap = pos[b] - pos[a]
            if gap < min_gap:
                shift = (min_gap - gap) / 2
                pos[a] -= shift
                pos[b] += shift
                moved = True
        if not moved:
            break
    return pos


def before_after(ax, before, after, ylabel, ylim, min_gap, panel):
    labels = [str(f) for f in per_fov.fov] + ["mean"]
    after_all = list(after) + [np.mean(after)]
    ypos = spread(after_all, min_gap)
    for (b, a), y, lab in zip(zip(before, after), ypos, labels):
        color = BLUE if a > b else RED
        ax.plot([0, 1], [b, a], "-o", color=color, ms=5, lw=1.6)
        ax.text(1.06, y, lab, va="center", fontsize=7.5, color=color)
    ax.plot([0, 1], [np.mean(before), np.mean(after)], "-o", color="#111", ms=7, lw=2.6, zorder=5)
    ax.text(1.06, ypos[-1], "mean", va="center", fontsize=8, fontweight="bold")
    ax.set_xlim(-0.25, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["before", "after"], fontsize=9.5)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_ylim(*ylim)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.02, 1.06, panel, transform=ax.transAxes, fontsize=13, fontweight="bold",
            va="top", ha="right")


fig = plt.figure(figsize=(16.5, 4.6), dpi=200)
gs = fig.add_gridspec(1, 5, width_ratios=[1.1, 1.1, 1, 1, 1], wspace=0.42)

before_after(fig.add_subplot(gs[0, 0]), per_fov.before_purity, per_fov.after_purity,
             "myeloid marker purity", (0.75, 0.88), 0.005, "(a)")
before_after(fig.add_subplot(gs[0, 1]), per_fov.before_lyve1_pct_lacking_bam,
             per_fov.after_lyve1_pct_lacking_bam,
             "Lyve1-positive cells lacking Mrc1 and Cd163 (%)", (0, 100), 2.6, "(b)")


def fmt_adv(v):
    if np.isinf(v):
        return "none\n(no move\nbetween cells)"
    if v == 0:
        return "1x\n(no adv.)"
    return {3.0: "20x", 1.5: "4.5x", 0.5: "1.65x"}.get(v, f"{np.exp(v):.3g}x")


panels = [
    ("radius_um", "candidate radius (µm)", lambda v: f"{v:g}", canon["radius_um"], False),
    ("home_adv", "current-cell advantage", fmt_adv, canon["home_adv"], True),
    ("decay_um", "spatial decay length (µm)", lambda v: "no\npenalty" if np.isinf(v) else f"{v:g}",
     canon["decay_um"], True),
]
for i, (param, xlabel, fmt, cval, categorical) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i + 2])
    s = sweep[sweep.param == param]
    vals = s["value"].to_numpy(dtype=float)
    xs = np.arange(len(s)) if categorical else vals
    ax.plot(xs, s.delta_purity, "-o", color="#333", ms=5, lw=1.6, zorder=3)
    ax.axhline(0, ls="--", lw=1.2, color=GREY, zorder=1)
    j = int(np.where(np.isclose(vals, cval))[0][0])
    ax.plot([xs[j]], [s.delta_purity.iloc[j]], "o", ms=10, mfc="none", mec=BLUE, mew=2, zorder=4)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_xticks(xs)
    ax.set_xticklabels([fmt(v) for v in vals], fontsize=8)
    ax.set_ylim(-0.012, 0.13)
    ax.spines[["top", "right"]].set_visible(False)
    if i == 0:
        ax.set_ylabel("change in marker purity after reassignment (FOV 514)", fontsize=9)
        ax.text(-0.03, 1.06, "(c)", transform=ax.transAxes, fontsize=13, fontweight="bold",
                va="top", ha="right")
        ax.text(0.97, -0.009, "no change", fontsize=7.5, color="#666", ha="right",
                transform=ax.get_yaxis_transform(), va="center")
    else:
        ax.tick_params(axis="y", labelsize=8)
    if param == "decay_um":
        last = s.delta_purity.iloc[-1]
        ax.annotate("no spatial\nconstraint", xy=(xs[-1], last), xytext=(xs[-1] - 1.2, last + 0.01),
                    fontsize=7.5, color=RED, ha="center",
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1))

fig.text(0.62, -0.09, "circled point: configuration used for the reported result",
         ha="center", fontsize=8, color="#444")
fig.savefig(OUT, bbox_inches="tight", dpi=200)
plt.close(fig)
print("Saved:", OUT)
