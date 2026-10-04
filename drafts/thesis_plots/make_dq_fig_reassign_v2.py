"""Figure 20 (version 2) -- transcript reassignment does not materially improve marker purity.

Same layout as make_dq_fig_reassign.py, but every value is read from the saved
outputs of agents/segmentation/reassign_v2.py instead of being typed in:
(a) marker purity before and after reassignment per FOV (main_per_fov.csv), with
    the mean of the five FOVs; (b) change in purity on FOV 514 under
    one-parameter-at-a-time sweeps (sweep_fov514.csv). Canonical configuration:
    candidate radius 14.4 um, decay length 3.6 um, current-cell advantage
    1.65-fold, two rounds. The previous version and its output are left untouched.

Run: conda run -n thesis_research python thesis_plots/make_dq_fig_reassign_v2.py
"""
import json
import pathlib

import matplotlib as mpl
import numpy as np
import pandas as pd

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES = pathlib.Path("D:/thesis-research/agents/outputs/segmentation_reassign_v2")
OUT = RES / "dq_fig_reassign_v2d.png"   # v2/v2b/v2c.png: earlier renders, labels since corrected
GREY, BLUE, RED = "#9A9A9A", "#0072B2", "#D62728"

if OUT.exists():
    raise FileExistsError(f"{OUT} exists; refusing to overwrite")
per_fov = pd.read_csv(RES / "main_per_fov.csv")
sweep = pd.read_csv(RES / "sweep_fov514.csv")
canon = json.load(open(RES / "main_pooled.json"))["config"]

fig = plt.figure(figsize=(13.4, 4.6), dpi=200)
gs = fig.add_gridspec(1, 4, width_ratios=[1.15, 1, 1, 1], wspace=0.34)

# ---- (a) before / after per field ----------------------------------------------
axA = fig.add_subplot(gs[0, 0])
for f, b, a in zip(per_fov.fov, per_fov.before_purity, per_fov.after_purity):
    color = BLUE if a > b else RED
    axA.plot([0, 1], [b, a], "-o", color=color, ms=5, lw=1.6)
    axA.text(1.06, a, str(f), va="center", fontsize=7.5, color=color)
mb, ma = per_fov.before_purity.mean(), per_fov.after_purity.mean()
axA.plot([0, 1], [mb, ma], "-o", color="#111", ms=7, lw=2.6, zorder=5)
axA.text(1.06, ma - 0.004, "mean", va="center", fontsize=8, fontweight="bold")
axA.set_xlim(-0.25, 1.45)
axA.set_xticks([0, 1])
axA.set_xticklabels(["before", "after"], fontsize=9.5)
axA.set_ylabel("myeloid marker purity")
axA.set_ylim(0.75, 0.87)
axA.spines[["top", "right"]].set_visible(False)
axA.text(-0.02, 1.06, "(a)", transform=axA.transAxes, fontsize=13, fontweight="bold",
         va="top", ha="right")

# ---- (b) parameter sweeps on FOV 514 --------------------------------------------
def fmt_adv(v):
    if np.isinf(v):
        return "none\n(no move\nbetween cells)"
    if v == 0:
        return "1x\n(no adv.)"
    return {3.0: "20x", 1.5: "4.5x", 0.5: "1.65x"}.get(v, f"{np.exp(v):.3g}x")   # exp(advantage)


panels = [
    ("radius_um", "candidate radius (µm)", lambda v: f"{v:g}", canon["radius_um"], False),
    ("home_adv", "current-cell advantage", fmt_adv, canon["home_adv"], True),
    ("decay_um", "spatial decay length (µm)", lambda v: "no penalty" if np.isinf(v) else f"{v:g}",
     canon["decay_um"], True),
]
for i, (param, xlabel, fmt, cval, categorical) in enumerate(panels):
    ax = fig.add_subplot(gs[0, i + 1])
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
    ax.set_ylim(-0.012, 0.115)
    ax.spines[["top", "right"]].set_visible(False)
    if i == 0:
        ax.set_ylabel("change in marker purity after reassignment (FOV 514)", fontsize=9)
        ax.text(-0.03, 1.06, "(b)", transform=ax.transAxes, fontsize=13, fontweight="bold",
                va="top", ha="right")
        ax.text(0.97, -0.009, "no change", fontsize=7.5, color="#666", ha="right",
                transform=ax.get_yaxis_transform(), va="center")
    else:
        ax.tick_params(axis="y", labelsize=8)
    if param == "decay_um":
        last = s.delta_purity.iloc[-1]
        ax.annotate("no spatial\nconstraint", xy=(xs[-1], last), xytext=(xs[-1] - 1.2, last + 0.008),
                    fontsize=7.5, color=RED, ha="center",
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1))

fig.text(0.5, -0.09, "circled point: configuration used for the reported result",
         ha="center", fontsize=8, color="#444")
fig.savefig(OUT, bbox_inches="tight", dpi=200)
plt.close(fig)
print("Saved:", OUT)
print(f"panel (a) mean of per-FOV purity: {mb:.3f} -> {ma:.3f}")
