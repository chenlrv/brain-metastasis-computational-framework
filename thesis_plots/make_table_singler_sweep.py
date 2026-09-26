"""Table 2 -- SingleR tumor candidates under increasingly strict score thresholds.

Shows why a single score threshold cannot replace the supervised refinement:
every candidate in the two sham-injected control slices is a false positive by
construction, and the floor at which the controls first return no candidates
retains almost none of the candidates in tumor-bearing tissue.

For each slice, cells whose SingleR label is "Tumor" are counted
  * with no further filter (the raw label), and
  * after the candidate filters -- a margin of more than 0.08 between
    score_tumor and the next-best reference score, best-class consistency, and
    score_tumor above each floor in FLOORS.
These are the filters of identify_tumor_cells.py (_get_tumor_candidates_ids),
with the score_tumor floor varied.

Saves: thesis_plots/singler_sweep.csv (per slice) and prints the pooled table.

Run: conda run -n thesis_research python thesis_plots/make_table_singler_sweep.py
"""
import pandas as pd

ROOT = "D:/thesis-research"
ANNOT = ROOT + "/outputs/cell_annotation/{slide}/05/{s}/slice_{s}_final_cell_annotations_refined_tabula_brain_tumor_avinoam.csv"
OUT_CSV = ROOT + "/thesis_plots/singler_sweep.csv"

SLICES = {1: ("L321", "tumor-bearing"), 2: ("L321", "tumor-bearing"),
          3: ("L321", "control"), 4: ("L34", "control"),
          5: ("L34", "tumor-bearing"), 6: ("L34", "tumor-bearing")}
FLOORS = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
DELTA = 0.08


def main():
    rows = []
    for s, (slide, kind) in SLICES.items():
        df = pd.read_csv(ANNOT.format(slide=slide, s=s))
        t = df[df["predicted_cell_type"] == "Tumor"].copy()
        t["next_best"] = t[["score_brain_struct", "score_brain_immune"]].max(axis=1)
        t["delta"] = (t["score_tumor"] - t["next_best"]).abs()
        passes = t[(t["delta"] > DELTA) & (t["score_tumor"] > t["next_best"])]
        row = {"slice": s, "slide": slide, "type": kind, "n_cells": len(df),
               "raw_tumor_label": len(t)}
        for f in FLOORS:
            row[f"floor_{f:.2f}"] = int((passes["score_tumor"] > f).sum())
        rows.append(row)

    per_slice = pd.DataFrame(rows)
    per_slice.to_csv(OUT_CSV, index=False)

    tumor = per_slice[per_slice["type"] == "tumor-bearing"]
    ctrl = per_slice[per_slice["type"] == "control"]
    base = tumor["floor_0.20"].sum()
    print("floor   tumor-bearing   control   % of the 0.20 pool")
    print(f"raw     {tumor['raw_tumor_label'].sum():>13,}   {ctrl['raw_tumor_label'].sum():>7,}")
    for f in FLOORS:
        c = f"floor_{f:.2f}"
        print(f"{f:.2f}    {tumor[c].sum():>13,}   {ctrl[c].sum():>7,}   "
              f"{100 * tumor[c].sum() / base:6.1f}")
    print("\nSaved:", OUT_CSV)


if __name__ == "__main__":
    main()
