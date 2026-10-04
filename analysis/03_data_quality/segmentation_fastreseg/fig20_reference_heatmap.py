"""Figure 20 (thesis): FastReseg's reference for the anomaly genes, slice-1 version.

Runs fig20_heatmap_base.py unchanged except that the top bar shows the share of the
QC-passed non-tumor cells of the whole of slice 1 in each cell group, as assigned
by FastReseg in the all-FOV run (eval_qcpassed/cell_groups_nontumor.csv from
07_evaluate_slice1.py), and that the output goes to the slice-1 folder.
The heatmap scores come from figures/score_matrix.csv (04_export_reference_scores.R); the
reference is the same in both runs. Refuses to overwrite.

Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/fig20_reference_heatmap.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
S1 = pathlib.Path(PROJECT_ROOT_STR + "/agents/outputs/segmentation_fastreseg_slice1")

src = (HERE / "fig20_heatmap_base.py").read_text(encoding="utf-8")
replacements = [
    ('OUT = FIG / "fig_A7_reference_heatmap.png"',
     f'OUT = pathlib.Path("{(S1 / "figures" / "fig20_reference_heatmap_slice1.png").as_posix()}")'),
    ('types = pd.read_csv(FIG / "cell_types_5fov.csv")',
     f'types = pd.read_csv("{(S1 / "eval_qcpassed" / "cell_groups_nontumor.csv").as_posix()}")'),
    ('share = types.loc[~types.tumor, "type"].value_counts(normalize=True)',
     'share = types["group"].value_counts(normalize=True)'),
]
for old, new in replacements:
    if src.count(old) != 1:
        raise ValueError(f"fig20_heatmap_base.py changed; cannot find: {old}")
    start = src.index(old)
    end = src.index("\n", start) if old.startswith("OUT") else start + len(old)
    src = src[:start] + new + src[end:]
(S1 / "figures").mkdir(exist_ok=True)
exec(compile(src, str(HERE / "fig20_heatmap_base.py") + "[slice1]", "exec"), {"__name__": "__main__"})
