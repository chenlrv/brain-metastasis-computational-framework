"""Rerun the classifier-comparison outputs with the random forest at library defaults.

The random forest in figure_3_model_comparison.py and figure_4_spatial_refinement.py
now uses the scikit-learn default of 100 trees (previously 300). This runner
regenerates every output that involves the random forest -- Figure 8, Figures
9-14 and Table 3 -- into a new folder, thesis_plots/rerun_rf_default_trees/,
leaving the earlier outputs in thesis_plots/ untouched. Each script is executed
unchanged except for its output folder. Refuses to run if the new folder exists.

Run: conda run -n thesis_research python thesis_plots/rerun_rf_default_trees.py
"""
import os

for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
OUT_SUB = "rerun_rf_default_trees"
OUT_DIR = HERE / OUT_SUB
sys.path.insert(0, str(HERE))

# script -> its output-folder line, redirected into OUT_DIR
SCRIPTS = {
    "figure_3_model_comparison.py": 'OUT_DIR = BASE_DIR / "thesis_plots"',
    "figure_4_spatial_refinement.py": 'OUT_DIR = BASE_DIR / "thesis_plots"',
    "control_specificity_validation.py": 'OUT_DIR = BASE_DIR / "thesis_plots" / "control_specificity"',
}


def redirected_source(script, line):
    src = (HERE / script).read_text(encoding="utf-8")
    if src.count(line) != 1:
        raise ValueError(f"{script}: expected one occurrence of {line!r}")
    new = line.replace('"thesis_plots"', f'"thesis_plots" / "{OUT_SUB}"')
    return src.replace(line, new)


def main():
    if OUT_DIR.exists():
        raise FileExistsError(f"{OUT_DIR} already exists; refusing to overwrite")
    OUT_DIR.mkdir()
    for script, line in SCRIPTS.items():
        print(f"\n===== {script}", flush=True)
        code = compile(redirected_source(script, line), str(HERE / script), "exec")
        exec(code, {"__name__": "__main__", "__file__": str(HERE / script)})
    print(f"\nOutputs in {OUT_DIR}")


if __name__ == "__main__":
    main()
