"""Rerun the Chapter-3 detection analysis on the final tumor calls, as a new version.

table05_detection_strength.py (Table 5) and fig16_probe_detection.py
(Figure 16) define non-tumor cells from resources/cache/with_tumor_prediction/,
whose pred_tumor_XGBoost column predates the final Stage-3 model (20,679 tumor
cells instead of 20,873). This runner executes each script with only two things
redirected: the tumor-prediction cache, to the with_tumor_prediction_final/
version written by write_final_tumor_calls.py, and every output path, to
thesis_plots/rerun_final_tumor_calls/. The scripts themselves, the original
cache and the original outputs are left untouched, and an existing output in
the new folder is never overwritten.

Run: conda run -n thesis_research python analysis/03_data_quality/probe_detection/run_on_final_tumor_calls.py
"""
from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
import os

for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(var, "1")

import pathlib

HERE = pathlib.Path(__file__).resolve().parent
OUT_SUB = "rerun_final_tumor_calls"
OUT_DIR = pathlib.Path(PROJECT_ROOT_STR) / "thesis_plots" / OUT_SUB
CACHE_OLD = "/resources/cache/with_tumor_prediction/"
CACHE_NEW = "/resources/cache/with_tumor_prediction_final/"

# script -> output files it writes (each redirected into OUT_DIR)
SCRIPTS = {
    "table05_detection_strength.py": [
        "detection_reliability_all6.png", "detection_reliability_all6.csv"],
    "fig16_probe_detection.py": [
        "detection_reliability_all6.csv",  # read as input: the rerun table above
        "acceptance_bar_all6.csv", "dq_fig1_detection.png"],
}


def redirected_source(script, outputs):
    src = (HERE / script).read_text(encoding="utf-8")
    if src.count(CACHE_OLD) != 1:
        raise ValueError(f"{script}: expected one tumor-prediction cache path")
    src = src.replace(CACHE_OLD, CACHE_NEW)
    for name in outputs:
        old, new = f"/thesis_plots/{name}", f"/thesis_plots/{OUT_SUB}/{name}"
        if src.count(old) != 1:
            raise ValueError(f"{script}: expected one occurrence of {old}")
        src = src.replace(old, new)
    return src


def main():
    OUT_DIR.mkdir(exist_ok=True)
    written = {n for outs in SCRIPTS.values() for n in outs[-2:]} | {"detection_reliability_all6.csv"}
    clash = [n for n in written if (OUT_DIR / n).exists()]
    if clash:
        raise FileExistsError(f"{OUT_DIR} already holds {clash}; refusing to overwrite")
    for script, outputs in SCRIPTS.items():
        print(f"\n===== {script}", flush=True)
        code = compile(redirected_source(script, outputs), str(HERE / script), "exec")
        exec(code, {"__name__": "__main__", "__file__": str(HERE / script)})
    print(f"\nOutputs in {OUT_DIR}")


if __name__ == "__main__":
    main()
