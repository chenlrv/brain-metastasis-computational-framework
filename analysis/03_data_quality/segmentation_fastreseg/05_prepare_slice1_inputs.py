"""FastReseg inputs for ALL FOVs of slice 1 (257 FOVs).

Streams the raw L321 transcript table and writes one transcript file per slice-1
FOV, restricted to the 958 panel genes (negative and system-control probes
removed), with a transcript ID per row. The count matrix and cluster labels are
the same as for the five-FOV run (02_prepare_reference_inputs.py) and are read from there, not
copied. Writes only to agents/outputs/segmentation_fastreseg_slice1/inputs/,
which must not exist.

Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/05_prepare_slice1_inputs.py
"""
from thesis_research.config import PROJECT_ROOT_STR, L321_TX_FILE  # noqa: E402
import pathlib
import time

import anndata as ad
import pandas as pd

ROOT = pathlib.Path(PROJECT_ROOT_STR)
TX = pathlib.Path(L321_TX_FILE)
FIVE = ROOT / "agents/outputs/segmentation_fastreseg/inputs"
OUT = ROOT / "agents/outputs/segmentation_fastreseg_slice1/inputs"
COLS = ["fov", "cell_ID", "cell", "x_global_px", "y_global_px", "z", "target"]


def main():
    if OUT.exists():
        raise FileExistsError(f"{OUT} already exists; refusing to overwrite")
    OUT.mkdir(parents=True)
    a = ad.read_h5ad(ROOT / "resources/cache/with_tumor_prediction_final/slice_1_adata.h5ad", backed="r")
    fovs = sorted(a.obs["fov"].astype(int).unique())
    a.file.close()
    genes = set(pd.read_csv(FIVE / "genes.csv")["gene"])
    fov_set = set(fovs)

    t0, n_in, n_kept, written = time.time(), 0, 0, set()
    counters = {f: 0 for f in fovs}
    for chunk in pd.read_csv(TX, usecols=COLS, chunksize=2_000_000):
        n_in += len(chunk)
        chunk = chunk[chunk["fov"].isin(fov_set) & chunk["target"].isin(genes)]
        for fov, sub in chunk.groupby("fov"):
            sub = sub.copy()
            start = counters[fov]
            sub["transcript_id"] = [f"t_{fov}_{i}" for i in range(start, start + len(sub))]
            counters[fov] += len(sub)
            sub[["transcript_id", "cell", "cell_ID", "x_global_px", "y_global_px", "z", "target"]] \
                .to_csv(OUT / f"tx_fov{fov}.csv", mode="a", index=False, header=fov not in written)
            written.add(fov)
        n_kept += len(chunk)
        print(f"  read {n_in:,} transcripts, kept {n_kept:,} ({time.time() - t0:.0f}s)", flush=True)

    missing = [f for f in fovs if f not in written]
    if missing:
        raise ValueError(f"no transcripts for FOVs {missing}")
    info = pd.DataFrame({"file_path": [str(OUT / f"tx_fov{f}.csv").replace("\\", "/") for f in fovs],
                         "slide": 1, "fov": fovs, "stage_X": 0.0, "stage_Y": 0.0})
    info.to_csv(OUT / "transDF_fileInfo.csv", index=False)
    pd.Series({"counts_dir": str(FIVE), "n_fovs": len(fovs), "n_transcripts": n_kept}).to_csv(OUT / "README.csv")
    print(f"{len(fovs)} FOVs, {n_kept:,} panel-gene transcripts written to {OUT}")


if __name__ == "__main__":
    main()
