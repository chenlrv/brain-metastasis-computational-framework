"""Write a new version of the tumor-prediction cache holding the final Stage-3 calls.

The pred_tumor_XGBoost column of resources/cache/with_tumor_prediction/ was
written on 2026-05-21, before Stage 3 moved to the xgboost 3.2.0 defaults, and
holds 20,679 tumor cells. The thesis reports the refit calls (20,873). This
script writes a copy of each slice to resources/cache/with_tumor_prediction_final/
in which pred_tumor_XGBoost holds the calls of fig15_final_tumor_calls.py (same
model, same candidates, same threshold); the earlier values are kept alongside as
pred_tumor_XGBoost_20260521. The original cache is only read, never modified,
and an existing output file is never overwritten.

Run: conda run -n thesis_research python analysis/03_data_quality/write_final_tumor_calls.py
"""
import pathlib
import sys

import anndata as ad

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "02_tumor_identification"))
from fig15_final_tumor_calls import (  # noqa: E402
    PRED_CACHE, PRED_COL, SLIDE_CACHE, SLICES_TO_SHOW, fit_current_xgboost, load_slice_refit,
)

OUT_CACHE = PRED_CACHE.parent / "with_tumor_prediction_final"
OLD_COL = PRED_COL + "_20260521"


def main():
    OUT_CACHE.mkdir(exist_ok=True)
    existing = [p.name for p in OUT_CACHE.glob("slice_*_adata.h5ad")]
    if existing:
        raise FileExistsError(f"{OUT_CACHE} already holds {existing}; refusing to overwrite")

    model, joint_var_names = fit_current_xgboost()
    total_old = total_new = 0
    for slice_id in SLICES_TO_SHOW:
        _, _, is_kept, _ = load_slice_refit(slice_id, model, joint_var_names)

        adata = ad.read_h5ad(PRED_CACHE / f"slice_{slice_id}_adata.h5ad")
        # load_slice_refit scores the cells in the order of the main slice cache;
        # the calls can only be copied across if both caches list the same cells.
        ref = ad.read_h5ad(SLIDE_CACHE / f"slice_{slice_id}_adata.h5ad", backed="r")
        same_order = len(is_kept) == adata.n_obs and (ref.obs_names == adata.obs_names).all()
        ref.file.close()
        if not same_order:
            raise ValueError(f"slice {slice_id}: cell order differs between caches")

        old = adata.obs[PRED_COL].to_numpy(dtype=bool)
        adata.obs[OLD_COL] = old
        adata.obs[PRED_COL] = is_kept
        out = OUT_CACHE / f"slice_{slice_id}_adata.h5ad"
        if out.exists():
            raise FileExistsError(out)
        adata.write_h5ad(out)

        total_old += int(old.sum())
        total_new += int(is_kept.sum())
        print(f"slice {slice_id}: {int(old.sum()):,} -> {int(is_kept.sum()):,} tumor "
              f"({int((old & ~is_kept).sum())} dropped, {int((~old & is_kept).sum())} added)",
              flush=True)
    print(f"total: {total_old:,} -> {total_new:,}; written to {OUT_CACHE}")


if __name__ == "__main__":
    main()
