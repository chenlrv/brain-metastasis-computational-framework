# Build report
Source: `D:\thesis-research`

55 thesis files, 7 data files, 277 draft files.

## `thesis_research/__init__.py` -> `thesis_research/__init__.py`

unchanged

## `thesis_research/pipeline/__init__.py` -> `thesis_research/pipeline/__init__.py`

unchanged

## `thesis_research/pipeline/run_pipeline.py` -> `thesis_research/pipeline/run_pipeline.py`

unchanged

## `thesis_research/pipeline/position_plots.py` -> `thesis_research/pipeline/position_plots.py`

unchanged

## `thesis_research/pipeline/cell_qc_plots.py` -> `thesis_research/pipeline/cell_qc_plots.py`

unchanged

## `thesis_research/pipeline/filters.py` -> `thesis_research/pipeline/filters.py`

unchanged

## `thesis_research/pipeline/fov_qc_plots.py` -> `thesis_research/pipeline/fov_qc_plots.py`

unchanged

## `thesis_research/pipeline/sample_slice.py` -> `thesis_research/pipeline/sample_slice.py`

unchanged

## `thesis_research/pipeline/utils.py` -> `thesis_research/pipeline/utils.py`

unchanged

## `thesis_research/pipeline/cell_type_annotation/__init__.py` -> `thesis_research/pipeline/cell_type_annotation/__init__.py`

unchanged

## `thesis_research/pipeline/cell_type_annotation/tumor_cells/__init__.py` -> `thesis_research/pipeline/cell_type_annotation/tumor_cells/__init__.py`

unchanged

## `thesis_research/pipeline/cell_type_annotation/tumor_cells/identify_tumor_cells.py` -> `thesis_research/pipeline/cell_type_annotation/tumor_cells/identify_tumor_cells.py`

```diff
@@ -0,0 +1 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -20 +21 @@
-BASE_DIR = r"D:/thesis-research/"
+BASE_DIR = PROJECT_ROOT_STR + r"/"
@@ -309 +310 @@
-    adata = ad.read_h5ad(fr"D:\thesis-research\resources\cache\slice_{slice_id}_adata.h5ad")
+    adata = ad.read_h5ad(fr"{PROJECT_ROOT_STR}/resources/cache/slice_{slice_id}_adata.h5ad")
@@ -311 +312 @@
-        fr"D:/thesis-research/outputs/cell_annotation/{slide_id}/05/{slice_id}/slice_{slice_id}_final_cell_annotations_refined_tabula_brain_tumor_avinoam.csv"
+        fr"{PROJECT_ROOT_STR}/outputs/cell_annotation/{slide_id}/05/{slice_id}/slice_{slice_id}_final_cell_annotations_refined_tabula_brain_tumor_avinoam.csv"
```

## `thesis_research/pipeline/cell_type_annotation/tumor_cells/refine_annotation_classifiers.py` -> `thesis_research/pipeline/cell_type_annotation/tumor_cells/refine_annotation_classifiers.py`

```diff
@@ -0,0 +1 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -30 +31 @@
-BASE_DIR = r"D:/thesis-research/"
+BASE_DIR = PROJECT_ROOT_STR + r"/"
@@ -193 +194 @@
-        # (thesis_plots/xgb_default_sensitivity.py) found no hyperparameter
+        # (analysis/02_tumor_identification/xgboost_default_sensitivity.py) found no hyperparameter
@@ -1036 +1037 @@
-    # xgboost 3.2.0 library defaults -- see thesis_plots/xgb_default_sensitivity.py.
+    # xgboost 3.2.0 library defaults -- see analysis/02_tumor_identification/xgboost_default_sensitivity.py.
@@ -1281 +1282 @@
-        df = pd.read_csv(rf"D:/thesis-research/resources/cosmx/{slide_id}/{slide_id}_metadata_file.csv")
+        df = pd.read_csv(rf"{PROJECT_ROOT_STR}/resources/cosmx/{slide_id}/{slide_id}_metadata_file.csv")
@@ -1310 +1311 @@
-        # plt.savefig(rf"D:/thesis-research/resources/cosmx/{slide_id}/{slide_id}_NC_ratio_kde.png", dpi=150)
+        # plt.savefig(rf"{PROJECT_ROOT_STR}/resources/cosmx/{slide_id}/{slide_id}_NC_ratio_kde.png", dpi=150)
```

## `thesis_research/pipeline/cell_type_annotation/tumor_cells/classifiers.py` -> `thesis_research/pipeline/cell_type_annotation/tumor_cells/classifiers.py`

```diff
@@ -0,0 +1 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -28 +29 @@
-BASE_DIR = r"D:/thesis-research/"
+BASE_DIR = PROJECT_ROOT_STR + r"/"
@@ -196 +197 @@
-        # (thesis_plots/xgb_default_sensitivity.py) found no hyperparameter
+        # (analysis/02_tumor_identification/xgboost_default_sensitivity.py) found no hyperparameter
```

## `thesis_research/utils/__init__.py` -> `thesis_research/utils/__init__.py`

unchanged

## `thesis_research/utils/columns.py` -> `thesis_research/utils/columns.py`

unchanged

## `thesis_research/utils/constants.py` -> `thesis_research/utils/constants.py`

```diff
@@ -3 +3 @@
-PROJECT_ROOT = Path(__file__).resolve().parents[2]
+from thesis_research.config import PROJECT_ROOT  # set by THESIS_PROJECT_ROOT
```

## `thesis_research/utils/entity_type.py` -> `thesis_research/utils/entity_type.py`

unchanged

## `thesis_plots/make_qc_metrics_overview_fig.py` -> `analysis/01_quality_control/fig02_qc_metrics_overview.py`

```diff
@@ -10 +10 @@
-    conda run -n thesis_research python thesis_plots/make_qc_metrics_overview_fig.py
+    conda run -n thesis_research python analysis/01_quality_control/fig02_qc_metrics_overview.py
@@ -14,0 +15 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -24 +25 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/make_qc_count_threshold_fig.py` -> `analysis/01_quality_control/fig03_count_threshold.py`

```diff
@@ -10 +10 @@
-Run: conda run -n thesis_research python thesis_plots/make_qc_count_threshold_fig.py
+Run: conda run -n thesis_research python analysis/01_quality_control/fig03_count_threshold.py
@@ -11,0 +12 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -19 +20 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `outputs/cell_annotation/annotate.R` -> `analysis/02_tumor_identification/01_singler_annotate.R`

```diff
@@ -0,0 +1,2 @@
+PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
+if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
@@ -18 +20 @@
-data_dir <- "D:/thesis-research/outputs/cell_annotation/D122_Reference_Avinoam/GSE103548_GeneCount_raw.tsv.gz"
+data_dir <- paste0(PROJECT_ROOT, "/outputs/cell_annotation/D122_Reference_Avinoam/GSE103548_GeneCount_raw.tsv.gz")
@@ -42 +44 @@
-setwd("D:/thesis-research/outputs/cell_annotation/L321/05")
+setwd(paste0(PROJECT_ROOT, "/outputs/cell_annotation/L321/05"))
```

## `outputs/cell_annotation/convert_annotations_to_df.R` -> `analysis/02_tumor_identification/02_singler_scores_to_table.R`

unchanged

## `thesis_plots/make_fig_singler_sets.py` -> `analysis/02_tumor_identification/fig04_07_singler_sets.py`

```diff
@@ -19 +19 @@
-Run: conda run -n thesis_research python thesis_plots/make_fig_singler_sets.py
+Run: conda run -n thesis_research python analysis/02_tumor_identification/fig04_07_singler_sets.py
@@ -20,0 +21 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -29 +30 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/make_table_singler_sweep.py` -> `analysis/02_tumor_identification/table02_singler_threshold_sweep.py`

```diff
@@ -18 +18 @@
-Run: conda run -n thesis_research python thesis_plots/make_table_singler_sweep.py
+Run: conda run -n thesis_research python analysis/02_tumor_identification/table02_singler_threshold_sweep.py
@@ -19,0 +20 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -22 +23 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/figure_3_model_comparison.py` -> `analysis/02_tumor_identification/fig08_classifier_comparison.py`

```diff
@@ -23,0 +24 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -50 +51 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
@@ -163 +164 @@
-    # (xgb_default_sensitivity.py) showed every hand-set value we previously
+    # (xgboost_default_sensitivity.py) showed every hand-set value we previously
```

## `thesis_plots/xgb_default_sensitivity.py` -> `analysis/02_tumor_identification/xgboost_default_sensitivity.py`

```diff
@@ -36,0 +37 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -57 +58 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/stage1_threshold_sensitivity.py` -> `analysis/02_tumor_identification/stage1_threshold_sensitivity.py`

```diff
@@ -31,0 +32 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -43 +44 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/stage1_threshold_heldout.py` -> `analysis/02_tumor_identification/stage1_threshold_heldout.py`

```diff
@@ -36 +36 @@
-Run: conda run -n thesis_research python thesis_plots/stage1_threshold_heldout.py
+Run: conda run -n thesis_research python analysis/02_tumor_identification/stage1_threshold_heldout.py
@@ -37,0 +38 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -54 +55 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/stage1_anchor_floor_heldout.py` -> `analysis/02_tumor_identification/stage1_anchor_floor_heldout.py`

```diff
@@ -11 +11 @@
-Run: conda run -n thesis_research python thesis_plots/stage1_anchor_floor_heldout.py
+Run: conda run -n thesis_research python analysis/02_tumor_identification/stage1_anchor_floor_heldout.py
```

## `thesis_plots/figure_4_spatial_refinement.py` -> `analysis/02_tumor_identification/fig09_14_spatial_refinement.py`

```diff
@@ -50,0 +51 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -74 +75 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
@@ -194 +195 @@
-    # balancing is deliberately not applied -- see figure_3_model_comparison.py.
+    # balancing is deliberately not applied -- see fig08_classifier_comparison.py.
@@ -204 +205 @@
-    # xgboost 3.2.0 library defaults -- see xgb_default_sensitivity.py.
+    # xgboost 3.2.0 library defaults -- see xgboost_default_sensitivity.py.
@@ -208 +209 @@
-    # scikit-learn defaults (100 trees), as in figure_3_model_comparison.py.
+    # scikit-learn defaults (100 trees), as in fig08_classifier_comparison.py.
```

## `thesis_plots/control_specificity_validation.py` -> `analysis/02_tumor_identification/table03_heldout_specificity.py`

```diff
@@ -10 +10 @@
-     and seed as figure_3_model_comparison.py). Every look-alike is scored by a
+     and seed as fig08_classifier_comparison.py). Every look-alike is scored by a
@@ -19 +19 @@
-figure_4_spatial_refinement.py, so the decision rules (P(tumor) > 0.5; the
+fig09_14_spatial_refinement.py, so the decision rules (P(tumor) > 0.5; the
@@ -27 +27 @@
-Run: conda run -n thesis_research python thesis_plots/control_specificity_validation.py
+Run: conda run -n thesis_research python analysis/02_tumor_identification/table03_heldout_specificity.py
@@ -28,0 +29 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -45 +46 @@
-from figure_4_spatial_refinement import fit_classifiers, score_candidates  # noqa: E402
+from fig09_14_spatial_refinement import fit_classifiers, score_candidates  # noqa: E402
@@ -51 +52 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/final_xgboost_refinement.py` -> `analysis/02_tumor_identification/fig15_final_tumor_calls.py`

```diff
@@ -33,0 +34 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -47 +48 @@
-from figure_4_spatial_refinement import _to_dense, build_reference_pool  # noqa: E402
+from fig09_14_spatial_refinement import _to_dense, build_reference_pool  # noqa: E402
@@ -53 +54 @@
-BASE_DIR = pathlib.Path(r"D:\thesis-research")
+BASE_DIR = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/update_tumor_prediction_cache.py` -> `analysis/03_data_quality/write_final_tumor_calls.py`

```diff
@@ -7 +7 @@
-in which pred_tumor_XGBoost holds the calls of final_xgboost_refinement.py (same
+in which pred_tumor_XGBoost holds the calls of fig15_final_tumor_calls.py (same
@@ -12 +12 @@
-Run: conda run -n thesis_research python thesis_plots/update_tumor_prediction_cache.py
+Run: conda run -n thesis_research python analysis/03_data_quality/write_final_tumor_calls.py
@@ -19,2 +19,2 @@
-sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
-from final_xgboost_refinement import (  # noqa: E402
+sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "02_tumor_identification"))
+from fig15_final_tumor_calls import (  # noqa: E402
```

## `thesis_plots/make_detection_reliability_6slice.py` -> `analysis/03_data_quality/probe_detection/table05_detection_strength.py`

```diff
@@ -12 +12 @@
-Run: conda run -n thesis_research python thesis_plots/make_detection_reliability_6slice.py
+Run: conda run -n thesis_research python analysis/03_data_quality/probe_detection/table05_detection_strength.py
@@ -13,0 +14 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -24 +25 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/make_dq_fig1_detection.py` -> `analysis/03_data_quality/probe_detection/fig16_probe_detection.py`

```diff
@@ -16 +16 @@
-Run: conda run -n thesis_research python thesis_plots/make_dq_fig1_detection.py
+Run: conda run -n thesis_research python analysis/03_data_quality/probe_detection/fig16_probe_detection.py
@@ -17,0 +18 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -28 +29 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/rerun_chapter3_final_calls.py` -> `analysis/03_data_quality/probe_detection/run_on_final_tumor_calls.py`

```diff
@@ -3 +3 @@
-make_detection_reliability_6slice.py (Table 4) and make_dq_fig1_detection.py
+table05_detection_strength.py (Table 5) and fig16_probe_detection.py
@@ -8 +8 @@
-version written by update_tumor_prediction_cache.py, and every output path, to
+version written by write_final_tumor_calls.py, and every output path, to
@@ -13 +13 @@
-Run: conda run -n thesis_research python thesis_plots/rerun_chapter3_final_calls.py
+Run: conda run -n thesis_research python analysis/03_data_quality/probe_detection/run_on_final_tumor_calls.py
@@ -14,0 +15 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -24 +25 @@
-OUT_DIR = HERE / OUT_SUB
+OUT_DIR = pathlib.Path(PROJECT_ROOT_STR) / "thesis_plots" / OUT_SUB
@@ -30 +31 @@
-    "make_detection_reliability_6slice.py": [
+    "table05_detection_strength.py": [
@@ -32 +33 @@
-    "make_dq_fig1_detection.py": [
+    "fig16_probe_detection.py": [
```

## `thesis_plots/make_dq_fig_reporter.py` -> `analysis/03_data_quality/fig17_tdtomato_prevalence.py`

```diff
@@ -12 +12 @@
-Run: conda run -n thesis_research python thesis_plots/make_dq_fig_reporter.py
+Run: conda run -n thesis_research python analysis/03_data_quality/fig17_tdtomato_prevalence.py
@@ -13,0 +14 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -26 +27 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/make_dq_fig_lyve1.py` -> `analysis/03_data_quality/fig18_lyve1_bam_markers.py`

```diff
@@ -14 +14 @@
-Run: conda run -n thesis_research python thesis_plots/make_dq_fig_lyve1.py
+Run: conda run -n thesis_research python analysis/03_data_quality/fig18_lyve1_bam_markers.py
@@ -15,0 +16 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -25 +26 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `score_genes/write_decontx_clusters.py` -> `analysis/03_data_quality/contamination_decontx/write_cluster_labels.py`

```diff
@@ -14 +14 @@
-Run: conda run -n thesis_research python score_genes/write_decontx_clusters.py
+Run: conda run -n thesis_research python analysis/03_data_quality/contamination_decontx/write_cluster_labels.py
@@ -26 +26 @@
-import run_decontx_correct as m  # noqa: E402
+import run_decontx as m  # noqa: E402
```

## `score_genes/run_decontx_correct.py` -> `analysis/03_data_quality/contamination_decontx/run_decontx.py`

```diff
@@ -5 +5 @@
-  1. python run_decontx_correct.py export
+  1. python run_decontx.py export
@@ -8,2 +8,2 @@
-        Rscript score_genes/run_decontx.R  D:/thesis-research/resources/cache/decontx
-        (or open run_decontx.R in RStudio, set WORKDIR_PARENT, and source() it)
+        Rscript analysis/03_data_quality/contamination_decontx/decontx_model.R  $THESIS_PROJECT_ROOT/resources/cache/decontx
+        (or open decontx_model.R in RStudio, set WORKDIR_PARENT, and source() it)
@@ -11 +11 @@
-  3. python run_decontx_correct.py assemble
+  3. python run_decontx.py assemble
@@ -14 +14 @@
-  python run_decontx_correct.py          # 'auto': export -> Rscript subprocess -> assemble
+  python run_decontx.py          # 'auto': export -> Rscript subprocess -> assemble
@@ -18 +18 @@
-(written by write_decontx_clusters.py) are passed as z. Every gene is kept --
+(written by write_cluster_labels.py) are passed as z. Every gene is kept --
@@ -21,0 +22 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -38 +39 @@
-    f"slice_{i}": f"D:/thesis-research/resources/cache/with_tumor_prediction/slice_{i}_adata.h5ad"
+    f"slice_{i}": f"{PROJECT_ROOT_STR}/resources/cache/with_tumor_prediction/slice_{i}_adata.h5ad"
@@ -41,2 +42,2 @@
-OUT_DIR = "D:/thesis-research/resources/cache/decontx"
-R_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_decontx.R")
+OUT_DIR = PROJECT_ROOT_STR + "/resources/cache/decontx"
+R_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "decontx_model.R")
@@ -215,3 +216,3 @@
-              f"  Rscript score_genes/run_decontx.R {OUT_DIR}\n"
-              f"  (or open score_genes/run_decontx.R in RStudio and source it)\n"
-              f"Then: python score_genes/run_decontx_correct.py assemble")
+              f"  Rscript analysis/03_data_quality/contamination_decontx/decontx_model.R {OUT_DIR}\n"
+              f"  (or open analysis/03_data_quality/contamination_decontx/decontx_model.R in RStudio and source it)\n"
+              f"Then: python analysis/03_data_quality/contamination_decontx/run_decontx.py assemble")
```

## `score_genes/run_decontx.R` -> `analysis/03_data_quality/contamination_decontx/decontx_model.R`

```diff
@@ -4 +4 @@
-#   Terminal: Rscript run_decontx.R [parent_dir]
+#   Terminal: Rscript decontx_model.R [parent_dir]
@@ -8,0 +9,2 @@
+PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
+if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
@@ -15 +17 @@
-WORKDIR_PARENT <- "D:/thesis-research/resources/cache/decontx"
+WORKDIR_PARENT <- paste0(PROJECT_ROOT, "/resources/cache/decontx")
```

## `thesis_plots/make_dq_fig_decontx.py` -> `analysis/03_data_quality/contamination_decontx/fig19_decontx_before_after.py`

```diff
@@ -12 +12 @@
-Run: conda run -n thesis_research python thesis_plots/make_dq_fig_decontx.py
+Run: conda run -n thesis_research python analysis/03_data_quality/contamination_decontx/fig19_decontx_before_after.py
@@ -13,0 +14 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -25 +26 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```

## `thesis_plots/decontx_partition_sensitivity.py` -> `analysis/03_data_quality/contamination_decontx/partition_sensitivity.py`

```diff
@@ -8 +8 @@
-               (write_decontx_clusters.py, seed 0), with k changed
+               (write_cluster_labels.py, seed 0), with k changed
@@ -16 +16 @@
-Requires the DecontX export of slice 1 (run_decontx_correct.py export), i.e.
+Requires the DecontX export of slice 1 (run_decontx.py export), i.e.
@@ -24 +24 @@
-  conda run -n thesis_research python thesis_plots/decontx_partition_sensitivity.py
+  conda run -n thesis_research python analysis/03_data_quality/contamination_decontx/partition_sensitivity.py
@@ -25,0 +26 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -38,4 +39,4 @@
-ROOT = "D:/thesis-research"
-sys.path.insert(0, os.path.join(ROOT, "score_genes"))
-import run_decontx_correct as m  # noqa: E402
-import write_decontx_clusters as w  # noqa: E402
+ROOT = PROJECT_ROOT_STR
+sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
+import run_decontx as m  # noqa: E402
+import write_cluster_labels as w  # noqa: E402
@@ -46 +47 @@
-CLUSTERS = os.path.join(ROOT, "score_genes", "decontx_clusters")
+CLUSTERS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "decontx_clusters")
@@ -53 +54 @@
-    """Same preprocessing and seed as write_decontx_clusters.py, for each k."""
+    """Same preprocessing and seed as write_cluster_labels.py, for each k."""
@@ -65 +66 @@
-    """One parent dir per partition, so run_decontx.R processes only that run."""
+    """One parent dir per partition, so decontx_model.R processes only that run."""
@@ -79 +80 @@
-                                "score_genes/run_decontx_correct.py export first")
+                                "analysis/03_data_quality/contamination_decontx/run_decontx.py export first")
```

## `agents/segmentation/03_filter_tx.py` -> `analysis/03_data_quality/segmentation_fastreseg/01_extract_reference_fovs.py`

```diff
@@ -11 +11 @@
-  conda run -n thesis_research python agents/segmentation/03_filter_tx.py
+  conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/01_extract_reference_fovs.py
@@ -12,0 +13 @@
+from thesis_research.config import PROJECT_ROOT_STR, L321_TX_FILE  # noqa: E402
@@ -17,3 +18,2 @@
-TX = (r"D:\20251214_CosMx_ReuvenStein\20251214_CosMx_ReuvenStein.tar\Analysis"
-      r"\L321__1__31_12_2025_12_32_59_204\flatFiles\L321\L321_tx_file.csv")
-OUT = r"D:\thesis-research\agents\outputs\segmentation\tx_by_fov"
+TX = L321_TX_FILE
+OUT = PROJECT_ROOT_STR + r"/agents/outputs/segmentation/tx_by_fov"
```

## `agents/segmentation/fastreseg/01_prep_inputs.py` -> `analysis/03_data_quality/segmentation_fastreseg/02_prepare_reference_inputs.py`

```diff
@@ -17 +17 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/01_prep_inputs.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/02_prepare_reference_inputs.py
@@ -18,0 +19 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -27 +28 @@
-ROOT = pathlib.Path(r"D:/thesis-research")
+ROOT = pathlib.Path(PROJECT_ROOT_STR)
```

## `agents/segmentation/fastreseg/02_run_fastreseg.R` -> `analysis/03_data_quality/segmentation_fastreseg/03_run_fastreseg_reference_fovs.R`

```diff
@@ -3 +3 @@
-# Inputs come from 01_prep_inputs.py. Only data-specific settings are given:
+# Inputs come from 02_prepare_reference_inputs.py. Only data-specific settings are given:
@@ -14 +14,3 @@
-# Run: Rscript agents/segmentation/fastreseg/02_run_fastreseg.R
+# Run: Rscript analysis/03_data_quality/segmentation_fastreseg/03_run_fastreseg_reference_fovs.R
+PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
+if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
@@ -16 +18 @@
-.libPaths(c("D:/R-libs/fastreseg", .libPaths()))
+.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))
@@ -23 +25 @@
-base <- "D:/thesis-research/agents/outputs/segmentation_fastreseg"
+base <- paste0(PROJECT_ROOT, "/agents/outputs/segmentation_fastreseg")
```

## `agents/segmentation/fastreseg/04_export_scores.R` -> `analysis/03_data_quality/segmentation_fastreseg/04_export_reference_scores.R`

```diff
@@ -3 +3 @@
-# 02_run_fastreseg.R, for plotting. Writes only to
+# 03_run_fastreseg_reference_fovs.R, for plotting. Writes only to
@@ -5 +5,3 @@
-.libPaths(c("D:/R-libs/fastreseg", .libPaths()))
+PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
+if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
+.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))
@@ -7 +9 @@
-base <- "D:/thesis-research/agents/outputs/segmentation_fastreseg"
+base <- paste0(PROJECT_ROOT, "/agents/outputs/segmentation_fastreseg")
```

## `agents/segmentation/fastreseg/slice1_01_prep_inputs.py` -> `analysis/03_data_quality/segmentation_fastreseg/05_prepare_slice1_inputs.py`

```diff
@@ -6 +6 @@
-the same as for the five-FOV run (01_prep_inputs.py) and are read from there, not
+the same as for the five-FOV run (02_prepare_reference_inputs.py) and are read from there, not
@@ -10 +10 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/slice1_01_prep_inputs.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/05_prepare_slice1_inputs.py
@@ -11,0 +12 @@
+from thesis_research.config import PROJECT_ROOT_STR, L321_TX_FILE  # noqa: E402
@@ -18,3 +19,2 @@
-ROOT = pathlib.Path("D:/thesis-research")
-TX = pathlib.Path("D:/20251214_CosMx_ReuvenStein/20251214_CosMx_ReuvenStein.tar/Analysis/"
-                  "L321__1__31_12_2025_12_32_59_204/flatFiles/L321/L321_tx_file.csv")
+ROOT = pathlib.Path(PROJECT_ROOT_STR)
+TX = pathlib.Path(L321_TX_FILE)
```

## `agents/segmentation/fastreseg/slice1_02_run_fastreseg.R` -> `analysis/03_data_quality/segmentation_fastreseg/06_run_fastreseg_slice1.R`

```diff
@@ -2,2 +2,2 @@
-# Identical to 02_run_fastreseg.R except for the input/output folders: counts and
-# clusters are those of the five-FOV run; transcripts come from slice1_01_prep_inputs.py.
+# Identical to 03_run_fastreseg_reference_fovs.R except for the input/output folders: counts and
+# clusters are those of the five-FOV run; transcripts come from 05_prepare_slice1_inputs.py.
@@ -5 +5 @@
-# Inputs come from 01_prep_inputs.py. Only data-specific settings are given:
+# Inputs come from 02_prepare_reference_inputs.py. Only data-specific settings are given:
@@ -16 +16,3 @@
-# Run: Rscript agents/segmentation/fastreseg/02_run_fastreseg.R
+# Run: Rscript analysis/03_data_quality/segmentation_fastreseg/06_run_fastreseg_slice1.R
+PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
+if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
@@ -18 +20 @@
-.libPaths(c("D:/R-libs/fastreseg", .libPaths()))
+.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))
@@ -25,2 +27,2 @@
-base <- "D:/thesis-research/agents/outputs/segmentation_fastreseg_slice1"
-counts_dir <- "D:/thesis-research/agents/outputs/segmentation_fastreseg/inputs"
+base <- paste0(PROJECT_ROOT, "/agents/outputs/segmentation_fastreseg_slice1")
+counts_dir <- paste0(PROJECT_ROOT, "/agents/outputs/segmentation_fastreseg/inputs")
```

## `agents/segmentation/fastreseg/slice1_04_evaluate_qcpassed.py` -> `analysis/03_data_quality/segmentation_fastreseg/07_evaluate_slice1.py`

```diff
@@ -9,2 +9,2 @@
-(slice1_01_prep_inputs.py). "After" = FastReseg's updated tables
-(slice1_02_run_fastreseg.R; transcripts it trimmed are absent). Tumor cells are
+(05_prepare_slice1_inputs.py). "After" = FastReseg's updated tables
+(06_run_fastreseg_slice1.R; transcripts it trimmed are absent). Tumor cells are
@@ -25,0 +26 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -32 +33 @@
-ROOT = pathlib.Path("D:/thesis-research")
+ROOT = pathlib.Path(PROJECT_ROOT_STR)
```

## `agents/segmentation/fastreseg/slice1_05_sparsity.py` -> `analysis/03_data_quality/segmentation_fastreseg/08_sparsity_slice1.py`

```diff
@@ -7 +7 @@
-Inputs: the slice-1 transcript files of slice1_01_prep_inputs.py (vendor cell
+Inputs: the slice-1 transcript files of 05_prepare_slice1_inputs.py (vendor cell
@@ -11 +11 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/slice1_05_sparsity.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/08_sparsity_slice1.py
@@ -12,0 +13 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -19 +20 @@
-ROOT = pathlib.Path("D:/thesis-research")
+ROOT = pathlib.Path(PROJECT_ROOT_STR)
```

## `agents/segmentation/fastreseg/slice1_06_figure20.py` -> `analysis/03_data_quality/segmentation_fastreseg/fig20_reference_heatmap.py`

```diff
@@ -3 +3 @@
-Runs 06_figure_A2.py unchanged except that the top bar shows the share of the
+Runs fig20_heatmap_base.py unchanged except that the top bar shows the share of the
@@ -6,2 +6,2 @@
-slice1_04_evaluate_qcpassed.py), and that the output goes to the slice-1 folder.
-The heatmap scores come from figures/score_matrix.csv (04_export_scores.R); the
+07_evaluate_slice1.py), and that the output goes to the slice-1 folder.
+The heatmap scores come from figures/score_matrix.csv (04_export_reference_scores.R); the
@@ -10 +10 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/slice1_06_figure20.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/fig20_reference_heatmap.py
@@ -11,0 +12 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -15 +16 @@
-S1 = pathlib.Path("D:/thesis-research/agents/outputs/segmentation_fastreseg_slice1")
+S1 = pathlib.Path(PROJECT_ROOT_STR + "/agents/outputs/segmentation_fastreseg_slice1")
@@ -17 +18 @@
-src = (HERE / "06_figure_A2.py").read_text(encoding="utf-8")
+src = (HERE / "fig20_heatmap_base.py").read_text(encoding="utf-8")
@@ -28 +29 @@
-        raise ValueError(f"06_figure_A2.py changed; cannot find: {old}")
+        raise ValueError(f"fig20_heatmap_base.py changed; cannot find: {old}")
@@ -33 +34 @@
-exec(compile(src, str(HERE / "06_figure_A2.py") + "[slice1]", "exec"), {"__name__": "__main__"})
+exec(compile(src, str(HERE / "fig20_heatmap_base.py") + "[slice1]", "exec"), {"__name__": "__main__"})
```

## `agents/segmentation/fastreseg/06_figure_A2.py` -> `analysis/03_data_quality/segmentation_fastreseg/fig20_heatmap_base.py`

```diff
@@ -9 +9 @@
-misplaced. Inputs from 04_export_scores.R and 05_figures.py (score_matrix.csv,
+misplaced. Inputs from 04_export_reference_scores.R and 05_figures.py (score_matrix.csv,
@@ -12 +12 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/06_figure_A2.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/fig20_heatmap_base.py
@@ -13,0 +14 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -24 +25 @@
-FIG = pathlib.Path("D:/thesis-research/agents/outputs/segmentation_fastreseg/figures")
+FIG = pathlib.Path(PROJECT_ROOT_STR + "/agents/outputs/segmentation_fastreseg/figures")
```

## `agents/segmentation/fastreseg/slice1_07_figS1_groups.py` -> `analysis/03_data_quality/segmentation_fastreseg/figS1_cell_groups.py`

```diff
@@ -6 +6 @@
-slice1_04_evaluate_qcpassed.py), the same as the top bar of Figure 20; tumor cells
+07_evaluate_slice1.py), the same as the top bar of Figure 20; tumor cells
@@ -12 +12 @@
-Run: conda run -n thesis_research python agents/segmentation/fastreseg/slice1_07_figS1_groups.py
+Run: conda run -n thesis_research python analysis/03_data_quality/segmentation_fastreseg/figS1_cell_groups.py
@@ -13,0 +14 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -24 +25 @@
-ROOT = pathlib.Path("D:/thesis-research")
+ROOT = pathlib.Path(PROJECT_ROOT_STR)
```

## `thesis_plots/make_dq_fig_tumor_spatial.py` -> `analysis/03_data_quality/fig21_tumor_spatial.py`

```diff
@@ -12 +12 @@
-final_xgboost_refinement.py (the source of the final-calls figure), so the two
+fig15_final_tumor_calls.py (the source of the final-calls figure), so the two
@@ -18 +18 @@
-  conda run -n thesis_research python thesis_plots/make_dq_fig_tumor_spatial.py
+  conda run -n thesis_research python analysis/03_data_quality/fig21_tumor_spatial.py
@@ -19,0 +20 @@
+from thesis_research.config import PROJECT_ROOT_STR  # noqa: E402
@@ -26 +27,4 @@
-from final_xgboost_refinement import fit_current_xgboost, load_slice_refit
+import pathlib  # noqa: E402
+import sys  # noqa: E402
+sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "02_tumor_identification"))
+from fig15_final_tumor_calls import fit_current_xgboost, load_slice_refit
@@ -28 +32 @@
-ROOT = "D:/thesis-research"
+ROOT = PROJECT_ROOT_STR
```
