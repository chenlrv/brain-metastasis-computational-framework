**Table S1 | Analysis scripts underlying each figure, table, and analysis step**

| Display item | Analysis | Script | Primary output |
| --- | --- | --- | --- |
| *Chapter 2 — Quality control* |  |  |  |
| Figure 1 | Segmented-cell area and FOV layout per slide | `thesis_research/pipeline/position_plots.py` (`plot_cells_positions_with_area`, `plot_fov_positions_all_slices`), driven by `thesis_research/pipeline/run_pipeline.py` | `cell_positions.png`, `fov_positions_sliced.png`<sup>a</sup> |
| Figure 2 | Per-cell QC metrics across six tissue slices | `thesis_plots/make_qc_metrics_overview_fig.py` | `qc_metrics_overview.png` |
| Figure 3 | Per-cell transcript-count distributions and the low-count threshold | `thesis_plots/make_qc_count_threshold_fig.py`, applying the cutoff rule of `thesis_research/pipeline/cell_qc_plots.py` (`_low_count_flag`) | `qc_count_threshold_hist.png` |
| Table 1 | Slices profiled and cell yield after low-transcript QC | `cell_qc_plots.py` (counts); `thesis_plots/make_nature_tables.py` (layout) | `tables/thesis_tables.html`<sup>b</sup> |
| *Chapter 2 — Tumor-cell identification* |  |  |  |
| — | Reference-based annotation with SingleR | `outputs/cell_annotation/annotate.R`, `convert_annotations_to_df.R` | SingleR score CSVs |
| Figure 4 | Cells labeled Tumor by SingleR, before filtering | `thesis_plots/make_fig_singler_sets.py` | `fig_set_raw_labels.png`; `fig_singler_sets.csv` |
| Figure 5 | Tumor candidates across the six slices | `thesis_plots/make_fig_singler_sets.py` | `fig_set_candidates.png` |
| Figure 6 | Tumor anchors, the positive training class | `thesis_plots/make_fig_singler_sets.py` | `fig_set_anchors.png` |
| Figure 7 | Healthy look-alikes, the negative training class | `thesis_plots/make_fig_singler_sets.py` | `fig_set_lookalikes.png` |
| Table 2 | SingleR candidates under increasingly strict score thresholds | `thesis_plots/make_table_singler_sweep.py` | `singler_sweep.csv` |
| Figure 8 | Classifier comparison on the joint reference pool | `thesis_plots/figure_3_model_comparison.py`<sup>c</sup> | `figure_3_model_comparison.png/.csv` |
| — | XGBoost default-parameter sensitivity analysis | `thesis_plots/xgb_default_sensitivity.py` | `xgb_hparam_sensitivity/cv_metrics.csv`, `slice_calls.csv`, `agreement.csv` |
| Table 3 | False-positive tumor calls on control-tissue look-alikes withheld from training (out-of-fold and leave-one-control-slice-out) | `thesis_plots/control_specificity_validation.py` | `control_specificity/oof_false_positives.csv`, `loso_false_positives.csv` |
| Figures 9–14 | Candidates retained and rejected by each classifier, one figure per slice | `thesis_plots/figure_4_spatial_refinement.py`<sup>c</sup> | `figure_4_spatial_refinement_slice{1..6}.png`; `figure_4_spatial_refinement.csv` |
| Figure 15 | Final XGBoost tumor calls, all six slices | `thesis_plots/final_xgboost_refinement.py` | `final_xgboost_refinement.png/.csv` |
| — | Stage-1 threshold calibration sweeps | `thesis_plots/stage1_threshold_sensitivity.py` | `stage1_sensitivity/sweep_a_candidate_floor.csv`, `sweep_b_delta_margin.csv`, `sweep_c_reference_floor.csv` |
| Table 4 | Parameter settings of the tumor-refinement pipeline | `…/tumor_cells/identify_tumor_cells.py`, `refine_annotation_classifiers.py` (values); `make_nature_tables.py` (layout) | `tables/thesis_tables.html`<sup>b</sup> |
| *Chapter 3 — Data-quality assessment* |  |  |  |
| — | Final tumor calls (Figure 15) written to the per-slice data, defining the non-tumor cells analyzed in this chapter | `thesis_plots/update_tumor_prediction_cache.py` | `resources/cache/with_tumor_prediction_final/slice_{1..6}_adata.h5ad` |
| Figure 16 | Probe detection against the per-slice acceptance threshold | `thesis_plots/make_dq_fig1_detection.py`, run on the final tumor calls by `thesis_plots/rerun_chapter3_final_calls.py`<sup>f</sup> | `rerun_final_tumor_calls/dq_fig1_detection.png`; `acceptance_bar_all6.csv` |
| Table 5 | Detection strength of the eight custom add-on probes | `thesis_plots/make_detection_reliability_6slice.py`, run on the final tumor calls by `thesis_plots/rerun_chapter3_final_calls.py`<sup>f</sup> | `rerun_final_tumor_calls/detection_reliability_all6.csv` |
| Figure 17 | tdTomato prevalence against pan-myeloid transcripts, and the control-to-tumor contrast | `thesis_plots/make_dq_fig_reporter.py` | `dq_fig_reporter.png`; `reporter_prevalence_all6.csv` |
| Figure 18 | Lyve1-positive cells against the canonical BAM markers, in composition and in space | `thesis_plots/make_dq_fig_lyve1.py` | `dq_fig_lyve1.png` |
| — | Cluster partitions supplied to DecontX | `score_genes/write_decontx_clusters.py`<sup>d</sup> | `score_genes/decontx_clusters/slice_{1..6}_clusters.csv` |
| — | Contamination estimation and correction (DecontX) | `score_genes/run_decontx_correct.py` (export and assembly); `score_genes/run_decontx.R` (DecontX)<sup>e</sup> | `slice_{1..6}_decontx.h5ad`; `decontx_contamination.csv` |
| Figure 19 | Probe positivity before and after DecontX correction | `thesis_plots/make_dq_fig_decontx.py` | `dq_fig_decontx.png`; `decontx_before_after.csv` |
| — | DecontX partition sensitivity (k-means k = 10, 25, 50 and a five-cluster Leiden partition, slice 1) | `thesis_plots/decontx_partition_sensitivity.py` | `decontx_partition_sensitivity.csv` |
| — | Segmentation-error correction with FastReseg, all 257 FOVs of slice 1 | `agents/segmentation/fastreseg/slice1_01_prep_inputs.py`, `slice1_02_run_fastreseg.R`<sup>g</sup> | `agents/outputs/segmentation_fastreseg_slice1/fastreseg_out/` |
| — | Read-outs before and after correction (transcripts removed and moved, Lyve1 and tdTomato positivity) and per-cell sparsity, QC-passed non-tumor cells | `agents/segmentation/fastreseg/slice1_04_evaluate_qcpassed.py`, `slice1_05_sparsity.py` | `…/segmentation_fastreseg_slice1/eval_qcpassed/summary.json`, `sparsity.json` |
| Figure 20 | FastReseg reference scores for the anomaly genes across the 12 cell groups | `agents/segmentation/fastreseg/04_export_scores.R`<sup>g</sup>, `slice1_06_figure20.py`<sup>h</sup> | `…/segmentation_fastreseg_slice1/figures/fig20_reference_heatmap_slice1.png` |
| Figure 21 | Spatial distribution of the tumor-cell calls, all six slices | `thesis_plots/make_dq_fig_tumor_spatial.py` | `dq_fig_tumor_spatial.png` |

Every figure, table, and analysis step in this thesis and the script that produced it. Script paths are relative to the repository root. Running the scripts in the order given in the repository README regenerates each output from the processed data.

<sup>a</sup>Written to `outputs/<run_id>/<sample_id>/`, where `run_id` is a fresh identifier generated on each pipeline run (`run_pipeline.py`), so these files have no stable path.

<sup>b</sup>Values are hardcoded in `make_nature_tables.py` rather than read from the analysis outputs.

<sup>c</sup>The script filename numbering predates the thesis figure numbering and does not match it.

<sup>d</sup>The partitions are a DecontX input, not an output; they are provided in the repository so that the correction can be reproduced exactly.

<sup>e</sup>DecontX is run from R (celda). The Python driver exports the count matrix and cluster labels, calls the R script, and assembles the corrected matrix; corrected counts are rounded to integers on assembly.

<sup>f</sup>The runner executes the listed script unchanged, redirecting only its tumor-call input to `with_tumor_prediction_final/` and its outputs to `thesis_plots/rerun_final_tumor_calls/`.

<sup>g</sup>The count matrix and the InSituType cell groups from which FastReseg builds its reference are prepared by `agents/segmentation/fastreseg/01_prep_inputs.py`; the reference profiles exported by `04_export_scores.R` come from the run of `02_run_fastreseg.R` on the same inputs and are identical to those of the slice-1 run.

<sup>h</sup>Executes `06_figure_A2.py` with only the cell-group shares (from `eval_qcpassed/cell_groups_nontumor.csv`) and the output path changed.

BAM, border-associated macrophage; FOV, field of view; QC, quality control.
