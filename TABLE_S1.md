**Table S1 | Analysis scripts underlying each figure, table, and analysis step**

| Display item | Analysis | Script | Primary output |
| --- | --- | --- | --- |
| *Chapter 2 — Quality control* |  |  |  |
| Figure 1 | Segmented-cell area and FOV layout per slide | `thesis_research/pipeline/position_plots.py` (`plot_cells_positions_with_area`, `plot_fov_positions_all_slices`), run by `thesis_research/pipeline/run_pipeline.py` | `cell_positions.png`, `fov_positions_sliced.png`<sup>a</sup> |
| Figure 2 | Per-cell QC metrics across six tissue slices | `analysis/01_quality_control/fig02_qc_metrics_overview.py` | `thesis_plots/qc_metrics_overview.png` |
| Figure 3 | Per-cell transcript-count distributions and the low-count threshold | `analysis/01_quality_control/fig03_count_threshold.py`, applying the cutoff rule of `thesis_research/pipeline/cell_qc_plots.py` (`_low_count_flag`) | `thesis_plots/qc_count_threshold_hist.png` |
| Table 1 | Slices profiled and cell yield after low-transcript QC | `thesis_research/pipeline/cell_qc_plots.py` (`run_cell_qc`), run by `run_pipeline.py`<sup>b</sup> | per-slice cell counts |
| *Chapter 2 — Tumor-cell identification* |  |  |  |
| — | Reference-based annotation with SingleR | `analysis/02_tumor_identification/01_singler_annotate.R`, `02_singler_scores_to_table.R`<sup>c</sup> | `outputs/cell_annotation/<slide>/05/<slice>/slice_<n>_final_cell_annotations_…csv` |
| Figures 4–7 | SingleR tumor labels, tumor candidates, tumor anchors and healthy look-alikes | `analysis/02_tumor_identification/fig04_07_singler_sets.py` | `thesis_plots/fig_set_raw_labels.png`, `fig_set_candidates.png`, `fig_set_anchors.png`, `fig_set_lookalikes.png`; `fig_singler_sets.csv` |
| Table 2 | SingleR candidates under increasingly strict score thresholds | `analysis/02_tumor_identification/table02_singler_threshold_sweep.py` | `thesis_plots/singler_sweep.csv` |
| — | Classifier training and refinement (Stages 2 and 3) | `thesis_research/pipeline/cell_type_annotation/tumor_cells/identify_tumor_cells.py`, `refine_annotation_classifiers.py` | reference pools and tumor calls |
| Figure 8 | Classifier comparison on the joint reference pool | `analysis/02_tumor_identification/fig08_classifier_comparison.py` | `thesis_plots/figure_3_model_comparison.png/.csv`<sup>d</sup> |
| — | XGBoost default-parameter sensitivity analysis | `analysis/02_tumor_identification/xgboost_default_sensitivity.py` | `thesis_plots/xgb_hparam_sensitivity/cv_metrics.csv`, `slice_calls.csv`, `agreement.csv` |
| Table 3 | False-positive tumor calls on control-tissue look-alikes withheld from training | `analysis/02_tumor_identification/table03_heldout_specificity.py` | `thesis_plots/control_specificity/oof_false_positives.csv`, `loso_false_positives.csv` |
| Figures 9–14 | Candidates retained and rejected by each classifier, one figure per slice | `analysis/02_tumor_identification/fig09_14_spatial_refinement.py` | `thesis_plots/figure_4_spatial_refinement_slice{1..6}.png`; `figure_4_spatial_refinement.csv`<sup>d</sup> |
| Figure 15 | Final XGBoost tumor calls, all six slices | `analysis/02_tumor_identification/fig15_final_tumor_calls.py` | `thesis_plots/final_xgboost_refinement.png/.csv`<sup>d</sup> |
| — | Stage-1 threshold calibration sweeps | `analysis/02_tumor_identification/stage1_threshold_sensitivity.py` | `thesis_plots/stage1_sensitivity/sweep_a_candidate_floor.csv`, `sweep_b_delta_margin.csv`, `sweep_c_reference_floor.csv` |
| Table 4 | Parameter settings of the tumor-refinement pipeline | values as defined in `identify_tumor_cells.py` and `refine_annotation_classifiers.py`<sup>e</sup> | — |
| *Chapter 3 — Data-quality assessment* |  |  |  |
| — | Final tumor calls (Figure 15) written to the per-slice data, defining the non-tumor cells analyzed in this chapter | `analysis/03_data_quality/write_final_tumor_calls.py` | `resources/cache/with_tumor_prediction_final/slice_{1..6}_adata.h5ad` |
| Figure 16 | Probe detection against the per-slice acceptance threshold | `analysis/03_data_quality/probe_detection/fig16_probe_detection.py`, run on the final tumor calls by `probe_detection/run_on_final_tumor_calls.py`<sup>f</sup> | `thesis_plots/rerun_final_tumor_calls/dq_fig1_detection.png`; `acceptance_bar_all6.csv` |
| Table 5 | Detection strength of the eight custom add-on probes | `analysis/03_data_quality/probe_detection/table05_detection_strength.py`, run on the final tumor calls by `probe_detection/run_on_final_tumor_calls.py`<sup>f</sup> | `thesis_plots/rerun_final_tumor_calls/detection_reliability_all6.csv` |
| Figure 17 | tdTomato prevalence against pan-myeloid transcripts, and the control-to-tumor contrast | `analysis/03_data_quality/fig17_tdtomato_prevalence.py` | `thesis_plots/dq_fig_reporter.png`; `reporter_prevalence_all6.csv` |
| Figure 18 | Lyve1-positive cells against the canonical BAM markers, in composition and in space | `analysis/03_data_quality/fig18_lyve1_bam_markers.py` | `thesis_plots/dq_fig_lyve1.png` |
| — | Cluster partitions supplied to DecontX | `analysis/03_data_quality/contamination_decontx/write_cluster_labels.py`<sup>g</sup> | `contamination_decontx/decontx_clusters/slice_{1..6}_clusters.csv` |
| — | Contamination estimation and correction (DecontX) | `analysis/03_data_quality/contamination_decontx/run_decontx.py` (export and assembly); `decontx_model.R` (DecontX)<sup>h</sup> | `resources/cache/decontx/slice_{1..6}_decontx.h5ad`; `decontx_contamination.csv` |
| Figure 19 | Probe positivity before and after DecontX correction | `analysis/03_data_quality/contamination_decontx/fig19_decontx_before_after.py` | `thesis_plots/dq_fig_decontx.png`; `decontx_before_after.csv` |
| — | DecontX partition sensitivity (k-means k = 10, 25, 50 and a five-cluster Leiden partition, slice 1) | `analysis/03_data_quality/contamination_decontx/partition_sensitivity.py` | `thesis_plots/decontx_partition_sensitivity.csv` |
| — | FastReseg reference profiles (count matrix, InSituType cell groups and score matrix) | `analysis/03_data_quality/segmentation_fastreseg/01_extract_reference_fovs.py` to `04_export_reference_scores.R`<sup>i</sup> | `agents/outputs/segmentation_fastreseg/inputs/`, `figures/score_matrix.csv` |
| — | Segmentation-error correction with FastReseg, all 257 FOVs of slice 1 | `analysis/03_data_quality/segmentation_fastreseg/05_prepare_slice1_inputs.py`, `06_run_fastreseg_slice1.R` | `agents/outputs/segmentation_fastreseg_slice1/fastreseg_out/` |
| — | Read-outs before and after correction and per-cell sparsity, QC-passed non-tumor cells | `analysis/03_data_quality/segmentation_fastreseg/07_evaluate_slice1.py`, `08_sparsity_slice1.py` | `agents/outputs/segmentation_fastreseg_slice1/eval_qcpassed/summary.json`, `sparsity.json` |
| Figure 20 | FastReseg reference scores for the anomaly genes across the 12 cell groups | `analysis/03_data_quality/segmentation_fastreseg/fig20_reference_heatmap.py`<sup>j</sup> | `agents/outputs/segmentation_fastreseg_slice1/figures/fig20_reference_heatmap_slice1.png` |
| Figure 21 | Spatial distribution of the tumor-cell calls, all six slices | `analysis/03_data_quality/fig21_tumor_spatial.py` | `thesis_plots/dq_fig_tumor_spatial.png` |

Every figure, table, and analysis step in this thesis and the script that produced it. Script paths are relative to the repository root; output paths are relative to the project data folder set by `THESIS_PROJECT_ROOT` (README, section 2). Running the scripts in the order given in the README regenerates each output from the processed data.

<sup>a</sup>Written to `outputs/<run_id>/<sample_id>/`, where `run_id` is a fresh identifier generated on each pipeline run, so these files have no stable path.

<sup>b</sup>The per-slice cell counts are computed during the pipeline's QC step; `fig03_count_threshold.py` re-derives the same cutoff from the raw vendor metadata and checks its per-slice counts against Table 1.

<sup>c</sup>Both scripts process one slice per run; the slice is set inside each script.

<sup>d</sup>Output file names predate the thesis figure numbering and do not match it.

<sup>e</sup>Table 4 documents parameter values; it is not produced by a script.

<sup>f</sup>The runner executes the listed script unchanged, redirecting only its tumor-call input to `with_tumor_prediction_final/` and its outputs to `thesis_plots/rerun_final_tumor_calls/`.

<sup>g</sup>The partitions are a DecontX input, not an output; they are included in the repository so that the correction can be reproduced exactly.

<sup>h</sup>DecontX is run from R (celda). The Python driver exports the count matrix and cluster labels, calls the R script, and assembles the corrected matrix; corrected counts are rounded to integers on assembly.

<sup>i</sup>The reference profiles are estimated from all QC-passed cells of slice 1, grouped by the vendor's InSituType cell typing, and are identical for the slice-1 run. Steps 01–04 run FastReseg on five FOVs to obtain them and export their score matrix for Figure 20.

<sup>j</sup>Executes `fig20_heatmap_base.py` with only the cell-group shares (from `eval_qcpassed/cell_groups_nontumor.csv`) and the output path changed.

BAM, border-associated macrophage; FOV, field of view; QC, quality control.
