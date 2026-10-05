# Supervised Tumor-Cell Identification and Data-Quality Assessment in Spatial Transcriptomics of Brain Metastasis

Analysis code for the MSc thesis of Chen Arviv (Tel Aviv University, Shmunis School of
Biomedicine and Cancer Research / Sagol School of Neuroscience), under the joint
supervision of Prof. Tal Pupko and Prof. Reuven Stein.

The study profiles a D122 lung-carcinoma brain-metastasis model in
*Ms4a3<sup>cre</sup>:R26<sup>tdT</sup>:Cx3cr1<sup>Gfp</sup>* reporter mice on the
NanoString CosMx platform: six sagittal slices across two slides, 926,318 segmented
cells, 958 genes including eight custom add-on probes.

**Every figure, table, and analysis step in the thesis is mapped to the script that
produced it in [Supplementary Table S1](TABLE_S1.md).**

## Repository layout

```
thesis_research/                 shared library: ingest and QC pipeline, SingleR candidate
                                 filters, classifiers, path configuration (config.py)
analysis/
  01_quality_control/            Figures 2-3
  02_tumor_identification/       SingleR annotation, Figures 4-15, Tables 2-3
  03_data_quality/               Figures 16-21, Table 5
    probe_detection/             probe signal against the negative-control background
    contamination_decontx/       contamination correction (DecontX)
    segmentation_fastreseg/      segmentation-error correction (FastReseg)
drafts/                          exploratory code and draft text, not part of the thesis
tools/                           how this repository was assembled (see section 7)
TABLE_S1.md                      script behind every thesis display item
```

Scripts are named after what they produce (`fig18_lyve1_bam_markers.py` makes
Figure 18); numbered prefixes give the run order within a folder.

## 1. Environment

```bash
conda create -n thesis_research python=3.12.13
conda activate thesis_research
pip install -r requirements.txt
pip install -e .            # makes the thesis_research package importable everywhere
```

`requirements.txt` pins the exact Python package versions the thesis results were
generated with.

R scripts (SingleR, DecontX, FastReseg) were run with R 4.5.2 and these package versions:

| Package | Version |
|---|---|
| SingleR | 2.12.0 |
| celda (DecontX) | 1.26.0 |
| FastReseg | 1.1.2 |
| TabulaMurisSenisData | 1.16.0 |
| SingleCellExperiment | 1.32.0 |
| SummarizedExperiment | 1.40.0 |
| scuttle | 1.20.0 |
| scater | 1.38.0 |
| celldex | 1.20.0 |
| Seurat | 5.4.0 |
| Matrix | 1.7.4 |
| BiocParallel | 1.44.0 |
| data.table | 1.18.2.1 |
| ggplot2 | 4.0.2 |
| ggrepel | 0.9.6 |

FastReseg may be installed into a separate R library; set `FASTRESEG_RLIB` to that
folder and the FastReseg scripts add it to `.libPaths()`.

Fixed random seeds are set inside the scripts that use stochastic methods; the
classifier comparison and the held-out specificity analysis use seed 42 (Tables 3 and 4
of the thesis).

## 2. Data layout and configuration

Raw data, caches and results are **not** in this repository. They live in one project
data folder, laid out as:

```
<project data folder>/
  resources/
    cosmx/<sample_id>/                  vendor export per slide (L321, L34)
      <sample_id>_exprMat_file.csv
      <sample_id>_metadata_file.csv
      <sample_id>_fov_positions_file.csv
      <sample_id>_fov_slices.csv
      <sample_id>_slice_types.csv
    cache/                              per-slice AnnData written by the pipeline
      sample_<sample_id>_adata.h5ad, slice_{1..6}_adata.h5ad
      with_tumor_prediction/            earlier tumor calls
      with_tumor_prediction_final/      final tumor calls (written in Chapter 3, step 1)
      decontx/
  outputs/cell_annotation/              SingleR references and score tables
    D122_Reference_Avinoam/             GSE103548 LLC1 bulk reference
  thesis_plots/                         figures and result tables
  agents/outputs/                       FastReseg inputs and results
```

Every script reads its inputs from, and writes its outputs to, this folder. Point the
scripts at it with an environment variable; there is no default, so no script runs
until it is set:

```bash
export THESIS_PROJECT_ROOT=/path/to/project-data-folder
export THESIS_L321_TX_FILE=/path/to/L321_tx_file.csv      # FastReseg steps only
export FASTRESEG_RLIB=/path/to/R/library                  # optional, see section 1
```

Some scripts refuse to overwrite existing outputs; others overwrite them. To
regenerate results without touching existing ones, point `THESIS_PROJECT_ROOT` at a
copy of the data folder.

## 3. Running the analysis

Run every command from the repository root. Chapters are ordered; each step depends on
the outputs of the ones before it.

### Chapter 2 — Quality control

```bash
python -m thesis_research.pipeline.run_pipeline                    # ingest, slice split, QC; Figure 1, Table 1 counts
python analysis/01_quality_control/fig02_qc_metrics_overview.py    # Figure 2
python analysis/01_quality_control/fig03_count_threshold.py        # Figure 3
```

The pipeline reads the vendor export, splits each slide into its three sagittal slices
by FOV coordinates, writes per-slice AnnData to `resources/cache/`, and applies the
low-transcript filter (per slice, `threshold = max(20, P5)` of the per-cell total-count
distribution) in `thesis_research/pipeline/cell_qc_plots.py`.

### Chapter 2 — Tumor-cell identification

```bash
Rscript analysis/02_tumor_identification/01_singler_annotate.R          # once per slice
Rscript analysis/02_tumor_identification/02_singler_scores_to_table.R   # once per slice
python analysis/02_tumor_identification/fig04_07_singler_sets.py        # Figures 4-7
python analysis/02_tumor_identification/table02_singler_threshold_sweep.py  # Table 2
python -m thesis_research.pipeline.cell_type_annotation.tumor_cells.refine_annotation_classifiers
python analysis/02_tumor_identification/fig08_classifier_comparison.py  # Figure 8
python analysis/02_tumor_identification/xgboost_default_sensitivity.py  # XGBoost defaults check
python analysis/02_tumor_identification/table03_heldout_specificity.py  # Table 3
python analysis/02_tumor_identification/fig09_14_spatial_refinement.py  # Figures 9-14
python analysis/02_tumor_identification/fig15_final_tumor_calls.py      # Figure 15
python analysis/02_tumor_identification/stage1_threshold_sensitivity.py # threshold counts and anchor purity
python analysis/02_tumor_identification/stage1_threshold_heldout.py    # threshold sensitivity, retrained (candidate floor x margin)
python analysis/02_tumor_identification/stage1_anchor_floor_heldout.py # threshold sensitivity, retrained (anchor floor)
```

SingleR is run against a three-part reference (brain structural, brain immune, and an
LLC1 carcinoma reference from GEO accession GSE103548) and returns per-cell
`score_brain_struct`, `score_brain_immune`, `score_tumor` and a predicted label. Tumor
candidates are the cells passing the filters in Table 4 of the thesis; the candidate
filters and reference-class construction are in `identify_tumor_cells.py`
(`_get_tumor_candidates_ids`, `_get_tumor_ref_ids`, `_get_healthy_ref_ids`).
XGBoost (library defaults) is refit on the joint reference pool and applied to all
candidates; cells with *P*(tumor) > 0.5 are tumor.

### Chapter 3 — Data-quality assessment

```bash
python analysis/03_data_quality/write_final_tumor_calls.py              # final calls -> with_tumor_prediction_final/
python analysis/03_data_quality/probe_detection/run_on_final_tumor_calls.py  # Table 5, Figure 16
python analysis/03_data_quality/fig17_tdtomato_prevalence.py            # Figure 17
python analysis/03_data_quality/fig18_lyve1_bam_markers.py              # Figure 18
```

`run_on_final_tumor_calls.py` runs `table05_detection_strength.py` and
`fig16_probe_detection.py` unchanged, with only their tumor-call input redirected to
`with_tumor_prediction_final/` and their outputs to `thesis_plots/rerun_final_tumor_calls/`.

Contamination correction with DecontX. All six slices use the same k-means partition
(30 truncated-SVD components, k = 25) as DecontX's population groups; the labels are
included in `contamination_decontx/decontx_clusters/`:

```bash
python analysis/03_data_quality/contamination_decontx/run_decontx.py export
Rscript analysis/03_data_quality/contamination_decontx/decontx_model.R $THESIS_PROJECT_ROOT/resources/cache/decontx
python analysis/03_data_quality/contamination_decontx/run_decontx.py assemble
python analysis/03_data_quality/contamination_decontx/fig19_decontx_before_after.py  # Figure 19
python analysis/03_data_quality/contamination_decontx/partition_sensitivity.py       # partition sensitivity, slice 1
```

Segmentation-error correction with FastReseg on all 257 FOVs of slice 1. Its reference
profiles come from all QC-passed cells of slice 1, grouped by the vendor's InSituType
cell typing; steps 01-04 obtain them from a run on five FOVs, whose reference is
identical to that of the slice-1 run, and export the score matrix for Figure 20:

```bash
python  analysis/03_data_quality/segmentation_fastreseg/01_extract_reference_fovs.py
python  analysis/03_data_quality/segmentation_fastreseg/02_prepare_reference_inputs.py
Rscript analysis/03_data_quality/segmentation_fastreseg/03_run_fastreseg_reference_fovs.R
Rscript analysis/03_data_quality/segmentation_fastreseg/04_export_reference_scores.R
python  analysis/03_data_quality/segmentation_fastreseg/05_prepare_slice1_inputs.py
Rscript analysis/03_data_quality/segmentation_fastreseg/06_run_fastreseg_slice1.R
python  analysis/03_data_quality/segmentation_fastreseg/07_evaluate_slice1.py   # read-outs
python  analysis/03_data_quality/segmentation_fastreseg/08_sparsity_slice1.py   # sparsity
python  analysis/03_data_quality/segmentation_fastreseg/fig20_reference_heatmap.py  # Figure 20
python  analysis/03_data_quality/fig21_tumor_spatial.py                          # Figure 21
```

## 4. Verifying a run

After the QC pipeline, the headline numbers should match Table 1 of the thesis:

| Quantity | Expected |
|---|---|
| Segmented cells before QC | 926,318 |
| Removed by low-transcript filter | 80,211 |
| Retained | 846,107 (91.34%) |
| FOVs | 1,205 |
| Per-slice retention | 88.86 / 91.88 / 92.74 / 74.49 / 93.07 / 95.20 % |

The classifier comparison should give out-of-fold ROC-AUC above 0.99 for all five
classifiers, and XGBoost accuracy 0.970, precision 0.985, recall 0.963 and F1 0.974,
with tumor as the positive class (`figure_3_model_comparison.csv`).

The final calls should retain zero tumor cells in the two control slices (3 and 4) out
of 343 and 287 candidates, and 20,873 tumor cells in the four tumor-bearing slices
(4,267 / 6,272 / 5,039 / 5,295). These zeros are in-sample: the control candidates are
the classifier's negative training class. On look-alikes withheld from training
(Table 3), XGBoost calls 13 of the 630 tumor out-of-fold, and 30 of 630 when each
control slice is withheld in turn.

## 5. Known limitations

- **No R lock file.** R package versions are listed in section 1, but the R
  environment is not restored automatically.
- **Pipeline figure paths.** `run_pipeline.py` generates a fresh identifier on every
  run, so Figure 1 is written to a new `outputs/<run_id>/<sample_id>/` folder each time.
- **SingleR scripts run one slice at a time.** The slice is set inside
  `01_singler_annotate.R` and `02_singler_scores_to_table.R` and must be changed for
  each slice.
- **Tables are laid out by hand.** Tables 1-5 were formatted in the thesis from the
  outputs listed in Table S1; no script in `analysis/` renders them. (An earlier
  rendering script, `drafts/thesis_plots/make_nature_tables.py`, uses the old table
  numbering and is kept for provenance only.)
- **Two tumor-call caches.** `with_tumor_prediction/` holds earlier calls (20,679 tumor
  cells); the Chapter 3 analyses use `with_tumor_prediction_final/` (20,873, the calls
  reported in the thesis) through `run_on_final_tumor_calls.py`. Run directly,
  `table05_detection_strength.py` and `fig16_probe_detection.py` reproduce the earlier
  non-tumor count (n = 825,428 instead of 825,234).
- **Output file names** (`figure_3_model_comparison.png`, `figure_4_spatial_refinement_*`,
  `dq_fig_*`) predate the thesis numbering; Table S1 maps them.

## 6. Data and code availability

Processed per-slice single-cell data, the assembled SingleR reference objects and the
raw CosMx exports (transcript tables and segmentation masks) are available from the
author on request. The cluster labels supplied to DecontX are included in
`analysis/03_data_quality/contamination_decontx/decontx_clusters/`.

All animal procedures, tumor implantation, tissue processing, and histological
preparation were performed by Mr. Avinoam Ratzabi, Tel Aviv University, in accordance
with institutional and ethical guidelines.

## 7. How this repository was assembled

The analyses were developed in a separate working repository.
`tools/build_from_working_repo.py` copied the code behind every thesis display item
into this structure, and all other code and draft text into `drafts/` with its
original paths. The copied thesis code differs from the working versions only in file
names, in reading the project folder from `THESIS_PROJECT_ROOT` instead of a fixed
path, and in imports between scripts that now sit in different folders; every changed
line is listed in `tools/build_report.md`.
