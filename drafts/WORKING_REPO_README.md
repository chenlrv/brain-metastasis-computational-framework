# Supervised Tumor-Cell Identification and Data-Quality Assessment in Spatial Transcriptomics of Brain Metastasis

Analysis code for the MSc thesis of Chen Arviv (Tel Aviv University, Shmunis School of
Biomedicine and Cancer Research / Sagol School of Neuroscience), under the joint
supervision of Prof. Tal Pupko and Prof. Reuven Stein.

The study profiles a D122 lung-carcinoma brain-metastasis model in
*Ms4a3<sup>cre</sup>:R26<sup>tdT</sup>:Cx3cr1<sup>Gfp</sup>* reporter mice on the
NanoString CosMx platform: six sagittal slices across two slides, 926,318 segmented
cells, 958 genes including eight custom add-on probes.

**Every figure, table, and analysis step in the thesis is mapped to the script that
produced it in [Supplementary Table S1](thesis_plots/tables/thesis_table_s1.md).**
That table is the authoritative map; this README gives the order in which the scripts
run and what each stage needs.

---

## 1. Environment

Analyses are implemented in Python and R.

```bash
conda create -n thesis_research python=3.12.13
conda activate thesis_research
pip install -r requirements.txt
```

`requirements.txt` pins the exact Python package versions the thesis results were
generated with, for the scripts listed in Table S1.

All Python entry points are then run through that environment:

```bash
conda run -n thesis_research python <script>
```

R scripts (SingleR, decontX) were run with R 4.5.2 and these package versions:

| Package | Version |
|---|---|
| SingleR | 2.12.0 |
| celda (DecontX) | 1.26.0 |
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

FastReseg 1.1.2 (segmentation-error correction, Figure 20) was built with the same R
4.5.2 and installed with its dependencies into a separate library,
`D:/R-libs/fastreseg`, which the FastReseg R scripts add to `.libPaths()`.

Fixed random seeds are set inside the scripts that use stochastic methods; the
classifier comparison and the held-out specificity analysis use seed 42 (Tables 3 and 4
of the thesis).

## 2. Data layout

Raw CosMx exports and derived caches are **not** in this repository. Place them as:

```
resources/
  cosmx/<sample_id>/            # vendor export per slide (L321, L34)
    <sample_id>_exprMat_file.csv
    <sample_id>_metadata_file.csv
    <sample_id>_fov_positions_file.csv
    <sample_id>_fov_slices.csv
    <sample_id>_slice_types.csv
  cache/                        # per-slice AnnData written by the pipeline
    slice_{1..6}*.h5ad
    with_tumor_prediction/      # earlier tumor calls (kept unmodified)
    with_tumor_prediction_final/  # final tumor calls, written in Stage 3
    decontx/slice_{1..6}_decontx.h5ad
outputs/
  cell_annotation/              # SingleR references and score tables
    D122_Reference_Avinoam/     # GSE103548 LLC1 bulk reference
```

The FastReseg analysis (Stage 4) also reads the raw slide-L321 transcript table
(`L321_tx_file.csv`, part of the vendor's flat-file export), whose location is set in
`agents/segmentation/fastreseg/slice1_01_prep_inputs.py`.

Paths are resolved from `thesis_research/utils/constants.py`, which derives everything
from the repository root. Some standalone plotting scripts still contain absolute
`D:/thesis-research/...` paths — see "Known gaps".

## 3. Running the analysis

Stages are ordered; each depends on the outputs of the ones before it.

### Stage 0 — Ingest, QC, and slice assignment

```bash
conda run -n thesis_research python -c "from thesis_research.pipeline.run_pipeline import run_pipeline; run_pipeline(run_exploration=True, run_qc=True)"
```

Reads the vendor export, splits each slide into its three sagittal slices by FOV
coordinates, writes per-slice AnnData to `resources/cache/`, and applies the
low-transcript filter — per slice, `threshold = max(20, P5)` of the per-cell total-count
distribution — in `thesis_research/pipeline/cell_qc_plots.py`
(`run_cell_qc`, `_low_count_flag`).

Produces **Figure 1** (`position_plots.py`) and the counts behind **Table 1**.
**Figure 2** (per-cell QC metrics) and **Figure 3** (the six-panel count-distribution
grid) are rendered separately; the latter re-derives the same cutoff from the raw
vendor metadata and checks its per-slice cell counts against Table 1:

```bash
conda run -n thesis_research python thesis_plots/make_qc_metrics_overview_fig.py    # Figure 2
conda run -n thesis_research python thesis_plots/make_qc_count_threshold_fig.py     # Figure 3
```

### Stage 1 — Reference-based annotation (SingleR)

```bash
Rscript outputs/cell_annotation/annotate.R
Rscript outputs/cell_annotation/convert_annotations_to_df.R
```

Builds the three-part reference — brain structural, brain immune, and an LLC1
carcinoma reference from GEO accession GSE103548 — and returns per-cell
`score_brain_struct`, `score_brain_immune`, `score_tumor`, plus the predicted label.
Tumor candidates are the cells passing the three filters in Table 4 of the thesis.
The raw labels, the three derived sets (Figures 4–7) and the score-threshold sweep
(Table 2) are drawn from these score tables:

```bash
conda run -n thesis_research python thesis_plots/make_fig_singler_sets.py      # Figures 4-7
conda run -n thesis_research python thesis_plots/make_table_singler_sweep.py   # Table 2
```

### Stage 2 — Supervised refinement of tumor candidates

```bash
conda run -n thesis_research python -m thesis_research.pipeline.cell_type_annotation.tumor_cells.refine_annotation_classifiers
```

Trains five classifiers (logistic regression, LogReg+PCA, LogReg+KNN on PCA, random
forest, XGBoost) on a reference pool whose negative class is drawn from tumor calls in
the sham-injected control slices. Reference-class construction lives in
`identify_tumor_cells.py` (`_get_tumor_ref_ids`, `_get_healthy_ref_ids`,
`_get_tumor_candidates_ids`).

### Stage 3 — Final tumor calls

XGBoost (library defaults) is refit on the joint reference pool and applied to all
candidates across the six slices; cells with *P*(tumor) > 0.5 are retained as refined
tumor. The runnable entry point is:

```bash
conda run -n thesis_research python thesis_plots/final_xgboost_refinement.py   # Figure 15
```

It writes the per-slice retained/rejected counts to `final_xgboost_refinement.csv`.
`thesis_plots/make_dq_fig_tumor_spatial.py` (Figure 21) refits the same model, so the
two figures show the same 20,873 tumor cells.

The final calls are then written to a new version of the per-slice cache, which the
Chapter 3 analyses use to exclude tumor cells (the earlier cache is left unmodified):

```bash
conda run -n thesis_research python thesis_plots/update_tumor_prediction_cache.py
```

### Stage 4 — Data-quality assessment

Probe detection and the lineage reporters:

```bash
conda run -n thesis_research python thesis_plots/rerun_chapter3_final_calls.py   # Table 5 values, Figure 16
conda run -n thesis_research python thesis_plots/make_dq_fig_reporter.py         # Figure 17
conda run -n thesis_research python thesis_plots/make_dq_fig_lyve1.py            # Figure 18
```

`rerun_chapter3_final_calls.py` runs `make_detection_reliability_6slice.py` and
`make_dq_fig1_detection.py` unchanged, with only their tumor-call input redirected to
`with_tumor_prediction_final/` and their outputs to `thesis_plots/rerun_final_tumor_calls/`.

Contamination correction with DecontX. All six slices use the same k-means partition
(30 truncated-SVD components, k = 25) as DecontX's population groups; the labels are
committed in `score_genes/decontx_clusters/` and the export step copies them into place:

```bash
conda run -n thesis_research python score_genes/run_decontx_correct.py export
Rscript score_genes/run_decontx.R resources/cache/decontx
conda run -n thesis_research python score_genes/run_decontx_correct.py assemble
conda run -n thesis_research python thesis_plots/make_dq_fig_decontx.py                # Figure 19
conda run -n thesis_research python thesis_plots/decontx_partition_sensitivity.py      # partition sensitivity, slice 1
```

Segmentation-error correction with FastReseg on all 257 fields of view of slice 1.
FastReseg's reference profiles come from the count matrix and the vendor's InSituType
cell groups prepared for an initial five-FOV run, so that run comes first; its
reference is identical to that of the slice-1 run, and its score matrix is exported
for Figure 20:

```bash
conda run -n thesis_research python agents/segmentation/fastreseg/01_prep_inputs.py
Rscript agents/segmentation/fastreseg/02_run_fastreseg.R
Rscript agents/segmentation/fastreseg/04_export_scores.R
conda run -n thesis_research python agents/segmentation/fastreseg/slice1_01_prep_inputs.py
Rscript agents/segmentation/fastreseg/slice1_02_run_fastreseg.R
conda run -n thesis_research python agents/segmentation/fastreseg/slice1_04_evaluate_qcpassed.py  # read-outs
conda run -n thesis_research python agents/segmentation/fastreseg/slice1_05_sparsity.py           # sparsity
conda run -n thesis_research python agents/segmentation/fastreseg/slice1_06_figure20.py           # Figure 20
```

Together these quantify each probe against the noise floor defined by the 11
negative-control probes, test the lineage reporters against their expected biology,
and test whether contamination and segmentation-error correction resolve the anomalies.

### Stage 5 — Thesis figures and tables

```bash
conda run -n thesis_research python thesis_plots/figure_3_model_comparison.py    # Figure 8
conda run -n thesis_research python thesis_plots/xgb_default_sensitivity.py      # XGBoost defaults check
conda run -n thesis_research python thesis_plots/control_specificity_validation.py # Table 3
conda run -n thesis_research python thesis_plots/figure_4_spatial_refinement.py  # Figures 9-14
conda run -n thesis_research python thesis_plots/stage1_threshold_sensitivity.py # threshold calibration
conda run -n thesis_research python thesis_plots/make_dq_fig_tumor_spatial.py    # Figure 21
conda run -n thesis_research python thesis_plots/make_nature_tables.py           # Tables 1 and 4
```

The Figures 9-14 and 21 scripts import from each other (`figure_4_spatial_refinement`,
`final_xgboost_refinement`) and from the `thesis_research` package, so run them from
the repository root with it on `PYTHONPATH`.

> **Note on numbering.** The `figure_N_*.py` filenames predate the thesis figure
> numbering and are offset from it. The mapping above and in Supplementary Table S1
> is authoritative; the filenames are not.

## 4. Repository layout

| Path | Contents |
|---|---|
| `thesis_research/pipeline/` | Ingest, QC, FOV/position plots, clustering, filters |
| `thesis_research/pipeline/cell_type_annotation/tumor_cells/` | SingleR filtering, classifier comparison, final tumor calls |
| `thesis_research/pipeline/cell_type_annotation/myeloid/` | **Superseded** reporter-*gated* myeloid typing (GFP/tdTomato → MDM vs resident) |
| `thesis_research/utils/` | Path constants, column names, entity types |
| `thesis_plots/` | Scripts generating every thesis figure and table |
| `thesis_plots/tables/` | Rendered Nature-style tables (HTML for Word, Markdown for drafts) |
| `outputs/cell_annotation/` | SingleR reference construction and score conversion (R) |
| `score_genes/` | Reporter-independent cell-type annotation: `score_genes` scoring with mirrored-FDR gating (`run_score_genes_*.py`), myeloid subtype gating (`myeloid_subtype_gate.py`, `myeloid_stage2*.py`), backbone classification (`classify_backbone.py`), final labels (`final_annotation.py`), plus decontX (`run_decontx.R`) and two-tier SingleR validation (`singler_two_tier.R`) |
| `agents/outputs/` | Intermediate validation reports and metrics |

Scripts prefixed with `_` at the repository root are exploratory one-offs kept for
provenance. **They are not part of the reproducible pipeline** and no thesis result
depends on them.

## 5. Known gaps

These are open items, listed so that the state of the repository is not overstated.
Each is also flagged in Supplementary Table S1.

- **No R lock file.** R package versions are listed in section 1 but there is no
  `renv.lock`, so the R environment is not restored automatically.
- **Pipeline outputs have no stable path.** `run_pipeline.py` generates a fresh
  `uuid4()` as `run_id` on every run, so Stage 0 figures land in a new
  `outputs/<run_id>/<sample_id>/` directory each time and cannot be cited or diffed.
  The pipeline should accept a fixed run identifier.
- **Tables 1 and 4** have their values hardcoded in `make_nature_tables.py` rather than
  read from the analysis outputs. The detection-strength table that script also renders
  holds values computed on the earlier tumor calls; thesis Table 5 is taken from
  `rerun_final_tumor_calls/detection_reliability_all6.csv` instead.
- **Two tumor-call caches.** `resources/cache/with_tumor_prediction/` holds the earlier
  tuned-XGBoost calls (20,679 tumor cells). The Chapter 3 scripts that read it are run
  through `rerun_chapter3_final_calls.py` against `with_tumor_prediction_final/`
  (20,873, the calls reported in the thesis); run directly, they reproduce the earlier
  non-tumor count (n = 825,428 instead of 825,234).
- **Absolute paths.** Several `thesis_plots/` scripts hardcode `D:/thesis-research/`
  instead of resolving from `constants.py`, so they will not run elsewhere unchanged.
  The FastReseg scripts additionally hardcode the raw transcript table and the
  `D:/R-libs/fastreseg` R library.
- **Figure numbering** in `figure_N_*.py` filenames does not match the thesis.
- **`myeloid/myeloid_typing.py` contradicts the thesis.** Its stage 1 gates on
  GFP/tdTomato, but the thesis concludes the reporters are unusable and states that
  annotation was reporter-independent. The reporter-independent framework used for the
  thesis lives in `score_genes/`. This module is kept for provenance and should be
  marked superseded in-file, or removed, so the two do not appear to disagree.
- **`score_genes/` is undocumented.** It holds the reporter-independent annotation
  framework in ~120 loosely-named scripts with no entry point or ordering. The
  cell-annotation chapter is not yet written, so no thesis display item depends on it
  today; once that chapter exists, those scripts need the same S1 treatment as the rest.

## 6. Verifying a run

After Stage 0 the headline numbers should match Table 1 of the thesis exactly:

| Quantity | Expected |
|---|---|
| Segmented cells before QC | 926,318 |
| Removed by low-transcript filter | 80,211 |
| Retained | 846,107 (91.34%) |
| FOVs | 1,205 |
| Per-slice retention | 88.86 / 91.88 / 92.74 / 74.49 / 93.07 / 95.20 % |

After Stage 2, out-of-fold ROC-AUC should exceed 0.99 for all five classifiers, and
XGBoost should reach accuracy 0.970, precision 0.985, recall 0.963, F1 0.974, with
tumor as the positive class (`figure_3_model_comparison.csv`).

After Stage 3, both control slices (3 and 4) should retain **zero** tumor cells out of
343 and 287 candidates respectively — the sharpest single check that the pipeline
behaves as designed — and the four tumor-bearing slices should retain 20,873 tumor
cells in total (4,267 / 6,272 / 5,039 / 5,295). These zeros are in-sample: the control
candidates are the classifier's negative training class. On look-alikes withheld from
training (Table 3), XGBoost should call 13 of the 630 tumor out-of-fold, and 30 of 630
when each control slice is withheld in turn.

## 7. Data and code availability

Processed per-slice single-cell data, the assembled SingleR reference objects and the
raw CosMx exports (transcript tables and segmentation masks) are available from the
author on request. The cluster labels supplied to DecontX are included in
`score_genes/decontx_clusters/`.

All animal procedures, tumor implantation, tissue processing, and histological
preparation were performed by Mr. Avinoam Ratzabi, Tel Aviv University, in accordance
with institutional and ethical guidelines.
