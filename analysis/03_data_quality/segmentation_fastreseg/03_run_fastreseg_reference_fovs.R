# FastReseg (Wu et al. 2025) on the five slice-1 FOVs, with library defaults.
#
# Inputs come from 02_prepare_reference_inputs.py. Only data-specific settings are given:
# the pixel size (0.12028 um, from the vendor's Area.um2 / Area), the column
# names of the transcript files, the cell ID of extracellular transcripts
# (CellId 0 in each FOV) and a fixed seed. Coordinates are global pixels, so
# the per-FOV offsets are zero and the y axis is not inverted. All other
# arguments are FastReseg 1.1.2 defaults.
#
# FastReseg and its dependencies are installed in a separate library on D:
# (D:/R-libs/fastreseg). Writes only to
# agents/outputs/segmentation_fastreseg/fastreseg_out/, which must not exist.
#
# Run: Rscript analysis/03_data_quality/segmentation_fastreseg/03_run_fastreseg_reference_fovs.R
PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")

.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))
suppressPackageStartupMessages({
  library(FastReseg)
  library(Matrix)
  library(data.table)
})

base <- paste0(PROJECT_ROOT, "/agents/outputs/segmentation_fastreseg")
inp <- file.path(base, "inputs")
out <- file.path(base, "fastreseg_out")
if (dir.exists(out)) stop(out, " already exists; refusing to overwrite")

counts <- as(readMM(file.path(inp, "counts.mtx")), "CsparseMatrix")
rownames(counts) <- fread(file.path(inp, "cells.csv"))$cell
colnames(counts) <- fread(file.path(inp, "genes.csv"))$gene
cl <- fread(file.path(inp, "clust.csv"))
clust <- setNames(cl$clust, cl$cell)
stopifnot(identical(names(clust), rownames(counts)))

fileInfo <- as.data.frame(fread(file.path(inp, "transDF_fileInfo.csv")))
cat("counts:", nrow(counts), "cells x", ncol(counts), "genes;",
    length(unique(clust)), "clusters;", nrow(fileInfo), "FOVs\n")

t0 <- Sys.time()
res <- fastReseg_full_pipeline(
  counts = counts,
  clust = clust,
  transDF_fileInfo = fileInfo,
  filepath_coln = "file_path",
  prefix_colns = NULL,                      # cell IDs are already unique (c_1_<fov>_<id>)
  fovOffset_colns = c("stage_X", "stage_Y"),
  pixel_size = 0.12028,
  transID_coln = "transcript_id",
  transGene_coln = "target",
  cellID_coln = "cell",
  spatLocs_colns = c("x_global_px", "y_global_px", "z"),
  invert_y = FALSE,
  extracellular_cellID = paste0("c_1_", fileInfo$fov, "_0"),
  path_to_output = out,
  seed_process = 42
)
cat("done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(res, file.path(out, "fastreseg_result.rds"))
cat("result components:", paste(names(res), collapse = ", "), "\n")
