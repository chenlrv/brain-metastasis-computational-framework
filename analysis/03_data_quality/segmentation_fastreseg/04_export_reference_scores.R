# Export FastReseg's transcript-score matrix (genes x cell types), computed with
# FastReseg's own scoreGenesInRef from the reference profiles of the run in
# 03_run_fastreseg_reference_fovs.R, for plotting. Writes only to
# automation/outputs/segmentation_fastreseg/figures/ (the file must not exist).
PROJECT_ROOT <- Sys.getenv("THESIS_PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("Set THESIS_PROJECT_ROOT to the project data root (see README)")
.libPaths(c(Sys.getenv("FASTRESEG_RLIB"), .libPaths()))
suppressMessages(library(FastReseg))
base <- paste0(PROJECT_ROOT, "/automation/outputs/segmentation_fastreseg")
out <- file.path(base, "figures")
dir.create(out, showWarnings = FALSE)
f <- file.path(out, "score_matrix.csv")
if (file.exists(f)) stop(f, " exists; refusing to overwrite")
r <- readRDS(file.path(base, "fastreseg_out", "fastreseg_result.rds"))
S <- as.matrix(scoreGenesInRef(genes = rownames(r$refProfiles), ref_profiles = as.matrix(r$refProfiles)))
write.csv(S, f)
cat("wrote", f, ":", nrow(S), "genes x", ncol(S), "cell types\n")
