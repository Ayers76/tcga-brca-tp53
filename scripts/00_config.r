DEBUG <- TRUE # TRUE: 200 samples；FALSE: all
N_DEBUG <- 200
SEED <- 76
mode_tag <- if (DEBUG) "DEBUG" else "FULL"

N_CORES <- 8

ENDPOINT <- "OS"  # "OS" or "PFI" (PFI needs data/TCGA-CDR-SupplementalTableS1.xlsx)

data_dir <- "data"
res_dir <- file.path("results", "summaries")
fig_dir <- file.path("results", "figures")
for (d in c(data_dir, res_dir, fig_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

tagged <- function(dir, name, ext) file.path(dir, paste0(name, "_", mode_tag, ".", ext))

subsample <- function(dds) {
  if (!DEBUG) return(dds)
  set.seed(SEED)
  dds[, sample(ncol(dds), N_DEBUG)]
}