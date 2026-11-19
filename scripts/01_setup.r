# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install packages
cran_pkgs <- c(
  "BiocManager",
  "tidyverse",
  "survival",
  "survminer",
  "broom",
  "readxl",
  "msigdbr"
)
bioc_pkgs <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "BiocParallel",
  "fgsea",
  "clusterProfiler",
  "org.Hs.eg.db"
)

missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)


# Download data (skip if already prepared)
se_file <- "data/brca_se.rds"
if (file.exists(se_file)) {
  message(se_file, " already exists. Skipping download.")
} else {
  library(TCGAbiolinks)
  library(SummarizedExperiment)

  query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  GDCdownload(query, files.per.chunk = 50)
  se <- GDCprepare(query)
  dir.create("data", showWarnings = FALSE)
  saveRDS(se, se_file)
}
