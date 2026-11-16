# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "TCGAbiolinks",
  "DESeq2",
  "SummarizedExperiment"
))

# Install CRAN packages
install.packages(c(
  "tidyverse",
  "survminer"
))


# Download data
library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = "Primary Tumor"
)
GDCdownload(query, files.per.chunk = 50)
se <- GDCprepare(query)
dir.create("data", showWarnings = FALSE)
saveRDS(se, "data/brca_se.rds")