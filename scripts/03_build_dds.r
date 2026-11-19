source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(SummarizedExperiment)
})

# Load data
se <- readRDS(file.path(data_dir, "brca_se.rds"))
meta <- readRDS(file.path(data_dir, "meta.rds")) |> filter(!is.na(subtype))

# Define categorical variables
meta$TP53 <- factor(meta$TP53, levels = c("wt", "mut"))
meta$subtype <- factor(meta$subtype)

# Extract protein-coding gene expression data
se2 <- se[rowData(se)$gene_type == "protein_coding", meta$barcode]
counts <- assay(se2, "unstranded")

# Model and Filter
dds <- DESeqDataSetFromMatrix(counts, colData = meta,
  design = ~ subtype + TP53 + subtype:TP53)
rowData(dds) <- rowData(se2)   # keep gene annotation (gene_name, gene_type, ...)
dds <- dds[rowSums(counts(dds) >= 10) >= 50, ]

# Save the FULL dds; downstream scripts call subsample() according to DEBUG
saveRDS(dds, file.path(data_dir, "dds.rds"))
