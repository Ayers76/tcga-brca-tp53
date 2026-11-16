library(SummarizedExperiment)
library(TCGAbiolinks)

se <- readRDS("data/brca_se.rds")
dim(se)
assayNames(se)
counts <- assay(se, "unstranded")
counts[1:5, 1:3]
head(rowData(se)$gene_name)
names(colData(se))