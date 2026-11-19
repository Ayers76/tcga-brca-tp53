source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

dds <- readRDS(file.path(data_dir, "dds.rds")) |> subsample()

# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# PCA plot colored by TP53 mutation status
pca_tp53 <- plotPCA(vsd, intgroup = "TP53")
ggsave(tagged(fig_dir, "PCA_TP53", "png"), plot = pca_tp53, width = 8, height = 6)

# PCA plot colored by molecular subtype
pca_subtype <- plotPCA(vsd, intgroup = "subtype")
ggsave(tagged(fig_dir, "PCA_Subtype", "png"), plot = pca_subtype, width = 8, height = 6)
