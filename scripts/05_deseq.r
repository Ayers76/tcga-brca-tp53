source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
})

dds <- readRDS(file.path(data_dir, "dds.rds")) |> subsample()
register(SnowParam(N_CORES))

# Model 1: ~ TP53
dds1 <- dds; design(dds1) <- ~ TP53
dds1 <- DESeq(dds1, parallel = TRUE)
res1 <- results(dds1, contrast = c("TP53", "mut", "wt"), alpha = 0.05)

# Model 2: ~ subtype + TP53 + subtype:TP53 (design from 03)
dds2 <- DESeq(dds, parallel = TRUE)
res2 <- results(dds2, contrast = c("TP53", "mut", "wt"), alpha = 0.05)

to_df <- function(res, dds) {
  df <- as.data.frame(res)
  df$id <- rownames(df)
  df$gene <- rowData(dds)[df$id, "gene_name"]
  df[, c("id", "gene", setdiff(names(df), c("id", "gene")))]
}

write.csv(to_df(res1, dds1), tagged(res_dir, "res1", "csv"), row.names = FALSE)
saveRDS(dds1, tagged(data_dir, "dds1", "rds"))

write.csv(to_df(res2, dds2), tagged(res_dir, "res2", "csv"), row.names = FALSE)
saveRDS(dds2, tagged(data_dir, "dds2", "rds"))
