source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(ggplot2)
})

r1 <- read.csv(tagged(res_dir, "res1", "csv"))
r2 <- read.csv(tagged(res_dir, "res2", "csv"))
dds2 <- readRDS(tagged(data_dir, "dds2", "rds"))

mark_sig <- function(df) df |>
  mutate(sig = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
r1 <- mark_sig(r1); r2 <- mark_sig(r2)

cat("Mode:", mode_tag, "\n")
cat("Significant genes (Model 1):", sum(r1$sig), "\n")
cat("Significant genes (Model 2):", sum(r2$sig), "\n")
cat("Overlap:", length(intersect(r1$id[r1$sig], r2$id[r2$sig])), "\n")

both <- inner_join(r1, r2, by = c("id", "gene"), suffix = c("_m1", "_m2"))
p_cmp <- ggplot(both, aes(log2FoldChange_m1, log2FoldChange_m2)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(linetype = "dashed") +
  labs(x = "log2FC, ~ TP53", y = "log2FC, ~ subtype + TP53 + subtype:TP53")
ggsave(tagged(fig_dir, "lfc_compare", "png"), p_cmp, width = 5, height = 5)

p_vol <- ggplot(r2, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8)
ggsave(tagged(fig_dir, "volcano", "png"), p_vol, width = 6, height = 5)

for (g in c("MDM2", "CDKN1A", "MKI67")) {
  id <- r2$id[r2$gene == g][1]
  if (is.na(id)) { message("Gene not found after filtering: ", g); next }
  png(tagged(fig_dir, paste0("counts_", g), "png"))
  plotCounts(dds2, gene = id, intgroup = "TP53", main = g)
  dev.off()
}

# plotCounts(dds2, gene = r2$id[r2$gene == "MDM2"], intgroup = "TP53")