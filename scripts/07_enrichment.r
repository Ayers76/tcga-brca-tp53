source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
  library(fgsea); library(msigdbr)
  library(clusterProfiler); library(org.Hs.eg.db)
})

r2 <- read.csv(tagged(res_dir, "res2", "csv"))

ranks_df <- r2 |>
  filter(!is.na(stat), !is.na(gene), gene != "") |>
  group_by(gene) |>
  slice_max(abs(stat), n = 1, with_ties = FALSE) |>
  ungroup()
ranks <- setNames(ranks_df$stat, ranks_df$gene)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

set.seed(SEED)
gsea <- fgsea(pathways, ranks, minSize = 15, maxSize = 500) |> arrange(padj)

gsea_out <- gsea |> mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(gsea_out, tagged(res_dir, "gsea_hallmark", "csv"), row.names = FALSE)

top <- gsea |> filter(padj < 0.05) |> slice_max(abs(NES), n = 20)
p_gsea <- ggplot(top, aes(NES, reorder(pathway, NES), fill = NES > 0)) +
  geom_col() +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                    labels = c(`FALSE` = "Higher in wt", `TRUE` = "Higher in mut")) +
  labs(x = "Normalized enrichment score", y = NULL, fill = NULL)
ggsave(tagged(fig_dir, "gsea_hallmark", "png"), p_gsea, width = 8, height = 6)

universe <- unique(na.omit(r2$gene[!is.na(r2$padj)]))
up   <- r2 |> filter(!is.na(padj), padj < 0.05, log2FoldChange >  1) |> pull(gene) |> unique()
down <- r2 |> filter(!is.na(padj), padj < 0.05, log2FoldChange < -1) |> pull(gene) |> unique()

run_go <- function(genes, label) {
  go <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 universe = universe, pAdjustMethod = "BH")
  if (is.null(go) || nrow(as.data.frame(go)) == 0) {
    message(label, ": no significant GO terms"); return(invisible(NULL))
  }
  write.csv(as.data.frame(go), tagged(res_dir, paste0("go_", label), "csv"), row.names = FALSE)
  ggsave(tagged(fig_dir, paste0("go_", label), "png"),
         dotplot(go, showCategory = 15), width = 8, height = 6)
}
run_go(up, "up"); run_go(down, "down")