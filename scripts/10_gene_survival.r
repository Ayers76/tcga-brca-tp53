source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(survival); library(broom)
})

dds_file <- file.path(data_dir, "dds.rds")
dds  <- readRDS(dds_file) |> subsample()
r2   <- read.csv(tagged(res_dir, "res2", "csv"))
surv <- readRDS(file.path(data_dir, "surv.rds")) |>
  filter(!is.na(age), !is.na(stage), !is.na(subtype))

# Cache vst; rebuild if dds.rds is newer than the cache
vsd_file <- tagged(data_dir, "vsd", "rds")
if (file.exists(vsd_file) && file.mtime(vsd_file) > file.mtime(dds_file)) {
  vsd <- readRDS(vsd_file)
} else {
  vsd <- vst(dds, blind = FALSE)
  saveRDS(vsd, vsd_file)
}

top <- r2 |> filter(!is.na(padj)) |> arrange(padj) |> slice_head(n = 20)
patients <- substr(colnames(vsd), 1, 12)

gene_cox <- lapply(seq_len(nrow(top)), function(i) {
  df <- data.frame(patient = patients,
                   z = as.numeric(scale(assay(vsd)[top$id[i], ]))) |>
    inner_join(surv, by = "patient")
  fit <- coxph(Surv(time, event) ~ z + TP53 + age + stage + strata(subtype), data = df)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "z") |>
    mutate(id = top$id[i], gene = top$gene[i])
}) |>
  bind_rows() |>
  mutate(padj = p.adjust(p.value, method = "BH")) |>
  arrange(p.value)

write.csv(gene_cox, tagged(res_dir, paste0("gene_cox_", ENDPOINT), "csv"), row.names = FALSE)
print(gene_cox |> select(gene, estimate, conf.low, conf.high, p.value, padj))
