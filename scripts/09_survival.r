source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(dplyr); library(survival); library(survminer); library(broom)
})

surv <- readRDS(file.path(data_dir, "surv.rds"))

km <- survfit(Surv(time, event) ~ TP53, data = surv)
p_km <- ggsurvplot(km, data = surv, pval = TRUE, risk.table = TRUE,
                   xlab = "Years", ylab = paste(ENDPOINT, "probability"),
                   legend.labs = c("wt", "mut"), legend.title = "TP53")
png(file.path(fig_dir, paste0("km_tp53_", ENDPOINT, ".png")), width = 1800, height = 1600, res = 220)
print(p_km); dev.off()

d <- surv |> filter(!is.na(age), !is.na(stage), !is.na(subtype))
cat("Multivariable sample:", nrow(d), "| Events:", sum(d$event), "\n")

cox_uni   <- coxph(Surv(time, event) ~ TP53, data = d)
cox_multi <- coxph(Surv(time, event) ~ TP53 + age + stage + subtype, data = d)

zph <- cox.zph(cox_multi)
print(zph)
png(file.path(fig_dir, paste0("schoenfeld_", ENDPOINT, ".png")), width = 2000, height = 1600, res = 200)
print(ggcoxzph(zph)); dev.off()

cox_strat <- coxph(Surv(time, event) ~ TP53 + age + stage + strata(subtype), data = d)

cox_tab <- bind_rows(
  univariate   = tidy(cox_uni,   exponentiate = TRUE, conf.int = TRUE),
  multivariate = tidy(cox_multi, exponentiate = TRUE, conf.int = TRUE),
  stratified   = tidy(cox_strat, exponentiate = TRUE, conf.int = TRUE),
  .id = "model"
)
write.csv(cox_tab, file.path(res_dir, paste0("cox_", ENDPOINT, ".csv")), row.names = FALSE)
print(cox_tab |> filter(term == "TP53mut"))