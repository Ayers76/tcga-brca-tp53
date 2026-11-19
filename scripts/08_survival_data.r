source("scripts/00_config.r")
suppressPackageStartupMessages({ library(dplyr) })

meta <- readRDS(file.path(data_dir, "meta.rds"))

collapse_stage <- function(s) {
  s <- sub("^Stage ", "", s)
  s <- sub("[ABC]$", "", s)
  ifelse(s %in% c("I", "II", "III", "IV"), s, NA)
}

surv <- meta |>
  transmute(patient,
    TP53    = factor(TP53, levels = c("wt", "mut")),
    subtype = relevel(factor(subtype), ref = "LumA"),
    age     = age_at_index,
    stage   = factor(collapse_stage(ajcc_pathologic_stage), levels = c("I", "II", "III", "IV")),
    OS      = as.integer(vital_status == "Dead"),
    OS.time = ifelse(OS == 1, days_to_death, days_to_last_follow_up)
  )

if (ENDPOINT == "PFI") {
  cdr_file <- file.path(data_dir, "TCGA-CDR-SupplementalTableS1.xlsx")
  if (!file.exists(cdr_file)) stop("File not found: ", cdr_file, ". Download it, or set ENDPOINT to \"OS\" in 00_config.r")
  cdr <- readxl::read_excel(cdr_file, sheet = 1) |>
    filter(type == "BRCA") |>
    transmute(patient  = bcr_patient_barcode,
              PFI      = as.integer(PFI),
              PFI.time = as.numeric(PFI.time))
  surv <- left_join(surv, cdr, by = "patient")
}

surv <- surv |>
  mutate(time  = .data[[paste0(ENDPOINT, ".time")]] / 365.25,
         event = .data[[ENDPOINT]]) |>
  filter(!is.na(time), time > 0, !is.na(event))

saveRDS(surv, file.path(data_dir, "surv.rds"))

cat(ENDPOINT, "| Patients:", nrow(surv), "| Events:", sum(surv$event), "\n")
print(table(surv$TP53, surv$event, dnn = c("TP53", "event")))