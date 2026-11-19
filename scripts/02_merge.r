source("scripts/00_config.r")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
})

# Query masked somatic mutation data from the TCGA-BRCA project
q_maf <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

# Download mutation data only if not already downloaded
# (GDCdownload saves to <directory>/TCGA-BRCA/Simple_Nucleotide_Variation/...)
maf_dir <- file.path(data_dir, "TCGA-BRCA", "Simple_Nucleotide_Variation")
if (!dir.exists(maf_dir)) {
  GDCdownload(q_maf, method = "api", directory = data_dir, files.per.chunk = 10)
} else {
  message("Mutation data already exists. Skipping download.")
}

# Prepare mutation data from the specified directory
maf <- GDCprepare(q_maf, directory = data_dir)

# Extract the 12-character patient identifier from the sample barcode
# Multiple samples from the same patient share the same patient identifier
maf <- maf |> mutate(patient = substr(Tumor_Sample_Barcode, 1, 12))

# Define mutation classifications considered silent or non-protein-altering
silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA")

# Identify patients with non-silent mutations in the TP53 gene
tp53_mut <- maf |>
  filter(Hugo_Symbol == "TP53", !Variant_Classification %in% silent) |>
  distinct(patient) |>
  pull(patient)

# Get the unique patients represented in the mutation dataset
sequenced <- unique(maf$patient)

# Load the previously processed BRCA SummarizedExperiment object
se <- readRDS(file.path(data_dir, "brca_se.rds"))

# Construct patient-level metadata
meta <- as.data.frame(colData(se)) |>
  mutate(patient = substr(barcode, 1, 12)) |>
  filter(patient %in% sequenced) |>             # 没测突变的病人去掉
  distinct(patient, .keep_all = TRUE) |>        # 一个病人只留一个样本
  mutate(TP53 = ifelse(patient %in% tp53_mut, "mut", "wt"),
         subtype = paper_BRCA_Subtype_PAM50)

# Contingency table of TP53 mutation status and molecular subtype
print(table(meta$TP53, meta$subtype, useNA = "ifany"))

saveRDS(meta, file.path(data_dir, "meta.rds"))
