q_maf <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  access        = "open"
)

GDCdownload(q_maf)
maf <- GDCprepare(q_maf)

library(dplyr)
maf <- maf |> mutate(patient = substr(Tumor_Sample_Barcode, 1, 12))

silent <- c("Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA")
tp53_mut <- maf |>
  filter(Hugo_Symbol == "TP53",
  !Variant_Classification %in% silent
  ) |>
  distinct(patient) |> 
  pull(patient)

sequenced <- unique(maf$patient)