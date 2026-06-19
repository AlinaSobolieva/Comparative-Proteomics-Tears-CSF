# ==============================================================================
# Bioinformatics Project
# Comparative Proteomics of Human Biofluids: Tears, Cerebrospinal Fluid and Plasma
#
# Authors: A. Akrapova, A. Sobolieva, O. Vlasiuk, Kh. Malysheva, H. Svitina
# ==============================================================================
# Data Processing and Quality Assessment
# ==============================================================================

# Load libraries
library(dplyr)
library(tidyr)
library(gt)

# Ensure data directory exists
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# Tears DDA --------------------------------------------------------------------
if (file.exists("protein.tsv")) {
  protein <- read.delim("protein.tsv")
  write.csv(
    protein,
    "data/Tears_DDA_protein_level_matrix.csv",
    row.names = FALSE
  )
} else {
  warning("protein.tsv not found in project directory.")
}

# CSF DIA ----------------------------------------------------------------------
if (file.exists("CSF_DIA_Quant.tsv")) {
  csf <- read.delim("CSF_DIA_Quant.tsv")
  
  csf_protein <- csf %>%
    mutate(
      Protein = sub(";.*", "", PG.ProteinAccessions)
    ) %>%
    group_by(
      Protein,
      R.Replicate
    ) %>%
    summarise(
      Quantity = sum(PG.Quantity, na.rm = TRUE),
      .groups = "drop"
    )
  
  csf_matrix <- csf_protein %>%
    pivot_wider(
      names_from = R.Replicate,
      values_from = Quantity
    )
  
  write.csv(
    csf_matrix,
    "data/CSF_protein_level_matrix.csv",
    row.names = FALSE
  )
} else {
  warning("CSF_DIA_Quant.tsv not found in project directory.")
}

# Depleted Plasma DIA ----------------------------------------------------------
if (file.exists("Depleted_plasma_report.tsv")) {
  plasma <- read.delim("Depleted_plasma_report.tsv")
  
  plasma_protein <- plasma %>%
    mutate(
      Protein = sub(";.*", "", PG.ProteinAccessions)
    ) %>%
    group_by(
      Protein,
      R.Replicate
    ) %>%
    summarise(
      Quantity = sum(PG.Quantity, na.rm = TRUE),
      .groups = "drop"
    )
  
  plasma_matrix <- plasma_protein %>%
    pivot_wider(
      names_from = R.Replicate,
      values_from = Quantity
    )
  
  write.csv(
    plasma_matrix,
    "data/Plasma_protein_level_matrix.csv",
    row.names = FALSE
  )
} else {
  warning("Depleted_plasma_report.tsv not found in project directory.")
}

# Tears DIA --------------------------------------------------------------------
if (file.exists("Report_training_data.pg_matrix.tsv")) {
  tears <- read.delim("Report_training_data.pg_matrix.tsv")
  
  write.csv(
    tears,
    "data/Tears_DIA_protein_level_matrix.csv",
    row.names = FALSE
  )
} else {
  warning("Report_training_data.pg_matrix.tsv not found in project directory.")
}

# Dataset Summary --------------------------------------------------------------
dataset_summary <- data.frame(
  Dataset = c(
    "Tears_DIA",
    "Tears_DDA",
    "CSF_DIA",
    "Depleted_plasma_DIA"
  ),
  Proteins = c(
    1478,
    323,
    1409,
    614
  ),
  Replicates = c(
    8,
    8,
    222,
    289
  ),
  Missing_values_percent = c(
    35.90,
    3.85,
    5.72,
    4.25
  ),
  Protein_level_format = c(
    "Already protein-level",
    "Already protein-level",
    "Converted from raw report",
    "Converted from raw report"
  )
)

tab <- gt(dataset_summary)

# Try saving as image if webshot2 package is available
tryCatch({
  gtsave(tab, "results/Dataset_summary.png")
}, error = function(e) {
  warning("gtsave failed. Make sure webshot2 package is installed to save tables as images.")
})
