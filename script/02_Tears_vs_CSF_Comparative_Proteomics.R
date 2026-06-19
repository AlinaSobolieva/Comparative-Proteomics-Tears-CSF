# ==============================================================================
# Comparative Proteomics of Tears and CSF
#
# Author: Alina Sobolieva
#
# Objective:
# To compare the proteomic profiles of tear fluid and cerebrospinal fluid 
# and evaluate the overlap between the two biofluids.
# ==============================================================================

# 1. Setup & Directories -------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggVennDiagram)
library(pheatmap)

# Ensure output directory exists (prevents crashes if results/ doesn't exist)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# 2. Load Data -----------------------------------------------------------------
tears <- read.csv("data/Tears_DIA_protein_level_matrix.csv")
csf   <- read.csv("data/CSF_protein_level_matrix.csv")

# Clean protein group names (extract first UniProt ID)
tears$Protein <- sub(";.*", "", tears$Protein.Group)

# 3. Protein Overlap Analysis --------------------------------------------------
tears_proteins <- unique(tears$Protein)
csf_proteins   <- unique(csf$Protein)

shared_proteins <- intersect(tears_proteins, csf_proteins)
tears_specific  <- setdiff(tears_proteins, csf_proteins)
csf_specific    <- setdiff(csf_proteins, tears_proteins)

# Save lists to CSV
write.csv(data.frame(Protein = shared_proteins), "results/shared_proteins.csv", row.names = FALSE)
write.csv(data.frame(Protein = tears_specific),  "results/tears_specific_proteins.csv", row.names = FALSE)
write.csv(data.frame(Protein = csf_specific),    "results/csf_specific_proteins.csv", row.names = FALSE)

# 4. Venn Diagram --------------------------------------------------------------
venn_list <- list(
  Tears = tears_proteins,
  CSF   = csf_proteins
)

venn_plot <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "#E9D5FF", high = "#6B21A8") +
  labs(
    title = "Protein Overlap Between Tears and CSF",
    fill  = "Proteins"
  ) +
  theme_minimal()

ggsave("results/Project2_Venn_Diagram.png", venn_plot, width = 8, height = 7, dpi = 600)

# 5. Concordance Plot (Fixed Alignment Bug) ------------------------------------

# Calculate means independently, keeping the Protein column as a key
tears_mean_df <- tears %>%
  filter(Protein %in% shared_proteins) %>%
  # Select sample columns dynamically (starts with Users document path or Pool_DIA)
  select(Protein, starts_with("C..Users"), starts_with("Pool_DIA")) %>% 
  rowwise() %>%
  mutate(Tears_Mean = mean(c_across(-Protein), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Protein, Tears = Tears_Mean)

csf_mean_df <- csf %>%
  filter(Protein %in% shared_proteins) %>%
  # Select everything except the Protein column
  rowwise() %>%
  mutate(CSF_Mean = mean(c_across(-Protein), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Protein, CSF = CSF_Mean)

# Join the tables by Protein to guarantee perfect row alignment
concordance_df <- inner_join(tears_mean_df, csf_mean_df, by = "Protein") %>%
  filter(Tears > 0, CSF > 0)

correlation <- cor(
  log10(concordance_df$Tears),
  log10(concordance_df$CSF),
  use = "complete.obs"
)

concordance_plot <- ggplot(concordance_df, aes(x = log10(Tears), y = log10(CSF))) +
  geom_point(color = "#7E22CE", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#4C1D95", se = TRUE) +
  labs(
    title    = "Concordance of Shared Proteins",
    subtitle = paste("Pearson r =", round(correlation, 3)),
    x        = "Mean abundance in Tears (log10)",
    y        = "Mean abundance in CSF (log10)"
  ) +
  theme_minimal()

ggsave("results/Project2_Concordance_Plot.png", concordance_plot, width = 8, height = 6, dpi = 600)

# 6. Missingness Heatmap (Fixed Alignment Bug) ---------------------------------

# Extract presence flags keeping the Protein column
tears_presence_df <- tears %>%
  filter(Protein %in% shared_proteins) %>%
  select(Protein, starts_with("C..Users"), starts_with("Pool_DIA")) %>%
  mutate(across(-Protein, ~ ifelse(. > 0, 1, 0)))

csf_presence_df <- csf %>%
  filter(Protein %in% shared_proteins) %>%
  select(Protein, everything()) %>%
  mutate(across(-Protein, ~ ifelse(. > 0, 1, 0)))

# Join presence dataframes together by Protein to ensure alignment
combined_presence_df <- inner_join(tears_presence_df, csf_presence_df, by = "Protein")

# Convert to matrix and format row names
presence_matrix <- as.matrix(combined_presence_df %>% select(-Protein))
rownames(presence_matrix) <- combined_presence_df$Protein

png("results/Project2_Missingness_Heatmap.png", width = 3000, height = 2000, res = 300)
pheatmap(
  presence_matrix,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = c("#E9D5FF", "#6B21A8"),
  main          = "Missingness of Shared Proteins"
)
dev.off()

# 7. Top 20 Proteins -----------------------------------------------------------
top20_tears <- tears %>%
  select(Protein, starts_with("C..Users"), starts_with("Pool_DIA")) %>%
  rowwise() %>%
  mutate(Mean_Abundance = mean(c_across(-Protein), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Protein, Mean_Abundance) %>%
  arrange(desc(Mean_Abundance)) %>%
  slice(1:20)

top20_csf <- csf %>%
  rowwise() %>%
  mutate(Mean_Abundance = mean(c_across(-Protein), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Protein, Mean_Abundance) %>%
  arrange(desc(Mean_Abundance)) %>%
  slice(1:20)

write.csv(top20_tears, "results/Top20_Tears.csv", row.names = FALSE)
write.csv(top20_csf,   "results/Top20_CSF.csv",   row.names = FALSE)

# ==================================================
# Principal Component Analysis (PCA)
#
# Objective:
# Evaluate global proteomic differences between
# Tears and CSF using shared proteins.
# ==================================================

# Select shared proteins
tears_shared <- tears %>%
  filter(Protein %in% shared_proteins)

csf_shared <- csf %>%
  filter(Protein %in% shared_proteins)

# Expression matrices
tears_expr <- tears_shared[, 7:26]
csf_expr   <- csf_shared[, -1]

rownames(tears_expr) <- tears_shared$Protein
rownames(csf_expr)   <- csf_shared$Protein

# Combine datasets
combined_expr <- cbind(tears_expr, csf_expr)

# Log transformation and missing value imputation
combined_expr <- log10(combined_expr + 1)
combined_expr[is.na(combined_expr)] <- 0

# PCA
pca <- prcomp(
  t(combined_expr),
  scale. = TRUE
)

# Variance explained
explained_var <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2),
  1
)

# Sample annotation
group <- c(
  rep("Tears", ncol(tears_expr)),
  rep("CSF", ncol(csf_expr))
)

# PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = group
)

# PCA plot
pca_plot <- ggplot(
  pca_df,
  aes(PC1, PC2, color = Group)
) +
  geom_point(
    size = 3,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      Tears = "#8B5CF6",
      CSF   = "#2563EB"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA of Shared Proteins",
    x = paste0("PC1 (", explained_var[1], "%)"),
    y = paste0("PC2 (", explained_var[2], "%)")
  )

ggsave(
  "results/Project2_PCA.png",
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# 9. Summary -------------------------------------------------------------------
results_summary <- data.frame(
  Shared_Proteins     = length(shared_proteins),
  Tears_Specific      = length(tears_specific),
  CSF_Specific        = length(csf_specific),
  Pearson_Correlation = round(correlation, 3)
)

write.csv(results_summary, "results/Project2_Summary.csv", row.names = FALSE)
