library(zellkonverter)
library(sccomp)
library(SingleCellExperiment)
library(tidyverse)
library(ggrepel)
library(ggthemes)

HUMAN_H5AD <- "/path/to/paired_blood_brain_cohort1.h5ad"
OUTPUT_DIR  <- "figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

adata_blbr <- readH5AD(HUMAN_H5AD, verbose = TRUE, use.dimnames = FALSE, reader = "R")

unique(adata_blbr$cell_type.v2)

adata_blbr_common <- adata_blbr[, adata_blbr$cell_type.v2 %in%
  c("CD16.mono", "CD14.mono", "CD8.T", "CD4.T", "CD16.NK", "B.cells")]

adata_blbr_common$pt.tissue      <- paste(adata_blbr_common$pt, adata_blbr_common$tissue)
adata_blbr_common$pt.tissue.side <- paste(adata_blbr_common$pt, adata_blbr_common$tissue, adata_blbr_common$side)
adata_blbr_common$cell_type.v2   <- droplevels(adata_blbr_common$cell_type.v2)

supp_table_1 <- as.data.frame(adata_blbr_common@colData) %>%
  count(pt, cell_type.v2, name = "n_cells") %>%
  pivot_wider(names_from = cell_type.v2, values_from = n_cells, values_fill = 0) %>%
  arrange(pt)

write.csv(supp_table_1, file.path(OUTPUT_DIR, "supp_table_1_cell_counts.csv"), row.names = FALSE)


sccomp_result <-
  adata_blbr_common |>
  sccomp_estimate(
    formula_composition          = ~tissue + chemistry,
    .sample                      = pt.tissue.side,
    .cell_group                  = cell_type.v2,
    bimodal_mean_variability_association = TRUE,
    percent_false_positive       = 5,
    cores                        = 1
  ) |>
  sccomp_remove_outliers(cores = 1) |>
  sccomp_test()

sccomp_result

sccomp_result |>
  sccomp_proportional_fold_change(
    formula_composition = ~tissue + chemistry,
    from = "brain",
    to   = "blood"
  ) |>
  select(cell_group, statement)

p_boxplot_tissue <- sccomp_result |>
  sccomp_boxplot(
    factor                  = "tissue",
    significance_threshold  = 0.01,
    remove_unwanted_effects = FALSE
  ) +
  theme_few(base_size = 25) +
  facet_wrap(~cell_type.v2, nrow = 3, ncol = 3)

ggsave(file.path(OUTPUT_DIR, "sccomp_boxplot_tissue.pdf"), p_boxplot_tissue, width = 14, height = 10)
