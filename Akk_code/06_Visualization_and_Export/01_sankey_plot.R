#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
})

root_dir <- ""./data/"投稿定图/cohort/species_survival_splits"
tables_dir <- file.path(root_dir, "tables", "two_cohort_targeted")
figures_dir <- file.path(root_dir, "figures", "two_cohort_targeted")

group_palette <- c("G1" = "#A35C72", "G2" = "#3C779F", "G3" = "#86D2B2", "G4" = "#81017E")

best_group_df <- fread(file.path(tables_dir, "pooled_best_os_group_effects.csv"))

binary_top10 <- best_group_df %>%
  filter(
    domain == "pooled",
    model_type == "adjusted_sex_age_bmi",
    split_label == "binary",
    is.finite(p_value),
    p_value < 0.05
  ) %>%
  arrange(desc(mean_relative_abundance), p_value, feature_short) %>%
  slice_head(n = 10) %>%
  select(feature, feature_short, mean_relative_abundance) %>%
  distinct(feature, .keep_all = TRUE)

sankey_df <- best_group_df %>%
  filter(
    domain == "pooled",
    model_type == "adjusted_sex_age_bmi",
    split_label %in% c("binary", "tertile")
  ) %>%
  select(feature, split_label, best_group) %>%
  distinct(feature, split_label, .keep_all = TRUE) %>%
  pivot_wider(names_from = split_label, values_from = best_group) %>%
  inner_join(binary_top10, by = "feature") %>%
  transmute(
    Species = feature_short,
    Binary = binary,
    Tertile = tertile,
    mean_relative_abundance = mean_relative_abundance
  ) %>%
  filter(!is.na(Binary), !is.na(Tertile)) %>%
  arrange(desc(mean_relative_abundance), Species)

fwrite(sankey_df, file.path(tables_dir, "pooled_top20_protective_group_sankey.csv"))

sankey_df$Species <- factor(sankey_df$Species, levels = rev(sankey_df$Species))
sankey_fill <- c(setNames(rep("#D9D9D9", nrow(sankey_df)), levels(sankey_df$Species)), group_palette)

p <- ggplot(
  sankey_df,
  aes(axis1 = Species, axis2 = Binary, axis3 = Tertile, y = 1)
) +
  ggalluvial::geom_alluvium(fill = "#CFCFCF", color = "white", alpha = 0.75, width = 0.22) +
  ggalluvial::geom_stratum(aes(fill = after_stat(stratum)), color = "white", width = 0.28) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Species", "Binary", "Tertile"),
    expand = c(0.08, 0.03)
  ) +
  scale_fill_manual(values = sankey_fill, guide = "none") +
  labs(
    title = "pooled | Significant species protective-group flow across OS splits",
    subtitle = "Top-10 abundant species among binary split p < 0.05; tertile labels follow the same species set",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

ggsave(file.path(figures_dir, "pooled_top20_protective_group_sankey.pdf"), p, width = 13, height = 10.5)
ggsave(file.path(figures_dir, "pooled_top20_protective_group_sankey.png"), p, width = 13, height = 10.5, dpi = 320)

message(sprintf("Regenerated sankey with %d species", nrow(sankey_df)))
