#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

base_dir <- ""./data/"投稿定图/cohort/structured_analysis/enterosignature_genus_all_samples"
tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")
species_summary_path <- file.path(
  base_dir,
  "maaslin3_species_by_enterosignature_response_no_covariates",
  "tables",
  "maaslin3_species_filtered_species_all_enterosignatures.csv"
)

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

raw_es_levels <- c(
  "ES1_Segatella",
  "ES2_Blautia",
  "ES3_Enterococcus",
  "ES4_Bacteroides",
  "ES5_Phocaeicola"
)

display_es_labels <- c(
  ES1_Segatella = "ES1 Segatella",
  ES2_Blautia = "ES2 Blautia",
  ES3_Enterococcus = "ES3 Enterococcus",
  ES4_Bacteroides = "ES4 Bacteroides",
  ES5_Phocaeicola = "ES5 Phocaeicola"
)

es_colors <- c(
  ES1_Segatella = "#B75B4D",
  ES2_Blautia = "#4C78A8",
  ES3_Enterococcus = "#6A994E",
  ES4_Bacteroides = "#7A6FAC",
  ES5_Phocaeicola = "#D08C36"
)

abbreviate_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  vapply(strsplit(x, " ", fixed = TRUE), function(parts) {
    parts <- parts[nzchar(parts)]
    if (length(parts) < 2) {
      return(parts[1])
    }
    paste(c(paste0(substr(parts[1], 1, 1), "."), parts[-1]), collapse = " ")
  }, character(1))
}

top_n <- 8L

species_summary <- fread(species_summary_path) %>%
  filter(
    enterosignature %in% raw_es_levels,
    is.finite(mean_relative_abundance),
    mean_relative_abundance > 0
  ) %>%
  mutate(
    enterosignature = factor(enterosignature, levels = raw_es_levels),
    enterosignature_label = display_es_labels[as.character(enterosignature)],
    species_label = abbreviate_species(species)
  )

dominant_species <- species_summary %>%
  group_by(enterosignature) %>%
  arrange(desc(mean_relative_abundance), desc(prevalence), species, .by_group = TRUE) %>%
  mutate(rank_in_enterosignature = row_number()) %>%
  filter(rank_in_enterosignature <= top_n) %>%
  ungroup() %>%
  mutate(
    enterosignature_label = factor(enterosignature_label, levels = display_es_labels),
    species_facet_label = paste(enterosignature, sprintf("%02d", top_n + 1L - rank_in_enterosignature), species_label, sep = "___"),
    species_facet_label = factor(species_facet_label, levels = unique(species_facet_label)),
    mean_relative_abundance_percent = 100 * mean_relative_abundance,
    prevalence_percent = 100 * prevalence
  )

dominant_species_table <- dominant_species %>%
  transmute(
    enterosignature = as.character(enterosignature),
    enterosignature_label = as.character(enterosignature_label),
    rank_in_enterosignature,
    species,
    species_label,
    prevalence,
    prevalence_percent,
    mean_relative_abundance,
    mean_relative_abundance_percent,
    mean_read_count
  ) %>%
  arrange(enterosignature, rank_in_enterosignature)

fwrite(
  dominant_species_table,
  file.path(tables_dir, "enterosignature_dominant_species_top8.csv")
)

axis_labels <- setNames(
  sub("^.*___.*___", "", levels(dominant_species$species_facet_label)),
  levels(dominant_species$species_facet_label)
)

plot_theme <- theme_classic(base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.28, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 6.7, face = "italic"),
    axis.text.x = element_text(size = 6.7),
    axis.title.x = element_text(size = 7.2, face = "bold", margin = margin(t = 4)),
    strip.background = element_rect(fill = "grey95", color = "grey55", linewidth = 0.28),
    strip.text = element_text(size = 7.4, face = "bold", color = "black"),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0),
    plot.subtitle = element_text(size = 7.0, color = "grey35", hjust = 0),
    plot.margin = margin(4, 6, 4, 4),
    panel.spacing.x = unit(0.8, "lines"),
    legend.position = "none"
  )

p <- ggplot(dominant_species, aes(x = mean_relative_abundance, y = species_facet_label, fill = enterosignature)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.18) +
  geom_text(
    aes(label = percent(mean_relative_abundance, accuracy = 0.1)),
    hjust = -0.12,
    size = 2.05,
    color = "grey15"
  ) +
  facet_wrap(~ enterosignature_label, nrow = 1, scales = "free_y") +
  scale_y_discrete(labels = axis_labels) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.23))
  ) +
  scale_fill_manual(values = es_colors) +
  labs(
    title = "Dominant species by enterosignature",
    subtitle = "Top 8 species ranked by mean relative abundance within each ES subtype",
    x = "Mean relative abundance",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  plot_theme

ggsave(
  file.path(figures_dir, "enterosignature_dominant_species_top8.pdf"),
  p,
  width = 12.0,
  height = 4.6,
  device = "pdf"
)

ggsave(
  file.path(figures_dir, "enterosignature_dominant_species_top8.png"),
  p,
  width = 12.0,
  height = 4.6,
  dpi = 360
)

message("Wrote dominant species table and figures:")
message(file.path(tables_dir, "enterosignature_dominant_species_top8.csv"))
message(file.path(figures_dir, "enterosignature_dominant_species_top8.pdf"))
message(file.path(figures_dir, "enterosignature_dominant_species_top8.png"))
