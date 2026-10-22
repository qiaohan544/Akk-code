#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(magick)
})

base_dir <- ""./data/"投稿定图/cohort/structured_analysis/enterosignature_genus_all_samples"
figures_dir <- file.path(base_dir, "figures")
tables_dir <- file.path(base_dir, "tables")
species_tables_dir <- file.path(base_dir, "maaslin3_species_by_enterosignature_response_no_covariates", "tables")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

fig2_path <- file.path(base_dir, "Fig2.pdf")
species_summary_path <- file.path(species_tables_dir, "maaslin3_species_filtered_species_all_enterosignatures.csv")
assignment_path <- file.path(tables_dir, "enterosignature_sample_assignments.csv")

raw_es_levels <- c(
  "ES1_Segatella",
  "ES2_Blautia",
  "ES3_Enterococcus",
  "ES4_Bacteroides",
  "ES5_Phocaeicola"
)

display_es_levels <- c(
  "ES_Segatella",
  "ES_Blautia",
  "ES_Enterococcus",
  "ES_Bacteroides",
  "ES_Phocaeicola"
)

enterosignature_name_map <- c(
  ES1_Segatella = "ES_Segatella",
  ES2_Blautia = "ES_Blautia",
  ES3_Enterococcus = "ES_Enterococcus",
  ES4_Bacteroides = "ES_Bacteroides",
  ES5_Phocaeicola = "ES_Phocaeicola"
)

short_es_labels <- c(
  ES_Segatella = "Segatella",
  ES_Blautia = "Blautia",
  ES_Enterococcus = "Enterococcus",
  ES_Bacteroides = "Bacteroides",
  ES_Phocaeicola = "Phocaeicola"
)

response_colors <- c(
  Non_Responder = "#AFC0D6",
  Responder = "#D87563"
)

abbreviate_species <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  vapply(strsplit(x, " ", fixed = TRUE), function(parts) {
    parts <- parts[nzchar(parts)]
    if (length(parts) < 2) return(parts[1])
    paste(c(paste0(substr(parts[1], 1, 1), "."), parts[-1]), collapse = " ")
  }, character(1))
}

compute_cluster_order <- function(df) {
  z_matrix <- df %>%
    select(species, enterosignature, global_z_log10_relative_abundance) %>%
    as.data.table()
  z_matrix <- data.table::dcast(z_matrix, species ~ enterosignature, value.var = "global_z_log10_relative_abundance")

  species_names <- z_matrix$species
  z_values <- as.matrix(z_matrix[, -1, with = FALSE])
  rownames(z_values) <- species_names

  row_order <- rownames(z_values)
  col_order <- colnames(z_values)

  if (nrow(z_values) >= 2) {
    row_order <- rownames(z_values)[hclust(dist(z_values), method = "ward.D2")$order]
  }
  if (ncol(z_values) >= 2) {
    col_order <- colnames(z_values)[hclust(dist(t(z_values)), method = "ward.D2")$order]
  }

  list(row_order = row_order, col_order = col_order)
}

select_top_species <- function(candidate_df, n_top = 20L) {
  candidate_df %>%
    filter(mean_relative_abundance > 0) %>%
    group_by(species) %>%
    summarise(
      max_mean_relative_abundance = max(mean_relative_abundance, na.rm = TRUE),
      max_marker_score = max(marker_score, na.rm = TRUE),
      max_specificity = max(specificity, na.rm = TRUE),
      marker_enterosignature = as.character(enterosignature[which.max(mean_relative_abundance)]),
      .groups = "drop"
    ) %>%
    arrange(desc(max_mean_relative_abundance), desc(max_marker_score), desc(max_specificity), species) %>%
    slice_head(n = n_top) %>%
    select(species, marker_enterosignature)
}

species_summary <- fread(species_summary_path) %>%
  mutate(
    enterosignature = recode(enterosignature, !!!enterosignature_name_map),
    enterosignature = factor(enterosignature, levels = display_es_levels)
  )

all_species <- sort(unique(species_summary$species))
complete_grid <- CJ(species = all_species, enterosignature = display_es_levels) %>%
  mutate(enterosignature = factor(enterosignature, levels = display_es_levels))

abundance_df <- complete_grid %>%
  left_join(
    species_summary %>% select(enterosignature, species, prevalence, mean_relative_abundance),
    by = c("enterosignature", "species")
  ) %>%
  mutate(
    prevalence = ifelse(is.na(prevalence), 0, prevalence),
    mean_relative_abundance = ifelse(is.na(mean_relative_abundance), 0, mean_relative_abundance)
  ) %>%
  group_by(species) %>%
  mutate(
    species_total_abundance = sum(mean_relative_abundance, na.rm = TRUE),
    specificity = ifelse(species_total_abundance > 0, mean_relative_abundance / species_total_abundance, 0)
  ) %>%
  ungroup()

pseudocount <- min(abundance_df$mean_relative_abundance[abundance_df$mean_relative_abundance > 0], na.rm = TRUE) / 2
abundance_df <- abundance_df %>%
  mutate(
    log10_relative_abundance = log10(mean_relative_abundance + pseudocount),
    marker_score = mean_relative_abundance * specificity
  )

selected_species <- select_top_species(abundance_df, n_top = 20L)
selected_species_order <- selected_species$species

selected_abundance_df <- abundance_df %>%
  filter(species %in% selected_species_order) %>%
  mutate(global_z_log10_relative_abundance = as.numeric(scale(log10_relative_abundance)))

cluster_orders <- selected_abundance_df %>%
  compute_cluster_order()

clustered_species_order <- cluster_orders$row_order
clustered_es_order <- cluster_orders$col_order

panel_e_df <- selected_abundance_df %>%
  left_join(selected_species, by = "species") %>%
  mutate(
    marker_flag = as.character(enterosignature) == marker_enterosignature,
    species_label = abbreviate_species(species),
    species_label = factor(species_label, levels = rev(abbreviate_species(clustered_species_order))),
    enterosignature = factor(as.character(enterosignature), levels = clustered_es_order)
  )

panel_e_table <- panel_e_df %>%
  transmute(
    species,
    species_label = as.character(species_label),
    enterosignature = as.character(enterosignature),
    marker_enterosignature,
    marker_flag,
    prevalence,
    mean_relative_abundance,
    log10_relative_abundance,
    global_z_log10_relative_abundance,
    specificity
  ) %>%
  arrange(match(species, clustered_species_order), enterosignature)
fwrite(panel_e_table, file.path(tables_dir, "figure_E_species_global_z_relative_abundance.csv"))

z_limit <- max(abs(panel_e_df$global_z_log10_relative_abundance), na.rm = TRUE)
z_limit <- max(1.5, min(3.5, z_limit))

shared_panel_theme <- theme_classic(base_size = 7) +
  theme(
    axis.line = element_line(linewidth = 0.28, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 9.0, hjust = 0),
    plot.subtitle = element_text(size = 6.2, color = "grey35", hjust = 0),
    plot.margin = margin(2, 3, 2, 2)
  )

panel_e <- ggplot(panel_e_df, aes(x = enterosignature, y = species_label)) +
  geom_point(
    aes(size = mean_relative_abundance, fill = global_z_log10_relative_abundance),
    shape = 21,
    color = "grey30",
    stroke = 0.22,
    alpha = 0.96
  ) +
  geom_point(
    data = panel_e_df %>% filter(marker_flag),
    shape = 21,
    aes(size = mean_relative_abundance),
    stroke = 0.45,
    color = "black",
    fill = NA
  ) +
  scale_x_discrete(drop = FALSE, labels = short_es_labels[clustered_es_order]) +
  scale_size_continuous(
    trans = "sqrt",
    range = c(0.9, 6.5),
    breaks = c(1e-04, 1e-03, 1e-02, 1e-01),
    labels = label_scientific(digits = 1),
    name = "Mean relative\nabundance"
  ) +
  scale_fill_gradient2(
    low = "#2F6F9F",
    mid = "#F7F7F2",
    high = "#B24A48",
    midpoint = 0,
    limits = c(-z_limit, z_limit),
    oob = squish,
    name = "Global z-score\nlog10 mean RA"
  ) +
  guides(
    size = guide_legend(order = 1, override.aes = list(fill = "white")),
    fill = guide_colorbar(order = 2, barheight = unit(1.0, "cm"), barwidth = unit(0.18, "cm"))
  ) +
  labs(
    title = "Characteristic species across enterosignatures",
    subtitle = "Top 20 species; bubble size shows mean relative abundance and color shows global z-score",
    x = NULL,
    y = NULL
  ) +
  shared_panel_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.0),
    axis.text.y = element_text(size = 5.7, face = "italic"),
    legend.position = "right",
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.1)
  )

assignment_df <- fread(assignment_path) %>%
  filter(SPT == "Baseline", Response %in% c("Responder", "Non_Responder"), enterosignature %in% raw_es_levels) %>%
  mutate(
    enterosignature = recode(enterosignature, !!!enterosignature_name_map),
    enterosignature = factor(enterosignature, levels = display_es_levels),
    Response = factor(Response, levels = c("Non_Responder", "Responder"))
  )

response_df <- assignment_df %>%
  count(enterosignature, Response, name = "n") %>%
  group_by(enterosignature) %>%
  mutate(total_n = sum(n), fraction = n / total_n) %>%
  ungroup()

response_table <- response_df %>%
  mutate(percent = 100 * fraction) %>%
  arrange(enterosignature, Response)
fwrite(response_table, file.path(tables_dir, "figure_G_enterosignature_response_composition.csv"))

panel_g <- ggplot(response_df, aes(x = fraction, y = enterosignature, fill = Response)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.25) +
  geom_text(
    aes(label = ifelse(fraction >= 0.12, percent(fraction, accuracy = 1), "")),
    position = position_stack(vjust = 0.5),
    size = 2.0,
    color = "black"
  ) +
  geom_text(
    data = response_df %>% distinct(enterosignature, total_n),
    aes(x = 1.04, y = enterosignature, label = paste0("n=", total_n)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.05,
    color = "grey25"
  ) +
  scale_y_discrete(drop = FALSE, labels = short_es_labels) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.20), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = response_colors, labels = c(Non_Responder = "NR", Responder = "R"), name = NULL) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Response composition by enterosignature",
    subtitle = "Baseline samples with evaluable response and enterosignature",
    x = "Fraction of samples",
    y = NULL,
    tag = "G"
  ) +
  shared_panel_theme +
  theme(
    axis.text.x = element_text(size = 6.0),
    axis.text.y = element_text(size = 6.6, face = "bold"),
    axis.title.x = element_text(size = 6.8),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 6.1),
    legend.key.size = unit(0.26, "cm"),
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(6, 3, 2, 2)
  )

panel_e_pdf <- file.path(figures_dir, "figure_E_species_global_z_relative_abundance.pdf")
panel_e_png <- file.path(figures_dir, "figure_E_species_global_z_relative_abundance.png")
panel_g_pdf <- file.path(figures_dir, "figure_G_enterosignature_response_composition.pdf")
panel_g_png <- file.path(figures_dir, "figure_G_enterosignature_response_composition.png")
ggsave(panel_e_pdf, panel_e, width = 5.2, height = 3.8, device = "pdf")
ggsave(panel_e_png, panel_e, width = 5.2, height = 3.8, dpi = 600)
ggsave(panel_g_pdf, panel_g, width = 4.1, height = 2.2, device = "pdf")
ggsave(panel_g_png, panel_g, width = 4.1, height = 2.2, dpi = 600)

fig2_prefix <- file.path(tempdir(), "fig2_panel_eg_base")
system2("pdftoppm", c("-png", "-r", "300", fig2_path, fig2_prefix), stdout = TRUE, stderr = TRUE)
fig2_png <- paste0(fig2_prefix, "-1.png")
if (!file.exists(fig2_png)) {
  stop("Failed to rasterize Fig2.pdf with pdftoppm: ", fig2_png)
}

fig2 <- image_read(fig2_png)
panel_e_image <- image_read(panel_e_png)
panel_g_image <- image_read(panel_g_png)

panel_e_width <- 930
panel_e_height <- round(panel_e_width * image_info(panel_e_image)$height / image_info(panel_e_image)$width)
panel_e_resized <- image_resize(panel_e_image, paste0(panel_e_width, "x", panel_e_height, "!"))

panel_g_width <- 820
panel_g_height <- round(panel_g_width * image_info(panel_g_image)$height / image_info(panel_g_image)$width)
panel_g_resized <- image_resize(panel_g_image, paste0(panel_g_width, "x", panel_g_height, "!"))

fig2_with_eg <- image_composite(
  fig2,
  panel_e_resized,
  operator = "over",
  offset = "+1515+970"
)
fig2_with_eg <- image_composite(
  fig2_with_eg,
  panel_g_resized,
  operator = "over",
  offset = "+1580+2420"
)

out_png <- file.path(base_dir, "Fig2_with_EG.png")
out_pdf <- file.path(base_dir, "Fig2_with_EG.pdf")
legacy_png <- file.path(base_dir, "Fig2_with_E.png")
legacy_pdf <- file.path(base_dir, "Fig2_with_E.pdf")
image_write(fig2_with_eg, out_png, density = "300x300")
image_write(fig2_with_eg, out_pdf, format = "pdf", density = "300x300")
image_write(fig2_with_eg, legacy_png, density = "300x300")
image_write(fig2_with_eg, legacy_pdf, format = "pdf", density = "300x300")

message("Wrote: ", panel_e_pdf)
message("Wrote: ", panel_e_png)
message("Wrote: ", file.path(tables_dir, "figure_E_species_global_z_relative_abundance.csv"))
message("Wrote: ", panel_g_pdf)
message("Wrote: ", panel_g_png)
message("Wrote: ", file.path(tables_dir, "figure_G_enterosignature_response_composition.csv"))
message("Wrote: ", out_png)
message("Wrote: ", out_pdf)
