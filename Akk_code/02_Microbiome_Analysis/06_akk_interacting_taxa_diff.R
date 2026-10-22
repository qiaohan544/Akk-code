#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, scipen = 999)

min_prev_each <- 0.25
min_mean_abundance <- 0.0005
sig_p_cutoff <- 0.05
max_taxa <- 20L

cohort_dir <- ""./data/"投稿定图/cohort"
out_dir <- file.path(cohort_dir, "structured_analysis", "akk_t2_t3_diff_sparcc")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("T2", "T3")

shorten_species_label <- function(x) {
  parts <- strsplit(as.character(x), "\\s+")
  vapply(parts, function(p) {
    if (length(p) < 2) {
      return(paste(p, collapse = " "))
    }
    paste0(substr(p[1], 1, 1), ". ", p[2])
  }, character(1))
}

select_network_taxa <- function(df) {
  shared_taxa_rank <- df %>%
    group_by(taxon) %>%
    summarise(
      n_groups = n_distinct(akk_group),
      any_sig = any(p_value <= sig_p_cutoff, na.rm = TRUE),
      min_p = min(p_value, na.rm = TRUE),
      max_abs_rho = max(abs(sparcc_rho), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_groups == 2) %>%
    arrange(desc(any_sig), min_p, desc(max_abs_rho), taxon)

  selected_taxa <- shared_taxa_rank %>%
    filter(any_sig) %>%
    slice_head(n = max_taxa) %>%
    pull(taxon)

  if (length(selected_taxa) == 0) {
    selected_taxa <- shared_taxa_rank %>%
      slice_head(n = max_taxa) %>%
      pull(taxon)
  }

  selected_taxa
}

plot_shifts <- function(level_name) {
  in_file <- file.path(tables_dir, sprintf("%s_akk_sparcc_correlations_t2_t3.csv", level_name))
  out_pdf <- file.path(figures_dir, sprintf("%s_largest_sparcc_shifts_t3_t2.pdf", level_name))
  out_png <- file.path(figures_dir, sprintf("%s_largest_sparcc_shifts_t3_t2.png", level_name))
  out_table <- file.path(tables_dir, sprintf("%s_network_selected_sparcc_shifts_t3_minus_t2.csv", level_name))

  df <- fread(in_file, data.table = FALSE)

  filtered_df <- df %>%
    filter(
      is.finite(sparcc_rho),
      is.finite(p_value),
      prevalence >= min_prev_each,
      mean_relative_abundance >= min_mean_abundance,
      akk_group %in% group_levels
    )

  selected_taxa <- select_network_taxa(filtered_df)
  if (length(selected_taxa) == 0) {
    stop(sprintf("No taxa available for %s shifts plot.", level_name))
  }

  sel_df <- filtered_df %>%
    filter(taxon %in% selected_taxa)

  wide <- as.data.table(sel_df)[, .(
    taxon,
    akk_group,
    sparcc_rho,
    p_value,
    fdr,
    prevalence,
    mean_relative_abundance
  )]

  wide_dt <- dcast(
    wide,
    taxon ~ akk_group,
    value.var = c("sparcc_rho", "p_value", "fdr", "prevalence", "mean_relative_abundance")
  )

  req_cols <- c("sparcc_rho_T2", "sparcc_rho_T3")
  if (!all(req_cols %in% colnames(wide_dt))) {
    stop(sprintf("Missing T2/T3 SparCC columns for %s.", level_name))
  }

  shifts <- as_tibble(wide_dt) %>%
    mutate(
      delta_rho_t3_minus_t2 = sparcc_rho_T3 - sparcc_rho_T2,
      shift_direction = ifelse(delta_rho_t3_minus_t2 >= 0, "T3 > T2", "T3 < T2"),
      p_min = pmin(p_value_T2, p_value_T3, na.rm = TRUE),
      fdr_min = pmin(fdr_T2, fdr_T3, na.rm = TRUE)
    ) %>%
    arrange(desc(delta_rho_t3_minus_t2), desc(abs(delta_rho_t3_minus_t2)), taxon)

  if (level_name == "species") {
    shifts <- shifts %>% mutate(taxon_display = shorten_species_label(taxon))
  } else {
    shifts <- shifts %>% mutate(taxon_display = taxon)
  }

  shifts <- shifts %>%
    mutate(
      taxon_plot = factor(taxon, levels = rev(taxon)),
      label_txt = sprintf("%.3f", delta_rho_t3_minus_t2)
    )

  y_label_map <- shifts$taxon_display
  names(y_label_map) <- shifts$taxon

  p <- ggplot(shifts, aes(x = delta_rho_t3_minus_t2, y = taxon_plot, fill = shift_direction)) +
    geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "#7A7A7A") +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_text(
      aes(
        label = label_txt,
        hjust = ifelse(delta_rho_t3_minus_t2 >= 0, -0.05, 1.05)
      ),
      size = 2.2,
      color = "#1A1A1A"
    ) +
    scale_fill_manual(values = c("T3 > T2" = "#2E6F9E", "T3 < T2" = "#B24C45")) +
    scale_y_discrete(labels = y_label_map) +
    labs(
      title = sprintf("%s-level SparCC shifts among network-selected taxa", tools::toTitleCase(level_name)),
      subtitle = "Shift direction defined as Delta rho = rho(T3) - rho(T2); bars sorted from highest to lowest Delta rho",
      x = expression(Delta * "SparCC (T3 - T2)"),
      y = NULL,
      fill = "Direction"
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      axis.text.y = element_text(size = 6.8),
      legend.position = "top",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 7.2, color = "#3A3A3A"),
      plot.margin = margin(6, 26, 6, 8)
    ) +
    coord_cartesian(clip = "off")

  ggsave(out_pdf, p, width = 7.6, height = 5.3, units = "in")
  ggsave(out_png, p, width = 7.6, height = 5.3, units = "in", dpi = 600)

  fwrite(as.data.table(shifts), out_table)

  list(
    output_pdf = out_pdf,
    output_png = out_png,
    output_table = out_table,
    n_taxa = nrow(shifts),
    max_shift = max(shifts$delta_rho_t3_minus_t2, na.rm = TRUE),
    min_shift = min(shifts$delta_rho_t3_minus_t2, na.rm = TRUE)
  )
}

species_out <- plot_shifts("species")
genus_out <- plot_shifts("genus")

message("Rendered SparCC shift plots from current network-selected taxa:")
message(" - ", species_out$output_pdf, " (n=", species_out$n_taxa, ", range=[", sprintf("%.3f", species_out$min_shift), ", ", sprintf("%.3f", species_out$max_shift), "])")
message(" - ", genus_out$output_pdf, " (n=", genus_out$n_taxa, ", range=[", sprintf("%.3f", genus_out$min_shift), ", ", sprintf("%.3f", genus_out$max_shift), "])")
