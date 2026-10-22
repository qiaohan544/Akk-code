#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(SpiecEasi)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

cohort_dir <- ""./data/"投稿定图/cohort"
rds_path <- file.path(cohort_dir, "phyloseq_object.rds")
out_dir <- file.path(cohort_dir, "structured_analysis", "akk_t2_t3_diff_sparcc")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")

target_projects <- c("PRJNA751792", "PRJNA762360")
depth_cutoff <- 20000L
akk_cutoff_low <- 0.0002313
akk_cutoff_high <- 0.005841

sparcc_iter <- 20L
sparcc_inner_iter <- 10L
boot_reps <- 999L
boot_ncpus <- 8L
set.seed(20260502)

strict_cfg <- list(
  species = list(min_prev_each = 0.25, min_mean = 0.0005),
  genus = list(min_prev_each = 0.25, min_mean = 0.0005)
)

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

clean_chr <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA")] <- NA_character_
  x_chr
}

extract_sample_data <- function(obj) {
  sample_df <- attr(obj, "sam_data", exact = TRUE)
  class(sample_df) <- "data.frame"
  rn <- attr(sample_df, "row.names", exact = TRUE)
  if (!is.null(rn)) {
    rownames(sample_df) <- rn
  }
  sample_df
}

extract_matrix_like <- function(obj, attr_name) {
  x <- attr(obj, attr_name, exact = TRUE)
  x_dimnames <- attr(x, "dimnames", exact = TRUE)
  suppressWarnings(class(x) <- NULL)
  if (!is.null(x_dimnames)) {
    dimnames(x) <- x_dimnames
  }
  x
}

aggregate_taxa_counts <- function(otu_mat, tax_df, level_col, sample_ids) {
  tax_level <- clean_chr(tax_df[[level_col]])
  tax_level[is.na(tax_level) | tax_level %in% c("?", "unknown", "Unknown", "unclassified")] <- NA_character_
  keep <- tax_df$feature %in% rownames(otu_mat) & !is.na(tax_level)
  rowsum(otu_mat[tax_df$feature[keep], sample_ids, drop = FALSE], group = tax_level[keep], reorder = TRUE)
}

vec_to_symmat <- function(x, labels, diag_value = 1) {
  n <- length(labels)
  mat <- matrix(diag_value, nrow = n, ncol = n, dimnames = list(labels, labels))
  mat[upper.tri(mat)] <- x
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  mat
}

to_sym_from_sparcc <- function(x, labels, diag_value) {
  if (is.matrix(x)) {
    out <- x
    rownames(out) <- rownames(out) %||% labels
    colnames(out) <- colnames(out) %||% labels
    return(out[labels, labels, drop = FALSE])
  }
  vec_to_symmat(x, labels, diag_value = diag_value)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

run_differential_t2_t3 <- function(taxa_counts, group_df, level_name, min_prev_each, min_mean, anchor_taxon) {
  taxa_counts <- taxa_counts[, group_df$sample_id, drop = FALSE]
  rel <- sweep(taxa_counts, 2, colSums(taxa_counts), "/")
  rel[!is.finite(rel)] <- 0

t2_ids <- group_df$sample_id[group_df$AKK_tertile == "T2"]
t3_ids <- group_df$sample_id[group_df$AKK_tertile == "T3"]

  prev_t2 <- rowMeans(taxa_counts[, t2_ids, drop = FALSE] > 0)
  prev_t3 <- rowMeans(taxa_counts[, t3_ids, drop = FALSE] > 0)
  mean_rel <- rowMeans(rel)

  keep <- prev_t2 >= min_prev_each & prev_t3 >= min_prev_each & mean_rel >= min_mean
  if (anchor_taxon %in% rownames(rel)) {
    keep[anchor_taxon] <- FALSE
  }

  taxa_kept <- rownames(rel)[keep]
  if (length(taxa_kept) == 0) {
    return(data.frame())
  }

  diff_rows <- lapply(taxa_kept, function(tx) {
    x2 <- rel[tx, t2_ids]
    x3 <- rel[tx, t3_ids]
    wt <- tryCatch(wilcox.test(x3, x2, exact = FALSE), error = function(e) NULL)
    data.frame(
      tax_level = level_name,
      taxon = tx,
      n_t2 = length(x2),
      n_t3 = length(x3),
      prevalence_t2 = prev_t2[tx],
      prevalence_t3 = prev_t3[tx],
      mean_rel_t2 = mean(x2),
      mean_rel_t3 = mean(x3),
      median_rel_t2 = median(x2),
      median_rel_t3 = median(x3),
      log2fc_t3_vs_t2 = log2((median(x3) + 1e-6) / (median(x2) + 1e-6)),
      delta_prevalence = prev_t3[tx] - prev_t2[tx],
      p_value = if (!is.null(wt)) wt$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(diff_rows) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr, p_value, desc(abs(log2fc_t3_vs_t2)), taxon)
  out
}

run_sparcc_akk_anchor <- function(taxa_counts, group_df, group_name, level_name, min_prev_each, min_mean, anchor_taxon) {
  sample_ids <- group_df$sample_id[group_df$AKK_tertile == group_name]
  mat_counts <- t(taxa_counts[, sample_ids, drop = FALSE])

  rel <- sweep(mat_counts, 1, rowSums(mat_counts), "/")
  rel[!is.finite(rel)] <- 0
  prevalence <- colMeans(mat_counts > 0)
  mean_rel <- colMeans(rel)

  keep <- prevalence >= min_prev_each & mean_rel >= min_mean
  if (anchor_taxon %in% names(keep)) {
    keep[anchor_taxon] <- TRUE
  }

  filtered_counts <- mat_counts[, keep, drop = FALSE]
  sparcc_input <- filtered_counts + 1

  if (!(anchor_taxon %in% colnames(sparcc_input))) {
    stop(sprintf("Anchor taxon '%s' not retained for %s-%s SparCC.", anchor_taxon, level_name, group_name))
  }

  sparcc_fit <- sparcc(sparcc_input, iter = sparcc_iter, inner_iter = sparcc_inner_iter)
  boot_fit <- sparccboot(
    sparcc_input,
    R = boot_reps,
    sparcc.params = list(iter = sparcc_iter, inner_iter = sparcc_inner_iter),
    ncpus = boot_ncpus
  )
  p_fit <- pval.sparccboot(boot_fit, sided = "both")
  min_empirical_p <- 1 / (boot_reps + 1)

  labels <- colnames(sparcc_input)
  cor_mat <- to_sym_from_sparcc(sparcc_fit$Cor, labels, diag_value = 1)
  p_mat <- to_sym_from_sparcc(p_fit$pvals, labels, diag_value = NA_real_)

  out <- data.frame(
    tax_level = level_name,
    akk_group = group_name,
    sparcc_input = "reads_plus1",
    taxon = setdiff(labels, anchor_taxon),
    n_samples = length(sample_ids),
    prevalence = prevalence[setdiff(labels, anchor_taxon)],
    mean_relative_abundance = mean_rel[setdiff(labels, anchor_taxon)],
    sparcc_rho = cor_mat[anchor_taxon, setdiff(labels, anchor_taxon)],
    p_value = p_mat[anchor_taxon, setdiff(labels, anchor_taxon)],
    stringsAsFactors = FALSE
  )
  out$p_value <- ifelse(is.na(out$p_value), NA_real_, pmax(out$p_value, min_empirical_p))
  out$fdr <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_value, -abs(out$sparcc_rho), out$taxon), ]

  filtered <- data.frame(
    tax_level = level_name,
    akk_group = group_name,
    sparcc_input = "reads_plus1",
    taxon = labels,
    prevalence = prevalence[labels],
    mean_relative_abundance = mean_rel[labels],
    stringsAsFactors = FALSE
  )

  list(result = out, filtered_taxa = filtered)
}

message("Loading phyloseq-like object")
obj <- readRDS(rds_path)
sample_df <- extract_sample_data(obj)
tax_mat <- extract_matrix_like(obj, "tax_table")
otu_mat <- extract_matrix_like(obj, "otu_table")

sample_df$sample_id <- rownames(sample_df)
sample_df$project <- clean_chr(sample_df$project)
sample_df$SPT <- clean_chr(sample_df$SPT)

bacterial_taxa <- rownames(tax_mat)[tax_mat[, "Kingdom"] == "Bacteria"]
otu_mat <- otu_mat[intersect(rownames(otu_mat), bacterial_taxa), , drop = FALSE]
sample_df$read_depth_bacteria <- colSums(otu_mat[, sample_df$sample_id, drop = FALSE])

base_df <- sample_df %>%
  filter(
    Host == "Human",
    AssayType == "WGS",
    project %in% target_projects,
    SPT == "Baseline",
    read_depth_bacteria >= depth_cutoff
  ) %>%
  select(sample_id, project)

sample_ids <- intersect(base_df$sample_id, colnames(otu_mat))
otu_subset <- otu_mat[, sample_ids, drop = FALSE]
tax_df <- data.frame(feature = rownames(tax_mat), tax_mat, check.names = FALSE, stringsAsFactors = FALSE)

species_counts <- aggregate_taxa_counts(otu_subset, tax_df, "Species", sample_ids)
genus_counts <- aggregate_taxa_counts(otu_subset, tax_df, "Genus", sample_ids)

if (!("Akkermansia muciniphila" %in% rownames(species_counts))) {
  stop("Akkermansia muciniphila not found in species-level matrix.")
}
if (!("Akkermansia" %in% rownames(genus_counts))) {
  stop("Akkermansia not found in genus-level matrix.")
}

species_rel <- sweep(species_counts, 2, colSums(species_counts), "/")
species_rel[!is.finite(species_rel)] <- 0

group_df <- base_df %>%
  mutate(
    AKK_rel_abundance = as.numeric(species_rel["Akkermansia muciniphila", sample_id]),
    AKK_tertile = case_when(
      AKK_rel_abundance <= akk_cutoff_low ~ "T1",
      AKK_rel_abundance <= akk_cutoff_high ~ "T2",
      AKK_rel_abundance > akk_cutoff_high ~ "T3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(AKK_tertile %in% c("T2", "T3"))

if (nrow(group_df) == 0) {
  stop("No T2/T3 samples after filtering.")
}

group_summary <- group_df %>%
  group_by(project, AKK_tertile) %>%
  summarise(n = n(), .groups = "drop")
message(sprintf("T2/T3 sample counts: %s", paste(sprintf("%s-%s=%d", group_summary$project, group_summary$AKK_tertile, group_summary$n), collapse = "; ")))

fwrite(group_df, file.path(tables_dir, "akk_t2_t3_sample_assignments.csv"))
fwrite(group_summary, file.path(tables_dir, "akk_t2_t3_group_counts_by_project.csv"))

run_level_pipeline <- function(level_name, counts_mat, anchor_taxon, cfg) {
  diff_df <- run_differential_t2_t3(
    taxa_counts = counts_mat,
    group_df = group_df,
    level_name = level_name,
    min_prev_each = cfg$min_prev_each,
    min_mean = cfg$min_mean,
    anchor_taxon = anchor_taxon
  )

  t2_fit <- run_sparcc_akk_anchor(
    taxa_counts = counts_mat,
    group_df = group_df,
    group_name = "T2",
    level_name = level_name,
    min_prev_each = cfg$min_prev_each,
    min_mean = cfg$min_mean,
    anchor_taxon = anchor_taxon
  )
  t3_fit <- run_sparcc_akk_anchor(
    taxa_counts = counts_mat,
    group_df = group_df,
    group_name = "T3",
    level_name = level_name,
    min_prev_each = cfg$min_prev_each,
    min_mean = cfg$min_mean,
    anchor_taxon = anchor_taxon
  )

  sparcc_df <- bind_rows(t2_fit$result, t3_fit$result)
  filtered_df <- bind_rows(t2_fit$filtered_taxa, t3_fit$filtered_taxa)

  sparcc_comp <- full_join(
    t2_fit$result %>% select(taxon, t2_rho = sparcc_rho, t2_p = p_value, t2_fdr = fdr),
    t3_fit$result %>% select(taxon, t3_rho = sparcc_rho, t3_p = p_value, t3_fdr = fdr),
    by = "taxon"
  ) %>%
    mutate(delta_rho_t3_minus_t2 = t3_rho - t2_rho) %>%
    arrange(desc(abs(delta_rho_t3_minus_t2)))

  fwrite(diff_df, file.path(tables_dir, sprintf("%s_t2_vs_t3_differential_taxa.csv", level_name)))
  fwrite(sparcc_df, file.path(tables_dir, sprintf("%s_akk_sparcc_correlations_t2_t3.csv", level_name)))
  fwrite(t2_fit$result, file.path(tables_dir, sprintf("%s_akk_sparcc_correlations_t2.csv", level_name)))
  fwrite(t3_fit$result, file.path(tables_dir, sprintf("%s_akk_sparcc_correlations_t3.csv", level_name)))
  fwrite(filtered_df, file.path(tables_dir, sprintf("%s_sparcc_filtered_taxa_t2_t3.csv", level_name)))
  fwrite(sparcc_comp, file.path(tables_dir, sprintf("%s_akk_sparcc_correlation_comparison_t2_vs_t3.csv", level_name)))

  top_sparcc <- sparcc_df %>%
    filter(is.finite(sparcc_rho), is.finite(p_value)) %>%
    group_by(akk_group) %>%
    slice_head(n = 12) %>%
    ungroup() %>%
    mutate(taxon = factor(taxon, levels = rev(unique(taxon[order(akk_group, sparcc_rho)]))))

  p_sparcc <- ggplot(top_sparcc, aes(x = sparcc_rho, y = taxon, color = akk_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey60") +
    geom_point(aes(size = -log10(pmax(p_value, .Machine$double.xmin))), alpha = 0.9) +
    facet_wrap(~ akk_group, scales = "free_y") +
    scale_color_manual(values = c("T2" = "#3C779F", "T3" = "#D97A6A")) +
    labs(
      title = sprintf("%s-level SparCC with AKK anchor in T2 vs T3", tools::toTitleCase(level_name)),
      subtitle = sprintf("Strict filter: prevalence >= %.2f in each group, mean relative abundance >= %.4f", cfg$min_prev_each, cfg$min_mean),
      x = "SparCC correlation",
      y = NULL,
      size = expression(-log[10](p))
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(figures_dir, sprintf("%s_akk_sparcc_top_correlations_t2_t3.pdf", level_name)), p_sparcc, width = 7.4, height = 4.8)
  ggsave(file.path(figures_dir, sprintf("%s_akk_sparcc_top_correlations_t2_t3.png", level_name)), p_sparcc, width = 7.4, height = 4.8, dpi = 600)

  volcano_df <- diff_df %>%
    filter(is.finite(log2fc_t3_vs_t2), is.finite(p_value)) %>%
    mutate(sig = ifelse(fdr <= 0.05, "FDR<=0.05", "NS"))

  p_diff <- ggplot(volcano_df, aes(x = log2fc_t3_vs_t2, y = -log10(pmax(p_value, .Machine$double.xmin)), color = sig)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_manual(values = c("NS" = "#9AA1A8", "FDR<=0.05" = "#A35C72")) +
    labs(
      title = sprintf("%s-level differential taxa: T3 vs T2", tools::toTitleCase(level_name)),
      subtitle = sprintf("Strict filter: prevalence >= %.2f in each group, mean relative abundance >= %.4f", cfg$min_prev_each, cfg$min_mean),
      x = expression(log[2]("FC (T3/T2) by median relative abundance")),
      y = expression(-log[10](p)),
      color = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(figures_dir, sprintf("%s_t2_vs_t3_differential_volcano.pdf", level_name)), p_diff, width = 6.8, height = 4.8)
  ggsave(file.path(figures_dir, sprintf("%s_t2_vs_t3_differential_volcano.png", level_name)), p_diff, width = 6.8, height = 4.8, dpi = 600)

  list(
    diff_df = diff_df,
    sparcc_df = sparcc_df,
    top_diff = diff_df %>% slice_head(n = 10),
    top_t2 = t2_fit$result %>% slice_head(n = 8),
    top_t3 = t3_fit$result %>% slice_head(n = 8)
  )
}

species_out <- run_level_pipeline(
  level_name = "species",
  counts_mat = species_counts,
  anchor_taxon = "Akkermansia muciniphila",
  cfg = strict_cfg$species
)

genus_out <- run_level_pipeline(
  level_name = "genus",
  counts_mat = genus_counts,
  anchor_taxon = "Akkermansia",
  cfg = strict_cfg$genus
)

report_lines <- c(
  "# T2 vs T3 differential taxa and within-group SparCC (AKK tertile)",
  "",
  "## Cohort and grouping",
  sprintf("- Cohorts: %s.", paste(target_projects, collapse = ", ")),
  "- Samples: baseline Human WGS with bacterial reads >= 20000.",
  sprintf("- Fixed AKK tertile cutoffs: low=%.7f, high=%.6f.", akk_cutoff_low, akk_cutoff_high),
  "- Analysis groups: T2 (middle AKK) vs T3 (high AKK).",
  "",
  "## Strict filtering",
  sprintf("- Species level: prevalence >= %.2f in both T2 and T3; mean relative abundance >= %.4f.", strict_cfg$species$min_prev_each, strict_cfg$species$min_mean),
  sprintf("- Genus level: prevalence >= %.2f in both T2 and T3; mean relative abundance >= %.4f.", strict_cfg$genus$min_prev_each, strict_cfg$genus$min_mean),
  sprintf("- SparCC settings: iter=%d, inner_iter=%d, bootstrap R=%d, ncpus=%d.", sparcc_iter, sparcc_inner_iter, boot_reps, boot_ncpus),
  "- SparCC input matrix: filtered reads counts with +1 pseudocount (not relative-abundance matrix).",
  "",
  "## Differential taxa (T3 vs T2)",
  sprintf("- Species: tested taxa=%d; FDR<=0.05 taxa=%d.", nrow(species_out$diff_df), sum(species_out$diff_df$fdr <= 0.05, na.rm = TRUE)),
  sprintf("- Genus: tested taxa=%d; FDR<=0.05 taxa=%d.", nrow(genus_out$diff_df), sum(genus_out$diff_df$fdr <= 0.05, na.rm = TRUE)),
  "",
  "## Top species differential taxa",
  if (nrow(species_out$top_diff) > 0) paste(sprintf("- %s: log2FC=%.3f, p=%.3g, FDR=%.3g.", species_out$top_diff$taxon, species_out$top_diff$log2fc_t3_vs_t2, species_out$top_diff$p_value, species_out$top_diff$fdr), collapse = "\n") else "- No tested taxa.",
  "",
  "## Top genus differential taxa",
  if (nrow(genus_out$top_diff) > 0) paste(sprintf("- %s: log2FC=%.3f, p=%.3g, FDR=%.3g.", genus_out$top_diff$taxon, genus_out$top_diff$log2fc_t3_vs_t2, genus_out$top_diff$p_value, genus_out$top_diff$fdr), collapse = "\n") else "- No tested taxa.",
  "",
  "## Key files",
  "- tables/species_t2_vs_t3_differential_taxa.csv",
  "- tables/genus_t2_vs_t3_differential_taxa.csv",
  "- tables/species_akk_sparcc_correlations_t2.csv",
  "- tables/species_akk_sparcc_correlations_t3.csv",
  "- tables/genus_akk_sparcc_correlations_t2.csv",
  "- tables/genus_akk_sparcc_correlations_t3.csv",
  "- figures/species_t2_vs_t3_differential_volcano.pdf",
  "- figures/genus_t2_vs_t3_differential_volcano.pdf",
  "- figures/species_akk_sparcc_top_correlations_t2_t3.pdf",
  "- figures/genus_akk_sparcc_top_correlations_t2_t3.pdf"
)

writeLines(report_lines, file.path(out_dir, "akk_t2_t3_diff_sparcc_report.md"))

pdf(file.path(out_dir, "akk_t2_t3_diff_sparcc_report.pdf"), width = 11, height = 8.5, onefile = TRUE)
grid.newpage()
grid.text(
  paste(report_lines, collapse = "\n"),
  x = unit(0.04, "npc"),
  y = unit(0.96, "npc"),
  just = c("left", "top"),
  gp = gpar(fontsize = 8.2, fontfamily = "sans")
)
dev.off()

message("Completed T2 vs T3 differential + SparCC analysis: ", out_dir)