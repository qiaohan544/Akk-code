#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(splines)
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- ""./data/"投稿定图/cohort"
analysis_dir <- file.path(root_dir, "species_survival_splits")
tables_dir <- file.path(analysis_dir, "tables", "two_cohort_targeted")
figures_dir <- file.path(analysis_dir, "figures", "two_cohort_targeted")
rds_path <- file.path(root_dir, "phyloseq_object.rds")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

target_projects <- c("PRJNA751792", "PRJNA762360")
depth_cutoff <- 20000L
pseudocount <- 1e-6
spline_df <- 3L
top_curve_n <- 8L

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

coerce_numeric <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", x_chr, fixed = TRUE)))
}

coerce_event <- function(x) {
  x_chr <- trimws(tolower(as.character(x)))
  x_chr[x_chr == ""] <- NA_character_
  out <- rep(NA_integer_, length(x_chr))
  out[x_chr %in% c("1", "dead", "deceased", "death", "yes", "true", "event")] <- 1L
  out[x_chr %in% c("0", "alive", "living", "censored", "no", "false")] <- 0L
  unresolved <- is.na(out) & !is.na(x_chr)
  if (any(unresolved)) {
    x_num <- suppressWarnings(as.numeric(x_chr[unresolved]))
    out[unresolved] <- ifelse(is.na(x_num), NA_integer_, ifelse(x_num > 0, 1L, 0L))
  }
  out
}

safe_p <- function(anova_obj) {
  p_col <- grep("^Pr|^P\\(", colnames(anova_obj), value = TRUE)
  if (length(p_col) == 0 || nrow(anova_obj) < 2) {
    return(NA_real_)
  }
  as.numeric(anova_obj[2, p_col[1]])
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  names(sort(table(x), decreasing = TRUE))[1]
}

make_prediction_curve <- function(fit, model_df, species_name) {
  x_grid <- seq(
    stats::quantile(model_df$abundance_log10, 0.02, na.rm = TRUE, names = FALSE),
    stats::quantile(model_df$abundance_log10, 0.98, na.rm = TRUE, names = FALSE),
    length.out = 120
  )
  ref_x <- stats::median(model_df$abundance_log10, na.rm = TRUE)

  template <- data.frame(
    os_time = stats::median(model_df$os_time, na.rm = TRUE),
    death_event = 0L,
    abundance_log10 = x_grid,
    project = factor(mode_value(as.character(model_df$project)), levels = levels(model_df$project)),
    sex = factor(mode_value(as.character(model_df$sex)), levels = levels(model_df$sex)),
    age_years = stats::median(model_df$age_years, na.rm = TRUE),
    bmi = stats::median(model_df$bmi, na.rm = TRUE)
  )
  ref_row <- template[1, , drop = FALSE]
  ref_row$abundance_log10 <- ref_x

  design_grid <- model.matrix(delete.response(terms(fit)), template)
  design_ref <- model.matrix(delete.response(terms(fit)), ref_row)
  beta <- stats::coef(fit)
  keep_cols <- names(beta)
  design_grid <- design_grid[, keep_cols, drop = FALSE]
  design_ref <- design_ref[, keep_cols, drop = FALSE]
  contrast <- sweep(design_grid, 2, as.numeric(design_ref[1, ]), "-")
  var_beta <- stats::vcov(fit)[keep_cols, keep_cols, drop = FALSE]
  log_hr <- as.numeric(contrast %*% beta)
  se_log_hr <- sqrt(rowSums((contrast %*% var_beta) * contrast))

  data.frame(
    feature_short = species_name,
    abundance_log10 = x_grid,
    abundance_plus_pseudocount = 10^x_grid,
    abundance = pmax(10^x_grid - pseudocount, 0),
    hr = exp(log_hr),
    conf_low = exp(log_hr - 1.96 * se_log_hr),
    conf_high = exp(log_hr + 1.96 * se_log_hr),
    stringsAsFactors = FALSE
  )
}

message("Reading phyloseq object")
obj <- readRDS(rds_path)
sample_df <- extract_sample_data(obj)
tax_mat <- extract_matrix_like(obj, "tax_table")
otu_mat <- extract_matrix_like(obj, "otu_table")

sample_df$sample_id <- rownames(sample_df)
bacterial_taxa <- rownames(tax_mat)[tax_mat[, "Kingdom"] == "Bacteria"]
otu_mat <- otu_mat[intersect(rownames(otu_mat), bacterial_taxa), , drop = FALSE]
sample_df$read_depth_bacteria <- colSums(otu_mat[, sample_df$sample_id, drop = FALSE])

sample_df$project <- trimws(as.character(sample_df$project))
sample_df$os_time <- if ("OS_months(汇总)" %in% names(sample_df)) {
  coerce_numeric(sample_df[["OS_months(汇总)"]])
} else {
  coerce_numeric(sample_df$Overall_survival_days)
}
sample_df$death_event <- coerce_event(sample_df$Death_Event)
sample_df$sex <- factor(trimws(as.character(sample_df$Host_sex)), levels = c("Female", "Male"))
sample_df$age_years <- coerce_numeric(if ("Host_age" %in% names(sample_df)) sample_df$Host_age else sample_df$Age)
sample_df$bmi <- coerce_numeric(sample_df$BMI)

keep_samples <- sample_df$Host == "Human" &
  sample_df$AssayType == "WGS" &
  sample_df$SPT == "Baseline" &
  sample_df$project %in% target_projects &
  sample_df$read_depth_bacteria >= depth_cutoff &
  is.finite(sample_df$os_time) &
  sample_df$os_time > 0 &
  sample_df$death_event %in% c(0L, 1L) &
  !is.na(sample_df$sex) &
  is.finite(sample_df$age_years) &
  is.finite(sample_df$bmi)

sample_df <- sample_df[keep_samples, , drop = FALSE]
sample_df$project <- factor(sample_df$project, levels = target_projects)
otu_mat <- otu_mat[, sample_df$sample_id, drop = FALSE]
rel_abund <- sweep(otu_mat, 2, colSums(otu_mat), "/")

species_filter_path <- file.path(tables_dir, "pooled_species_filter_summary.csv")
if (!file.exists(species_filter_path)) {
  stop("Missing pooled_species_filter_summary.csv; run targeted two-cohort analysis first.")
}

species_filter_df <- fread(species_filter_path) %>%
  filter(domain == "pooled", feature %in% rownames(rel_abund)) %>%
  arrange(desc(mean_relative_abundance), feature_short)

if (nrow(species_filter_df) == 0) {
  stop("No pooled species remained for nonlinear Cox analysis.")
}

result_rows <- list()
curve_fits <- list()

message(sprintf("Running adjusted nonlinear Cox models for %d species", nrow(species_filter_df)))
for (row_idx in seq_len(nrow(species_filter_df))) {
  feature_id <- species_filter_df$feature[row_idx]
  feature_short <- species_filter_df$feature_short[row_idx]
  abundance <- as.numeric(rel_abund[feature_id, sample_df$sample_id])
  model_df <- sample_df %>%
    transmute(
      os_time = os_time,
      death_event = death_event,
      project = project,
      sex = sex,
      age_years = age_years,
      bmi = bmi,
      abundance = abundance,
      abundance_log10 = log10(abundance + pseudocount)
    ) %>%
    filter(is.finite(abundance_log10))

  if (nrow(model_df) < 60 || sum(model_df$death_event) < 20 || length(unique(model_df$abundance_log10)) < 6) {
    next
  }

  null_fit <- tryCatch(
    coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi, data = model_df, ties = "efron"),
    error = function(e) NULL
  )
  linear_fit <- tryCatch(
    coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi + abundance_log10, data = model_df, ties = "efron"),
    error = function(e) NULL
  )
  spline_fit <- tryCatch(
    coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi + ns(abundance_log10, df = spline_df), data = model_df, ties = "efron"),
    error = function(e) NULL
  )

  if (is.null(null_fit) || is.null(linear_fit) || is.null(spline_fit)) {
    next
  }

  linear_summary <- summary(linear_fit)
  overall_lrt <- tryCatch(anova(null_fit, spline_fit, test = "LRT"), error = function(e) NULL)
  nonlinear_lrt <- tryCatch(anova(linear_fit, spline_fit, test = "LRT"), error = function(e) NULL)
  overall_p <- if (!is.null(overall_lrt)) safe_p(overall_lrt) else NA_real_
  nonlinear_p <- if (!is.null(nonlinear_lrt)) safe_p(nonlinear_lrt) else NA_real_

  linear_term <- linear_summary$coefficients["abundance_log10", , drop = FALSE]
  linear_conf <- linear_summary$conf.int["abundance_log10", , drop = FALSE]

  result_rows[[length(result_rows) + 1L]] <- data.frame(
    feature = feature_id,
    feature_short = feature_short,
    phylum = species_filter_df$phylum[row_idx],
    prevalence = species_filter_df$prevalence[row_idx],
    mean_relative_abundance = species_filter_df$mean_relative_abundance[row_idx],
    n = nrow(model_df),
    events = sum(model_df$death_event),
    abundance_min = min(model_df$abundance),
    abundance_median = median(model_df$abundance),
    abundance_max = max(model_df$abundance),
    linear_hr_per_log10 = unname(linear_conf[, "exp(coef)"]),
    linear_conf_low = unname(linear_conf[, "lower .95"]),
    linear_conf_high = unname(linear_conf[, "upper .95"]),
    linear_p = unname(linear_term[, "Pr(>|z|)"]),
    spline_overall_p = overall_p,
    nonlinear_p = nonlinear_p,
    aic_null = AIC(null_fit),
    aic_linear = AIC(linear_fit),
    aic_spline = AIC(spline_fit),
    delta_aic_spline_vs_linear = AIC(spline_fit) - AIC(linear_fit),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  curve_fits[[feature_short]] <- list(fit = spline_fit, data = model_df)
}

summary_df <- bind_rows(result_rows) %>%
  mutate(
    linear_fdr = p.adjust(linear_p, method = "BH"),
    spline_overall_fdr = p.adjust(spline_overall_p, method = "BH"),
    nonlinear_fdr = p.adjust(nonlinear_p, method = "BH")
  ) %>%
  arrange(nonlinear_p, spline_overall_p, linear_p, feature_short)

fwrite(summary_df, file.path(tables_dir, "pooled_os_microbiome_nonlinear_cox_summary.csv"))

if (nrow(summary_df) == 0) {
  stop("No nonlinear Cox models were successfully fitted.")
}

plot_summary <- summary_df %>%
  slice_min(order_by = nonlinear_p, n = min(25L, nrow(summary_df)), with_ties = FALSE) %>%
  mutate(
    feature_short = factor(feature_short, levels = rev(feature_short)),
    neg_log10_nonlinear_p = -log10(pmax(nonlinear_p, .Machine$double.xmin)),
    direction = ifelse(delta_aic_spline_vs_linear < 0, "Spline better", "Linear better/equal")
  )

nonlinear_bar <- ggplot(plot_summary, aes(x = neg_log10_nonlinear_p, y = feature_short, fill = direction)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
  scale_fill_manual(values = c("Spline better" = "#A35C72", "Linear better/equal" = "#9AA3A8")) +
  labs(
    title = "Nonlinear Cox screen for OS-associated microbial abundance features",
    subtitle = sprintf("Adjusted for cohort, sex, age, BMI; ns(log10 abundance + %.0e), df=%d", pseudocount, spline_df),
    x = expression(-log[10]("nonlinear LRT p")),
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(figures_dir, "pooled_os_microbiome_nonlinear_cox_barplot.pdf"), nonlinear_bar, width = 9.5, height = 7.2)
ggsave(file.path(figures_dir, "pooled_os_microbiome_nonlinear_cox_barplot.png"), nonlinear_bar, width = 9.5, height = 7.2, dpi = 320)

linear_nonlinear_df <- summary_df %>%
  mutate(
    feature_class = case_when(
      linear_fdr <= 0.05 & nonlinear_fdr <= 0.05 ~ "Linear + Nonlinear",
      linear_fdr <= 0.05 ~ "Linear only",
      nonlinear_fdr <= 0.05 ~ "Nonlinear only",
      TRUE ~ "Not significant"
    ),
    linear_log2_hr = log2(linear_hr_per_log10),
    neg_log10_linear_p = -log10(pmax(linear_p, .Machine$double.xmin)),
    spline_gain = -delta_aic_spline_vs_linear,
    neg_log10_nonlinear_p = -log10(pmax(nonlinear_p, .Machine$double.xmin)),
    is_akk = feature_short == "Akkermansia muciniphila"
  )

linear_nonlinear_palette <- c(
  "Not significant" = "#BEBEBE",
  "Linear only" = "#3C779F",
  "Nonlinear only" = "#A35C72",
  "Linear + Nonlinear" = "#4D9D59"
)

p_linear <- ggplot(linear_nonlinear_df, aes(x = linear_log2_hr, y = neg_log10_linear_p, color = feature_class)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(alpha = 0.8, size = 2) +
  geom_point(data = linear_nonlinear_df %>% filter(is_akk), shape = 21, fill = "#F7F7F7", color = "black", size = 3.2, stroke = 0.9) +
  geom_text(
    data = linear_nonlinear_df %>% filter(is_akk),
    aes(label = "AKK"),
    nudge_x = 0.08,
    nudge_y = 0.15,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_color_manual(values = linear_nonlinear_palette) +
  labs(
    title = "Linear association summary",
    subtitle = "Adjusted Cox: project + sex + age + BMI + log10(abundance + 1e-6)",
    x = expression(log[2]("HR per 1-unit log10 abundance")),
    y = expression(-log[10]("linear p")),
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

p_nonlinear <- ggplot(linear_nonlinear_df, aes(x = spline_gain, y = neg_log10_nonlinear_p, color = feature_class)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(alpha = 0.8, size = 2) +
  geom_point(data = linear_nonlinear_df %>% filter(is_akk), shape = 21, fill = "#F7F7F7", color = "black", size = 3.2, stroke = 0.9) +
  geom_text(
    data = linear_nonlinear_df %>% filter(is_akk),
    aes(label = "AKK"),
    nudge_x = 0.15,
    nudge_y = 0.15,
    size = 3,
    color = "black",
    check_overlap = TRUE
  ) +
  scale_color_manual(values = linear_nonlinear_palette) +
  labs(
    title = "Nonlinear association summary",
    subtitle = "LRT: spline(df=3) vs linear model",
    x = "Spline improvement: -Delta AIC (spline - linear)",
    y = expression(-log[10]("nonlinear p")),
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

linear_nonlinear_panel <- p_linear + p_nonlinear +
  plot_annotation(
    title = "OS-associated microbial abundance features: linear vs nonlinear summary",
    subtitle = "AKK is highlighted to bridge split-based and spline-based findings"
  )

fwrite(
  linear_nonlinear_df %>%
    select(feature, feature_short, phylum, linear_hr_per_log10, linear_p, linear_fdr, nonlinear_p, nonlinear_fdr, delta_aic_spline_vs_linear, feature_class),
  file.path(tables_dir, "pooled_os_microbiome_linear_nonlinear_summary.csv")
)

ggsave(file.path(figures_dir, "pooled_os_microbiome_linear_nonlinear_summary_panel.pdf"), linear_nonlinear_panel, width = 13, height = 6.2)
ggsave(file.path(figures_dir, "pooled_os_microbiome_linear_nonlinear_summary_panel.png"), linear_nonlinear_panel, width = 13, height = 6.2, dpi = 320)

make_curve_plot <- function(curve_df, title_text, subtitle_text = NULL) {
  ggplot(curve_df, aes(x = abundance_plus_pseudocount, y = hr)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = "#A35C72", alpha = 0.18) +
    geom_line(color = "#A35C72", linewidth = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
    scale_x_continuous(trans = "log10", labels = label_number(accuracy = 0.0001)) +
    scale_y_continuous(trans = "log10") +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = expression("Relative abundance" + 10^-6),
      y = "Adjusted hazard ratio vs median abundance"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
}

akk_name <- "Akkermansia muciniphila"
if (akk_name %in% names(curve_fits)) {
  akk_curve <- make_prediction_curve(curve_fits[[akk_name]]$fit, curve_fits[[akk_name]]$data, akk_name)
  akk_row <- summary_df %>% filter(feature_short == akk_name) %>% slice_head(n = 1)
  fwrite(akk_curve, file.path(tables_dir, "pooled_akk_os_nonlinear_cox_curve.csv"))
  akk_plot <- make_curve_plot(
    akk_curve,
    title_text = "AKK abundance and OS risk: adjusted nonlinear Cox model",
    subtitle_text = sprintf(
      "Spline overall p=%.3g; nonlinear p=%.3g; linear HR/log10=%.2f",
      akk_row$spline_overall_p,
      akk_row$nonlinear_p,
      akk_row$linear_hr_per_log10
    )
  )
  ggsave(file.path(figures_dir, "pooled_akk_os_nonlinear_cox_curve.pdf"), akk_plot, width = 7.5, height = 5.2)
  ggsave(file.path(figures_dir, "pooled_akk_os_nonlinear_cox_curve.png"), akk_plot, width = 7.5, height = 5.2, dpi = 320)
}

top_curve_species <- summary_df %>%
  filter(feature_short %in% names(curve_fits)) %>%
  slice_min(order_by = nonlinear_p, n = min(top_curve_n, nrow(.)), with_ties = FALSE) %>%
  pull(feature_short)

curve_df <- bind_rows(lapply(top_curve_species, function(species_name) {
  make_prediction_curve(curve_fits[[species_name]]$fit, curve_fits[[species_name]]$data, species_name)
}))

if (nrow(curve_df) > 0) {
  curve_panel <- ggplot(curve_df, aes(x = abundance_plus_pseudocount, y = hr)) +
    geom_ribbon(aes(ymin = conf_low, ymax = conf_high), fill = "#3C779F", alpha = 0.15) +
    geom_line(color = "#3C779F", linewidth = 0.75) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey45") +
    facet_wrap(~feature_short, scales = "free_x", ncol = 2) +
    scale_x_continuous(trans = "log10", labels = label_number(accuracy = 0.0001)) +
    scale_y_continuous(trans = "log10") +
    labs(
      title = sprintf("Top %d nonlinear microbiome-OS Cox curves", length(top_curve_species)),
      subtitle = "Adjusted HR relative to each species' median abundance",
      x = expression("Relative abundance" + 10^-6),
      y = "Adjusted hazard ratio"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

  fwrite(curve_df, file.path(tables_dir, "pooled_top_nonlinear_os_cox_curves.csv"))
  ggsave(file.path(figures_dir, "pooled_top_nonlinear_os_cox_curves.pdf"), curve_panel, width = 11, height = 9)
  ggsave(file.path(figures_dir, "pooled_top_nonlinear_os_cox_curves.png"), curve_panel, width = 11, height = 9, dpi = 320)
}

message("Finished nonlinear microbiome-OS Cox analysis")