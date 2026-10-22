#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

options(stringsAsFactors = FALSE, scipen = 999)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("^--file=", args)][1])
source_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

plot_df <- read.csv(file.path(source_dir, "figure_plot_data.csv"), check.names = FALSE)
sample_df <- read.csv(gzfile(file.path(source_dir, "nonlinear_cox_sample_level_source.csv.gz")), check.names = FALSE)

draw_barplot <- function(out_path, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") pdf(out_path, width = 9.5, height = 7.2, onefile = FALSE) else png(out_path, width = 9.5, height = 7.2, units = "in", res = 320)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  df <- plot_df[order(plot_df$neg_log10_nonlinear_p), ]
  cols <- ifelse(df$direction == "Spline better", "#A35C72", "#9AA3A8")
  par(mar = c(5, 12, 4, 2))
  barplot(
    df$neg_log10_nonlinear_p,
    names.arg = df$feature_short,
    horiz = TRUE,
    las = 1,
    col = cols,
    border = NA,
    xlab = "-log10(nonlinear LRT p)",
    main = "Nonlinear Cox screen for OS-associated microbial abundance features",
    cex.names = 0.72
  )
  abline(v = -log10(0.05), lty = 2, col = "grey45")
  legend("bottomright", legend = c("Spline better", "Linear better/equal"), fill = c("#A35C72", "#9AA3A8"), bty = "n")
}

safe_p <- function(anova_obj) {
  p_col <- grep("^Pr|^P\\(", colnames(anova_obj), value = TRUE)
  if (length(p_col) == 0 || nrow(anova_obj) < 2) return(NA_real_)
  as.numeric(anova_obj[2, p_col[1]])
}

fit_feature <- function(df) {
  df$project <- factor(df$project, levels = c("PRJNA751792", "PRJNA762360"))
  df$sex <- factor(df$sex, levels = c("Female", "Male"))
  null_fit <- coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi, data = df, ties = "efron")
  linear_fit <- coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi + abundance_log10, data = df, ties = "efron")
  spline_fit <- coxph(Surv(os_time, death_event) ~ project + sex + age_years + bmi + ns(abundance_log10, df = 3), data = df, ties = "efron")
  linear_summary <- summary(linear_fit)
  data.frame(
    feature = df$feature[1],
    feature_short = df$feature_short[1],
    linear_p = unname(linear_summary$coefficients["abundance_log10", "Pr(>|z|)"]),
    spline_overall_p = safe_p(anova(null_fit, spline_fit, test = "LRT")),
    nonlinear_p = safe_p(anova(linear_fit, spline_fit, test = "LRT")),
    delta_aic_spline_vs_linear = AIC(spline_fit) - AIC(linear_fit)
  )
}

draw_barplot(file.path(source_dir, "reproduced_pooled_os_microbiome_nonlinear_cox_barplot.pdf"), "pdf")
draw_barplot(file.path(source_dir, "reproduced_pooled_os_microbiome_nonlinear_cox_barplot.png"), "png")

features <- unique(plot_df$feature)
rows <- lapply(features, function(f) fit_feature(sample_df[sample_df$feature == f, , drop = FALSE]))
recomputed <- do.call(rbind, rows)
write.csv(recomputed, file.path(source_dir, "recomputed_top25_nonlinear_cox_statistics.csv"), row.names = FALSE)
message("Wrote reproduced figure and recomputed top-25 statistics to: ", source_dir)
