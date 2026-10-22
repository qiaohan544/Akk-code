#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(splines)
})

options(stringsAsFactors = FALSE, scipen = 999)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("^--file=", args)][1])
source_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

curve_df <- read.csv(file.path(source_dir, "figure_plot_data.csv"), check.names = FALSE)
sample_df <- read.csv(gzfile(file.path(source_dir, "top_curve_species_sample_level_source.csv.gz")), check.names = FALSE)

draw_curves <- function(out_path, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") pdf(out_path, width = 11, height = 9, onefile = FALSE) else png(out_path, width = 11, height = 9, units = "in", res = 320)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  species <- unique(curve_df$feature_short)
  par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))
  for (sp in species) {
    df <- curve_df[curve_df$feature_short == sp, ]
    df <- df[order(df$abundance_plus_pseudocount), ]
    x <- df$abundance_plus_pseudocount
    plot(
      x, df$hr,
      type = "n",
      log = "xy",
      xlab = "Relative abundance + 1e-6",
      ylab = "Adjusted hazard ratio",
      main = sp,
      las = 1,
      bty = "l",
      ylim = range(c(df$conf_low, df$conf_high), na.rm = TRUE)
    )
    polygon(c(x, rev(x)), c(df$conf_low, rev(df$conf_high)), col = adjustcolor("#3C779F", alpha.f = 0.15), border = NA)
    lines(x, df$hr, col = "#3C779F", lwd = 2)
    abline(h = 1, lty = 2, col = "grey45")
  }
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

draw_curves(file.path(source_dir, "reproduced_pooled_top_nonlinear_os_cox_curves.pdf"), "pdf")
draw_curves(file.path(source_dir, "reproduced_pooled_top_nonlinear_os_cox_curves.png"), "png")

features <- unique(sample_df$feature)
rows <- lapply(features, function(f) fit_feature(sample_df[sample_df$feature == f, , drop = FALSE]))
recomputed <- do.call(rbind, rows)
write.csv(recomputed, file.path(source_dir, "recomputed_top_curve_species_nonlinear_cox_statistics.csv"), row.names = FALSE)
message("Wrote reproduced curves and recomputed top-curve statistics to: ", source_dir)
