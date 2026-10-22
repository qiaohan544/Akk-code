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

palette <- c(
  "Not significant" = "#BEBEBE",
  "Linear only" = "#3C779F",
  "Nonlinear only" = "#A35C72",
  "Linear + Nonlinear" = "#4D9D59"
)

draw_panel <- function(out_path, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") pdf(out_path, width = 13, height = 6.2, onefile = FALSE) else png(out_path, width = 13, height = 6.2, units = "in", res = 320)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  layout(matrix(c(1, 2, 3), 1), widths = c(1, 1, 0.35))
  cols <- palette[plot_df$feature_class]

  par(mar = c(5, 5, 4, 1))
  plot(plot_df$linear_log2_hr, plot_df$neg_log10_linear_p, pch = 16, col = cols, xlab = "log2(HR per 1-unit log10 abundance)", ylab = "-log10(linear p)", main = "Linear association summary", las = 1, bty = "l")
  abline(h = -log10(0.05), lty = 2, col = "grey50")
  abline(v = 0, lty = 3, col = "grey60")
  akk <- plot_df$feature_short == "Akkermansia muciniphila"
  points(plot_df$linear_log2_hr[akk], plot_df$neg_log10_linear_p[akk], pch = 21, bg = "#F7F7F7", col = "black", cex = 1.3)
  text(plot_df$linear_log2_hr[akk], plot_df$neg_log10_linear_p[akk], labels = "AKK", pos = 4, cex = 0.8)

  plot(plot_df$spline_gain, plot_df$neg_log10_nonlinear_p, pch = 16, col = cols, xlab = "Spline improvement: -Delta AIC (spline - linear)", ylab = "-log10(nonlinear p)", main = "Nonlinear association summary", las = 1, bty = "l")
  abline(h = -log10(0.05), lty = 2, col = "grey50")
  abline(v = 0, lty = 3, col = "grey60")
  points(plot_df$spline_gain[akk], plot_df$neg_log10_nonlinear_p[akk], pch = 21, bg = "#F7F7F7", col = "black", cex = 1.3)
  text(plot_df$spline_gain[akk], plot_df$neg_log10_nonlinear_p[akk], labels = "AKK", pos = 4, cex = 0.8)

  par(mar = c(5, 0.5, 4, 1))
  plot.new()
  legend("topleft", legend = names(palette), col = palette, pch = 16, bty = "n", title = "Feature class")
  mtext("OS-associated microbial abundance features: linear vs nonlinear summary", outer = TRUE, line = -1.4, font = 2, cex = 1.2)
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
  linear_conf <- linear_summary$conf.int["abundance_log10", , drop = FALSE]
  linear_term <- linear_summary$coefficients["abundance_log10", , drop = FALSE]
  data.frame(
    feature = df$feature[1],
    feature_short = df$feature_short[1],
    linear_hr_per_log10 = unname(linear_conf[, "exp(coef)"]),
    linear_p = unname(linear_term[, "Pr(>|z|)"]),
    spline_overall_p = safe_p(anova(null_fit, spline_fit, test = "LRT")),
    nonlinear_p = safe_p(anova(linear_fit, spline_fit, test = "LRT")),
    delta_aic_spline_vs_linear = AIC(spline_fit) - AIC(linear_fit)
  )
}

draw_panel(file.path(source_dir, "reproduced_pooled_os_microbiome_linear_nonlinear_summary_panel.pdf"), "pdf")
draw_panel(file.path(source_dir, "reproduced_pooled_os_microbiome_linear_nonlinear_summary_panel.png"), "png")

features <- unique(sample_df$feature)
rows <- lapply(features, function(f) fit_feature(sample_df[sample_df$feature == f, , drop = FALSE]))
recomputed <- do.call(rbind, rows)
recomputed$linear_fdr <- p.adjust(recomputed$linear_p, method = "BH")
recomputed$nonlinear_fdr <- p.adjust(recomputed$nonlinear_p, method = "BH")
write.csv(recomputed, file.path(source_dir, "recomputed_nonlinear_cox_statistics.csv"), row.names = FALSE)
message("Wrote reproduced figure and recomputed statistics to: ", source_dir)
