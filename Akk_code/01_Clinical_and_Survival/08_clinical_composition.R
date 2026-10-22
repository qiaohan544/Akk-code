#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("^--file=", args)][1])
source_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

plot_df <- read.csv(file.path(source_dir, "figure_plot_data.csv"), check.names = FALSE)
sample_df <- read.csv(gzfile(file.path(source_dir, "akk_tertile_composition_sample_source.csv.gz")), check.names = FALSE)
sample_df$group <- factor(sample_df$group, levels = c("G1", "G2", "G3"))

vars <- c("os_12_month", "response", "antibiotics", "sex", "age_group", "bmi_group")
titles <- c(os_12_month = "OS_12_month", response = "Response", antibiotics = "Antibiotics", sex = "Sex", age_group = "Age group", bmi_group = "BMI group")
palette <- c("#A35C72", "#F1E4E8", "#3C779F", "#E1E8F1", "#86D2B2", "#CAE4D8", "#81017E", "#A0B3BC")

recompute_one <- function(variable) {
  df <- sample_df[!is.na(sample_df$group) & !is.na(sample_df[[variable]]) & sample_df[[variable]] != "", ]
  tab <- table(df$group, df[[variable]])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.frame(variable = variable, p_value = NA_real_, method = "insufficient_levels"))
  }
  fit <- suppressWarnings(chisq.test(tab))
  if (!is.finite(fit$p.value) || any(fit$expected < 5)) {
    fit <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 5000))
  }
  data.frame(variable = variable, p_value = fit$p.value, method = fit$method)
}

draw_panel <- function(out_path, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") pdf(out_path, width = 13, height = 8.5, onefile = FALSE) else png(out_path, width = 13, height = 8.5, units = "in", res = 320)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  par(mfrow = c(2, 3), mar = c(4, 4, 4, 1))
  for (v in vars) {
    sub <- plot_df[plot_df$variable == v, ]
    vals <- unique(sub$value)
    mat <- matrix(0, nrow = length(vals), ncol = 3, dimnames = list(vals, c("G1", "G2", "G3")))
    for (g in colnames(mat)) {
      gd <- sub[sub$group == g, ]
      denom <- sum(gd$Freq)
      if (denom > 0) mat[as.character(gd$value), g] <- gd$Freq / denom
    }
    barplot(mat, beside = FALSE, col = palette[seq_len(nrow(mat))], ylim = c(0, 1), ylab = "Composition", main = titles[v], las = 1)
    legend("bottom", legend = rownames(mat), fill = palette[seq_len(nrow(mat))], bty = "n", cex = 0.75)
  }
}

stats <- do.call(rbind, lapply(vars, recompute_one))
write.csv(stats, file.path(source_dir, "recomputed_chisq_results.csv"), row.names = FALSE)
draw_panel(file.path(source_dir, "reproduced_pooled_akk_tertile_composition_panel.pdf"), "pdf")
draw_panel(file.path(source_dir, "reproduced_pooled_akk_tertile_composition_panel.png"), "png")
message("Wrote reproduced composition panel and recomputed chi-square results to: ", source_dir)
