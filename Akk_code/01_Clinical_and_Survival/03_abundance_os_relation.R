#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", script_args[grep(file_arg, script_args)][1])
source_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

plot_data_path <- file.path(source_dir, "figure_plot_data.csv")
assignment_path <- file.path(source_dir, "statistical_source_group_assignments_binary_tertile.csv.gz")
split_summary_path <- file.path(source_dir, "logrank_split_summary_binary_tertile.csv")

plot_df <- read.csv(plot_data_path, check.names = FALSE)
assign_df <- read.csv(gzfile(assignment_path), check.names = FALSE)
split_summary_df <- read.csv(split_summary_path, check.names = FALSE)

group_cols <- c(G1 = "#A35C72", G2 = "#3C779F", G3 = "#86D2B2", G4 = "#81017E")
split_levels <- c("binary", "tertile")

scale_size <- function(x, min_cex = 0.7, max_cex = 3.0) {
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep((min_cex + max_cex) / 2, length(x)))
  }
  min_cex + (x - rng[1]) / diff(rng) * (max_cex - min_cex)
}

draw_one_split <- function(df, split_name) {
  split_df <- df[df$split_label == split_name, , drop = FALSE]
  if (nrow(split_df) == 0) {
    plot.new()
    title(main = split_name)
    return(invisible(NULL))
  }

  x <- split_df$hr
  y <- split_df$neg_log10_adjusted_p
  cex <- scale_size(split_df$mean_relative_abundance)
  cols <- group_cols[as.character(split_df$best_group)]
  cols[is.na(cols)] <- "grey40"
  x_lim <- range(c(1, x[is.finite(x)]), na.rm = TRUE)
  y_lim <- range(c(0, -log10(0.05), y[is.finite(y)]), na.rm = TRUE)
  x_pad <- diff(x_lim) * 0.08
  y_pad <- diff(y_lim) * 0.10

  plot(
    x, y,
    pch = 16,
    cex = cex,
    col = adjustcolor(cols, alpha.f = 0.85),
    xlim = c(max(0, x_lim[1] - x_pad), x_lim[2] + x_pad),
    ylim = c(max(0, y_lim[1] - y_pad), y_lim[2] + y_pad),
    xlab = "Hazard ratio of best OS group vs all other groups",
    ylab = "-log10(adjusted p value)",
    main = if (split_name == "binary") "Binary split" else "Tertile split",
    las = 1,
    bty = "l"
  )
  abline(v = 1, lty = 2, col = "grey65")
  abline(h = -log10(0.05), lty = 3, col = "grey70")

  label_df <- split_df[split_df$is_labeled_in_figure %in% c(TRUE, "TRUE", "True", "true"), , drop = FALSE]
  if (nrow(label_df) > 0) {
    text(label_df$hr, label_df$neg_log10_adjusted_p, labels = label_df$feature_short, pos = 3, cex = 0.55)
  }
}

draw_bubble_figure <- function(out_path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(out_path, width = 15, height = 6.5, onefile = FALSE)
  } else {
    png(out_path, width = 15, height = 6.5, units = "in", res = 320)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.32))
  par(mar = c(5, 5, 4, 1))
  draw_one_split(plot_df, "binary")
  draw_one_split(plot_df, "tertile")
  par(mar = c(5, 0.5, 4, 1))
  plot.new()
  legend(
    "topleft",
    legend = names(group_cols),
    col = group_cols,
    pch = 16,
    pt.cex = 1.5,
    bty = "n",
    title = "Group"
  )
  legend(
    "left",
    legend = c("low", "medium", "high"),
    pch = 16,
    pt.cex = c(0.7, 1.8, 3.0),
    col = "grey35",
    bty = "n",
    title = "Mean relative abundance"
  )
  mtext("pooled | Best OS group bubble plot adjusted for sex, age, BMI", outer = TRUE, line = -1.5, font = 2, cex = 1.2)
}

draw_bubble_figure(file.path(source_dir, "reproduced_pooled_best_os_group_bubble_adjusted_sex_age_bmi.pdf"), "pdf")
draw_bubble_figure(file.path(source_dir, "reproduced_pooled_best_os_group_bubble_adjusted_sex_age_bmi.png"), "png")

cox_one_best_group <- function(df, best_group) {
  model_df <- df[!is.na(df$sex) & is.finite(df$age_years) & is.finite(df$bmi), , drop = FALSE]
  if (nrow(model_df) == 0) {
    return(data.frame(hr = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, model_n = 0L, model_events = 0L))
  }
  model_df$focus_group <- factor(ifelse(model_df$group == best_group, best_group, "Other"), levels = c("Other", best_group))
  fit <- tryCatch(
    coxph(Surv(os_days, death_event) ~ focus_group + sex + age_years + bmi, data = model_df, ties = "efron"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(hr = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, model_n = nrow(model_df), model_events = sum(model_df$death_event)))
  }
  fit_summary <- summary(fit)
  coef_name <- grep("^focus_group", rownames(fit_summary$coefficients), value = TRUE)[1]
  if (is.na(coef_name)) {
    return(data.frame(hr = NA_real_, p_value = NA_real_, conf_low = NA_real_, conf_high = NA_real_, model_n = nrow(model_df), model_events = sum(model_df$death_event)))
  }
  data.frame(
    hr = unname(fit_summary$coefficients[coef_name, "exp(coef)"]),
    p_value = unname(fit_summary$coefficients[coef_name, "Pr(>|z|)"]),
    conf_low = unname(fit_summary$conf.int[coef_name, "lower .95"]),
    conf_high = unname(fit_summary$conf.int[coef_name, "upper .95"]),
    model_n = nrow(model_df),
    model_events = sum(model_df$death_event)
  )
}

logrank_one_split <- function(df) {
  fit <- tryCatch(
    survdiff(Surv(os_days, death_event) ~ group, data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(logrank_chisq = NA_real_, logrank_p = NA_real_))
  }
  data.frame(
    logrank_chisq = unname(fit$chisq),
    logrank_p = pchisq(unname(fit$chisq), df = length(fit$n) - 1L, lower.tail = FALSE)
  )
}

keys <- unique(plot_df[, c("feature", "feature_short", "split_label", "best_group")])
cox_rows <- vector("list", nrow(keys))
logrank_rows <- vector("list", nrow(keys))

for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]
  current <- assign_df[assign_df$feature == key$feature & assign_df$split_label == key$split_label, , drop = FALSE]
  cox_rows[[i]] <- cbind(key, cox_one_best_group(current, key$best_group))
  logrank_rows[[i]] <- cbind(key[1:3], logrank_one_split(current))
}

recomputed_cox <- do.call(rbind, cox_rows)
recomputed_logrank <- do.call(rbind, logrank_rows)
recomputed_cox$best_group_fdr <- ave(
  recomputed_cox$p_value,
  recomputed_cox$split_label,
  FUN = function(x) p.adjust(x, method = "BH")
)

write.csv(recomputed_cox, file.path(source_dir, "recomputed_adjusted_best_group_cox.csv"), row.names = FALSE)
write.csv(recomputed_logrank, file.path(source_dir, "recomputed_logrank_binary_tertile.csv"), row.names = FALSE)

message("Wrote reproduced figure and recomputed statistics to: ", source_dir)
