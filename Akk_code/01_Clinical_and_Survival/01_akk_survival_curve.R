#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
})

options(stringsAsFactors = FALSE, scipen = 999)
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("^--file=", args)][1])
source_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else getwd()

curve_df <- read.csv(file.path(source_dir, "figure_plot_data.csv"), check.names = FALSE)
assign_df <- read.csv(gzfile(file.path(source_dir, "km_sample_level_source.csv.gz")), check.names = FALSE)
assign_df$group <- factor(assign_df$group, levels = c("G1", "G2", "G3"))
group_cols <- c(G1 = "#A35C72", G2 = "#3C779F", G3 = "#86D2B2")

draw_km <- function(out_path, type = c("pdf", "png")) {
  type <- match.arg(type)
  if (type == "pdf") pdf(out_path, width = 8, height = 5.5, onefile = FALSE) else png(out_path, width = 8, height = 5.5, units = "in", res = 320)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  par(mar = c(5, 5, 4, 2))
  plot(NA, xlim = range(curve_df$time), ylim = c(0, 1), xlab = "Overall survival (months)", ylab = "Survival probability", main = "pooled | AKK tertile overall survival", las = 1, bty = "l")
  for (g in names(group_cols)) {
    df <- curve_df[curve_df$strata == g, ]
    lines(df$time, df$surv, type = "s", lwd = 2, col = group_cols[g])
  }
  counts <- aggregate(death_event ~ group, assign_df, function(x) c(n = length(x), events = sum(x)))
  labels <- vapply(seq_len(nrow(counts)), function(i) sprintf("%s (n=%d, events=%d)", counts$group[i], counts$death_event[i, "n"], counts$death_event[i, "events"]), character(1))
  legend("bottomleft", legend = labels, col = group_cols[as.character(counts$group)], lwd = 2, bty = "n")
}

fit <- survdiff(Surv(os_days, death_event) ~ group, data = assign_df)
recomputed <- data.frame(logrank_chisq = unname(fit$chisq), df = length(fit$n) - 1L, logrank_p = pchisq(unname(fit$chisq), df = length(fit$n) - 1L, lower.tail = FALSE))
write.csv(recomputed, file.path(source_dir, "recomputed_logrank_from_script.csv"), row.names = FALSE)
draw_km(file.path(source_dir, "reproduced_pooled_akk_tertile_km.pdf"), "pdf")
draw_km(file.path(source_dir, "reproduced_pooled_akk_tertile_km.png"), "png")
message("Wrote reproduced KM plot and recomputed log-rank to: ", source_dir)
