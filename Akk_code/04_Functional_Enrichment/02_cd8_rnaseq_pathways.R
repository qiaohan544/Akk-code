suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
table_dir <- file.path("rnaseq_results", "tables")
fig_root <- file.path("rnaseq_results", "figures")
fig_dir <- file.path(fig_root, "GSEA_enrichment_plots")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gmt_file <- "m5.go.bp.v2026.1.Mm.symbols.gmt"
de_file <- file.path(table_dir, "DESeq2_Isobutyrate_vs_NC_all_genes.csv")
gsea_file <- file.path(table_dir, "GSEA_GO_BP_all_pathways.csv")

selected_ids <- c(
  "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
  "GOBP_RESPONSE_TO_INTERLEUKIN_2",
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
  "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION"
)

de <- read.csv(de_file, check.names = FALSE)
rank_df <- de[!is.na(de$symbol) & nzchar(de$symbol) & is.finite(de$stat),
              c("symbol", "stat")]
rank_df <- rank_df[order(abs(rank_df$stat), decreasing = TRUE), ]
rank_df <- rank_df[!duplicated(rank_df$symbol), ]
ranks <- rank_df$stat
names(ranks) <- rank_df$symbol
ranks <- sort(ranks, decreasing = TRUE)

pathways <- gmtPathways(gmt_file)
missing_ids <- setdiff(selected_ids, names(pathways))
if (length(missing_ids)) stop("Pathways missing from GMT: ", paste(missing_ids, collapse = ", "))

gsea <- read.csv(gsea_file, check.names = FALSE)
gsea_selected <- gsea[match(selected_ids, gsea$pathway), ]
if (anyNA(gsea_selected$pathway)) stop("Selected pathways missing from GSEA results.")

pretty_path <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- tools::toTitleCase(tolower(gsub("_", " ", x)))
  x <- gsub("t Cell", "T-cell", x, fixed = TRUE)
  x <- gsub("T Cell", "T-cell", x, fixed = TRUE)
  x <- gsub("Interleukin 2", "interleukin-2", x, fixed = TRUE)
  x <- gsub("Leukocyte", "Leukocyte", x, fixed = TRUE)
  paste(strwrap(x, width = 58), collapse = "\n")
}

theme_gsea <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4, color = "black"),
        axis.ticks = element_line(linewidth = 0.35),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 9),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        plot.subtitle = element_text(size = 8.5, color = "#333333", hjust = 0),
        plot.margin = margin(3, 7, 3, 7))

calculate_running_es <- function(stats, genes, exponent = 1) {
  hits <- names(stats) %in% genes
  nh <- sum(hits)
  n <- length(stats)
  if (nh == 0 || nh == n) stop("Invalid gene-set overlap with ranked list.")
  hit_weights <- abs(stats)^exponent
  hit_weights[!hits] <- 0
  hit_weights <- hit_weights / sum(hit_weights)
  miss_weights <- rep(1 / (n - nh), n)
  running <- cumsum(ifelse(hits, hit_weights, -miss_weights))
  peak <- if (abs(max(running)) >= abs(min(running))) which.max(running) else which.min(running)
  list(running = running, hits = hits, hit_ranks = which(hits), peak = peak,
       ES = running[peak])
}

make_gsea_plot <- function(pathway_id) {
  stat_row <- gsea_selected[gsea_selected$pathway == pathway_id, ][1, ]
  calc <- calculate_running_es(ranks, intersect(pathways[[pathway_id]], names(ranks)))
  n <- length(ranks)
  curve_df <- data.frame(rank = seq_len(n), running_ES = calc$running)
  hit_df <- data.frame(rank = calc$hit_ranks)

  # Downsample only the ranked-metric panel to keep vector PDFs compact.
  metric_idx <- unique(round(seq(1, n, length.out = min(n, 3000))))
  metric_df <- data.frame(rank = metric_idx, stat = unname(ranks[metric_idx]))
  display_limits <- quantile(unname(ranks), c(0.01, 0.99), na.rm = TRUE)
  metric_df$display_stat <- pmax(display_limits[1], pmin(display_limits[2], metric_df$stat))
  metric_df$positive <- pmax(metric_df$display_stat, 0)
  metric_df$negative <- pmin(metric_df$display_stat, 0)

  p_label <- if (stat_row$pval < 0.001) format(stat_row$pval, scientific = TRUE, digits = 2) else sprintf("%.3f", stat_row$pval)
  subtitle <- sprintf("ES = %.3f     NES = %.2f     P = %s     FDR = %.3f     Gene-set size = %d",
                      stat_row$ES, stat_row$NES, p_label, stat_row$padj, stat_row$size)
  line_color <- if (stat_row$NES >= 0) "#B23A48" else "#3D6D9E"

  p_es <- ggplot(curve_df, aes(rank, running_ES)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#777777") +
    geom_vline(xintercept = calc$peak, linewidth = 0.35, linetype = "dashed",
               color = "#777777") +
    geom_line(linewidth = 0.75, color = line_color) +
    annotate("point", x = calc$peak, y = calc$ES, size = 1.8, color = line_color) +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    labs(x = NULL, y = "Running ES",
         title = pretty_path(pathway_id), subtitle = subtitle) +
    theme_gsea +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank())

  p_hits <- ggplot(hit_df, aes(rank, 0)) +
    geom_segment(aes(xend = rank, y = -0.48, yend = 0.48), linewidth = 0.22,
                 color = "#161616") +
    geom_vline(xintercept = calc$peak, linewidth = 0.35, linetype = "dashed",
               color = "#777777") +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.5, 0.5)) +
    labs(x = NULL, y = "Gene hits") +
    theme_void(base_family = "Helvetica") +
    theme(axis.title.y = element_text(size = 8, angle = 90, margin = margin(r = 6)),
          plot.margin = margin(0, 7, 0, 7))

  p_metric <- ggplot(metric_df, aes(rank)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "#777777") +
    geom_area(aes(y = positive), fill = "#D96B66", alpha = 0.80) +
    geom_area(aes(y = negative), fill = "#5B83B4", alpha = 0.80) +
    geom_vline(xintercept = calc$peak, linewidth = 0.35, linetype = "dashed",
               color = "#777777") +
    scale_x_continuous(limits = c(1, n), expand = c(0, 0),
                       labels = scales::label_number(big.mark = ",")) +
    labs(x = "Rank in ordered gene list", y = "Wald statistic") +
    theme_gsea

  plot_obj <- p_es / p_hits / p_metric +
    plot_layout(heights = c(3.0, 0.65, 1.35))

  if (stat_row$NES >= 0) {
    core_idx <- calc$hit_ranks[calc$hit_ranks <= calc$peak]
  } else {
    core_idx <- calc$hit_ranks[calc$hit_ranks >= calc$peak]
  }
  core_df <- data.frame(
    pathway = pathway_id,
    gene = names(ranks)[core_idx],
    rank = core_idx,
    rank_statistic = unname(ranks[core_idx]),
    ES = stat_row$ES,
    NES = stat_row$NES,
    pvalue = stat_row$pval,
    FDR = stat_row$padj
  )
  list(plot = plot_obj, core = core_df,
    summary = data.frame(pathway = pathway_id, ES = stat_row$ES,
                            NES = stat_row$NES, pvalue = stat_row$pval,
                            FDR = stat_row$padj,
                            gene_set_size = stat_row$size,
                            peak_rank = calc$peak,
                            core_gene_count = nrow(core_df)))
}

plot_results <- lapply(selected_ids, make_gsea_plot)
names(plot_results) <- selected_ids

for (id in selected_ids) {
  out_base <- file.path(fig_dir, id)
  ggsave(paste0(out_base, ".pdf"), plot_results[[id]]$plot,
         width = 178, height = 135, units = "mm", device = cairo_pdf)
  ggsave(paste0(out_base, ".png"), plot_results[[id]]$plot,
         width = 178, height = 135, units = "mm", dpi = 600, bg = "white")
}

combined_pdf <- file.path(fig_root, "Figure4_five_traditional_GSEA_enrichment_plots.pdf")
cairo_pdf(combined_pdf, width = 178 / 25.4, height = 135 / 25.4,
          family = "Helvetica", onefile = TRUE)
for (id in selected_ids) print(plot_results[[id]]$plot)
dev.off()

summary_df <- do.call(rbind, lapply(plot_results, `[[`, "summary"))
core_df <- do.call(rbind, lapply(plot_results, `[[`, "core"))
write.csv(summary_df, file.path(table_dir, "five_GSEA_enrichment_plot_summary.csv"),
          row.names = FALSE)
write.csv(core_df, file.path(table_dir, "five_GSEA_core_enrichment_genes.csv"),
          row.names = FALSE)

cat(sprintf("Traditional GSEA plots generated: %d\n", length(plot_results)))
cat(sprintf("Combined PDF pages: %d\n", length(plot_results)))
