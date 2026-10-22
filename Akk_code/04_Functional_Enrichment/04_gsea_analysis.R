suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(BiocParallel)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

set.seed(20260701)

root <- normalizePath(file.path(dirname(getwd()), "转录组"), mustWork = FALSE)
if (!file.exists(file.path(root, "mh.all.v2025.1.Mm.symbols.gmt"))) {
  root <- normalizePath(getwd(), mustWork = TRUE)
}

gmt_file <- file.path(root, "mh.all.v2025.1.Mm.symbols.gmt")
de_dir <- file.path(root, "immune_de_nature", "tables", "de")
out_dir <- file.path(root, "immune_de_nature", "hallmark_gsea")
table_dir <- file.path(out_dir, "tables")
pdf_dir <- file.path(out_dir, "pdf")
png_dir <- file.path(out_dir, "png")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

comparisons <- c("AB_vs_AE", "AB_vs_A")
comparison_labels <- c(AB_vs_AE = "AB vs AE", AB_vs_A = "AB vs A")
comparison_cols <- c("AB vs AE" = "#6A3D9A", "AB vs A" = "#D55E00")

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pathways <- strsplit(lines, "\t")
  names(pathways) <- vapply(pathways, `[`, character(1), 1)
  lapply(pathways, function(x) unique(x[-c(1, 2)]))
}

clean_pathway_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

pathways_raw <- read_gmt(gmt_file)
canonical_gene <- setNames(
  unlist(pathways_raw, use.names = FALSE),
  toupper(unlist(pathways_raw, use.names = FALSE))
)
canonical_gene <- canonical_gene[!duplicated(names(canonical_gene))]
pathways <- lapply(pathways_raw, function(x) unname(canonical_gene[toupper(x)]))

make_rank_vector <- function(comp) {
  de <- read_csv(file.path(de_dir, paste0("DESeq2_", comp, "_all.csv")), show_col_types = FALSE) %>%
    filter(!is.na(stat), !is.na(gene)) %>%
    mutate(
      gene_upper = toupper(gene),
      canonical = unname(canonical_gene[gene_upper])
    ) %>%
    filter(!is.na(canonical)) %>%
    group_by(canonical) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranks <- de$stat
  names(ranks) <- de$canonical
  sort(ranks, decreasing = TRUE)
}

gsea_list <- lapply(seq_along(comparisons), function(i) {
  comp <- comparisons[i]
  set.seed(20260701 + i)
  ranks <- make_rank_vector(comp)
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0, BPPARAM = SerialParam()) %>%
    as_tibble() %>%
    mutate(
      comparison = comp,
      comparison_label = comparison_labels[comp],
      pathway_label = clean_pathway_name(pathway),
      leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")
    ) %>%
    arrange(padj, desc(abs(NES)))
  write_csv(res, file.path(table_dir, paste0("hallmark_fgsea_", comp, ".csv")))
  res
})

gsea_all <- bind_rows(gsea_list) %>%
  mutate(
    direction = ifelse(NES >= 0, "Enriched in numerator group", "Enriched in denominator group"),
    neg_log10_fdr = -log10(pmax(padj, .Machine$double.xmin))
  )
write_csv(gsea_all, file.path(table_dir, "hallmark_fgsea_all_comparisons.csv"))

metabolic_pathways <- c(
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_PEROXISOME",
  "HALLMARK_XENOBIOTIC_METABOLISM"
)

top_terms <- gsea_all %>%
  filter(padj < 0.25) %>%
  group_by(comparison) %>%
  arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

metabolic_terms <- gsea_all %>%
  filter(pathway %in% metabolic_pathways, padj < 0.25)

plot_terms <- bind_rows(top_terms, metabolic_terms)

if (nrow(plot_terms) == 0) {
  plot_terms <- gsea_all %>%
    group_by(comparison) %>%
    arrange(padj, desc(abs(NES)), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup()
}

plot_terms <- plot_terms %>%
  distinct(comparison_label, pathway, .keep_all = TRUE) %>%
  mutate(
    link_value = pmax(0.25, pmin(3.5, neg_log10_fdr)),
    p_label = ifelse(pval < 0.001, formatC(pval, format = "e", digits = 2), formatC(pval, format = "f", digits = 3)),
    fdr_label = ifelse(padj < 0.001, formatC(padj, format = "e", digits = 2), formatC(padj, format = "f", digits = 3))
  )

nes_limit <- max(abs(plot_terms$NES), na.rm = TRUE)
nes_limit <- ceiling(nes_limit * 10) / 10
nes_col_fun <- circlize::colorRamp2(
  c(-nes_limit, 0, nes_limit),
  c("#2F5597", "#F7F7F7", "#D55E00")
)

plot_terms <- plot_terms %>%
  mutate(link_col = circlize::add_transparency(nes_col_fun(NES), transparency = 0.12))

write_csv(plot_terms, file.path(table_dir, "hallmark_gsea_chord_plot_data.csv"))

pathway_order <- plot_terms %>%
  group_by(pathway_label) %>%
  summarise(best_padj = min(padj), max_abs_nes = max(abs(NES)), .groups = "drop") %>%
  arrange(best_padj, desc(max_abs_nes)) %>%
  pull(pathway_label)

grid_df <- plot_terms %>%
  transmute(from = comparison_label, to = pathway_label, value = link_value)

sector_order <- c(unname(comparison_labels[comparisons]), pathway_order)
grid_col <- c(
  comparison_cols,
  setNames(rep("#D9D9D9", length(pathway_order)), pathway_order)
)

truncate_text <- function(x, max_n = 27) {
  ifelse(nchar(x) > max_n, paste0(substr(x, 1, max_n - 3), "..."), x)
}

draw_pvalue_table <- function(df) {
  table_df <- df %>%
    arrange(factor(comparison_label, levels = unname(comparison_labels[comparisons])), desc(NES), padj) %>%
    mutate(
      pathway_short = truncate_text(pathway_label, 25),
      nes_label = sprintf("%.2f", NES)
    )

  x0 <- 0.745
  y_top <- 0.83
  row_h <- 0.026
  grid.text("Exact GSEA statistics", x = unit(x0, "npc"), y = unit(0.89, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 8.5, fontface = "bold"))
  grid.text("Pathway", x = unit(x0, "npc"), y = unit(y_top, "npc"), just = c("left", "top"), gp = gpar(fontsize = 6.2, fontface = "bold"))
  grid.text("NES", x = unit(0.895, "npc"), y = unit(y_top, "npc"), just = c("right", "top"), gp = gpar(fontsize = 6.2, fontface = "bold"))
  grid.text("P", x = unit(0.945, "npc"), y = unit(y_top, "npc"), just = c("right", "top"), gp = gpar(fontsize = 6.2, fontface = "bold"))
  grid.text("FDR", x = unit(0.995, "npc"), y = unit(y_top, "npc"), just = c("right", "top"), gp = gpar(fontsize = 6.2, fontface = "bold"))
  grid.lines(x = unit(c(x0, 0.995), "npc"), y = unit(c(y_top - 0.014, y_top - 0.014), "npc"), gp = gpar(col = "grey45", lwd = 0.5))

  current_y <- y_top - 0.027
  last_comp <- ""
  for (i in seq_len(nrow(table_df))) {
    row <- table_df[i, ]
    if (row$comparison_label != last_comp) {
      if (last_comp != "") current_y <- current_y - row_h * 0.30
      grid.text(row$comparison_label, x = unit(x0, "npc"), y = unit(current_y, "npc"),
                just = c("left", "top"), gp = gpar(fontsize = 6.2, fontface = "bold", col = comparison_cols[row$comparison_label]))
      current_y <- current_y - row_h * 0.78
      last_comp <- row$comparison_label
    }
    grid.text(row$pathway_short, x = unit(x0, "npc"), y = unit(current_y, "npc"),
              just = c("left", "top"), gp = gpar(fontsize = 5.5, col = "grey10"))
    grid.text(row$nes_label, x = unit(0.895, "npc"), y = unit(current_y, "npc"),
              just = c("right", "top"), gp = gpar(fontsize = 5.5, col = nes_col_fun(row$NES)))
    grid.text(row$p_label, x = unit(0.945, "npc"), y = unit(current_y, "npc"),
              just = c("right", "top"), gp = gpar(fontsize = 5.5, col = "grey10"))
    grid.text(row$fdr_label, x = unit(0.995, "npc"), y = unit(current_y, "npc"),
              just = c("right", "top"), gp = gpar(fontsize = 5.5, col = "grey10"))
    current_y <- current_y - row_h
  }
}

draw_chord <- function() {
  circos.clear()
  circos.par(
    start.degree = 88,
    canvas.xlim = c(-1.78, 2.20),
    canvas.ylim = c(-1.42, 1.42),
    gap.after = c(rep(7, length(comparisons) - 1), 18, rep(2.2, length(pathway_order) - 1), 18),
    track.margin = c(0.01, 0.01),
    points.overflow.warning = FALSE
  )

  chordDiagram(
    x = grid_df,
    order = sector_order,
    grid.col = grid_col,
    col = plot_terms$link_col,
    transparency = 0,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.16)
  )

  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      is_comp <- sector %in% unname(comparison_labels[comparisons])
      circos.text(
        x = mean(xlim),
        y = ylim[1] + 0.55,
        labels = sector,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = ifelse(is_comp, 0.82, 0.47),
        font = ifelse(is_comp, 2, 1)
      )
    }
  )

  lgd <- packLegend(
    Legend(labels = names(comparison_cols), legend_gp = gpar(fill = comparison_cols), title = "Comparison"),
    Legend(col_fun = nes_col_fun, title = "NES", at = c(-nes_limit, 0, nes_limit), labels = sprintf("%.1f", c(-nes_limit, 0, nes_limit))),
    Legend(labels = c("0.5", "1.5", "3.0"), type = "lines", legend_gp = gpar(lwd = c(1.2, 3.2, 6.0), col = "grey35"), title = expression(-log[10]~FDR)),
    direction = "vertical",
    gap = unit(3, "mm")
  )
  grid.text("Hallmark GSEA chord plot", x = unit(0.035, "npc"), y = unit(0.965, "npc"),
            just = c("left", "top"), gp = gpar(fontsize = 13, fontface = "bold"))
  grid.text("NES is mapped to a continuous link colour gradient; link width denotes -log10(FDR). Exact P and FDR values are listed at right.",
            x = unit(0.035, "npc"), y = unit(0.935, "npc"), just = c("left", "top"), gp = gpar(fontsize = 7.5, col = "grey30"))
  draw(lgd, x = unit(0.045, "npc"), y = unit(0.07, "npc"), just = c("left", "bottom"))
  draw_pvalue_table(plot_terms)
}

pdf(file.path(pdf_dir, "hallmark_gsea_AB_AE_AB_A_chord.pdf"), width = 12.4, height = 8.8, useDingbats = FALSE)
draw_chord()
dev.off()

ragg::agg_png(file.path(png_dir, "hallmark_gsea_AB_AE_AB_A_chord.png"), width = 12.4, height = 8.8, units = "in", res = 600, background = "white")
draw_chord()
dev.off()

summary_lines <- c(
  "Hallmark GSEA analysis",
  paste0("GMT: ", gmt_file),
  paste0("Comparisons: ", paste(comparisons, collapse = ", ")),
  "Ranking metric: DESeq2 Wald statistic, case-insensitive matched to mouse Hallmark symbols",
  "Chord inclusion: FDR < 0.25, up to top 8 pathways per comparison plus all significant Hallmark metabolic pathways",
  "Chord aesthetics: continuous NES colour gradient; link width = -log10(FDR); exact nominal P and FDR values printed in the right-side table",
  "",
  "Significant Hallmark counts (FDR < 0.25):",
  capture.output(print(gsea_all %>% group_by(comparison) %>% summarise(n_fdr_lt_0.25 = sum(padj < 0.25), .groups = "drop"))),
  "",
  "Output:",
  "pdf/hallmark_gsea_AB_AE_AB_A_chord.pdf",
  "png/hallmark_gsea_AB_AE_AB_A_chord.png"
)
writeLines(gsub("top 12", "top 8", summary_lines), file.path(out_dir, "hallmark_gsea_summary.txt"))

message("Done. Hallmark GSEA outputs written to: ", out_dir)
