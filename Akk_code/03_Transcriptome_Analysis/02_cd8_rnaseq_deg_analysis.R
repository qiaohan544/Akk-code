#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(grid)
})

base_dir <- ""./data/"异丁酸处理CD8转录组"
out_dir <- file.path(base_dir, "rnaseq_results")
count_dir <- file.path(out_dir, "counts")
table_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
gmt_file <- file.path(base_dir, "m5.go.bp.v2026.1.Mm.symbols.gmt")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("NC1", "NC4", "NC5", "I2", "I3", "I4")
group <- factor(c("NC", "NC", "NC", "Isobutyrate", "Isobutyrate", "Isobutyrate"),
                levels = c("NC", "Isobutyrate"))

read_fc <- function(strand) {
  files <- file.path(count_dir, sprintf("%s.s%d.gene_counts.txt", samples, strand))
  xs <- lapply(files, read.delim, check.names = FALSE, comment.char = "#")
  count_cols <- vapply(xs, function(x) names(x)[ncol(x)], character(1))
  mat <- do.call(cbind, Map(function(x, nm) x[[nm]], xs, count_cols))
  storage.mode(mat) <- "integer"
  rownames(mat) <- xs[[1]]$Geneid
  colnames(mat) <- samples
  list(counts = mat, annotation = xs[[1]][, 1:6])
}

assigned_fraction <- function(strand) {
  files <- file.path(count_dir, sprintf("%s.s%d.gene_counts.txt.summary", samples, strand))
  ans <- vapply(files, function(f) {
    x <- read.delim(f, check.names = FALSE, row.names = 1)
    x["Assigned", 1] / sum(x[, 1])
  }, numeric(1))
  setNames(ans, samples)
}

fractions <- rbind(s0 = assigned_fraction(0), s1 = assigned_fraction(1),
                   s2 = assigned_fraction(2))
mean_fraction <- rowMeans(fractions)
best_stranded <- if (mean_fraction["s1"] >= mean_fraction["s2"]) 1L else 2L
other_stranded <- if (best_stranded == 1L) 2L else 1L
use_strand <- if (mean_fraction[paste0("s", best_stranded)] >= 0.75 * mean_fraction["s0"] &&
                  mean_fraction[paste0("s", best_stranded)] >= 4 * mean_fraction[paste0("s", other_stranded)]) {
  best_stranded
} else 0L
write.csv(cbind(strand_mode = rownames(fractions), fractions,
                mean_assigned_fraction = mean_fraction),
          file.path(table_dir, "featureCounts_strand_assessment.csv"), row.names = FALSE)
writeLines(sprintf("Selected featureCounts strand mode: %d", use_strand),
           file.path(table_dir, "selected_strand_mode.txt"))

fc <- read_fc(use_strand)
counts <- fc$counts
gene_id_clean <- sub("\\..*$", "", rownames(counts))
symbol <- mapIds(org.Mm.eg.db, keys = gene_id_clean, keytype = "ENSEMBL",
                 column = "SYMBOL", multiVals = "first")
symbol <- unname(symbol)

coldata <- data.frame(row.names = samples, group = group)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds, quiet = TRUE)
res_raw <- results(dds, contrast = c("group", "Isobutyrate", "NC"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, contrast = c("group", "Isobutyrate", "NC"),
                        res = res_raw, type = "normal", quiet = TRUE)

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
res_df$ensembl_id <- sub("\\..*$", "", res_df$gene_id)
res_df$symbol <- symbol[match(res_df$ensembl_id, gene_id_clean)]
raw_df <- as.data.frame(res_raw)
res_df$stat <- raw_df[res_df$gene_id, "stat"]
res_df <- res_df[, c("gene_id", "ensembl_id", "symbol", "baseMean", "log2FoldChange",
                     "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
write.csv(res_df, file.path(table_dir, "DESeq2_Isobutyrate_vs_NC_all_genes.csv"), row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
norm_out <- data.frame(gene_id = rownames(norm_counts), symbol = res_df$symbol[match(rownames(norm_counts), res_df$gene_id)],
                       norm_counts, check.names = FALSE)
write.csv(norm_out, file.path(table_dir, "DESeq2_normalized_counts.csv"), row.names = FALSE)

# Preranked GSEA: DESeq2 Wald statistic, one representative per mouse gene symbol.
rank_df <- data.frame(symbol = res_df$symbol, stat = res_df$stat)
rank_df <- rank_df[!is.na(rank_df$symbol) & nzchar(rank_df$symbol) & is.finite(rank_df$stat), ]
rank_df <- rank_df[order(abs(rank_df$stat), decreasing = TRUE), ]
rank_df <- rank_df[!duplicated(rank_df$symbol), ]
ranks <- rank_df$stat
names(ranks) <- rank_df$symbol
ranks <- sort(ranks, decreasing = TRUE)

pathways <- gmtPathways(gmt_file)

# GO BP over-representation analysis (ORA) for the volcano-plot companion
# panels. Up- and down-regulated genes are tested separately against all genes
# entering the DESeq2 test. The supplied GO BP GMT is used as the ontology.
universe_symbols <- unique(res_df$symbol[!is.na(res_df$symbol) & nzchar(res_df$symbol)])
up_symbols <- unique(res_df$symbol[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1 &
                                   !is.na(res_df$symbol) & nzchar(res_df$symbol)])
down_symbols <- unique(res_df$symbol[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1 &
                                     !is.na(res_df$symbol) & nzchar(res_df$symbol)])

go_ora <- function(query, direction) {
  query <- intersect(unique(query), universe_symbols)
  ans <- lapply(names(pathways), function(term) {
    genes <- intersect(unique(pathways[[term]]), universe_symbols)
    m <- length(genes)
    k <- length(intersect(query, genes))
    if (m < 10 || m > 500 || k < 2) return(NULL)
    p <- phyper(k - 1, m, length(universe_symbols) - m, length(query), lower.tail = FALSE)
    data.frame(pathway = term, direction = direction, overlap = k,
               query_size = length(query), pathway_size = m,
               gene_ratio = k / length(query),
               fold_enrichment = (k / length(query)) / (m / length(universe_symbols)),
               pvalue = p,
               genes = paste(sort(intersect(query, genes)), collapse = ";"),
               stringsAsFactors = FALSE)
  })
  ans <- do.call(rbind, ans)
  if (is.null(ans) || nrow(ans) == 0) return(data.frame())
  ans$padj <- p.adjust(ans$pvalue, method = "BH")
  ans[order(ans$padj, ans$pvalue, -ans$fold_enrichment), ]
}

go_up <- go_ora(up_symbols, "Up in isobutyrate")
go_down <- go_ora(down_symbols, "Down in isobutyrate")
go_all <- rbind(go_up, go_down)
write.csv(go_up, file.path(table_dir, "GO_BP_ORA_upregulated_genes.csv"), row.names = FALSE)
write.csv(go_down, file.path(table_dir, "GO_BP_ORA_downregulated_genes.csv"), row.names = FALSE)
write.csv(go_all, file.path(table_dir, "GO_BP_ORA_all_results.csv"), row.names = FALSE)

set.seed(20260822)
fg <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15,
                      maxSize = 500, eps = 0,
                      BPPARAM = BiocParallel::SerialParam())
fg <- as.data.frame(fg)
fg <- fg[order(fg$padj, -abs(fg$NES)), ]
fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
write.csv(fg, file.path(table_dir, "GSEA_GO_BP_all_pathways.csv"), row.names = FALSE)

theme_nature <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(axis.title = element_text(color = "black", size = 9),
        axis.text = element_text(color = "black", size = 8),
        axis.line = element_line(linewidth = 0.45),
        axis.ticks = element_line(linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 10, hjust = 0),
        legend.title = element_blank(), legend.text = element_text(size = 8))

pretty_path <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("t Cell", "T Cell", x, fixed = TRUE)
  x <- gsub("Mrna", "mRNA", x, fixed = TRUE)
  x <- gsub("Dna", "DNA", x, fixed = TRUE)
  x <- gsub("Stat Protein", "STAT Protein", x, fixed = TRUE)
  vapply(x, function(s) paste(strwrap(s, width = 34), collapse = "\n"),
         FUN.VALUE = character(1), USE.NAMES = FALSE)
}

# Volcano plot.
vol <- res_df[is.finite(res_df$log2FoldChange) & is.finite(res_df$padj), ]
vol$status <- "Not significant"
vol$status[vol$padj < 0.05 & vol$log2FoldChange >= 1] <- "Up"
vol$status[vol$padj < 0.05 & vol$log2FoldChange <= -1] <- "Down"
vol$status <- factor(vol$status, levels = c("Not significant", "Down", "Up"))
vol$mlog10_raw <- -log10(pmax(vol$padj, 1e-300))
vol$mlog10 <- pmin(vol$mlog10_raw, 100)
cd8_label_genes <- c(
  "Ccl3", "Ccl4", "Ccl5", "Gzmb", "Prf1", "Crtam", "Xcl1", "Tnf",
  "Il2rb", "Il2ra", "Tcf7", "Eomes", "Tbx21", "Bcl6", "Ctla4", "Tigit",
  "Pdcd1", "Tox", "Cxcr3", "Ccr5", "Nkg7", "Cd69"
)
label_df <- vol[match(tolower(cd8_label_genes), tolower(vol$symbol)), , drop = FALSE]
label_df <- label_df[!is.na(label_df$symbol) & is.finite(label_df$padj), ]
write.csv(label_df[, c("symbol", "baseMean", "log2FoldChange", "stat", "pvalue", "padj")],
          file.path(table_dir, "CD8T_molecules_labeled_in_volcano.csv"), row.names = FALSE)

p_volcano <- ggplot(vol, aes(log2FoldChange, mlog10, color = status)) +
  geom_point(size = 0.7, alpha = 0.62, stroke = 0) +
  geom_point(data = label_df, aes(log2FoldChange, mlog10), inherit.aes = FALSE,
             shape = 21, size = 1.8, stroke = 0.35, color = "#202020", fill = "#F2C14E") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
  geom_text_repel(data = label_df, aes(label = symbol), size = 2.6,
                  color = "black", min.segment.length = 0, box.padding = 0.25,
                  point.padding = 0.22, max.overlaps = Inf, force = 1.5,
                  max.time = 3, seed = 20260822) +
  scale_color_manual(values = c("Not significant" = "#B8B8B8", "Down" = "#4C78A8", "Up" = "#D64F4F")) +
  scale_y_continuous(limits = c(0, 115), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = expression(log[2]~fold~change), y = expression(-log[10]~FDR~"(capped at 100)"),
       title = "Differential expression: isobutyrate versus control") +
  theme_nature + theme(legend.position = "top")

# GSEA bars styled after the supplied reference. UP and DOWN show eight ranked
# representative pathways; NES values are printed directly on the bars.
fg_sig <- fg[is.finite(fg$padj) & fg$padj < 0.05, ]
requested_gsea_up <- c(
  "GOBP_MRNA_PROCESSING",
  "GOBP_REGULATION_OF_TRANSLATION",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_RESPONSE_TO_TUMOR_CELL",
  "GOBP_PROTEIN_FOLDING",
  "GOBP_NUCLEAR_TRANSPORT",
  "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY",
  "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ORGANELLE"
)
requested_gsea_down <- c(
  "GOBP_POSITIVE_REGULATION_OF_TYROSINE_PHOSPHORYLATION_OF_STAT_PROTEIN",
  "GOBP_RETICULOPHAGY",
  "GOBP_GLYCERALDEHYDE_3_PHOSPHATE_METABOLIC_PROCESS",
  "GOBP_MITOCHONDRIAL_DNA_METABOLIC_PROCESS",
  "GOBP_NEGATIVE_REGULATION_OF_CYCLIN_DEPENDENT_PROTEIN_SERINE_THREONINE_KINASE_ACTIVITY",
  "GOBP_ALDEHYDE_BIOSYNTHETIC_PROCESS",
  "GOBP_GLYCOLYTIC_PROCESS_THROUGH_GLUCOSE_6_PHOSPHATE",
  "GOBP_INOSITOL_PHOSPHATE_METABOLIC_PROCESS"
)
select_requested_gsea <- function(ids) {
  idx <- match(ids, fg$pathway)
  if (anyNA(idx)) stop("Requested GSEA pathways missing: ", paste(ids[is.na(idx)], collapse = ", "))
  fg[idx, , drop = FALSE]
}
g_up <- select_requested_gsea(requested_gsea_up)
g_down <- select_requested_gsea(requested_gsea_down)
write.csv(rbind(transform(g_up, display_group = "UP"),
                transform(g_down, display_group = "Down")),
          file.path(table_dir, "GSEA_pathways_displayed_in_Figure1.csv"), row.names = FALSE)
make_gsea_bar <- function(z, fill, title) {
  if (nrow(z) == 0) {
    return(ggplot() + annotate("text", 0, 0, label = "No pathways at FDR < 0.05", size = 3) +
             xlim(-1, 1) + ylim(-1, 1) + theme_void(base_family = "Helvetica") + ggtitle(title) +
             theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)))
  }
  z$label <- pretty_path(z$pathway)
  z$score <- -log10(pmax(z$padj, 1e-300))
  z$label <- factor(z$label, levels = rev(z$label))
  label_size <- if (nrow(z) > 6) 6.0 else 7.2
  nes_size <- if (nrow(z) > 6) 2.05 else 2.35
  ggplot(z, aes(score, label)) +
    geom_col(width = 0.78, fill = fill) +
    geom_text(aes(label = sprintf("NES %.2f", NES)), hjust = 1.06,
              size = nes_size, fontface = "bold", color = "#222222") +
    labs(x = expression(-log[10]~GSEA~FDR), y = NULL, title = title) +
    theme_nature +
    theme(axis.text.y = element_text(size = label_size), axis.ticks.y = element_blank())
}
p_gsea_up <- make_gsea_bar(g_up, "#DFA7AF", "Enriched in isobutyrate")
p_gsea_down <- make_gsea_bar(g_down, "#AABCCB", "Enriched in control")
p_combo_gsea <- (p_gsea_up | p_gsea_down) / p_volcano +
  plot_layout(heights = c(1.20, 1.65)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(file.path(fig_dir, "Figure1_volcano_GSEA.pdf"), p_combo_gsea,
       width = 183, height = 232, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, "Figure1_volcano_GSEA.png"), p_combo_gsea,
       width = 183, height = 232, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "volcano_plot.pdf"), p_volcano,
       width = 100, height = 100, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, "GSEA_GO_BP_barplot.pdf"), p_gsea_up | p_gsea_down,
       width = 183, height = 96, units = "mm", device = cairo_pdf)

# Main figure requested by the user: GO BP ORA bars above the volcano plot.
make_go_bar <- function(z, fill, title, n = 5) {
  z <- z[is.finite(z$pvalue) & z$pvalue < 0.05, , drop = FALSE]
  z <- head(z[order(z$pvalue, z$padj, -z$fold_enrichment), , drop = FALSE], n)
  if (nrow(z) == 0) {
    return(ggplot() + annotate("text", 0, 0, label = "No GO terms at P < 0.05", size = 3) +
             xlim(-1, 1) + ylim(-1, 1) + theme_void(base_family = "Helvetica") + ggtitle(title) +
             theme(plot.title = element_text(face = "bold", size = 10, hjust = 0)))
  }
  z$label <- factor(pretty_path(z$pathway), levels = rev(pretty_path(z$pathway)))
  z$score <- -log10(pmax(z$pvalue, 1e-300))
  ggplot(z, aes(score, label)) +
    geom_col(width = 0.84, fill = fill) +
    labs(x = expression(-log[10]~italic(P)), y = NULL, title = title) +
    theme_nature +
    theme(axis.text.y = element_text(size = 7.2), axis.ticks.y = element_blank())
}

p_go_up <- make_go_bar(go_up, "#DFA7AF", "Upregulated genes")
p_go_down <- make_go_bar(go_down, "#AABCCB", "Downregulated genes")
p_combo_go <- (p_go_up | p_go_down) / p_volcano +
  plot_layout(heights = c(0.72, 1.65)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(file.path(fig_dir, "Supplementary_volcano_GO_enrichment.pdf"), p_combo_go,
       width = 183, height = 205, units = "mm", device = cairo_pdf)
ggsave(file.path(fig_dir, "Supplementary_volcano_GO_enrichment.png"), p_combo_go,
       width = 183, height = 205, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "GO_BP_ORA_barplot.pdf"), p_go_up | p_go_down,
       width = 183, height = 72, units = "mm", device = cairo_pdf)

# T-cell pathway heatmap: significant pathways only; genes are selected from leading edges.
tcell_pattern <- paste(c("(^|_)T_CELL(_|$)", "T_LYMPHOCYTE", "CD8_POSITIVE.*T_CELL",
                         "CD4_POSITIVE.*T_CELL", "REGULATORY_T_CELL", "MEMORY_T_CELL",
                         "NK_T_CELL", "THYMOCYTE", "LYMPHOCYTE_MEDIATED_IMMUNITY",
                         "LEUKOCYTE_MEDIATED_IMMUNITY", "RESPONSE_TO_TUMOR_CELL",
                         "CELL_KILLING", "CYTOTOXIC"), collapse = "|")
tfg <- fg_sig[grepl(tcell_pattern, fg_sig$pathway, perl = TRUE), ]
tfg <- tfg[order(tfg$padj, -abs(tfg$NES)), ]
write.csv(tfg, file.path(table_dir, "GSEA_Tcell_pathways_FDR_lt_0.05.csv"), row.names = FALSE)

make_empty_heatmap <- function() {
  cairo_pdf(file.path(fig_dir, "Figure2_Tcell_GSEA_heatmap.pdf"), width = 7.2, height = 4.5)
  grid.newpage()
  grid.text("No T-cell-related GO BP pathway reached GSEA FDR < 0.05",
            gp = gpar(fontfamily = "Helvetica", fontsize = 12))
  dev.off()
}

if (nrow(tfg) == 0) {
  make_empty_heatmap()
} else {
  selected_paths <- head(tfg, 8)
  leading_lists <- strsplit(selected_paths$leadingEdge, ";", fixed = TRUE)
  names(leading_lists) <- selected_paths$pathway
  gene_frequency <- sort(table(unlist(leading_lists)), decreasing = TRUE)
  candidate_genes <- names(gene_frequency)

  vsd <- vst(dds, blind = FALSE)
  vmat <- assay(vsd)
  vsym <- mapIds(org.Mm.eg.db, keys = sub("\\..*$", "", rownames(vmat)),
                 keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  keep <- !is.na(vsym) & !duplicated(vsym)
  vmat <- vmat[keep, , drop = FALSE]
  rownames(vmat) <- unname(vsym[keep])
  candidate_genes <- intersect(candidate_genes, rownames(vmat))

  stat_by_symbol <- ranks[candidate_genes]
  priority <- as.numeric(gene_frequency[candidate_genes]) * 10 + abs(stat_by_symbol)
  candidate_genes <- candidate_genes[order(priority, decreasing = TRUE)]
  heat_genes <- head(candidate_genes, 60)
  hmat <- vmat[heat_genes, samples, drop = FALSE]
  hmat <- t(scale(t(hmat)))
  hmat[!is.finite(hmat)] <- 0

  first_path <- vapply(heat_genes, function(g) {
    hit <- names(leading_lists)[vapply(leading_lists, function(x) g %in% x, logical(1))]
    if (length(hit)) hit[1] else "Other"
  }, character(1))
  split_labels <- pretty_path(first_path)

  membership <- do.call(rbind, lapply(names(leading_lists), function(p) {
    data.frame(pathway = p, gene = intersect(leading_lists[[p]], heat_genes))
  }))
  write.csv(membership, file.path(table_dir, "Tcell_heatmap_gene_pathway_membership.csv"), row.names = FALSE)

  ha <- HeatmapAnnotation(
    Group = group,
    col = list(Group = c(NC = "#4C78A8", Isobutyrate = "#D64F4F")),
    annotation_name_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
    annotation_legend_param = list(title_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
                                   labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8)))
  ht <- Heatmap(hmat, name = "Row z-score",
                col = colorRamp2(c(-2, 0, 2), c("#2C5A85", "#F7F7F7", "#B2182B")),
                top_annotation = ha, column_split = group, cluster_columns = TRUE,
                cluster_column_slices = FALSE, cluster_rows = TRUE,
                row_split = split_labels, show_row_dend = FALSE,
                row_names_gp = gpar(fontfamily = "Helvetica", fontsize = 6.5),
                column_names_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
                row_title_gp = gpar(fontfamily = "Helvetica", fontsize = 7),
                column_title_gp = gpar(fontfamily = "Helvetica", fontsize = 9, fontface = "bold"),
                heatmap_legend_param = list(title_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
                                            labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8)),
                column_title = "Leading-edge genes from significant T-cell-related GO BP pathways")
  heat_height <- max(6.5, min(11, 2.5 + 0.12 * nrow(hmat)))
  cairo_pdf(file.path(fig_dir, "Figure2_Tcell_GSEA_heatmap.pdf"), width = 8.3, height = heat_height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(fig_dir, "Figure2_Tcell_GSEA_heatmap.png"), width = 8.3, height = heat_height,
      units = "in", res = 600, type = "cairo", bg = "white")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# PCA provides a compact diagnostic alongside the requested main figures.
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group, label = name)) +
  geom_point(size = 3) + geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(NC = "#4C78A8", Isobutyrate = "#D64F4F")) +
  labs(x = sprintf("PC1 (%d%%)", percent_var[1]), y = sprintf("PC2 (%d%%)", percent_var[2]),
       title = "RNA-seq sample structure") + theme_nature + theme(legend.position = "top")
ggsave(file.path(fig_dir, "QC_PCA.pdf"), p_pca, width = 90, height = 80,
       units = "mm", device = cairo_pdf)

summary_lines <- c(
  sprintf("Selected featureCounts strand mode: %d", use_strand),
  sprintf("Genes tested: %d", nrow(res_df)),
  sprintf("DE genes (FDR < 0.05, |log2FC| >= 1): %d", sum(vol$status != "Not significant")),
  sprintf("  Up in isobutyrate: %d", sum(vol$status == "Up")),
  sprintf("  Down in isobutyrate: %d", sum(vol$status == "Down")),
  sprintf("Significant GO BP ORA terms for upregulated genes: %d", sum(go_up$padj < 0.05)),
  sprintf("Significant GO BP ORA terms for downregulated genes: %d", sum(go_down$padj < 0.05)),
  sprintf("Significant GO BP GSEA pathways (FDR < 0.05): %d", nrow(fg_sig)),
  sprintf("Significant T-cell-related pathways: %d", nrow(tfg))
)
writeLines(summary_lines, file.path(out_dir, "analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
