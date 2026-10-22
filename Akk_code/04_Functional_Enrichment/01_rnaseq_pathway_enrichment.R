suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ggpubr)
  library(rstatix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(RColorBrewer)
  library(viridis)
  library(patchwork)
  library(grid)
  library(ragg)
})

set.seed(20260505)
root <- normalizePath(file.path(dirname(getwd()), "转录组"), mustWork = FALSE)
if (!file.exists(file.path(root, "group.csv"))) {
  root <- normalizePath(getwd(), mustWork = TRUE)
}
out_dir <- file.path(root, "immune_de_nature")
dirs <- list(
  de = file.path(out_dir, "tables", "de"),
  immune = file.path(out_dir, "tables", "immune"),
  pdf = file.path(out_dir, "figures", "pdf"),
  png = file.path(out_dir, "figures", "png"),
  logs = file.path(out_dir, "logs")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

save_plot <- function(plot, name, width = 7.2, height = 5.6) {
  ggsave(file.path(dirs$pdf, paste0(name, ".pdf")), plot, width = width, height = height, device = "pdf")
  ggsave(file.path(dirs$png, paste0(name, ".png")), plot, width = width, height = height, dpi = 600, device = ragg::agg_png)
}

save_pheatmap <- function(mat, name, width = 6.5, height = 7, ...) {
  ph <- pheatmap(mat, silent = TRUE, ...)
  pdf(file.path(dirs$pdf, paste0(name, ".pdf")), width = width, height = height)
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
  ragg::agg_png(file.path(dirs$png, paste0(name, ".png")), width = width, height = height, units = "in", res = 600, background = "white")
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
}

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = unit(3.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.margin = margin(5.5, 6, 5.5, 5.5)
    )
}

group_levels <- c("A", "AE", "AB")
group_cols <- c(A = "#3C5488", AE = "#00A087", AB = "#E64B35")
deg_cols <- c(Up = "#E64B35", Down = "#3C5488", NS = "#BDBDBD")
comparisons <- list(c("AE", "A"), c("AB", "A"), c("AB", "AE"))
comparison_names <- c("AE_vs_A", "AB_vs_A", "AB_vs_AE")

meta <- read_csv(file.path(root, "group.csv"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = group_levels))
stopifnot(all(!is.na(meta$group)))

counts <- read.csv(file.path(root, "AKK.reads.symbol.csv"), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "numeric"
counts <- round(counts)
counts <- counts[, meta$sample, drop = FALSE]
counts <- counts[rowSums(counts) >= 10, , drop = FALSE]
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]

coldata <- data.frame(group = meta$group, row.names = meta$sample)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ group)
dds$group <- relevel(dds$group, ref = "A")
dds <- DESeq(dds, quiet = TRUE)
norm_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(norm_counts) %>% rownames_to_column("gene"), file.path(dirs$de, "normalized_counts.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
write.csv(as.data.frame(vst_mat) %>% rownames_to_column("gene"), file.path(dirs$de, "vst_expression.csv"), row.names = FALSE)

plot_volcano <- function(res_df, comp_name) {
  fdr <- 0.05
  lfc <- 1
  res_df <- res_df %>%
    mutate(
      plot_p = pmax(pvalue, .Machine$double.xmin, na.rm = TRUE),
      neg_log10_p = -log10(plot_p),
      regulation = case_when(
        padj < fdr & log2FoldChange >= lfc ~ "Up",
        padj < fdr & log2FoldChange <= -lfc ~ "Down",
        TRUE ~ "NS"
      ),
      regulation = factor(regulation, levels = c("Up", "Down", "NS"))
    )
  label_df <- res_df %>%
    filter(regulation != "NS") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 12)
  counts_text <- res_df %>% count(regulation) %>% tidyr::complete(regulation, fill = list(n = 0))
  subtitle <- paste0("Up: ", counts_text$n[counts_text$regulation == "Up"],
                     " | Down: ", counts_text$n[counts_text$regulation == "Down"],
                     " | FDR < 0.05, |log2FC| >= 1")
  ggplot(res_df, aes(log2FoldChange, neg_log10_p)) +
    geom_point(aes(color = regulation), size = 1.2, alpha = 0.68, stroke = 0) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", linewidth = 0.28, color = "grey35") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.28, color = "grey35") +
    geom_text_repel(data = label_df, aes(label = gene), size = 2.4, max.overlaps = 30,
                    min.segment.length = 0, segment.size = 0.22, box.padding = 0.25) +
    scale_color_manual(values = deg_cols, drop = FALSE) +
    labs(title = gsub("_", " ", comp_name), subtitle = subtitle,
         x = expression(log[2]~fold~change), y = expression(-log[10]~italic(P)), color = NULL) +
    theme_nature(8)
}

all_de <- list()
for (i in seq_along(comparisons)) {
  comp <- comparisons[[i]]
  comp_name <- comparison_names[i]
  res <- results(dds, contrast = c("group", comp[1], comp[2]), alpha = 0.05)
  shrink <- tryCatch(lfcShrink(dds, contrast = c("group", comp[1], comp[2]), res = res, type = "normal"), error = function(e) res)
  res_df <- as.data.frame(shrink) %>%
    rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(
      comparison = comp_name,
      regulation = case_when(
        padj < 0.05 & log2FoldChange >= 1 ~ "Up",
        padj < 0.05 & log2FoldChange <= -1 ~ "Down",
        TRUE ~ "NS"
      )
    )
  all_de[[comp_name]] <- res_df
  write.csv(res_df, file.path(dirs$de, paste0("DESeq2_", comp_name, "_all.csv")), row.names = FALSE)
  write.csv(filter(res_df, regulation != "NS"), file.path(dirs$de, paste0("DESeq2_", comp_name, "_significant.csv")), row.names = FALSE)
  write.csv(filter(res_df, regulation == "Up"), file.path(dirs$de, paste0("DESeq2_", comp_name, "_up.csv")), row.names = FALSE)
  write.csv(filter(res_df, regulation == "Down"), file.path(dirs$de, paste0("DESeq2_", comp_name, "_down.csv")), row.names = FALSE)
  save_plot(plot_volcano(res_df, comp_name), paste0("volcano_", comp_name), width = 6.2, height = 5.2)

  top_genes <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 50) %>%
    pull(gene)
  if (length(top_genes) >= 2) {
    hm <- vst_mat[top_genes, meta$sample, drop = FALSE]
    hm <- t(scale(t(hm)))
    hm[is.na(hm)] <- 0
    ann_col <- data.frame(Group = meta$group, row.names = meta$sample)
    ann_colors <- list(Group = group_cols)
    save_pheatmap(hm,
            paste0("heatmap_top50_", comp_name),
            annotation_col = ann_col,
            annotation_colors = ann_colors,
            color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
            breaks = seq(-2.5, 2.5, length.out = 102),
            cluster_cols = FALSE,
            cluster_rows = TRUE,
            show_rownames = TRUE,
            fontsize = 6.5,
            border_color = NA,
            width = 6.4,
            height = 7.2)
  }
}

all_de_df <- bind_rows(all_de)
write.csv(all_de_df, file.path(dirs$de, "DESeq2_all_comparisons.csv"), row.names = FALSE)

deg_summary <- all_de_df %>%
  filter(regulation != "NS") %>%
  count(comparison, regulation) %>%
  tidyr::complete(comparison = comparison_names, regulation = c("Up", "Down"), fill = list(n = 0))
write.csv(deg_summary, file.path(dirs$de, "DEG_count_summary.csv"), row.names = FALSE)

p_deg_bar <- ggplot(deg_summary, aes(comparison, n, fill = regulation)) +
  geom_col(width = 0.66, position = position_dodge(width = 0.72), color = "black", linewidth = 0.25) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72), vjust = -0.35, size = 2.6) +
  scale_fill_manual(values = deg_cols[c("Up", "Down")]) +
  labs(x = NULL, y = "Number of DEGs", fill = NULL, title = "Differentially expressed genes") +
  theme_nature(8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p_deg_bar, "bar_DEG_counts", width = 4.8, height = 4.1)

immune <- read.csv(file.path(root, "immucellai_mouse_result.csv"), check.names = FALSE)
immune <- immune %>% left_join(meta, by = "sample") %>% mutate(group = factor(group, levels = group_levels))
stopifnot(all(!is.na(immune$group)))
write.csv(immune, file.path(dirs$immune, "immucellai_mouse_with_group.csv"), row.names = FALSE)

immune_value_cols <- setdiff(colnames(immune), c("sample", "group"))
immune_long <- immune %>%
  pivot_longer(all_of(immune_value_cols), names_to = "cell_type", values_to = "score")
write.csv(immune_long, file.path(dirs$immune, "immucellai_mouse_long.csv"), row.names = FALSE)

immune_group_summary <- immune_long %>%
  group_by(cell_type, group) %>%
  summarise(mean = mean(score, na.rm = TRUE), sd = sd(score, na.rm = TRUE), sem = sd / sqrt(n()), n = n(), .groups = "drop")
write.csv(immune_group_summary, file.path(dirs$immune, "immucellai_group_summary.csv"), row.names = FALSE)

kw_stats <- immune_long %>%
  group_by(cell_type) %>%
  kruskal_test(score ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  arrange(p.adj)
write.csv(kw_stats, file.path(dirs$immune, "immune_kruskal_BH.csv"), row.names = FALSE)

p_to_symbol <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p <= 0.0001 ~ "****",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

pairwise_stats <- bind_rows(lapply(unique(immune_long$cell_type), function(cell) {
  cell_df <- filter(immune_long, cell_type == cell)
  bind_rows(lapply(comparisons, function(comp) {
    x <- cell_df$score[cell_df$group == comp[1]]
    y <- cell_df$score[cell_df$group == comp[2]]
    wt <- tryCatch(wilcox.test(x, y, exact = FALSE), error = function(e) NULL)
    data.frame(
      cell_type = cell,
      group1 = comp[1],
      group2 = comp[2],
      n1 = length(x),
      n2 = length(y),
      statistic = if (is.null(wt)) NA_real_ else unname(wt$statistic),
      p = if (is.null(wt)) NA_real_ else wt$p.value,
      stringsAsFactors = FALSE
    )
  }))
})) %>%
  group_by(cell_type) %>%
  mutate(p.adj = p.adjust(p, method = "BH"), p.adj = ifelse(is.na(p.adj), 1, p.adj), p.adj.signif = p_to_symbol(p.adj)) %>%
  ungroup() %>%
  arrange(cell_type, p.adj)
write.csv(pairwise_stats, file.path(dirs$immune, "immune_pairwise_wilcox_BH.csv"), row.names = FALSE)

top_cells <- kw_stats %>% slice_head(n = min(12, nrow(kw_stats))) %>% pull(cell_type)
immune_top <- immune_long %>% filter(cell_type %in% top_cells) %>% mutate(cell_type = factor(cell_type, levels = top_cells))
stat_top <- pairwise_stats %>%
  filter(cell_type %in% top_cells) %>%
  group_by(cell_type) %>%
  mutate(y.position = max(immune_top$score[immune_top$cell_type == unique(cell_type)], na.rm = TRUE) * (1.12 + 0.12 * row_number())) %>%
  ungroup()

p_immune_bar <- ggplot(immune_top, aes(group, score, fill = group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.66, color = "black", linewidth = 0.28, alpha = 0.72) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.28) +
  geom_point(aes(color = group), position = position_jitter(width = 0.08, height = 0), size = 1.65, alpha = 0.95) +
  stat_pvalue_manual(stat_top, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, size = 2.4, bracket.size = 0.22) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = NULL, y = "ImmuCellAI score", title = "Immune infiltration differences") +
  theme_nature(8) +
  theme(legend.position = "none")
save_plot(p_immune_bar, "bar_scatter_top_immune_cells_with_stats", width = 8.8, height = 7.0)

score_df <- immune_long %>% filter(cell_type == "Infiltration_score")
score_stats <- pairwise_stats %>%
  filter(cell_type == "Infiltration_score") %>%
  mutate(y.position = max(score_df$score, na.rm = TRUE) * (1.08 + 0.12 * row_number()))
p_score <- ggplot(score_df, aes(group, score, fill = group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6, color = "black", linewidth = 0.28, alpha = 0.72) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.16, linewidth = 0.28) +
  geom_point(aes(color = group), position = position_jitter(width = 0.08), size = 2.1) +
  stat_pvalue_manual(score_stats, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE, size = 3, bracket.size = 0.25) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  labs(x = NULL, y = "Infiltration score", title = "Global immune infiltration") +
  theme_nature(8) +
  theme(legend.position = "none")
save_plot(p_score, "bar_scatter_infiltration_score_with_stats", width = 3.4, height = 4.2)

bubble_df <- pairwise_stats %>%
  left_join(
    immune_group_summary %>% select(cell_type, group, mean) %>% pivot_wider(names_from = group, values_from = mean),
    by = "cell_type"
  )
bubble_df$mean_diff <- mapply(
  function(i, g1, g2) bubble_df[[g1]][i] - bubble_df[[g2]][i],
  seq_len(nrow(bubble_df)), bubble_df$group1, bubble_df$group2
)
bubble_df <- bubble_df %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    neg_log10_fdr = -log10(pmax(p.adj, .Machine$double.xmin)),
    cell_type = factor(cell_type, levels = rev(kw_stats$cell_type))
  )
write.csv(bubble_df, file.path(dirs$immune, "immune_pairwise_effect_bubble_data.csv"), row.names = FALSE)

p_bubble <- ggplot(bubble_df, aes(comparison, cell_type)) +
  geom_point(aes(size = neg_log10_fdr, fill = mean_diff), shape = 21, color = "black", stroke = 0.18, alpha = 0.9) +
  geom_text(aes(label = p.adj.signif), size = 2.2, color = "black") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, name = "Mean difference") +
  scale_size_continuous(range = c(1.2, 5.2), name = expression(-log[10]~FDR)) +
  labs(x = NULL, y = NULL, title = "Pairwise immune infiltration statistics") +
  theme_nature(7) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_bubble, "bubble_immune_pairwise_stats", width = 5.6, height = 8.8)

immune_mat <- immune %>% select(sample, all_of(immune_value_cols)) %>% column_to_rownames("sample") %>% as.matrix()
immune_mat_scaled <- t(scale(t(immune_mat)))
immune_mat_scaled[is.na(immune_mat_scaled)] <- 0
ann_col <- data.frame(Group = immune$group, row.names = immune$sample)
save_pheatmap(immune_mat_scaled,
              "heatmap_immune_infiltration_scores",
              annotation_col = ann_col,
              annotation_colors = list(Group = group_cols),
              color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(101),
              breaks = seq(-2.5, 2.5, length.out = 102),
              cluster_cols = FALSE,
              cluster_rows = TRUE,
              border_color = NA,
              fontsize = 6.8,
              width = 6.8,
              height = 7.4)

summary_lines <- c(
  paste0("Analysis date: ", Sys.Date()),
  paste0("Samples: ", nrow(meta)),
  paste0("Groups: ", paste(levels(meta$group), collapse = ", ")),
  paste0("Genes tested after filtering: ", nrow(counts)),
  paste0("DE comparisons: ", paste(comparison_names, collapse = ", ")),
  paste0("Immune cell features tested: ", length(immune_value_cols)),
  "DEG threshold: FDR < 0.05 and |log2FC| >= 1",
  "Immune statistics: Kruskal-Wallis across groups and pairwise Wilcoxon tests, BH adjusted"
)
writeLines(summary_lines, file.path(dirs$logs, "analysis_summary.txt"))

sessionInfo_text <- capture.output(sessionInfo())
writeLines(sessionInfo_text, file.path(dirs$logs, "sessionInfo.txt"))
cat("DONE grouped immune infiltration and differential expression analysis\n")
