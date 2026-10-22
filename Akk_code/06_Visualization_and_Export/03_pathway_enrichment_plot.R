suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(RColorBrewer)
  library(ragg)
  library(fgsea)
  library(stringr)
})

set.seed(20260505)
root <- normalizePath(file.path(dirname(getwd()), "转录组"), mustWork = FALSE)
if (!file.exists(file.path(root, "group.csv"))) root <- normalizePath(getwd(), mustWork = TRUE)
out_root <- file.path(root, "immune_de_nature")
out_dir <- file.path(out_root, "revised_ordered_figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, name, width, height) {
  ggsave(file.path(out_dir, "pdf", paste0(name, ".pdf")), plot, width = width, height = height, device = "pdf")
  ggsave(file.path(out_dir, "png", paste0(name, ".png")), plot, width = width, height = height, dpi = 600, device = ragg::agg_png)
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
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      plot.margin = margin(5, 6, 5, 5)
    )
}

group_levels <- c("A", "AE", "AB")
group_cols <- c(A = "#3C5488", AE = "#00A087", AB = "#E64B35")
deg_cols <- c(Up = "#E64B35", Down = "#3C5488", NS = "#BDBDBD")
comparisons <- c("AE_vs_A", "AB_vs_A", "AB_vs_AE")

meta <- read_csv(file.path(root, "group.csv"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = group_levels))

plot_volcano <- function(comp_name) {
  df <- read.csv(file.path(out_root, "tables", "de", paste0("DESeq2_", comp_name, "_all.csv")), check.names = FALSE)
  p_cut <- 0.05
  lfc_cut <- 0.5
  p_floor <- 1e-50
  fdr_floor <- 1e-50
  y_cap <- 25
  df <- df %>%
    mutate(
      pvalue_plot = ifelse(is.na(pvalue), 1, pmax(pvalue, p_floor)),
      neg_log10_p_raw = -log10(pvalue_plot),
      neg_log10_p = pmin(neg_log10_p_raw, y_cap),
      padj_plot = ifelse(is.na(padj), 1, pmax(padj, fdr_floor)),
      neg_log10_fdr_raw = -log10(padj_plot),
      capped = neg_log10_p_raw > y_cap,
      regulation = case_when(
        pvalue < p_cut & log2FoldChange >= lfc_cut ~ "Up",
        pvalue < p_cut & log2FoldChange <= -lfc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      regulation = factor(regulation, levels = c("Up", "Down", "NS"))
    )
  label_df <- df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    filter(regulation != "NS" | row_number() <= 10) %>%
    slice_head(n = 14)
  count_df <- df %>% count(regulation) %>% complete(regulation = c("Up", "Down", "NS"), fill = list(n = 0))
  subtitle <- paste0("Threshold: P < 0.05, |log2FC| >= 0.5; ",
                     "Up = ", count_df$n[count_df$regulation == "Up"],
                     ", Down = ", count_df$n[count_df$regulation == "Down"],
                     "; -log10(P) capped at ", y_cap)
  ggplot(df, aes(log2FoldChange, neg_log10_p)) +
    geom_point(aes(color = regulation), size = 1.05, alpha = 0.68, stroke = 0) +
    geom_point(data = filter(df, capped), shape = 24, size = 1.35, fill = "black", color = "black", alpha = 0.75) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.28, color = "grey35") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.28, color = "grey35") +
    geom_text_repel(data = label_df, aes(label = gene), size = 2.35, max.overlaps = 30,
                    min.segment.length = 0, segment.size = 0.22, box.padding = 0.25) +
    scale_color_manual(values = deg_cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = gsub("_", " ", comp_name), subtitle = subtitle,
         x = expression(log[2]~fold~change), y = expression(-log[10]~italic(P)), color = NULL) +
    theme_nature(8)
}

for (i in seq_along(comparisons)) {
  p <- plot_volcano(comparisons[i])
  save_plot(p, sprintf("%02d_volcano_%s_log2FC0.5_FDRcapped", i, comparisons[i]), 6.2, 5.2)
}

# Pairwise immune infiltration statistics ternary bubble plot
pw <- read.csv(file.path(out_root, "tables", "immune", "immune_pairwise_wilcox_BH.csv"), check.names = FALSE)
summary_immune <- read.csv(file.path(out_root, "tables", "immune", "immucellai_group_summary.csv"), check.names = FALSE)
wide_mean <- summary_immune %>% select(cell_type, group, mean) %>% pivot_wider(names_from = group, values_from = mean)
pw2 <- pw %>% left_join(wide_mean, by = "cell_type")
pw2$mean_diff <- mapply(function(i, g1, g2) pw2[[g1]][i] - pw2[[g2]][i], seq_len(nrow(pw2)), pw2$group1, pw2$group2)
pw2 <- pw2 %>%
  mutate(
    comparison = paste(group1, "vs", group2),
    signed_score = sign(mean_diff) * -log10(pmax(p.adj, .Machine$double.xmin)),
    p_label = case_when(p.adj <= 0.001 ~ "***", p.adj <= 0.01 ~ "**", p.adj <= 0.05 ~ "*", TRUE ~ ""),
    cell_type = factor(cell_type, levels = rev(unique(cell_type[order(tapply(p.adj, cell_type, min, na.rm = TRUE)[cell_type])])))
  )
write.csv(pw2, file.path(out_dir, "tables", "04_pairwise_immune_heatmap_data.csv"), row.names = FALSE)

ternary_df <- wide_mean %>%
  mutate(
    total = A + AE + AB,
    A_prop = ifelse(total > 0, A / total, 1 / 3),
    AE_prop = ifelse(total > 0, AE / total, 1 / 3),
    AB_prop = ifelse(total > 0, AB / total, 1 / 3),
    x = AB_prop + 0.5 * AE_prop,
    y = (sqrt(3) / 2) * AE_prop,
    dominant_group = factor(c("A", "AE", "AB")[max.col(cbind(A_prop, AE_prop, AB_prop), ties.method = "first")], levels = group_levels)
  ) %>%
  left_join(
    pw2 %>%
      group_by(cell_type) %>%
      summarise(
        min_pairwise_p = ifelse(all(is.na(p)), 1, min(p, na.rm = TRUE)),
        max_sig_score = ifelse(all(is.na(p)), 0, max(-log10(pmax(p, .Machine$double.xmin)), na.rm = TRUE)),
        sig_label = case_when(
          min_pairwise_p <= 0.001 ~ "***",
          min_pairwise_p <= 0.01 ~ "**",
          min_pairwise_p <= 0.05 ~ "*",
          TRUE ~ ""
        ),
        .groups = "drop"
      ),
    by = "cell_type"
  )
write.csv(ternary_df, file.path(out_dir, "tables", "04_immune_infiltration_ternary_bubble_data.csv"), row.names = FALSE)

triangle <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, sqrt(3) / 2, 0))
grid_lines <- bind_rows(lapply(seq(0.2, 0.8, by = 0.2), function(v) {
  data.frame(
    x = c(v, 1 - 0.5 * v, 0.5 * (1 - v), v, NA),
    y = c(0, (sqrt(3) / 2) * v, (sqrt(3) / 2) * (1 - v), 0, NA),
    line = v
  )
}))

p_immune_ternary <- ggplot() +
  geom_path(data = grid_lines, aes(x, y, group = line), color = "grey88", linewidth = 0.22) +
  geom_path(data = triangle, aes(x, y), linewidth = 0.45, color = "black") +
  geom_point(data = ternary_df, aes(x, y, size = max_sig_score, fill = dominant_group),
             shape = 21, color = "black", stroke = 0.22, alpha = 0.88) +
  ggrepel::geom_text_repel(data = filter(ternary_df, sig_label != ""), aes(x, y, label = cell_type),
                           size = 2.35, min.segment.length = 0, segment.size = 0.2,
                           box.padding = 0.24, max.overlaps = 50, color = "black") +
  geom_text(data = filter(ternary_df, sig_label != ""), aes(x, y, label = sig_label),
            size = 2.5, vjust = -0.7, color = "black") +
  geom_text(data = data.frame(x = c(-0.035, 1.035, 0.5), y = c(-0.045, -0.045, sqrt(3) / 2 + 0.045), label = c("A", "AB", "AE")),
            aes(x, y, label = label), fontface = "bold", size = 3) +
  scale_fill_manual(values = group_cols, name = "Highest mean") +
  scale_size_continuous(range = c(1.6, 6), name = expression(max~-log[10]~italic(P))) +
  coord_equal(xlim = c(-0.08, 1.08), ylim = c(-0.08, sqrt(3) / 2 + 0.09), clip = "off") +
  labs(title = "Immune infiltration ternary bubble plot", subtitle = "Position: relative mean abundance across A/AE/AB; stars: significant pairwise p-value", x = NULL, y = NULL) +
  theme_void(base_size = 8) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 7.5, hjust = 0),
    plot.margin = margin(5, 8, 5, 5)
  )
save_plot(p_immune_ternary, "04_immune_infiltration_ternary_bubble_plot", 5.8, 5.2)

# CD8 and antigen-presentation functional gene sets for AB vs AE GSEA and heatmap
cd8_sets <- list(
  CD8_CYTOTOXICITY = c("Cd8a", "Cd8b1", "Gzma", "Gzmb", "Gzmk", "Gzmm", "Prf1", "Nkg7", "Ctsw", "Cst7", "Ifng", "Fasl", "Ccl5", "Tbx21", "Eomes", "Lamp1", "Klrg1"),
  ANTIGEN_PRESENTATION_MHC_I = c("B2m", "H2-K1", "H2-D1", "H2-L", "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Psmb10", "Nlrc5", "Pdia3", "Calr", "Canx", "Erp57"),
  ANTIGEN_PRESENTATION_MHC_II = c("H2-Aa", "H2-Ab1", "H2-Eb1", "H2-DMa", "H2-DMb1", "Cd74", "Ciita", "Ctss", "Ifi30"),
  T_CELL_RECEPTOR_SIGNALING = c("Cd3d", "Cd3e", "Cd3g", "Cd247", "Lck", "Zap70", "Lat", "Lcp2", "Itk", "Themis", "Cd28", "Icos", "Il2ra"),
  CD8_EXHAUSTION_CHECKPOINT = c("Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tigit", "Entpd1", "Tox", "Tox2", "Nr4a1", "Nr4a2"),
  INTERFERON_GAMMA_RESPONSE = c("Ifng", "Ifngr1", "Ifngr2", "Stat1", "Irf1", "Cxcl9", "Cxcl10", "Cxcl11", "Gbp2", "Gbp5", "Isg15", "Ifit1", "Ifit3")
)

ab_ae <- read.csv(file.path(out_root, "tables", "de", "DESeq2_AB_vs_AE_all.csv"), check.names = FALSE)
rank_df <- ab_ae %>%
  filter(!is.na(stat), !is.na(gene)) %>%
  group_by(gene) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup()
ranks <- rank_df$stat
names(ranks) <- rank_df$gene
ranks <- sort(ranks, decreasing = TRUE)

rank_gene_map <- setNames(names(ranks), toupper(names(ranks)))
map_to_rank_genes <- function(genes) {
  mapped <- unname(rank_gene_map[toupper(genes)])
  unique(mapped[!is.na(mapped)])
}
present_sets <- lapply(cd8_sets, map_to_rank_genes)
present_sets <- present_sets[lengths(present_sets) >= 5]
fg <- fgsea(pathways = present_sets, stats = ranks, minSize = 5, maxSize = 500, eps = 0) %>%
  as.data.frame() %>%
  arrange(padj, desc(abs(NES))) %>%
  mutate(
    pathway_label = str_replace_all(pathway, "_", " "),
    pathway_label = factor(pathway_label, levels = rev(str_replace_all(pathway, "_", " "))),
    neg_log10_fdr = -log10(pmax(padj, .Machine$double.xmin)),
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  )
write.csv(fg, file.path(out_dir, "tables", "05_GSEA_AB_vs_AE_CD8_antigen_presentation_fgsea.csv"), row.names = FALSE)

p_gsea <- ggplot(fg, aes(NES, pathway_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey40") +
  geom_point(aes(size = size, color = neg_log10_fdr), alpha = 0.92) +
  scale_color_viridis_c(option = "magma", direction = 1, name = expression(-log[10]~FDR)) +
  scale_size_continuous(range = c(2.3, 6), name = "Genes") +
  labs(title = "AB vs AE GSEA: CD8 and antigen-presentation functions",
       subtitle = "Color scale shows -log10(FDR)",
       x = "Normalized enrichment score", y = NULL) +
  theme_nature(8)
save_plot(p_gsea, "05_GSEA_AB_vs_AE_CD8_antigen_presentation_dotplot", 6.2, 4.3)

# CD8 functional gene-set heatmap
vst <- read.csv(file.path(out_root, "tables", "de", "vst_expression.csv"), check.names = FALSE)
vst_mat <- as.matrix(vst[, meta$sample])
rownames(vst_mat) <- vst$gene
vst_gene_map <- setNames(rownames(vst_mat), toupper(rownames(vst_mat)))
map_to_vst_genes <- function(genes) {
  mapped <- unname(vst_gene_map[toupper(genes)])
  unique(mapped[!is.na(mapped)])
}
gene_set_map <- bind_rows(lapply(names(cd8_sets), function(nm) {
  data.frame(gene = map_to_vst_genes(cd8_sets[[nm]]), signature = nm)
})) %>%
  filter(!is.na(gene), nzchar(gene)) %>%
  distinct(gene, .keep_all = TRUE)
mat <- vst_mat[gene_set_map$gene, meta$sample, drop = FALSE]
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0
hm_df <- as.data.frame(mat_z) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "z") %>%
  left_join(meta, by = "sample") %>%
  left_join(gene_set_map, by = "gene") %>%
  mutate(
    sample = factor(sample, levels = meta$sample[order(meta$group)]),
    gene = factor(gene, levels = rev(gene_set_map$gene)),
    signature = factor(signature, levels = names(cd8_sets))
  )
write.csv(hm_df, file.path(out_dir, "tables", "06_CD8_functional_gene_set_heatmap_data.csv"), row.names = FALSE)

p_cd8_heat <- ggplot(hm_df, aes(sample, gene, fill = pmax(pmin(z, 2.5), -2.5))) +
  geom_tile(color = "white", linewidth = 0.12) +
  facet_grid(signature ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 0, name = "Z-score") +
  labs(title = "CD8-related functional gene-set expression", x = NULL, y = NULL) +
  theme_nature(6.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y = element_text(angle = 0, hjust = 0), panel.spacing.y = unit(1, "mm"))
save_plot(p_cd8_heat, "06_CD8_related_functional_gene_set_heatmap", 7.2, 8.8)

# Ordered file manifest
manifest <- data.frame(
  order = c(1, 2, 3, 4, 5, 6),
  figure = c(
    "01_volcano_AE_vs_A_log2FC0.5_FDRcapped",
    "02_volcano_AB_vs_A_log2FC0.5_FDRcapped",
    "03_volcano_AB_vs_AE_log2FC0.5_FDRcapped",
    "04_immune_infiltration_ternary_bubble_plot",
    "05_GSEA_AB_vs_AE_CD8_antigen_presentation_dotplot",
    "06_CD8_related_functional_gene_set_heatmap"
  )
)
write.csv(manifest, file.path(out_dir, "tables", "figure_order_manifest.csv"), row.names = FALSE)

writeLines(c(
  "Revised figure set generated for requested order.",
  "Volcano threshold: P < 0.05 and |log2FC| >= 0.5; -log10(P) capped at 25 for display.",
  "Ternary bubble plot significance and bubble size are based on pairwise p-values.",
  "GSEA dotplot (Figure 05) uses -log10(FDR) in color legend.",
  "GSEA package note: clusterProfiler/enrichplot unavailable because ggtree/ape failed against system CXXABI; fgsea was used for statistics with a clusterProfiler-style dotplot.",
  "Primary GSEA comparison: AB vs AE."
), file.path(out_dir, "README_revised_figures.txt"))

cat("DONE revised ordered figures\n")
