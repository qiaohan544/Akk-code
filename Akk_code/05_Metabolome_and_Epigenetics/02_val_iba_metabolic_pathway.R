suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

meta <- read_csv("../meta3.csv", show_col_types = FALSE)
kegg <- read_tsv("../kegg.all_unstratified.CPM.tsv", show_col_types = FALSE)
pathway <- read_tsv("../pathways.all_unstratified.CPM.tsv", show_col_types = FALSE)

meta_sub <- meta %>%
  filter(class == 2, time == 3, group %in% c("PD1+A+B", "PD1+ABX+Akk")) %>%
  mutate(group = factor(group, levels = c("PD1+A+B", "PD1+ABX+Akk"))) %>%
  arrange(group, ids)

sample_ids <- meta_sub$ids

feature_def <- tibble::tribble(
  ~feature_id, ~feature_label, ~feature_class, ~feature_type,
  "K00835", "K00835 valine-pyruvate aminotransferase", "Valine/Butyrate genes", "KO",
  "K00836", "K00836 diaminobutyrate transaminase", "Valine/Butyrate genes", "KO",
  "K07250", "K07250 4-aminobutyrate aminotransferase", "Valine/Butyrate genes", "KO",
  "K00019", "K00019 3-hydroxybutyrate dehydrogenase", "Valine/Butyrate genes", "KO",
  "K00020", "K00020 3-hydroxyisobutyrate dehydrogenase", "Valine/Butyrate genes", "KO",
  "K00929", "K00929 butyrate kinase", "Valine/Butyrate genes", "KO",
  "K08318", "K08318 4-hydroxybutyrate dehydrogenase", "Valine/Butyrate genes", "KO",
  "K00140", "K00140 methylmalonate-semialdehyde dehydrogenase", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00826", "K00826 branched-chain amino acid aminotransferase", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00167", "K00167 2-oxoisovalerate dehydrogenase E1 beta", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K11381", "K11381 2-oxoisovalerate dehydrogenase E1", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00186", "K00186 2-oxoisovalerate ferredoxin oxidoreductase alpha", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00187", "K00187 2-oxoisovalerate ferredoxin oxidoreductase beta", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00189", "K00189 2-oxoisovalerate ferredoxin oxidoreductase gamma", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00248", "K00248 butyryl-CoA dehydrogenase", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K00634", "K00634 phosphate butyryltransferase", "Valine degradation / isobutyrate synthesis genes", "KO",
  "K01186", "K01186 sialidase-1", "Mucin degradation genes", "KO",
  "K01190", "K01190 beta-galactosidase", "Mucin degradation genes", "KO",
  "K01205", "K01205 alpha-N-acetylglucosaminidase", "Mucin degradation genes", "KO",
  "K01206", "K01206 alpha-L-fucosidase", "Mucin degradation genes", "KO",
  "K01207", "K01207 beta-N-acetylhexosaminidase", "Mucin degradation genes", "KO",
  "K01443", "K01443 N-acetylglucosamine-6-phosphate deacetylase", "Mucin degradation genes", "KO",
  "K12373", "K12373 hexosaminidase", "Mucin degradation genes", "KO",
  "PWY0-42", "PWY0-42 2-methylcitrate cycle I", "Valine degradation / isobutyrate pathways", "Pathway",
  "PWY-5747", "PWY-5747 2-methylcitrate cycle II", "Valine degradation / isobutyrate pathways", "Pathway",
  "FUCCAT-PWY", "FUCCAT-PWY fucose degradation", "Mucin degradation pathways", "Pathway",
  "FUC-RHAMCAT-PWY", "FUC-RHAMCAT-PWY fucose/rhamnose degradation", "Mucin degradation pathways", "Pathway",
  "GLCMANNANAUT-PWY", "GLCMANNANAUT-PWY GlcNAc/ManNAc/Neu5Ac degradation", "Mucin degradation pathways", "Pathway",
  "P441-PWY", "P441-PWY N-acetylneuraminate degradation", "Mucin degradation pathways", "Pathway",
  "PWY-6317", "PWY-6317 D-galactose degradation I", "Mucin degradation pathways", "Pathway"
)

kegg2 <- kegg %>%
  mutate(feature_id = sub(":.*", "", GeneFamily))

pathway2 <- pathway %>%
  mutate(feature_id = sub(":.*", "", Family))

extract_long <- function(df, id_col, ids, feature_ids, value_name = "abundance") {
  out <- df %>%
    filter(feature_id %in% feature_ids) %>%
    select(all_of(c(id_col, "feature_id", ids)))
  melt(out, id.vars = c(id_col, "feature_id"), variable.name = "ids", value.name = value_name)
}

ko_long <- extract_long(kegg2, "GeneFamily", sample_ids, feature_def$feature_id[feature_def$feature_type == "KO"])
pw_long <- extract_long(pathway2, "Family", sample_ids, feature_def$feature_id[feature_def$feature_type == "Pathway"])

all_long <- bind_rows(
  ko_long %>% transmute(feature_id, feature_name = GeneFamily, ids, abundance),
  pw_long %>% transmute(feature_id, feature_name = Family, ids, abundance)
) %>%
  left_join(feature_def, by = "feature_id") %>%
  left_join(meta_sub, by = "ids")

target_feature_classes <- c(
  "Valine/Butyrate genes",
  "Valine degradation / isobutyrate pathways"
)

target_gene_classes <- c(
  "Valine/Butyrate genes"
)

valine_butyrate_long <- all_long %>%
  filter(feature_class %in% target_feature_classes)

valine_butyrate_gene_long <- all_long %>%
  filter(feature_class %in% target_gene_classes, feature_type == "KO")

analyze_feature <- function(df) {
  if (n_distinct(df$group) < 2) {
    return(NULL)
  }
  non_missing_abundance <- df$abundance[!is.na(df$abundance)]
  if (length(non_missing_abundance) == 0 || dplyr::n_distinct(non_missing_abundance) < 2) {
    return(NULL)
  }
  kw <- kruskal.test(abundance ~ group, data = df)
  pair <- tryCatch(
    compare_means(abundance ~ group, data = df, method = "wilcox.test", p.adjust.method = "BH"),
    error = function(e) NULL
  )
  list(kw = kw, pair = pair)
}

feature_stats <- lapply(split(all_long, all_long$feature_id), analyze_feature)
feature_stats <- feature_stats[!vapply(feature_stats, is.null, logical(1))]

kw_table <- bind_rows(lapply(names(feature_stats), function(fid) {
  tibble(
    feature_id = fid,
    KW_p = feature_stats[[fid]]$kw$p.value
  )
})) %>%
  left_join(feature_def, by = "feature_id") %>%
  arrange(feature_class, KW_p)

pairwise_table <- bind_rows(lapply(names(feature_stats), function(fid) {
  pair <- feature_stats[[fid]]$pair
  if (is.null(pair)) {
    return(NULL)
  }
  pair %>% mutate(feature_id = fid)
})) %>%
  left_join(feature_def, by = "feature_id") %>%
  left_join(kw_table %>% select(feature_id, KW_p), by = "feature_id") %>%
  arrange(feature_class, p.adj, p)

write_csv(kw_table, "class2_time3_valine_butyrate_mucin_kw.csv")
write_csv(pairwise_table, "class2_time3_valine_butyrate_mucin_pairwise.csv")

valine_butyrate_gene_kw <- kw_table %>%
  filter(feature_class %in% target_gene_classes, feature_type == "KO") %>%
  arrange(KW_p, feature_class, feature_id)

valine_butyrate_gene_pairwise <- pairwise_table %>%
  filter(feature_class %in% target_gene_classes, feature_type == "KO") %>%
  arrange(feature_id, p.adj, p)

write_csv(valine_butyrate_gene_kw, "class2_time3_valine_butyrate_genes_kw.csv")
write_csv(valine_butyrate_gene_pairwise, "class2_time3_valine_butyrate_genes_pairwise.csv")

sig_features <- kw_table %>%
  filter(KW_p < 0.05) %>%
  pull(feature_id)

if (length(sig_features) == 0) {
  sig_features <- kw_table %>% slice_head(n = min(8, n())) %>% pull(feature_id)
}

sig_long <- all_long %>%
  filter(feature_id %in% sig_features) %>%
  mutate(
    feature_label = factor(feature_label, levels = kw_table %>% filter(feature_id %in% sig_features) %>% pull(feature_label))
  )

boxplot_obj <- ggplot(sig_long, aes(x = group, y = abundance, fill = group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9, color = "#333333", linewidth = 0.25) +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.85, color = "#1F1F1F") +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 2.8) +
  scale_fill_manual(values = c("PD1+A+B" = "#2F6C8F", "PD1+ABX+Akk" = "#3F8C68")) +
  labs(x = NULL, y = "Relative abundance / CPM", title = "Class 2 time3: significant functional features") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F4F1EA", color = "#D0C8BA"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12)
  )

gene_plot_long <- valine_butyrate_gene_long %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = valine_butyrate_gene_kw$feature_label
    )
  )

gene_kw_labels <- valine_butyrate_gene_kw %>%
  mutate(
    feature_label = factor(feature_label, levels = valine_butyrate_gene_kw$feature_label),
    label = paste0("KW p=", ifelse(KW_p < 1e-4, formatC(KW_p, format = "e", digits = 2), signif(KW_p, 3))),
    x = 2,
    y = Inf
  )

gene_bar_scatter_obj <- ggplot(gene_plot_long, aes(x = group, y = abundance, fill = group)) +
  stat_summary(fun = mean, geom = "col", width = 0.62, color = "#333333", linewidth = 0.22, alpha = 0.88) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.25, color = "#333333") +
  geom_jitter(width = 0.10, size = 1.65, alpha = 0.82, color = "#1F1F1F") +
  geom_text(
    data = gene_kw_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 2.55
  ) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("PD1+A+B" = "#2F6C8F", "PD1+ABX+Akk" = "#3F8C68")) +
  labs(
    x = NULL,
    y = "Relative abundance / CPM",
    title = "Class 2 time3: valine and butyrate genes across treatment groups",
    subtitle = "Bars show mean +/- SE; points show individual samples"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F4F1EA", color = "#D0C8BA", linewidth = 0.25),
    strip.text = element_text(face = "bold", size = 7.8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, size = 8.2),
    axis.text.y = element_text(size = 8.2),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, color = "#4D4D4D")
  )

gene_panel_n <- n_distinct(gene_plot_long$feature_label)
gene_panel_ncol <- 3
gene_panel_nrow <- ceiling(gene_panel_n / gene_panel_ncol)

ggsave(
  "class2_time3_valine_butyrate_genes_bar_scatter_nature.pdf",
  gene_bar_scatter_obj,
  width = 11,
  height = 11,
  useDingbats = FALSE
)

heatmap_df <- valine_butyrate_long %>%
  select(feature_id, feature_label, feature_class, ids, abundance) %>%
  distinct()

heatmap_mat <- heatmap_df %>%
  select(feature_label, ids, abundance) %>%
  pivot_wider(names_from = ids, values_from = abundance) %>%
  as.data.frame()

rownames(heatmap_mat) <- heatmap_mat$feature_label
heatmap_mat$feature_label <- NULL
heatmap_mat <- as.matrix(heatmap_mat[, sample_ids, drop = FALSE])

heatmap_log <- log10(heatmap_mat + 1e-4)
heatmap_z <- t(scale(t(heatmap_log)))
heatmap_z[is.na(heatmap_z)] <- 0
heatmap_z[heatmap_z > 2.5] <- 2.5
heatmap_z[heatmap_z < -2.5] <- -2.5

row_info <- feature_def %>%
  filter(feature_label %in% rownames(heatmap_z)) %>%
  mutate(feature_label = factor(feature_label, levels = rownames(heatmap_z))) %>%
  arrange(feature_label)

row_split <- factor(
  row_info$feature_class,
  levels = c(
    "Valine/Butyrate genes",
    "Valine degradation / isobutyrate pathways"
  )
)
names(row_split) <- as.character(row_info$feature_label)

ha_top <- HeatmapAnnotation(
  Group = meta_sub$group,
  col = list(Group = c("PD1+A+B" = "#2F6C8F", "PD1+ABX+Akk" = "#3F8C68")),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  simple_anno_size = unit(4, "mm")
)

row_anno <- rowAnnotation(
  Class = row_split[rownames(heatmap_z)],
  col = list(Class = c(
    "Valine/Butyrate genes" = "#6C8EBF",
    "Valine degradation / isobutyrate pathways" = "#7FA35B"
  )),
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)

heatmap_cell_size <- unit(4.2, "mm")

valine_module_stats <- kw_table %>%
  filter(feature_class %in% c("Valine degradation / isobutyrate synthesis genes", "Valine degradation / isobutyrate pathways"))

valine_module_pairwise <- pairwise_table %>%
  filter(feature_class %in% c("Valine degradation / isobutyrate synthesis genes", "Valine degradation / isobutyrate pathways"))

write_csv(valine_module_stats, "class2_time3_valine_degradation_kw.csv")
write_csv(valine_module_pairwise, "class2_time3_valine_degradation_pairwise.csv")

valine_sig_features <- valine_module_stats %>%
  filter(KW_p < 0.05) %>%
  pull(feature_id)

if (length(valine_sig_features) == 0) {
  valine_sig_features <- valine_module_stats %>% slice_head(n = min(6, n())) %>% pull(feature_id)
}

if (length(valine_sig_features) > 1) {
  valine_sig_features <- head(valine_sig_features, -1)
}

valine_sig_long <- all_long %>%
  filter(feature_id %in% valine_sig_features) %>%
  mutate(
    feature_label = factor(
      feature_label,
      levels = valine_module_stats %>% filter(feature_id %in% valine_sig_features) %>% pull(feature_label)
    )
  )

valine_boxplot_obj <- ggplot(valine_sig_long, aes(x = group, y = abundance, fill = group)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.9, color = "#333333", linewidth = 0.25) +
  geom_jitter(width = 0.12, size = 1.7, alpha = 0.85, color = "#1F1F1F") +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 2.8) +
  scale_fill_manual(values = c("PD1+A+B" = "#2F6C8F", "PD1+ABX+Akk" = "#3F8C68")) +
  labs(x = NULL, y = "Relative abundance / CPM", title = "Class 2 time3: valine degradation-related features") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#F4F1EA", color = "#D0C8BA"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave(
  "class2_time3_valine_degradation_boxplots.pdf",
  valine_boxplot_obj,
  width = 11,
  height = max(5, ceiling(length(valine_sig_features) / 3) * 2.2),
  useDingbats = FALSE
)

pdf("class2_time3_valine_butyrate_mucin_heatmap_nature.pdf", width = 8.8, height = 5.2, useDingbats = FALSE)
ht <- Heatmap(
  heatmap_z,
  name = "Row z-score",
  top_annotation = ha_top,
  left_annotation = row_anno,
  row_split = row_split[rownames(heatmap_z)],
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.8),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 9.2, fontface = "bold"),
  column_title = "Class 2 time3 functional landscape: valine and butyrate genes/pathways",
  column_title_gp = gpar(fontsize = 10.8, fontface = "bold"),
  width = ncol(heatmap_z) * heatmap_cell_size,
  height = nrow(heatmap_z) * heatmap_cell_size,
  border = TRUE,
  col = colorRamp2(c(-2.5, 0, 2.5), c("#1F5A91", "#F7F3EC", "#A9342A")),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"), labels_gp = gpar(fontsize = 8))
)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(6, 8, 4, 8), "mm"))
grid.text(
  "Rows scaled independently after log10 transformation; samples ordered by treatment group.",
  x = unit(0.5, "npc"), y = unit(1.8, "mm"), gp = gpar(fontsize = 8, col = "#4D4D4D")
)
dev.off()

write_csv(
  cbind(feature = rownames(heatmap_z), as.data.frame(heatmap_z, check.names = FALSE)),
  "class2_time3_valine_butyrate_mucin_heatmap_matrix_zscore.csv"
)

message("Done:")
message(" - class2_time3_valine_butyrate_mucin_kw.csv")
message(" - class2_time3_valine_butyrate_mucin_pairwise.csv")
message(" - class2_time3_valine_butyrate_mucin_heatmap_nature.pdf")
message(" - class2_time3_valine_butyrate_mucin_heatmap_matrix_zscore.csv")
message(" - class2_time3_valine_butyrate_genes_kw.csv")
message(" - class2_time3_valine_butyrate_genes_pairwise.csv")
message(" - class2_time3_valine_butyrate_genes_bar_scatter_nature.pdf")
message(" - class2_time3_valine_degradation_kw.csv")
message(" - class2_time3_valine_degradation_pairwise.csv")
message(" - class2_time3_valine_degradation_boxplots.pdf")
