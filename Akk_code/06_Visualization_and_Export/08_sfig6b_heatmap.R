suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(ragg)
})

table_file <- "SFig6B plot data.xlsx"
out_dir <- "."

theme_heatmap_nature <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.28, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(size = base_size, angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = base_size - 0.2),
      axis.title = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.size = unit(3.6, "mm"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.3, colour = "grey30", hjust = 0),
      plot.margin = margin(5.5, 6, 5.5, 5.5),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

format_cell_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("\\bCD([0-9]) ", "CD\\1 ", x)
  x
}

# Read data
bubble_df <- read_excel(table_file)

# Rename columns to match the code
colnames(bubble_df) <- c("cell_type", "mean_diff", "comparison", "neg_log10_fdr", "sig_label")

bubble_df <- bubble_df %>%
  filter(!grepl("Infiltration_score|Infiltration score", cell_type, ignore.case = TRUE)) %>%
  mutate(
    mean_diff = as.numeric(mean_diff),
    neg_log10_fdr = as.numeric(neg_log10_fdr),
    comparison = factor(comparison, levels = c("A. m.+P. v. vs A. m.", "A. m.+P. v. vs A. m.+E. f.", "A. m.+E. f. vs A. m.")),
    sig_label = ifelse(is.na(sig_label), "", sig_label),
    is_significant = sig_label != ""
  ) %>%
  # 按每个细胞类型（行）进行 Z-score 标准化
  group_by(cell_type) %>%
  mutate(
    z_score = if(sd(mean_diff, na.rm = TRUE) == 0) 0 else (mean_diff - mean(mean_diff, na.rm = TRUE)) / sd(mean_diff, na.rm = TRUE)
  ) %>%
  ungroup()

cell_order <- bubble_df %>%
  group_by(cell_type) %>%
  summarise(
    target_z_score = sum(z_score[comparison == "A. m.+P. v. vs A. m."], na.rm = TRUE),
    any_significant = any(is_significant, na.rm = TRUE),
    max_sig_strength = max(neg_log10_fdr, na.rm = TRUE),
    max_abs_effect = max(abs(mean_diff), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 按照“A. m.+P. v. vs A. m.”组的 Z score 降序排列，取前 15 名
  arrange(desc(target_z_score)) %>%
  slice_head(n = 15) %>%
  # 再按照之前的显著性逻辑进行排序以便在图上呈现
  arrange(desc(any_significant), desc(max_sig_strength), desc(max_abs_effect), cell_type) %>%
  pull(cell_type)

# 只保留前15名的细胞类型进行绘图
bubble_df <- bubble_df %>%
  filter(cell_type %in% cell_order)

effect_limit <- max(abs(bubble_df$z_score), na.rm = TRUE)
effect_limit <- ceiling(effect_limit * 10) / 10

plot_df <- bubble_df %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(cell_order)),
    cell_label = factor(format_cell_label(as.character(cell_type)), levels = format_cell_label(rev(cell_order)))
  )

p <- ggplot(plot_df, aes(comparison, cell_label)) +
  geom_tile(aes(fill = z_score), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 2.8, colour = "black", fontface = "bold", vjust = 0.75) +
  scale_fill_gradient2(
    low = "#3C5488", mid = "white", high = "#E64B35",
    midpoint = 0, limits = c(-effect_limit, effect_limit), oob = scales::squish,
    name = "Z-score\n(Row scaled)"
  ) +
  guides(
    fill = guide_colorbar(order = 1, barheight = unit(24, "mm"), barwidth = unit(4, "mm"))
  ) +
  labs(
    title = "Pairwise immune infiltration statistics",
    subtitle = "Colour indicates row-scaled Z-score of mean difference, asterisks indicate significance"
  ) +
  coord_cartesian(clip = "off") +
  theme_heatmap_nature(7.5) +
  theme(legend.position = "right")

# Export source data
library(writexl)
write_xlsx(plot_df, file.path(out_dir, "Source_Data_heatmap.xlsx"))

ggsave(file.path(out_dir, "heatmap_immune_pairwise_stats.pdf"), p, width = 4.5, height = 8.9, device = "pdf", useDingbats = FALSE)
ggsave(file.path(out_dir, "heatmap_immune_pairwise_stats.png"), p, width = 4.5, height = 8.9, dpi = 600, device = ragg::agg_png, background = "white")

message("Heatmap plotted successfully!")
