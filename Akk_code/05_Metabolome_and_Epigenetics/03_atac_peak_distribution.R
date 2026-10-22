#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(
    paste(
      "Usage: make_atac_appendix_figures.R",
      "<annotated_report_tsv> <counted_dba_rds> <sample_sheet_tsv>",
      "<gene_set_xlsx> <outdir> [method=DESeq2]"
    )
  )
}

annotated_report_path <- args[[1]]
counted_dba_path <- args[[2]]
sample_sheet_path <- args[[3]]
gene_set_path <- args[[4]]
outdir <- args[[5]]
method <- if (length(args) >= 6) args[[6]] else "DESeq2"

userlib <- normalizePath("r_libs", mustWork = FALSE)
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(DiffBind)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(readxl)
  library(tidyr)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  out <- df
  list_cols <- vapply(out, is.list, logical(1))
  if (any(list_cols)) {
    for (col_name in names(out)[list_cols]) {
      out[[col_name]] <- vapply(
        out[[col_name]],
        function(x) {
          if (length(x) == 0 || all(is.na(x))) {
            return("")
          }
          paste(as.character(x), collapse = ";")
        },
        character(1)
      )
    }
  }
  write.table(
    out,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

clean_gene_tokens <- function(text) {
  if (is.na(text) || text == "") {
    return(character(0))
  }
  tokens <- trimws(unlist(strsplit(as.character(text), ";", fixed = TRUE)))
  tokens <- tokens[tokens != ""]
  tokens <- tokens[!grepl("^ENSG[0-9]+$", tokens)]
  unique(tokens)
}

build_peak_id <- function(df) {
  paste(df$chrom, df$start, df$end, sep = "_")
}

row_scale_matrix <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled
}

extract_ratio_numeric <- function(ratio_text) {
  nums <- strsplit(ratio_text, "/", fixed = TRUE)[[1]]
  as.numeric(nums[[1]]) / as.numeric(nums[[2]])
}

normalize_enrich_df <- function(enrich_obj) {
  out <- as.data.frame(enrich_obj)
  if (nrow(out) == 0) {
    out <- data.frame(
      ID = character(0),
      Description = character(0),
      GeneRatio = character(0),
      BgRatio = character(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      qvalue = numeric(0),
      geneID = character(0),
      Count = integer(0)
    )
  } else {
    out$GeneRatioNumeric <- vapply(
      out$GeneRatio,
      extract_ratio_numeric,
      numeric(1)
    )
  }
  out
}

plot_bar_terms <- function(
  df,
  title_text,
  fill_low,
  fill_high,
  out_prefix,
  fill_column = "p.adjust",
  fill_label = "Adjusted p"
) {
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }

  plot_height <- min(12, max(4.5, nrow(df) * 0.35))
  plot_df <- df %>%
    arrange(.data[[fill_column]], desc(Count)) %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(unique(Description))
      )
    )

  p <- ggplot(
    plot_df,
    aes(x = Count, y = Description, fill = .data[[fill_column]])
  ) +
    geom_col() +
    scale_fill_gradient(low = fill_low, high = fill_high) +
    theme_bw(base_size = 12) +
    labs(
      title = title_text,
      x = "Gene count",
      y = NULL,
      fill = fill_label
    )

  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = 10,
    height = plot_height
  )
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = 10,
    height = plot_height,
    dpi = 300
  )
}

plot_empty_message <- function(title_text, message_text, out_prefix) {
  p <- ggplot() +
    annotate("text", x = 1, y = 1, label = message_text, size = 5) +
    xlim(0, 2) +
    ylim(0, 2) +
    theme_void(base_size = 12) +
    labs(title = title_text)

  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = 9,
    height = 4.5
  )
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = 9,
    height = 4.5,
    dpi = 300
  )
}

get_method_columns <- function(method_name) {
  if (method_name == "DESeq2") {
    return(list(p = "DESeq2_p.value", fc = "DESeq2_Fold", fdr = "DESeq2_FDR"))
  }
  if (method_name == "edgeR") {
    return(list(p = "edgeR_p.value", fc = "edgeR_Fold", fdr = "edgeR_FDR"))
  }
  stop("method must be DESeq2 or edgeR")
}

annotated_df <- read.delim(
  annotated_report_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
sample_df <- read.delim(
  sample_sheet_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_df$Condition <- factor(
  sample_df$Condition,
  levels = c("Isobutyrate", "NC")
)
sample_df <- sample_df %>%
  arrange(Condition, Replicate, SampleID)

sample_order <- sample_df$SampleID
sample_condition_df <- sample_df %>%
  dplyr::select(SampleID, Condition) %>%
  tibble::column_to_rownames("SampleID")

method_cols <- get_method_columns(method)
p_col <- method_cols$p
fc_col <- method_cols$fc
fdr_col <- method_cols$fdr

annotated_df$peak_key <- build_peak_id(annotated_df)
annotated_df$promoter_gene_list <- lapply(
  annotated_df$promoter_gene_names,
  clean_gene_tokens
)

annotation_summary <- annotated_df %>%
  filter(region_annotation %in% c("Promoter", "Exon", "Intron", "Intergenic")) %>%
  count(region_annotation, name = "peak_count") %>%
  mutate(
    fraction = peak_count / sum(peak_count),
    percent = fraction * 100,
    label = sprintf("%s (%d, %.1f%%)", region_annotation, peak_count, percent),
    plot_label = sprintf("%s\n%d (%.1f%%)", region_annotation, peak_count, percent)
  )

write_tsv(
  annotation_summary,
  file.path(outdir, "peak_genomic_annotation_summary.tsv")
)

pie_plot <- ggplot(
  annotation_summary,
  aes(x = "", y = peak_count, fill = region_annotation)
) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  theme_void(base_size = 12) +
  labs(
    title = "Peak Genomic Annotation",
    fill = "Region"
  ) +
  geom_text(
    aes(label = plot_label),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      Promoter = "#d95f02",
      Exon = "#1b9e77",
      Intron = "#7570b3",
      Intergenic = "#66a61e"
    )
  )

ggsave(
  filename = file.path(outdir, "peak_genomic_annotation_pie.pdf"),
  plot = pie_plot,
  width = 7,
  height = 6
)
ggsave(
  filename = file.path(outdir, "peak_genomic_annotation_pie.png"),
  plot = pie_plot,
  width = 7,
  height = 6,
  dpi = 300
)

selected_peaks_df <- annotated_df %>%
  filter(
    !is.na(.data[[p_col]]),
    !is.na(.data[[fc_col]]),
    .data[[p_col]] <= 0.05,
    .data[[fc_col]] > 0,
    region_annotation == "Promoter"
  )

selected_peak_gene_records <- selected_peaks_df %>%
  dplyr::select(
    peak_id,
    peak_key,
    chrom,
    start,
    end,
    width,
    strand,
    all_of(fc_col),
    all_of(p_col),
    all_of(fdr_col),
    region_annotation,
    region_annotation_detail,
    promoter_bins,
    promoter_gene_names
  ) %>%
  mutate(promoter_gene_list = lapply(promoter_gene_names, clean_gene_tokens)) %>%
  tidyr::unnest_longer(promoter_gene_list, values_to = "gene_name") %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  distinct()

write_tsv(
  selected_peaks_df,
  file.path(outdir, sprintf("%s_selected_promoter_peaks.tsv", method))
)
write_tsv(
  selected_peak_gene_records,
  file.path(outdir, sprintf("%s_selected_promoter_peak_gene_records.tsv", method))
)

gene_set_df <- read_excel(gene_set_path, sheet = "merged_unique_genes") %>%
  transmute(
    gene_name_set = as.character(gene_name),
    gene_name_normalized = toupper(as.character(gene_name_normalized)),
    source_count = as.integer(source_count),
    source_files = as.character(source_files)
  ) %>%
  distinct()

selected_peak_gene_records <- selected_peak_gene_records %>%
  mutate(gene_name_normalized = toupper(gene_name))

overlap_gene_df <- selected_peak_gene_records %>%
  inner_join(gene_set_df, by = "gene_name_normalized") %>%
  mutate(gene_name = gene_name) %>%
  distinct(
    peak_id, peak_key, chrom, start, end, width, strand,
    gene_name, gene_name_normalized,
    source_count, source_files,
    .data[[fc_col]], .data[[p_col]], .data[[fdr_col]],
    region_annotation, region_annotation_detail, promoter_bins
  ) %>%
  arrange(.data[[p_col]], desc(.data[[fc_col]]), gene_name, chrom, start, end)

write_tsv(
  overlap_gene_df,
  file.path(outdir, sprintf("%s_promoter_gene_overlap_with_112_gene_set.tsv", method))
)

overlap_peak_df <- overlap_gene_df %>%
  group_by(peak_key, chrom, start, end, width, strand) %>%
  summarise(
    overlap_genes = paste(sort(unique(gene_name)), collapse = ";"),
    min_p_value = min(.data[[p_col]], na.rm = TRUE),
    max_logFC = max(.data[[fc_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_p_value, desc(max_logFC), chrom, start, end)

write_tsv(
  overlap_peak_df,
  file.path(outdir, sprintf("%s_promoter_overlap_peaks_for_heatmap.tsv", method))
)

counted_db <- readRDS(counted_dba_path)
count_df <- dba.peakset(counted_db, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
colnames(count_df)[1:3] <- c("chrom", "start", "end")
count_df$peak_key <- build_peak_id(count_df)

count_matrix_df <- count_df %>%
  dplyr::select(peak_key, all_of(sample_order))

heatmap_df <- overlap_peak_df %>%
  left_join(count_matrix_df, by = "peak_key") %>%
  mutate(row_label = paste0(overlap_genes, "\n", chrom, ":", start, "-", end))

heatmap_matrix <- as.matrix(heatmap_df[, sample_order, drop = FALSE])
rownames(heatmap_matrix) <- heatmap_df$row_label
heatmap_matrix <- log2(heatmap_matrix + 1)
heatmap_z <- row_scale_matrix(heatmap_matrix)

write_tsv(
  data.frame(peak_label = rownames(heatmap_matrix), heatmap_matrix, check.names = FALSE),
  file.path(outdir, sprintf("%s_overlap_peak_heatmap_log2_counts.tsv", method))
)
write_tsv(
  data.frame(peak_label = rownames(heatmap_z), heatmap_z, check.names = FALSE),
  file.path(outdir, sprintf("%s_overlap_peak_heatmap_row_zscore.tsv", method))
)

pdf(
  file.path(outdir, sprintf("%s_overlap_peak_heatmap.pdf", method)),
  width = 8,
  height = max(6, nrow(heatmap_z) * 0.22)
)
pheatmap(
  heatmap_z,
  color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_condition_df,
  fontsize_row = 7,
  fontsize_col = 10,
  border_color = NA,
  main = sprintf(
    "%s promoter-overlap peaks in 112-gene set",
    method
  )
)
dev.off()

selected_gene_symbols <- sort(unique(selected_peak_gene_records$gene_name))
background_gene_symbols <- sort(
  unique(unlist(lapply(annotated_df$promoter_gene_list, identity)))
)

selected_entrez <- bitr(
  selected_gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
background_entrez <- bitr(
  background_gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

selected_symbol_map <- selected_entrez %>%
  distinct(ENTREZID, .keep_all = TRUE)
background_symbol_map <- background_entrez %>%
  distinct(ENTREZID, .keep_all = TRUE)

ego <- enrichGO(
  gene = selected_symbol_map$ENTREZID,
  universe = background_symbol_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

enrich_df <- normalize_enrich_df(ego)

write_tsv(
  enrich_df,
  file.path(outdir, sprintf("%s_selected_promoter_gene_enrichGO_BP.tsv", method))
)

top10_df <- enrich_df %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(unique(Description))
    )
  )

if (nrow(top10_df) > 0) {
  top10_plot <- ggplot(
    top10_df,
    aes(x = Count, y = Description, fill = p.adjust)
  ) +
    geom_col() +
    scale_fill_gradient(low = "#d73027", high = "#fee090") +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf("Top 10 enriched pathways (%s, GO-BP)", method),
      x = "Gene count",
      y = NULL,
      fill = "Adjusted p"
    )

  ggsave(
    filename = file.path(outdir, sprintf("%s_enrichment_top10_barplot.pdf", method)),
    plot = top10_plot,
    width = 9,
    height = 5.5
  )
  ggsave(
    filename = file.path(outdir, sprintf("%s_enrichment_top10_barplot.png", method)),
    plot = top10_plot,
    width = 9,
    height = 5.5,
    dpi = 300
  )
}

broad_tcell_pattern <- paste(
  c(
    "\\bT[ -]?cell\\b",
    "\\bT lymphocyte\\b",
    "\\bthymocyte\\b",
    "\\bTCR\\b",
    "\\bTh1\\b",
    "\\bTh2\\b",
    "\\bTh17\\b",
    "\\bimmune\\b",
    "\\blymph(?:ocyte|oid)?\\b",
    "\\bcytokine\\b",
    "\\bchemokine\\b",
    "\\bantigen\\b",
    "\\bleukocyte\\b",
    "\\bcell adhesion\\b",
    "\\bCAM\\b",
    "\\bPD-1\\b",
    "\\bPD-L1\\b",
    "\\bNF-kappa\\b",
    "\\bJAK\\b",
    "\\bSTAT\\b",
    "\\bTGF\\b",
    "\\bMAPK\\b",
    "\\bPI3K\\b",
    "\\bmTOR\\b",
    "\\bNK\\b"
  ),
  collapse = "|"
)

tcell_df <- enrich_df %>%
  filter(grepl(broad_tcell_pattern, Description, ignore.case = TRUE)) %>%
  arrange(pvalue, desc(Count)) %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(unique(Description))
    )
  )

write_tsv(
  tcell_df,
  file.path(outdir, sprintf("%s_tcell_related_enrichment.tsv", method))
)

if (nrow(tcell_df) > 0) {
  plot_bar_terms(
    tcell_df,
    sprintf("Broad T-cell/immune related pathways (%s, GO-BP)", method),
    "#4575b4",
    "#e0f3f8",
    file.path(outdir, sprintf("%s_tcell_related_enrichment_barplot", method)),
    fill_column = "pvalue",
    fill_label = "p value"
  )

  tcell_top20_df <- tcell_df %>%
    arrange(pvalue, desc(Count)) %>%
    slice_head(n = 20)
  write_tsv(
    tcell_top20_df,
    file.path(outdir, sprintf("%s_tcell_related_enrichment_top20.tsv", method))
  )
  plot_bar_terms(
    tcell_top20_df,
    sprintf("Top 20 broad T-cell/immune pathways (%s, GO-BP)", method),
    "#2166ac",
    "#d1e5f0",
    file.path(outdir, sprintf("%s_tcell_related_enrichment_top20_barplot", method)),
    fill_column = "pvalue",
    fill_label = "p value"
  )
}

ekegg <- tryCatch(
  enrichKEGG(
    gene = selected_symbol_map$ENTREZID,
    universe = background_symbol_map$ENTREZID,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH"
  ),
  error = function(e) NULL
)

kegg_df <- normalize_enrich_df(ekegg)
write_tsv(
  kegg_df,
  file.path(outdir, sprintf("%s_selected_promoter_gene_enrichKEGG.tsv", method))
)

kegg_tcell_df <- kegg_df %>%
  filter(grepl(broad_tcell_pattern, Description, ignore.case = TRUE)) %>%
  arrange(pvalue, desc(Count))

write_tsv(
  kegg_tcell_df,
  file.path(outdir, sprintf("%s_tcell_related_kegg_enrichment.tsv", method))
)

if (nrow(kegg_tcell_df) > 0) {
  plot_bar_terms(
    kegg_tcell_df,
    sprintf("Broad T-cell/immune related pathways (%s, KEGG)", method),
    "#1b7837",
    "#d9f0d3",
    file.path(outdir, sprintf("%s_tcell_related_kegg_enrichment_barplot", method)),
    fill_column = "pvalue",
    fill_label = "p value"
  )

  kegg_tcell_top20_df <- kegg_tcell_df %>%
    filter(pvalue < 0.05) %>%
    arrange(pvalue, desc(Count)) %>%
    slice_head(n = 20)
  write_tsv(
    kegg_tcell_top20_df,
    file.path(outdir, sprintf("%s_tcell_related_kegg_enrichment_top20.tsv", method))
  )
  if (nrow(kegg_tcell_top20_df) > 0) {
    plot_bar_terms(
      kegg_tcell_top20_df,
      sprintf("Broad T-cell/immune KEGG pathways with p < 0.05 (%s)", method),
      "#238b45",
      "#c7e9c0",
      file.path(outdir, sprintf("%s_tcell_related_kegg_enrichment_top20_barplot", method)),
      fill_column = "pvalue",
      fill_label = "p value"
    )
  } else {
    plot_empty_message(
      sprintf("Broad T-cell/immune KEGG pathways with p < 0.05 (%s)", method),
      "No pathways passed p < 0.05",
      file.path(outdir, sprintf("%s_tcell_related_kegg_enrichment_top20_barplot", method))
    )
  }
}

summary_df <- data.frame(
  method = method,
  annotation_total_peaks = nrow(annotated_df),
  promoter_selected_peaks = nrow(selected_peaks_df),
  promoter_selected_gene_records = nrow(selected_peak_gene_records),
  promoter_selected_unique_genes = length(selected_gene_symbols),
  overlap_gene_set_size = nrow(gene_set_df),
  overlap_peak_count = nrow(overlap_peak_df),
  overlap_gene_count = length(unique(overlap_gene_df$gene_name)),
  enrichment_term_count = nrow(enrich_df),
  tcell_enrichment_term_count = nrow(tcell_df),
  kegg_enrichment_term_count = nrow(kegg_df),
  kegg_tcell_enrichment_term_count = nrow(kegg_tcell_df),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_df,
  file.path(outdir, sprintf("%s_appendix_analysis_summary.tsv", method))
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1) {
  script_path <- sub("^--file=", "", script_arg)
  file.copy(
    script_path,
    file.path(outdir, basename(script_path)),
    overwrite = TRUE
  )
}

cat("Appendix analysis completed in", outdir, "\n")
