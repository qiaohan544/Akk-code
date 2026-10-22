#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    paste(
      "Usage: make_selected_pathway_gsea_plots.R",
      "<ranking_tsv> <gmt_dir> <outdir>",
      "[method=DESeq2]"
    )
  )
}

ranking_path <- args[[1]]
gmt_dir <- args[[2]]
outdir <- args[[3]]
method <- if (length(args) >= 4) args[[4]] else "DESeq2"

userlib <- normalizePath("r_libs", mustWork = FALSE)
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(readr)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

target_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "LEE_DIFFERENTIATING_T_LYMPHOCYTE",
  "GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION",
  "GOBP_TRANSLATIONAL_INITIATION",
  "REACTOME_TRANSLATION"
)

display_title_map <- c(
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "Interferon Gamma Response",
  LEE_DIFFERENTIATING_T_LYMPHOCYTE = "Differentiating T Lymphocyte",
  GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION = "Epigenetic Regulation of Gene Expression",
  GOBP_TRANSLATIONAL_INITIATION = "Translational Initiation",
  REACTOME_TRANSLATION = "Reactome Translation"
)

write_tsv <- function(df, path) {
  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

save_plot_both <- function(plot_obj, file_prefix, width, height) {
  ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

ranking_df <- read.delim(
  ranking_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("gene_name", "score")
missing_cols <- setdiff(required_cols, colnames(ranking_df))
if (length(missing_cols) > 0) {
  stop(sprintf("Ranking table missing columns: %s", paste(missing_cols, collapse = ", ")))
}

gene_list <- ranking_df$score
names(gene_list) <- toupper(ranking_df$gene_name)
gene_list <- sort(gene_list, decreasing = TRUE)

gmt_files <- sort(list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE))
if (length(gmt_files) == 0) {
  stop("No GMT files found in gmt_dir.")
}

term2gene <- bind_rows(lapply(
  gmt_files,
  function(gmt_file) {
    collection_name <- sub("\\.gmt$", "", basename(gmt_file))
    read.gmt(gmt_file) %>%
      transmute(
        term = term,
        gene = toupper(gene),
        collection = collection_name
      )
  }
))

missing_terms <- setdiff(target_pathways, unique(term2gene$term))
if (length(missing_terms) > 0) {
  stop(
    sprintf(
      "Requested pathways not found in GMT files: %s",
      paste(missing_terms, collapse = ", ")
    )
  )
}

gsea_res <- GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene %>% dplyr::select(term, gene),
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_df <- as.data.frame(gsea_res@result) %>%
  left_join(
    term2gene %>%
      distinct(term, collection),
    by = c("ID" = "term")
  ) %>%
  arrange(pvalue, p.adjust, desc(abs(NES)))

selected_df <- gsea_df %>%
  filter(ID %in% target_pathways) %>%
  mutate(
    display_title = unname(display_title_map[ID]),
    subtitle_text = sprintf(
      "ES = %.3f    NES = %.2f    P = %.2g    FDR = %.3f    Gene-set size = %d",
      enrichmentScore, NES, pvalue, p.adjust, setSize
    )
  )

write_tsv(
  selected_df,
  file.path(outdir, sprintf("%s_selected_pathway_gsea_stats.tsv", method))
)

for (i in seq_len(nrow(selected_df))) {
  row <- selected_df[i, , drop = FALSE]
  pathway_id <- row$ID[[1]]
  plot_title <- paste0(row$display_title[[1]], "\n", row$subtitle_text[[1]])

  p <- gseaplot2(
    gsea_res,
    geneSetID = pathway_id,
    title = plot_title,
    color = "#b23a48",
    base_size = 11,
    rel_heights = c(1.5, 0.45, 1),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom = "line"
  )

  save_plot_both(
    p,
    file.path(outdir, sprintf("%s_%s_gsea_plot", method, pathway_id)),
    width = 8.5,
    height = 6.2
  )
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (length(script_path) == 1 && nzchar(script_path) && file.exists(script_path)) {
  file.copy(
    script_path,
    file.path(outdir, basename(script_path)),
    overwrite = TRUE
  )
}

message("wrote GSEA plots to: ", outdir)
