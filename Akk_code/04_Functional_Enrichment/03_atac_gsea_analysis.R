#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    paste(
      "Usage: run_atac_gsea.R",
      "<annotated_report_tsv> <gmt_dir> <outdir>",
      "[method=DESeq2] [min_size=10] [max_size=500]"
    )
  )
}

annotated_report_path <- args[[1]]
gmt_dir <- args[[2]]
outdir <- args[[3]]
method <- if (length(args) >= 4) args[[4]] else "DESeq2"
min_size <- if (length(args) >= 5) as.integer(args[[5]]) else 10L
max_size <- if (length(args) >= 6) as.integer(args[[6]]) else 500L

userlib <- normalizePath("r_libs", mustWork = FALSE)
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(fgsea)
  library(readr)
  library(tidyr)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(
    df,
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
  tokens <- tokens[!tolower(tokens) %in% c("na", "nan", "none")]
  unique(tokens)
}

get_method_columns <- function(method_name) {
  if (method_name == "DESeq2") {
    return(list(fc = "DESeq2_Fold", p = "DESeq2_p.value"))
  }
  if (method_name == "edgeR") {
    return(list(fc = "edgeR_Fold", p = "edgeR_p.value"))
  }
  stop("method must be DESeq2 or edgeR")
}

rank_gene_scores <- function(df, fc_col, p_col) {
  rows <- vector("list", nrow(df))
  keep_idx <- 0L

  for (i in seq_len(nrow(df))) {
    genes <- clean_gene_tokens(df$promoter_gene_names[[i]])
    if (length(genes) == 0) {
      next
    }

    fc <- suppressWarnings(as.numeric(df[[fc_col]][[i]]))
    pv <- suppressWarnings(as.numeric(df[[p_col]][[i]]))
    if (is.na(fc) || is.na(pv)) {
      next
    }

    pv <- max(pv, 1e-300)
    score <- fc * (-log10(pv))
    keep_idx <- keep_idx + 1L
    rows[[keep_idx]] <- data.frame(
      gene_name = toupper(genes),
      score = score,
      logFC = fc,
      pvalue = pv,
      stringsAsFactors = FALSE
    )
  }

  if (keep_idx == 0L) {
    stop("No promoter gene records available for ranking.")
  }

  rank_df <- bind_rows(rows[seq_len(keep_idx)]) %>%
    group_by(gene_name) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(score), gene_name)

  rank_vec <- rank_df$score
  names(rank_vec) <- rank_df$gene_name
  list(rank_df = rank_df, rank_vec = rank_vec)
}

annotated_df <- read.delim(
  annotated_report_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

method_cols <- get_method_columns(method)
ranking <- rank_gene_scores(
  annotated_df,
  fc_col = method_cols$fc,
  p_col = method_cols$p
)

write_tsv(
  ranking$rank_df,
  file.path(outdir, sprintf("%s_gene_ranking_for_gsea.tsv", method))
)

gmt_files <- sort(list.files(gmt_dir, pattern = "\\.gmt$", full.names = TRUE))
if (length(gmt_files) == 0) {
  stop("No GMT files found in gmt_dir.")
}

summary_rows <- list()
summary_idx <- 0L

for (gmt_file in gmt_files) {
  gmt_name <- basename(gmt_file)
  collection_name <- sub("\\.gmt$", "", gmt_name)

  term2gene <- read.gmt(gmt_file) %>%
    transmute(term = term, gene = toupper(gene))

  gsea_res <- tryCatch(
    GSEA(
      geneList = ranking$rank_vec,
      TERM2GENE = term2gene,
      minGSSize = min_size,
      maxGSSize = max_size,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(gsea_res)) {
    full_df <- data.frame()
  } else {
    full_df <- as.data.frame(gsea_res@result)
  }

  if (nrow(full_df) == 0) {
    full_df <- data.frame(
      ID = character(0),
      Description = character(0),
      setSize = integer(0),
      enrichmentScore = numeric(0),
      NES = numeric(0),
      pvalue = numeric(0),
      p.adjust = numeric(0),
      qvalue = numeric(0),
      rank = integer(0),
      leading_edge = character(0),
      core_enrichment = character(0)
    )
  } else {
    full_df <- full_df %>%
      arrange(p.adjust, pvalue, desc(abs(NES)))
  }

  fdr05_df <- full_df %>%
    filter(!is.na(p.adjust), p.adjust < 0.05)
  fdr25_df <- full_df %>%
    filter(!is.na(p.adjust), p.adjust < 0.25)

  write_tsv(
    full_df,
    file.path(outdir, sprintf("%s_%s_gsea_all.tsv", method, collection_name))
  )
  write_tsv(
    fdr05_df,
    file.path(outdir, sprintf("%s_%s_gsea_FDRlt0.05.tsv", method, collection_name))
  )
  write_tsv(
    fdr25_df,
    file.path(outdir, sprintf("%s_%s_gsea_FDRlt0.25.tsv", method, collection_name))
  )

  summary_idx <- summary_idx + 1L
  summary_rows[[summary_idx]] <- data.frame(
    collection = collection_name,
    total_terms = nrow(full_df),
    sig_terms_fdr_lt_0_05 = nrow(fdr05_df),
    sig_terms_fdr_lt_0_25 = nrow(fdr25_df),
    top_term = if (nrow(full_df) > 0) full_df$Description[[1]] else "",
    top_term_fdr = if (nrow(full_df) > 0) full_df$p.adjust[[1]] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
write_tsv(
  summary_df,
  file.path(outdir, sprintf("%s_gsea_collection_summary.tsv", method))
)

sig_all_fdr05 <- bind_rows(lapply(
  gmt_files,
  function(gmt_file) {
    collection_name <- sub("\\.gmt$", "", basename(gmt_file))
    path <- file.path(outdir, sprintf("%s_%s_gsea_FDRlt0.05.tsv", method, collection_name))
    df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(df) == 0) {
      return(df)
    }
    df$collection <- collection_name
    df
  }
))

sig_all_fdr25 <- bind_rows(lapply(
  gmt_files,
  function(gmt_file) {
    collection_name <- sub("\\.gmt$", "", basename(gmt_file))
    path <- file.path(outdir, sprintf("%s_%s_gsea_FDRlt0.25.tsv", method, collection_name))
    df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(df) == 0) {
      return(df)
    }
    df$collection <- collection_name
    df
  }
))

if (nrow(sig_all_fdr05) > 0) {
  sig_all_fdr05 <- sig_all_fdr05 %>%
    arrange(p.adjust, pvalue, desc(abs(NES)))
}
if (nrow(sig_all_fdr25) > 0) {
  sig_all_fdr25 <- sig_all_fdr25 %>%
    arrange(p.adjust, pvalue, desc(abs(NES)))
}

write_tsv(
  sig_all_fdr05,
  file.path(outdir, sprintf("%s_gsea_significant_pathways_FDRlt0.05.tsv", method))
)
write_tsv(
  sig_all_fdr25,
  file.path(outdir, sprintf("%s_gsea_significant_pathways_FDRlt0.25.tsv", method))
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

cat("GSEA completed in", outdir, "\n")
