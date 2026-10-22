suppressPackageStartupMessages({
  library(readxl)
  library(VennDiagram)
  library(grid)
})

base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)

atac_file <- file.path(base_dir, "ATACseq差异分析-p小于0.05&FC大于0-基因名去编号.xlsx")
rna_file <- file.path(base_dir, "转录组-p小于0.05&FC大于0.xlsx")
pathway_dir <- file.path(base_dir, "6通路tsv文件")

out_png <- file.path(base_dir, "ATAC_RNA_6通路_韦恩图.png")
out_pdf <- file.path(base_dir, "ATAC_RNA_6通路_韦恩图.pdf")
out_summary <- file.path(base_dir, "ATAC_RNA_6通路_韦恩图统计.tsv")

read_gene_column <- function(file, column, split_by = NULL) {
  df <- read_excel(file, col_types = "text")
  values <- trimws(df[[column]])
  values <- values[!is.na(values) & values != "" & toupper(values) != "NA"]

  if (!is.null(split_by)) {
    values <- unlist(strsplit(values, split_by, fixed = TRUE), use.names = FALSE)
    values <- trimws(values)
    values <- values[values != "" & toupper(values) != "NA"]
  }

  unique(toupper(values))
}

read_pathway_gene_sets <- function(folder) {
  files <- sort(list.files(folder, pattern = "\\.tsv$", full.names = TRUE))
  if (length(files) == 0) {
    stop("6通路tsv文件 目录下没有找到 tsv 文件。")
  }

  raw_records <- character()
  file_stats <- vector("list", length(files))

  for (i in seq_along(files)) {
    df <- read.delim(
      files[i],
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = ""
    )
    key_col <- names(df)[1]
    value_col <- names(df)[2]
    hit <- which(toupper(trimws(df[[key_col]])) == "GENE_SYMBOLS")
    if (length(hit) == 0) {
      genes <- character()
    } else {
      genes <- unlist(strsplit(df[hit[1], value_col], ",", fixed = TRUE), use.names = FALSE)
      genes <- trimws(genes)
      genes <- genes[genes != ""]
    }

    raw_records <- c(raw_records, genes)
    file_stats[[i]] <- data.frame(
      source_file = basename(files[i]),
      gene_count_raw = length(genes),
      gene_count_unique_case_insensitive = length(unique(toupper(genes))),
      stringsAsFactors = FALSE
    )
  }

  list(
    raw_records = raw_records,
    unique_genes = unique(toupper(raw_records)),
    file_stats = do.call(rbind, file_stats)
  )
}

atac_genes <- read_gene_column(atac_file, "promoter_gene_names", split_by = ";")
rna_genes <- read_gene_column(rna_file, "symbol")
pathway_info <- read_pathway_gene_sets(pathway_dir)
pathway_genes <- pathway_info$unique_genes

area1 <- length(atac_genes)
area2 <- length(rna_genes)
area3 <- length(pathway_genes)
pathway_raw_n <- length(pathway_info$raw_records)

n12 <- length(intersect(atac_genes, rna_genes))
n13 <- length(intersect(atac_genes, pathway_genes))
n23 <- length(intersect(rna_genes, pathway_genes))
n123 <- length(Reduce(intersect, list(atac_genes, rna_genes, pathway_genes)))

regions <- data.frame(
  region = c(
    "ATAC_only",
    "RNA_only",
    "Pathway_only",
    "ATAC_RNA_only",
    "ATAC_Pathway_only",
    "RNA_Pathway_only",
    "ATAC_RNA_Pathway"
  ),
  count = c(
    length(setdiff(atac_genes, union(rna_genes, pathway_genes))),
    length(setdiff(rna_genes, union(atac_genes, pathway_genes))),
    length(setdiff(pathway_genes, union(atac_genes, rna_genes))),
    length(setdiff(intersect(atac_genes, rna_genes), pathway_genes)),
    length(setdiff(intersect(atac_genes, pathway_genes), rna_genes)),
    length(setdiff(intersect(rna_genes, pathway_genes), atac_genes)),
    n123
  ),
  stringsAsFactors = FALSE
)

pairwise <- data.frame(
  comparison = c("ATAC_vs_RNA", "ATAC_vs_6Pathway", "RNA_vs_6Pathway"),
  overlap_count = c(n12, n13, n23),
  stringsAsFactors = FALSE
)

set_sizes <- data.frame(
  set_name = c(
    "ATAC_promoter_genes_unique_case_insensitive",
    "RNA_symbol_unique_case_insensitive",
    "Pathway_raw_records",
    "Pathway_unique_genes_case_insensitive"
  ),
  count = c(area1, area2, pathway_raw_n, area3),
  stringsAsFactors = FALSE
)

summary_long <- rbind(
  data.frame(
    section = "set_sizes",
    item = set_sizes$set_name,
    value = set_sizes$count,
    stringsAsFactors = FALSE
  ),
  data.frame(
    section = "pairwise_overlap",
    item = pairwise$comparison,
    value = pairwise$overlap_count,
    stringsAsFactors = FALSE
  ),
  data.frame(
    section = "venn_regions",
    item = regions$region,
    value = regions$count,
    stringsAsFactors = FALSE
  )
)

write.table(
  summary_long,
  file = out_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

draw_plot <- function() {
  grid.newpage()
  pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc")))
  venn.plot <- draw.triple.venn(
    area1 = area3,
    area2 = area2,
    area3 = area1,
    n12 = n23,
    n23 = n12,
    n13 = n13,
    n123 = n123,
    category = c(
      sprintf("6 pathway genes\n(%d unique)\n(%d records)", area3, pathway_raw_n),
      sprintf("RNA-seq (DEGs)\n(%d)", area2),
      sprintf("ATAC-seq (promoter genes)\n(%d)", area1)
    ),
    fill = c("#F4DEC8", "#CFE7C4", "#BFD5EA"),
    col = c("#F28C28", "#459B2E", "#2D73BE"),
    alpha = c(0.55, 0.55, 0.55),
    cex = 2.1,
    fontface = "bold",
    cat.cex = 1.4,
    cat.fontface = "bold",
    cat.col = c("#8A4D00", "#2B7A22", "#15317E"),
    cat.pos = c(145, 215, 0),
    cat.dist = c(0.04, 0.08, 0.08),
    cat.just = list(c(0, 1), c(1, 1), c(0.5, 0)),
    margin = 0.08,
    scaled = FALSE,
    lwd = 2,
    rotation.degree = 180
  )
  grid.draw(venn.plot)
  popViewport()
}

png(out_png, width = 2000, height = 2000, res = 220, bg = "white")
draw_plot()
dev.off()

pdf(out_pdf, width = 9, height = 9, bg = "white")
draw_plot()
dev.off()

message("ATAC unique genes: ", area1)
message("RNA unique genes: ", area2)
message("Pathway raw records: ", pathway_raw_n)
message("Pathway unique genes: ", area3)
message("ATAC ∩ RNA: ", n12)
message("ATAC ∩ Pathway: ", n13)
message("RNA ∩ Pathway: ", n23)
message("ATAC ∩ RNA ∩ Pathway: ", n123)
message("PNG: ", out_png)
message("PDF: ", out_pdf)
message("Summary: ", out_summary)
