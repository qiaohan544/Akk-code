#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(vegan)
})

options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- ""./data/"投稿定图/cohort"
rds_path <- file.path(root_dir, "phyloseq_object.rds")
out_dir <- file.path(root_dir, "structured_analysis", "enterosignature_genus_all_samples")
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

depth_cutoff <- 20000L
rank_k <- 5L
nmf_runs <- 20L
nmf_max_iter <- 600L
nmf_tol <- 1e-6
min_genus_prevalence <- 0.005
min_genus_mean <- 0.00001
pseudocount <- 1e-6

extract_sample_data <- function(obj) {
  sample_df <- attr(obj, "sam_data", exact = TRUE)
  class(sample_df) <- "data.frame"
  rn <- attr(sample_df, "row.names", exact = TRUE)
  if (!is.null(rn)) {
    rownames(sample_df) <- rn
  }
  sample_df
}

extract_matrix_like <- function(obj, attr_name) {
  x <- attr(obj, attr_name, exact = TRUE)
  x_dimnames <- attr(x, "dimnames", exact = TRUE)
  suppressWarnings(class(x) <- NULL)
  if (!is.null(x_dimnames)) {
    dimnames(x) <- x_dimnames
  }
  x
}

clean_chr <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA")] <- NA_character_
  x_chr
}

normalize_response <- function(x) {
  x_chr <- clean_chr(x)
  dplyr::case_when(
    x_chr %in% c("Responder", "R", "response", "Response") ~ "Responder",
    x_chr %in% c("Non_Responder", "Non-responder", "NR", "No response") ~ "Non_Responder",
    TRUE ~ NA_character_
  )
}

safe_glm <- function(formula, data) {
  tryCatch(
    suppressWarnings(glm(formula, data = data, family = binomial())),
    error = function(e) NULL
  )
}

extract_glm_term <- function(fit, term_name) {
  if (is.null(fit)) {
    return(NULL)
  }
  coef_df <- as.data.frame(summary(fit)$coefficients)
  coef_df$term <- rownames(coef_df)
  row <- coef_df[coef_df$term == term_name, , drop = FALSE]
  if (nrow(row) == 0) {
    return(NULL)
  }
  beta <- row$Estimate[1]
  se <- row$`Std. Error`[1]
  data.frame(
    term = term_name,
    beta = beta,
    odds_ratio = exp(beta),
    or_low = exp(beta - 1.96 * se),
    or_high = exp(beta + 1.96 * se),
    p_value = row$`Pr(>|z|)`[1],
    stringsAsFactors = FALSE
  )
}

kl_divergence <- function(x, reconstruction, eps = 1e-12) {
  reconstruction <- pmax(reconstruction, eps)
  positive <- x > 0
  sum(x[positive] * log(x[positive] / reconstruction[positive]) - x[positive] + reconstruction[positive]) +
    sum(reconstruction[!positive])
}

run_kl_nmf <- function(x, rank_k, nruns = 20L, max_iter = 600L, tol = 1e-6, seed = 20260502L) {
  set.seed(seed)
  x <- as.matrix(x)
  n_features <- nrow(x)
  n_samples <- ncol(x)
  eps <- 1e-12
  best <- NULL
  run_summary <- vector("list", nruns)

  for (run_idx in seq_len(nruns)) {
    w <- matrix(runif(n_features * rank_k, min = 0.001, max = 1), nrow = n_features, ncol = rank_k)
    h <- matrix(runif(rank_k * n_samples, min = 0.001, max = 1), nrow = rank_k, ncol = n_samples)
    w <- sweep(w, 2, colSums(w), "/")
    previous_obj <- Inf
    current_obj <- Inf

    for (iter_idx in seq_len(max_iter)) {
      reconstruction <- pmax(w %*% h, eps)
      h <- h * ((t(w) %*% (x / reconstruction)) / pmax(matrix(colSums(w), nrow = rank_k, ncol = n_samples), eps))
      h <- pmax(h, eps)

      reconstruction <- pmax(w %*% h, eps)
      w <- w * (((x / reconstruction) %*% t(h)) / pmax(matrix(rowSums(h), nrow = n_features, ncol = rank_k, byrow = TRUE), eps))
      w <- pmax(w, eps)

      w_colsum <- colSums(w)
      w <- sweep(w, 2, w_colsum, "/")
      h <- sweep(h, 1, w_colsum, "*")

      if (iter_idx %% 50L == 0L || iter_idx == max_iter) {
        current_obj <- kl_divergence(x, w %*% h, eps = eps)
        if (is.finite(previous_obj)) {
          rel_change <- abs(previous_obj - current_obj) / max(abs(previous_obj), eps)
          if (rel_change < tol) {
            break
          }
        }
        previous_obj <- current_obj
      }
    }

    reconstruction <- pmax(w %*% h, eps)
    rss <- sum((x - reconstruction)^2)
    evar <- 1 - rss / sum(x^2)
    cosine <- sum(x * reconstruction) / sqrt(sum(x^2) * sum(reconstruction^2))
    current_obj <- kl_divergence(x, reconstruction, eps = eps)

    run_summary[[run_idx]] <- data.frame(
      run = run_idx,
      objective = current_obj,
      explained_variance = evar,
      cosine_similarity = cosine,
      iterations = iter_idx,
      stringsAsFactors = FALSE
    )

    if (is.null(best) || evar > best$explained_variance) {
      best <- list(W = w, H = h, objective = current_obj, explained_variance = evar, cosine_similarity = cosine, run = run_idx, iterations = iter_idx)
    }
  }

  best$run_summary <- bind_rows(run_summary)
  best
}

label_enterosignatures <- function(w_norm) {
  top_genus <- apply(w_norm, 2, function(x) rownames(w_norm)[which.max(x)])
  labels <- paste0("ES", seq_along(top_genus), "_", gsub("[^A-Za-z0-9]+", "_", top_genus))
  labels
}

message("Reading phyloseq-like object")
obj <- readRDS(rds_path)
sample_df <- extract_sample_data(obj)
tax_mat <- extract_matrix_like(obj, "tax_table")
otu_mat <- extract_matrix_like(obj, "otu_table")

sample_df$sample_id <- rownames(sample_df)
sample_df$Response <- normalize_response(sample_df$Response)
sample_df$SPT <- clean_chr(sample_df$SPT)
sample_df$Sampling_Timepoint <- clean_chr(sample_df$Sampling_Timepoint)
sample_df$project <- clean_chr(sample_df$project)

bacterial_taxa <- rownames(tax_mat)[tax_mat[, "Kingdom"] == "Bacteria"]
otu_mat <- otu_mat[intersect(rownames(otu_mat), bacterial_taxa), sample_df$sample_id, drop = FALSE]
sample_df$read_depth_bacteria <- colSums(otu_mat[, sample_df$sample_id, drop = FALSE])

analysis_df <- sample_df %>%
  filter(Host == "Human", AssayType == "WGS", read_depth_bacteria >= depth_cutoff) %>%
  arrange(project, sample_id)

analysis_otu <- otu_mat[, analysis_df$sample_id, drop = FALSE]

tax_df <- data.frame(
  feature = rownames(tax_mat),
  genus = clean_chr(tax_mat[, "Genus"]),
  species = clean_chr(tax_mat[, "Species"]),
  stringsAsFactors = FALSE
)
tax_df$genus[is.na(tax_df$genus) | tax_df$genus %in% c("?", "unknown", "Unknown", "unclassified")] <- NA_character_
tax_df <- tax_df %>% filter(feature %in% rownames(analysis_otu), !is.na(genus))

genus_counts <- rowsum(analysis_otu[tax_df$feature, , drop = FALSE], group = tax_df$genus, reorder = TRUE)
genus_rel <- sweep(genus_counts, 2, colSums(genus_counts), "/")
genus_rel[!is.finite(genus_rel)] <- 0

genus_prevalence <- rowMeans(genus_rel > 0)
genus_mean <- rowMeans(genus_rel)
keep_genera <- rownames(genus_rel)[genus_prevalence >= min_genus_prevalence | genus_mean >= min_genus_mean]
nmf_input <- genus_rel[keep_genera, , drop = FALSE]
low_genera <- setdiff(rownames(genus_rel), keep_genera)
if (length(low_genera) > 0) {
  nmf_input <- rbind(nmf_input, Other_low_abundance = colSums(genus_rel[low_genera, , drop = FALSE]))
}
nmf_input <- nmf_input[rowSums(nmf_input) > 0, colSums(nmf_input) > 0, drop = FALSE]

write.table(nmf_input, file.path(tables_dir, "genus_relative_abundance_for_nmf.tsv"), sep = "\t", quote = FALSE, col.names = NA)

message(sprintf("Running KL-NMF enterosignatures on %d genera x %d samples", nrow(nmf_input), ncol(nmf_input)))
nmf_res <- run_kl_nmf(nmf_input, rank_k = rank_k, nruns = nmf_runs, max_iter = nmf_max_iter, tol = nmf_tol)
W <- nmf_res$W
H <- nmf_res$H
rownames(W) <- rownames(nmf_input)
colnames(W) <- paste0("ES", seq_len(rank_k))
rownames(H) <- colnames(W)
colnames(H) <- colnames(nmf_input)

W_norm <- sweep(W, 2, colSums(W), "/")
es_labels <- label_enterosignatures(W_norm)
colnames(W_norm) <- es_labels
rownames(H) <- es_labels

H_prop <- sweep(H, 2, colSums(H), "/")
H_prop[!is.finite(H_prop)] <- 0
rownames(H_prop) <- es_labels

primary_idx <- apply(H_prop, 2, which.max)
assignment_df <- analysis_df %>%
  select(sample_id, project, SPT, Sampling_Timepoint, Response, read_depth_bacteria, Antibiotics, Therapy, Host_disease) %>%
  mutate(
    enterosignature = factor(es_labels[primary_idx[sample_id]], levels = es_labels),
    enterosignature_score = H_prop[cbind(primary_idx[sample_id], match(sample_id, colnames(H_prop)))]
  )

es_score_df <- as.data.frame(t(H_prop))
es_score_df$sample_id <- rownames(es_score_df)
assignment_df <- assignment_df %>% left_join(es_score_df, by = "sample_id")

genus_meta <- data.frame(
  genus = rownames(genus_rel),
  prevalence = genus_prevalence[rownames(genus_rel)],
  mean_relative_abundance = genus_mean[rownames(genus_rel)],
  retained_for_nmf = rownames(genus_rel) %in% keep_genera,
  stringsAsFactors = FALSE
) %>% arrange(desc(mean_relative_abundance))

composition_df <- as.data.frame(W_norm) %>%
  mutate(genus = rownames(W_norm)) %>%
  pivot_longer(cols = all_of(es_labels), names_to = "enterosignature", values_to = "weight") %>%
  arrange(enterosignature, desc(weight), genus)

run_summary <- nmf_res$run_summary %>% arrange(desc(explained_variance))

fwrite(assignment_df, file.path(tables_dir, "enterosignature_sample_assignments.csv"))
fwrite(composition_df, file.path(tables_dir, "enterosignature_genus_composition_long.csv"))
fwrite(genus_meta, file.path(tables_dir, "genus_filter_summary.csv"))
fwrite(run_summary, file.path(tables_dir, "nmf_run_summary.csv"))
write.table(W_norm, file.path(tables_dir, "enterosignature_W_genus_weights.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(H_prop, file.path(tables_dir, "enterosignature_H_sample_weights.tsv"), sep = "\t", quote = FALSE, col.names = NA)

baseline_df <- assignment_df %>%
  filter(SPT == "Baseline", Response %in% c("Non_Responder", "Responder")) %>%
  mutate(response_num = ifelse(Response == "Responder", 1, 0))

response_counts <- baseline_df %>%
  count(enterosignature, Response, name = "n") %>%
  group_by(enterosignature) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

response_summary <- baseline_df %>%
  group_by(enterosignature) %>%
  summarise(
    n = n(),
    responders = sum(Response == "Responder"),
    nonresponders = sum(Response == "Non_Responder"),
    responder_rate = mean(Response == "Responder"),
    .groups = "drop"
  )

response_overall_test <- suppressWarnings(chisq.test(table(baseline_df$enterosignature, baseline_df$Response)))
response_test_df <- data.frame(
  test = "enterosignature_vs_response_chisq",
  n = nrow(baseline_df),
  statistic = unname(response_overall_test$statistic),
  df = unname(response_overall_test$parameter),
  p_value = response_overall_test$p.value,
  min_expected_n = min(response_overall_test$expected),
  low_expected_flag = any(response_overall_test$expected < 5),
  stringsAsFactors = FALSE
)

response_es_effects <- lapply(es_labels, function(es_name) {
  model_df <- baseline_df %>%
    mutate(es_z = as.numeric(scale(.data[[es_name]])))
  if (length(unique(model_df$response_num)) < 2 || sd(model_df$es_z, na.rm = TRUE) == 0) {
    return(NULL)
  }
  fit <- safe_glm(response_num ~ es_z + factor(project), model_df)
  term <- extract_glm_term(fit, "es_z")
  if (is.null(term)) {
    return(NULL)
  }
  term %>% mutate(enterosignature = es_name, n = nrow(model_df), responders = sum(model_df$response_num), nonresponders = sum(model_df$response_num == 0))
}) %>% bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

akk_genus <- genus_rel[rownames(genus_rel) == "Akkermansia", assignment_df$sample_id, drop = TRUE]
if (length(akk_genus) == 0) {
  akk_genus <- rep(0, nrow(assignment_df))
  names(akk_genus) <- assignment_df$sample_id
}
assignment_df$Akkermansia_rel_abundance <- as.numeric(akk_genus[assignment_df$sample_id])
baseline_df <- assignment_df %>%
  filter(SPT == "Baseline", Response %in% c("Non_Responder", "Responder")) %>%
  mutate(
    response_num = ifelse(Response == "Responder", 1, 0),
    akk_log = log10(Akkermansia_rel_abundance + pseudocount),
    akk_z = as.numeric(scale(akk_log)),
    akk_tertile = factor(paste0("T", ntile(rank(Akkermansia_rel_abundance, ties.method = "first"), 3)), levels = c("T1", "T2", "T3"))
  )

fwrite(assignment_df, file.path(tables_dir, "enterosignature_sample_assignments.csv"))

akk_interaction_effects <- lapply(es_labels, function(es_name) {
  model_df <- baseline_df %>% mutate(es_z = as.numeric(scale(.data[[es_name]])))
  if (length(unique(model_df$response_num)) < 2 || sd(model_df$es_z, na.rm = TRUE) == 0 || sd(model_df$akk_z, na.rm = TRUE) == 0) {
    return(NULL)
  }
  fit <- safe_glm(response_num ~ es_z * akk_z + factor(project), model_df)
  term <- extract_glm_term(fit, "es_z:akk_z")
  if (is.null(term)) {
    return(NULL)
  }
  term %>% mutate(enterosignature = es_name, n = nrow(model_df), responders = sum(model_df$response_num), nonresponders = sum(model_df$response_num == 0))
}) %>% bind_rows() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

akk_stratified_effects <- lapply(es_labels, function(es_name) {
  bind_rows(lapply(levels(baseline_df$akk_tertile), function(level_name) {
    model_df <- baseline_df %>% filter(akk_tertile == level_name) %>% mutate(es_z = as.numeric(scale(.data[[es_name]])))
    if (nrow(model_df) < 20 || length(unique(model_df$response_num)) < 2 || sd(model_df$es_z, na.rm = TRUE) == 0) {
      return(NULL)
    }
    fit <- safe_glm(response_num ~ es_z + factor(project), model_df)
    term <- extract_glm_term(fit, "es_z")
    if (is.null(term)) {
      return(NULL)
    }
    term %>% mutate(enterosignature = es_name, akk_tertile = level_name, n = nrow(model_df), responders = sum(model_df$response_num), nonresponders = sum(model_df$response_num == 0))
  }))
}) %>% bind_rows() %>%
  group_by(akk_tertile) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(akk_tertile, p_value)

alpha_df <- assignment_df %>%
  mutate(
    shannon = vegan::diversity(t(genus_rel[, sample_id, drop = FALSE]), index = "shannon"),
    simpson = vegan::diversity(t(genus_rel[, sample_id, drop = FALSE]), index = "simpson"),
    observed_genera = colSums(genus_rel[, sample_id, drop = FALSE] > 0)
  )

alpha_long <- alpha_df %>%
  select(sample_id, project, SPT, Response, enterosignature, shannon, simpson, observed_genera) %>%
  pivot_longer(cols = c(shannon, simpson, observed_genera), names_to = "metric", values_to = "value")

alpha_tests <- alpha_long %>%
  group_by(metric) %>%
  summarise(
    kruskal_p = suppressWarnings(kruskal.test(value ~ enterosignature)$p.value),
    .groups = "drop"
  )

clr_input <- t(genus_rel[, assignment_df$sample_id, drop = FALSE])
clr_mat <- log(clr_input + pseudocount)
clr_mat <- sweep(clr_mat, 1, rowMeans(clr_mat), "-")
aitchison_dist <- dist(clr_mat, method = "euclidean")
pcoa <- cmdscale(aitchison_dist, k = 2, eig = TRUE)
pcoa_df <- data.frame(
  sample_id = rownames(clr_mat),
  Aitchison_PCo1 = pcoa$points[, 1],
  Aitchison_PCo2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
) %>% left_join(assignment_df, by = "sample_id")
pcoa_var <- pcoa$eig[pcoa$eig > 0] / sum(pcoa$eig[pcoa$eig > 0])

adonis_df <- assignment_df %>% select(sample_id, enterosignature, project, SPT)
rownames(adonis_df) <- adonis_df$sample_id
adonis_df <- adonis_df[labels(aitchison_dist), , drop = FALSE]
adonis_fit <- vegan::adonis2(aitchison_dist ~ enterosignature + project + SPT, data = adonis_df, permutations = 999)
adonis_table <- as.data.frame(adonis_fit)
adonis_table$term <- rownames(adonis_table)

centroids <- clr_mat %>%
  as.data.frame() %>%
  mutate(sample_id = rownames(clr_mat)) %>%
  left_join(assignment_df %>% select(sample_id, enterosignature), by = "sample_id") %>%
  group_by(enterosignature) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

centroid_dist <- lapply(seq_len(nrow(assignment_df)), function(idx) {
  sample_id <- assignment_df$sample_id[idx]
  es <- assignment_df$enterosignature[idx]
  centroid <- centroids %>% filter(enterosignature == es) %>% select(-enterosignature)
  data.frame(
    sample_id = sample_id,
    enterosignature = es,
    aitchison_distance_to_es_centroid = sqrt(sum((as.numeric(clr_mat[sample_id, ]) - as.numeric(centroid[1, ]))^2)),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

fwrite(response_counts, file.path(tables_dir, "baseline_enterosignature_response_counts.csv"))
fwrite(response_summary, file.path(tables_dir, "baseline_enterosignature_response_summary.csv"))
fwrite(response_test_df, file.path(tables_dir, "baseline_enterosignature_response_chisq.csv"))
fwrite(response_es_effects, file.path(tables_dir, "baseline_enterosignature_response_logistic_effects.csv"))
fwrite(akk_interaction_effects, file.path(tables_dir, "baseline_enterosignature_akk_response_interactions.csv"))
fwrite(akk_stratified_effects, file.path(tables_dir, "baseline_enterosignature_response_by_akk_tertile.csv"))
fwrite(alpha_df, file.path(tables_dir, "enterosignature_alpha_diversity_by_sample.csv"))
fwrite(alpha_tests, file.path(tables_dir, "enterosignature_alpha_diversity_tests.csv"))
fwrite(pcoa_df, file.path(tables_dir, "enterosignature_aitchison_pcoa.csv"))
fwrite(adonis_table, file.path(tables_dir, "enterosignature_aitchison_permanova.csv"))
fwrite(centroid_dist, file.path(tables_dir, "enterosignature_aitchison_centroid_distances.csv"))

top_composition <- composition_df %>%
  group_by(enterosignature) %>%
  slice_max(order_by = weight, n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(genus = reorder(genus, weight))

p_comp <- ggplot(top_composition, aes(x = weight, y = genus, fill = enterosignature)) +
  geom_col(width = 0.8) +
  facet_wrap(~ enterosignature, scales = "free_y", ncol = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Genus composition of enterosignatures", x = "NMF genus weight", y = NULL, fill = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none", plot.title = element_text(face = "bold"))

p_response <- ggplot(response_counts, aes(x = enterosignature, y = prop, fill = Response)) +
  geom_col(width = 0.75, color = "white") +
  geom_text(aes(label = paste0(n, "\n", percent(prop, accuracy = 1))), position = position_stack(vjust = 0.5), size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Non_Responder" = "#B7C9E2", "Responder" = "#D97A6A")) +
  labs(title = "Baseline response composition by enterosignature", subtitle = sprintf("Chi-square p = %.3g", response_test_df$p_value[1]), x = NULL, y = "Response composition", fill = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

p_effect <- ggplot(response_es_effects, aes(x = odds_ratio, y = reorder(enterosignature, odds_ratio))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = or_low, xmax = or_high), width = 0.18, color = "grey40") +
  geom_point(size = 2.8, color = "#2C3E50") +
  scale_x_log10() +
  labs(title = "Baseline response association of ES abundance", x = "Odds ratio per 1 SD ES increase", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p_interaction <- ggplot(akk_interaction_effects, aes(x = odds_ratio, y = reorder(enterosignature, odds_ratio))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbar(aes(xmin = or_low, xmax = or_high), width = 0.18, color = "grey40") +
  geom_point(size = 2.8, color = "#7B3F8C") +
  scale_x_log10() +
  labs(title = "AKK-dependent ES-response interaction", x = "Interaction odds ratio", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p_alpha <- ggplot(alpha_long, aes(x = enterosignature, y = value, fill = enterosignature)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.7, show.legend = FALSE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(title = "Alpha diversity by enterosignature", x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

p_pcoa <- ggplot(pcoa_df, aes(x = Aitchison_PCo1, y = Aitchison_PCo2, color = enterosignature, shape = SPT)) +
  geom_point(alpha = 0.75, size = 1.8) +
  labs(
    title = "Aitchison distance PCoA by enterosignature",
    subtitle = sprintf("PCo1 %.1f%% | PCo2 %.1f%%", 100 * pcoa_var[1], 100 * pcoa_var[2]),
    x = "PCo1",
    y = "PCo2",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

centroid_plot_df <- centroid_dist %>% left_join(assignment_df %>% select(sample_id, SPT, project), by = "sample_id")
p_centroid <- ggplot(centroid_plot_df, aes(x = enterosignature, y = aitchison_distance_to_es_centroid, fill = enterosignature)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.7, show.legend = FALSE) +
  labs(title = "Within-ES Aitchison distance to centroid", x = NULL, y = "Aitchison distance") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold"))

ggsave(file.path(figures_dir, "enterosignature_genus_composition.pdf"), p_comp, width = 12, height = 8)
ggsave(file.path(figures_dir, "enterosignature_genus_composition.png"), p_comp, width = 12, height = 8, dpi = 300)
ggsave(file.path(figures_dir, "baseline_enterosignature_response_composition.pdf"), p_response, width = 8, height = 5.5)
ggsave(file.path(figures_dir, "baseline_enterosignature_response_composition.png"), p_response, width = 8, height = 5.5, dpi = 300)
ggsave(file.path(figures_dir, "baseline_enterosignature_response_effects.pdf"), p_effect, width = 7, height = 4.8)
ggsave(file.path(figures_dir, "baseline_enterosignature_response_effects.png"), p_effect, width = 7, height = 4.8, dpi = 300)
ggsave(file.path(figures_dir, "baseline_enterosignature_akk_interactions.pdf"), p_interaction, width = 7, height = 4.8)
ggsave(file.path(figures_dir, "baseline_enterosignature_akk_interactions.png"), p_interaction, width = 7, height = 4.8, dpi = 300)
ggsave(file.path(figures_dir, "enterosignature_alpha_diversity.pdf"), p_alpha, width = 11, height = 5.5)
ggsave(file.path(figures_dir, "enterosignature_alpha_diversity.png"), p_alpha, width = 11, height = 5.5, dpi = 300)
ggsave(file.path(figures_dir, "enterosignature_aitchison_pcoa.pdf"), p_pcoa, width = 8, height = 6)
ggsave(file.path(figures_dir, "enterosignature_aitchison_pcoa.png"), p_pcoa, width = 8, height = 6, dpi = 300)
ggsave(file.path(figures_dir, "enterosignature_aitchison_centroid_distance.pdf"), p_centroid, width = 8, height = 5.5)
ggsave(file.path(figures_dir, "enterosignature_aitchison_centroid_distance.png"), p_centroid, width = 8, height = 5.5, dpi = 300)

summary_panel <- (p_comp / (p_response | p_effect) / (p_alpha | p_pcoa)) +
  plot_layout(heights = c(1.25, 1, 1)) +
  plot_annotation(title = "Genus-level enterosignatures across baseline and non-baseline metagenomes")
ggsave(file.path(figures_dir, "enterosignature_summary_panel.pdf"), summary_panel, width = 15, height = 18)
ggsave(file.path(figures_dir, "enterosignature_summary_panel.png"), summary_panel, width = 15, height = 18, dpi = 300)

top_effect <- if (nrow(response_es_effects) > 0) response_es_effects[1, ] else NULL
top_interaction <- if (nrow(akk_interaction_effects) > 0) akk_interaction_effects[1, ] else NULL
report_lines <- c(
  "# Genus-level enterosignature analysis",
  "",
  "## Method",
  sprintf("- Input: all Human WGS samples with bacterial read depth >= %d, including baseline and non-baseline samples.", depth_cutoff),
  sprintf("- Taxonomic level: bacterial genera aggregated from the phyloseq-like object; %d genera/features were used for NMF after prevalence/mean-abundance filtering.", nrow(nmf_input)),
  sprintf("- Enterosignature algorithm: local R implementation of KL-divergence NMF following the enterosignature-paper workflow; rank = %d, runs = %d, max iterations = %d. Best run selected by explained variance.", rank_k, nmf_runs, nmf_max_iter),
  sprintf("- Best NMF run: run %d, explained variance %.4f, cosine similarity %.4f.", nmf_res$run, nmf_res$explained_variance, nmf_res$cosine_similarity),
  "",
  "## Baseline response association",
  sprintf("- Baseline samples with evaluable response: %d.", nrow(baseline_df)),
  sprintf("- Overall enterosignature vs Response chi-square p = %.3g.", response_test_df$p_value[1]),
  if (!is.null(top_effect)) sprintf("- Strongest continuous ES-response association: %s, OR %.2f (95%% CI %.2f-%.2f), p = %.3g.", top_effect$enterosignature, top_effect$odds_ratio, top_effect$or_low, top_effect$or_high, top_effect$p_value) else "- No estimable continuous ES-response association.",
  "",
  "## AKK-dependent response-related enterosignatures",
  if (!is.null(top_interaction)) sprintf("- Strongest ES x Akkermansia interaction: %s, interaction OR %.2f (95%% CI %.2f-%.2f), p = %.3g.", top_interaction$enterosignature, top_interaction$odds_ratio, top_interaction$or_low, top_interaction$or_high, top_interaction$p_value) else "- No estimable ES x Akkermansia interaction.",
  "- Stratified ES-response models by Akkermansia tertile are exported as a table for follow-up.",
  "",
  "## Diversity and Aitchison structure",
  paste0("- Alpha diversity Kruskal-Wallis p-values: ", paste(sprintf("%s=%.3g", alpha_tests$metric, alpha_tests$kruskal_p), collapse = "; "), "."),
  sprintf("- Aitchison PERMANOVA term for enterosignature: R2 = %.3f, p = %.3g.", adonis_table$R2[adonis_table$term == "enterosignature"], adonis_table$`Pr(>F)`[adonis_table$term == "enterosignature"]),
  "",
  "## Key output files",
  "- tables/enterosignature_sample_assignments.csv",
  "- tables/enterosignature_genus_composition_long.csv",
  "- tables/baseline_enterosignature_response_summary.csv",
  "- tables/baseline_enterosignature_akk_response_interactions.csv",
  "- figures/enterosignature_summary_panel.pdf",
  "- figures/enterosignature_aitchison_pcoa.pdf",
  "- figures/enterosignature_alpha_diversity.pdf"
)
writeLines(report_lines, file.path(out_dir, "enterosignature_analysis_report.md"))

message("Enterosignature genus-level analysis completed: ", out_dir)