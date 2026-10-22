#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(splines)
})

set.seed(20260627)

base_dir <- ""./data/"投稿定图/HGMA"
out_dir <- file.path(base_dir, "nature_pdf_figures")
risk_dir <- file.path(base_dir, "AKK_risk_curve_analysis")

matrix_path <- file.path(out_dir, "top30_species_high_tertile_disease_matrix_summary.csv")
disease_path <- file.path(out_dir, "cohort_disease_type_distribution_22_projects.csv")

matrix_df <- fread(matrix_path) %>%
  mutate(BioProject = as.character(BioProject))
disease_df <- fread(disease_path) %>%
  mutate(BioProject = as.character(BioProject))

red_categories <- c("High disease higher", "Mixed significant")

project_order <- matrix_df %>%
  mutate(red_signal = category %in% red_categories) %>%
  group_by(BioProject) %>%
  summarise(
    red_signal_count = sum(red_signal, na.rm = TRUE),
    n_control = first(n_control),
    n_disease = first(n_disease),
    .groups = "drop"
  ) %>%
  arrange(red_signal_count, BioProject) %>%
  pull(BioProject)

project_counts <- matrix_df %>%
  distinct(BioProject, n_control, n_disease) %>%
  mutate(
    project_rank = match(BioProject, project_order),
    y = length(project_order) - project_rank + 1,
    project_label = sprintf("%s (C=%d, D=%d)", BioProject, n_control, n_disease)
  ) %>%
  arrange(project_rank)

species_order <- matrix_df %>%
  mutate(red_signal = category %in% red_categories) %>%
  group_by(column_rank, species_name) %>%
  summarise(red_signal_count = sum(red_signal, na.rm = TRUE), .groups = "drop") %>%
  arrange(red_signal_count, species_name == "Akkermansia muciniphila", column_rank) %>%
  mutate(
    plot_col = row_number(),
    species_label = gsub("\\[|\\]", "", species_name)
  )

disease_order <- c("Healthy", sort(setdiff(unique(disease_df$Disease_specific), "Healthy")))
distinct_cols <- c(
  "#56B4E9", "#D55E00", "#009E73", "#0072B2", "#CC79A7",
  "#E69F00", "#7F7F7F", "#882255", "#44AA99", "#AA4499",
  "#117733", "#88CCEE", "#DDCC77", "#CC6677", "#332288",
  "#A6761D", "#E7298A"
)
disease_cols <- setNames(distinct_cols[seq_along(disease_order)], disease_order)
disease_cols["Healthy"] <- "#56B4E9"

disease_plot_df <- disease_df %>%
  filter(BioProject %in% project_order) %>%
  mutate(
    Disease_specific = factor(Disease_specific, levels = disease_order),
    BioProject = factor(BioProject, levels = project_order)
  ) %>%
  left_join(project_counts %>% select(BioProject, y), by = "BioProject") %>%
  arrange(y, Disease_specific) %>%
  group_by(BioProject, y) %>%
  mutate(
    x_end = cumsum(n),
    x_start = x_end - n
  ) %>%
  ungroup()

grid_df <- tidyr::expand_grid(
  plot_col = species_order$plot_col,
  BioProject = project_order
) %>%
  left_join(project_counts %>% select(BioProject, y), by = "BioProject")

red_df <- matrix_df %>%
  left_join(species_order %>% select(column_rank, plot_col), by = "column_rank") %>%
  mutate(
    red_signal = category %in% red_categories,
    project_rank = match(BioProject, project_order),
    y_center = length(project_order) - project_rank + 1,
    x_center = plot_col
  ) %>%
  filter(red_signal)

disease_plot <- ggplot(disease_plot_df, aes(fill = Disease_specific)) +
  geom_rect(
    aes(xmin = x_start, xmax = x_end, ymin = y - 0.39, ymax = y + 0.39),
    color = "white",
    linewidth = 0.20
  ) +
  scale_y_continuous(
    breaks = project_counts$y,
    labels = rep("", nrow(project_counts)),
    expand = expansion(add = 0.5)
  ) +
  scale_fill_manual(values = disease_cols, drop = FALSE, name = NULL) +
  labs(title = "Cohort disease composition", x = "Participants", y = NULL) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
    axis.text.x = element_text(size = 6.1, color = "black"),
    axis.title.x = element_text(size = 6.8),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.30),
    legend.position = "right",
    legend.justification = "left",
    legend.text = element_text(size = 5.6),
    legend.key.size = unit(3, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    plot.margin = margin(5, 2, 5, 5)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 0)))

matrix_plot <- ggplot() +
  geom_tile(
    data = grid_df,
    aes(x = plot_col, y = y),
    fill = "white",
    color = "#D9D9D9",
    linewidth = 0.18,
    width = 1,
    height = 1
  ) +
  geom_tile(
    data = red_df,
    aes(x = x_center, y = y_center),
    fill = "#D73027",
    color = "white",
    linewidth = 0.10,
    width = 1,
    height = 1,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = species_order$plot_col,
    labels = species_order$species_label,
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    breaks = project_counts$y,
    labels = project_counts$project_label,
    expand = expansion(add = 0.5)
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(
    title = "Top 30 species: high-tertile disease enrichment",
    subtitle = "Red cell: disease proportion significantly higher in the high tertile (within-species BH FDR < 0.05)",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 5.4, face = "italic", color = "black"),
    axis.text.y = element_text(size = 5.8, color = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "grey30", hjust = 0),
    plot.margin = margin(5, 5, 5, 1)
  )

top_panel <- matrix_plot + disease_plot +
  plot_layout(widths = c(1.85, 0.78), guides = "collect") &
  theme(legend.position = "right")

stacked_composition <- fread(file.path(out_dir, "top30_species_control_disease_stacked_plot_data_no_others.csv")) %>%
  mutate(
    group = factor(group, levels = c("Control", "Disease")),
    top_rank = as.integer(top_rank),
    display_species_name = gsub("\\[|\\]", "", display_species_name)
  ) %>%
  arrange(top_rank)

stack_levels <- stacked_composition %>%
  distinct(top_rank, display_species_name) %>%
  arrange(top_rank) %>%
  pull(display_species_name)

stacked_composition <- stacked_composition %>%
  mutate(display_species_name = factor(display_species_name, levels = rev(stack_levels)))

advanced_species_palette <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#7F3C8D", "#11A579",
  "#3969AC", "#F2B701", "#E73F74", "#80BA5A", "#E68310", "#008695",
  "#CF1C90", "#F97B72", "#4B4B8F", "#A5AA99", "#A6761D", "#666666"
)
stack_cols <- setNames(
  advanced_species_palette[seq_along(stack_levels)],
  stack_levels
)
stack_cols <- stack_cols[levels(stacked_composition$display_species_name)]

stack_ymax <- stacked_composition %>%
  group_by(group) %>%
  summarise(total = sum(mean_composition_percent), .groups = "drop") %>%
  summarise(ymax = max(ceiling(total / 5) * 5, 5)) %>%
  pull(ymax)

stacked_plot <- ggplot(stacked_composition, aes(x = group, y = mean_composition_percent, fill = display_species_name)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.16) +
  scale_fill_manual(values = stack_cols, breaks = rev(stack_levels), name = NULL) +
  scale_y_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, stack_ymax), clip = "off") +
  labs(
    title = "Top species composition",
    subtitle = "Mean relative abundance in control and disease groups",
    x = NULL,
    y = "Mean relative composition (%)"
  ) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.32, color = "black"),
    axis.ticks = element_line(linewidth = 0.28, color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.title.y = element_text(size = 6.8),
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "grey30", hjust = 0),
    legend.position = "right",
    legend.text = element_text(size = 4.4, lineheight = 0.82),
    legend.key.height = unit(0.18, "cm"),
    legend.key.width = unit(0.26, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5, 2, 5, 4)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = FALSE))

selected_projects <- matrix_df %>%
  filter(
    species_name %in% c("Akkermansia muciniphila", "Escherichia coli"),
    category %in% red_categories
  ) %>%
  distinct(BioProject) %>%
  mutate(project_rank = match(BioProject, project_order)) %>%
  arrange(project_rank) %>%
  pull(BioProject)

if (length(selected_projects) != 11) {
  warning(sprintf("Expected 11 selected projects but found %d.", length(selected_projects)))
}

selected_source <- matrix_df %>%
  filter(
    species_name %in% c("Akkermansia muciniphila", "Escherichia coli"),
    category %in% red_categories
  ) %>%
  arrange(match(BioProject, selected_projects), species_name) %>%
  select(BioProject, species_name, category, significant_comparisons, min_q_within_species, n_control, n_disease)

fwrite(selected_source, file.path(out_dir, "AKK_Ecoli_significant_11_cohorts_source_data.csv"))

control_col <- "#4DBBD5"
disease_col <- "#E64B35"
unadjusted_col <- "#3C5488"
adjusted_col <- "#00A087"
gray_col <- "#4A4A4A"

mode_value <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

fit_glm_generic <- function(df, x_col, adjusted = FALSE, include_project = FALSE) {
  df <- df %>% mutate(disease_binary = as.integer(disease_binary))
  if (nrow(df) < 20 || length(unique(df$disease_binary)) < 2 || length(unique(df[[x_col]])) < 4) {
    return(list(fit = NULL, status = "not_enough_complete_cases", formula = NA_character_, terms = "none"))
  }

  terms <- character()
  x_term <- if (length(unique(df[[x_col]])) >= 5 && nrow(df) >= 30) {
    sprintf("bs(%s, df = 3, degree = 3)", x_col)
  } else {
    x_col
  }
  if (adjusted) {
    if ("Age" %in% names(df) && sum(!is.na(df$Age)) >= 20 && length(unique(df$Age[!is.na(df$Age)])) > 1) {
      terms <- c(terms, "Age")
    }
    if ("BMI" %in% names(df) && sum(!is.na(df$BMI)) >= 20 && length(unique(df$BMI[!is.na(df$BMI)])) > 1) {
      terms <- c(terms, "BMI")
    }
    if ("Gender" %in% names(df) && length(unique(df$Gender[!is.na(df$Gender)])) > 1) {
      terms <- c(terms, "Gender")
      df$Gender <- factor(df$Gender)
    }
    if (include_project && "BioProject" %in% names(df) && length(unique(df$BioProject)) > 1) {
      terms <- c(terms, "BioProject")
      df$BioProject <- factor(df$BioProject)
    }
  }
  formula_text <- paste("disease_binary ~", paste(c(x_term, terms), collapse = " + "))

  tryCatch({
    fit <- suppressWarnings(glm(as.formula(formula_text), data = df, family = binomial()))
    list(fit = fit, status = "ok", formula = formula_text, terms = ifelse(length(terms) == 0, "none", paste(terms, collapse = "+")))
  }, error = function(e) {
    fallback <- paste("disease_binary ~", paste(c(x_col, terms), collapse = " + "))
    tryCatch({
      fit <- suppressWarnings(glm(as.formula(fallback), data = df, family = binomial()))
      list(fit = fit, status = "ok_linear_fallback", formula = fallback, terms = ifelse(length(terms) == 0, "none", paste(terms, collapse = "+")))
    }, error = function(e2) {
      list(fit = NULL, status = paste0("failed: ", conditionMessage(e2)), formula = fallback, terms = "none")
    })
  })
}

make_grid_generic <- function(df, x_col, n = 120) {
  rng <- range(df[[x_col]], na.rm = TRUE)
  if (isTRUE(all.equal(rng[1], rng[2]))) rng <- rng + c(-0.01, 0.01)
  seq(rng[1], rng[2], length.out = n)
}

predict_curve_generic <- function(fit, grid_df, x_col) {
  pred <- predict(fit, newdata = grid_df, type = "link", se.fit = TRUE)
  out <- data.frame(
    risk = plogis(pred$fit),
    risk_lower = plogis(pred$fit - 1.96 * pred$se.fit),
    risk_upper = plogis(pred$fit + 1.96 * pred$se.fit)
  )
  out[[x_col]] <- grid_df[[x_col]]
  out
}

marginal_predict_curve_generic <- function(fit, grid, covariate_df, x_col) {
  rows <- lapply(grid, function(x) {
    newdata <- covariate_df
    newdata[[x_col]] <- x
    pred <- predict(fit, newdata = newdata, type = "link", se.fit = TRUE)
    out <- data.frame(
      risk = mean(plogis(pred$fit), na.rm = TRUE),
      risk_lower = mean(plogis(pred$fit - 1.96 * pred$se.fit), na.rm = TRUE),
      risk_upper = mean(plogis(pred$fit + 1.96 * pred$se.fit), na.rm = TRUE)
    )
    out[[x_col]] <- x
    out
  })
  bind_rows(rows)
}

theme_risk <- theme_classic(base_size = 6.3, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.25, color = "black"),
    axis.ticks = element_line(linewidth = 0.22, color = "black"),
    axis.text = element_text(color = "black", size = 4.9),
    axis.title = element_text(size = 6.1),
    strip.background = element_rect(fill = "grey95", color = "grey55", linewidth = 0.22),
    strip.text = element_text(size = 5.1, face = "bold", lineheight = 0.86),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 5.3),
    legend.key.width = unit(0.36, "cm"),
    legend.key.height = unit(0.20, "cm"),
    panel.spacing = unit(0.48, "lines"),
    plot.title = element_text(face = "bold", size = 8.2, hjust = 0.5),
    plot.margin = margin(3, 4, 3, 4)
  )

make_selected_risk_panel <- function(samples, pred, x_col, x_label, plot_title) {
  pooled_label <- "Overall selected"
  panel_levels <- c(selected_projects, pooled_label)

  samples_sel <- samples %>%
    filter(BioProject %in% selected_projects, type %in% c("Control", "Disease")) %>%
    mutate(
      BioProject = as.character(BioProject),
      type = factor(type, levels = c("Control", "Disease")),
      disease_binary = as.integer(disease_binary),
      Age = suppressWarnings(as.numeric(Age)),
      BMI = suppressWarnings(as.numeric(BMI))
    )

  pred_sel <- pred %>%
    filter(BioProject %in% selected_projects) %>%
    mutate(
      BioProject = as.character(BioProject),
      model = recode(
        as.character(model),
        "unadjusted" = "Unadjusted",
        "adjusted_age_sex_BMI" = "Adjusted age/sex/BMI"
      )
    ) %>%
    filter(model %in% c("Unadjusted", "Adjusted age/sex/BMI"))

  pooled <- samples_sel %>%
    mutate(BioProject = pooled_label)
  pooled_fit_df <- samples_sel %>%
    mutate(BioProject = factor(BioProject, levels = selected_projects))
  pooled_unadj <- fit_glm_generic(pooled_fit_df, x_col, adjusted = FALSE)
  pooled_adj_df <- pooled_fit_df %>%
    filter(complete.cases(Age, BMI, Gender, BioProject))
  pooled_adj <- fit_glm_generic(pooled_adj_df, x_col, adjusted = TRUE, include_project = TRUE)

  pooled_pred <- list()
  if (!is.null(pooled_unadj$fit)) {
    grid <- data.frame(x = make_grid_generic(pooled, x_col, 140))
    names(grid) <- x_col
    pooled_pred[[length(pooled_pred) + 1]] <- predict_curve_generic(pooled_unadj$fit, grid, x_col) %>%
      mutate(BioProject = pooled_label, model = "Unadjusted")
  }
  if (!is.null(pooled_adj$fit)) {
    covariates <- pooled_adj_df %>%
      select(any_of(c("Age", "BMI", "Gender", "BioProject")))
    if ("Gender" %in% names(covariates)) covariates$Gender <- factor(covariates$Gender)
    if ("BioProject" %in% names(covariates)) covariates$BioProject <- factor(covariates$BioProject, levels = selected_projects)
    term_labels <- if (!is.null(pooled_adj$fit)) attr(terms(pooled_adj$fit), "term.labels") else character()
    if ("Gender" %in% term_labels && "Gender" %in% names(covariates)) {
      covariates$Gender <- factor(covariates$Gender, levels = levels(model.frame(pooled_adj$fit)$Gender))
    }
    if ("BioProject" %in% term_labels && "BioProject" %in% names(covariates)) {
      covariates$BioProject <- factor(covariates$BioProject, levels = levels(model.frame(pooled_adj$fit)$BioProject))
    }
    pooled_pred[[length(pooled_pred) + 1]] <- marginal_predict_curve_generic(
      pooled_adj$fit,
      make_grid_generic(pooled_adj_df, x_col, 140),
      covariates,
      x_col
    ) %>%
      mutate(BioProject = pooled_label, model = "Adjusted age/sex/BMI")
  }

  pred_plot <- bind_rows(pred_sel, bind_rows(pooled_pred)) %>%
    mutate(
      BioProject = factor(BioProject, levels = panel_levels),
      model = factor(model, levels = c("Unadjusted", "Adjusted age/sex/BMI"))
    )

  point_plot <- bind_rows(samples_sel, pooled) %>%
    mutate(BioProject = factor(BioProject, levels = panel_levels)) %>%
    group_by(BioProject, type) %>%
    mutate(.rand = runif(n())) %>%
    arrange(.rand, .by_group = TRUE) %>%
    slice_head(n = 350) %>%
    ungroup() %>%
    select(-.rand) %>%
    mutate(y_jitter = pmin(pmax(disease_binary + rnorm(n(), 0, 0.018), -0.045), 1.045))

  counts <- bind_rows(samples_sel, pooled) %>%
    mutate(BioProject = factor(BioProject, levels = panel_levels)) %>%
    count(BioProject, type) %>%
    pivot_wider(names_from = type, values_from = n, values_fill = 0) %>%
    mutate(
      panel_label = sprintf("%s\nC=%d, D=%d", BioProject, Control, Disease)
    )
  panel_labels <- setNames(counts$panel_label, as.character(counts$BioProject))
  x_limits <- range(c(point_plot[[x_col]], pred_plot[[x_col]]), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.025
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.01
  x_limits <- x_limits + c(-x_pad, x_pad)
  y_limits <- c(-0.06, 1.06)
  square_ratio <- diff(x_limits) / diff(y_limits)

  ggplot() +
    geom_ribbon(
      data = pred_plot,
      aes(x = .data[[x_col]], ymin = risk_lower, ymax = risk_upper, fill = model),
      alpha = 0.10,
      linewidth = 0
    ) +
    geom_line(
      data = pred_plot,
      aes(x = .data[[x_col]], y = risk, color = model),
      linewidth = 0.34
    ) +
    facet_wrap(~ BioProject, ncol = 3, labeller = as_labeller(panel_labels)) +
    scale_color_manual(
      values = c(
        Control = control_col,
        Disease = disease_col,
        Unadjusted = unadjusted_col,
        `Adjusted age/sex/BMI` = adjusted_col
      ),
      breaks = c("Unadjusted", "Adjusted age/sex/BMI", "Control", "Disease")
    ) +
    scale_fill_manual(
      values = c(Unadjusted = unadjusted_col, `Adjusted age/sex/BMI` = adjusted_col),
      breaks = c("Unadjusted", "Adjusted age/sex/BMI")
    ) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits, breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    coord_fixed(ratio = square_ratio, clip = "off") +
    labs(title = plot_title, x = x_label, y = "Predicted disease risk") +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 1.25, linewidth = 0.55), nrow = 1),
      fill = "none"
    ) +
    theme_risk
}

akk_samples <- fread(file.path(risk_dir, "AKK_risk_curve_input_samples.csv"))
akk_pred <- fread(file.path(risk_dir, "AKK_risk_curve_predictions.csv"))
ecoli_samples <- fread(file.path(out_dir, "Ecoli_risk_curve_input_samples_22_projects.csv"))
ecoli_pred <- fread(file.path(out_dir, "Ecoli_risk_curve_predictions_22_projects.csv"))

akk_risk <- make_selected_risk_panel(
  akk_samples,
  akk_pred,
  x_col = "log10_AKK",
  x_label = "log10(AKK abundance + pseudo-count)",
  plot_title = "Akkermansia muciniphila risk curves: 11 significant cohorts + pooled"
)

ecoli_risk <- make_selected_risk_panel(
  ecoli_samples,
  ecoli_pred,
  x_col = "log10_Ecoli",
  x_label = "log10(E. coli abundance + pseudo-count)",
  plot_title = "Escherichia coli risk curves: 11 significant cohorts + pooled"
)

bottom_panel <- akk_risk + ecoli_risk + plot_layout(widths = c(1, 1), guides = "collect") &
  theme(legend.position = "top")

top_combined_panel <- stacked_plot + top_panel +
  plot_layout(widths = c(0.70, 1.55), guides = "keep")

final_plot <- top_combined_panel / bottom_panel +
  plot_layout(heights = c(0.72, 1.65), guides = "keep")

pdf_out <- file.path(out_dir, "Fig_cohort_disease_top30_red_matrix_AKK_Ecoli_selected_risk_combined_NatureStyle.pdf")
png_out <- file.path(out_dir, "Fig_cohort_disease_top30_red_matrix_AKK_Ecoli_selected_risk_combined_NatureStyle.png")

ggsave(pdf_out, final_plot, width = 17.4, height = 20.6, device = "pdf", useDingbats = FALSE)
ggsave(png_out, final_plot, width = 17.4, height = 20.6, dpi = 600)

message("Selected cohorts: ", paste(selected_projects, collapse = ", "))
message("Wrote: ", pdf_out)
message("Wrote: ", png_out)
message("Wrote: ", file.path(out_dir, "AKK_Ecoli_significant_11_cohorts_source_data.csv"))
