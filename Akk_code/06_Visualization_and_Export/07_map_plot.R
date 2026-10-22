# Supplementary Figure 1: all panels use the same 1,198 samples.
# 使用方法：解压完整数据包，在 RStudio 打开本文件并运行全部代码。
# 仅使用 R 自带的 graphics/grDevices/utils/stats，无需安装第三方包。
# 输入为 data/ 下的样本注释原表、补注规则和地图轮廓。
# 地图对所有国家使用同一种展示方式，无单独国家的引线或强调标签。

find_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grepl("^--file=", a)]
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    f <- frames[[i]]$ofile
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(f)) return(dirname(normalizePath(f)))
  }
  getwd()
}
BASE_DIR <- find_script_dir()
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) BASE_DIR <- normalizePath(args[1], mustWork = TRUE)
DATA_DIR <- file.path(BASE_DIR, "data")
OUT_DIR <- file.path(BASE_DIR, "r_output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
EXCLUDED_SAMPLES <- c("SRR16554661")
EXPECTED_N <- 1198L

read_data <- function(name) {
  path <- file.path(DATA_DIR, name)
  if (!file.exists(path)) stop("Missing input: ", path,
    "\nUnzip the complete package and keep the data folder beside this script.")
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA"), fileEncoding = "UTF-8-BOM")
}
write_data <- function(x, name) {
  write.csv(x, file.path(OUT_DIR, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
is_missing <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "Missing", "Unknown", "N/A")
as_number <- function(x) suppressWarnings(as.numeric(as.character(x)))

# 1. 从原始样本注释表重新筛选，而不是从已汇总的柱高反推数据。
raw <- read_data("wgs_metadata_plot_data.csv")
baseline <- !is.na(raw$AssayType) & !is.na(raw$SPT) &
  toupper(trimws(raw$AssayType)) == "WGS" & tolower(trimws(raw$SPT)) == "baseline"
known_response <- raw$response %in% c("Responder", "Non_Responder")
keep <- baseline & !(raw$sample_id %in% EXCLUDED_SAMPLES) & known_response
keep[is.na(keep)] <- FALSE
x <- raw[keep, , drop = FALSE]
stopifnot(nrow(x) == EXPECTED_N, !anyDuplicated(x$sample_id),
          !any(x$sample_id %in% EXCLUDED_SAMPLES))
N <- nrow(x)

# 补注规则单独保存在数据表，保留原始字段，不改写原表。
rules <- read_data("annotation_rules.csv")
country_plot <- as.character(x$country)
country_plot[!is.na(country_plot) & country_plot == "US"] <- "USA"
cancer_plot <- as.character(x$cancer_type)
country_source <- rep("Supplied metadata", N)
cancer_source <- rep("Supplied metadata", N)
for (i in seq_len(nrow(rules))) {
  same_cohort <- x$project == rules$cohort[i]
  same_cohort[is.na(same_cohort)] <- FALSE
  if (!is_missing(rules$fill_country[i])) {
    j <- same_cohort & is_missing(country_plot)
    country_plot[j] <- rules$fill_country[i]
    country_source[j] <- rules$source[i]
  }
  if (!is_missing(rules$fill_cancer_type[i])) {
    j <- same_cohort & is_missing(cancer_plot)
    cancer_plot[j] <- rules$fill_cancer_type[i]
    cancer_source[j] <- rules$source[i]
  }
}
sex_plot <- as.character(x$sex)
sex_plot[is_missing(sex_plot)] <- "Missing"
d <- data.frame(
  sample_id = x$sample_id, cohort = x$project, assay_type = x$AssayType,
  baseline_annotation = x$SPT, patient_id = x$Patient_ID,
  country_original = x$country, country_plot = country_plot,
  sex = sex_plot, cancer_type_original = x$cancer_type, cancer_type_plot = cancer_plot,
  response = ifelse(x$response == "Responder", "R", "NR"),
  age_raw = x$Host_age, age_years = as_number(x$age_years),
  age_group = x$age_group, bmi_raw = x$BMI, bmi_num = as_number(x$bmi_num),
  age_uses_interval_midpoint = is.na(as_number(x$Host_age)) & !is.na(as_number(x$age_years)),
  country_annotation_source = country_source, cancer_annotation_source = cancer_source,
  stringsAsFactors = FALSE
)
stopifnot(!any(is_missing(d$country_plot)), !any(is_missing(d$cancer_type_plot)))
write_data(d, "sFigure1_plot_samples_1198.csv")
ex <- raw[baseline & !keep, c("sample_id", "project", "SPT", "response"), drop = FALSE]
ex$exclusion_reason <- ifelse(ex$sample_id %in% EXCLUDED_SAMPLES,
  "User-requested exclusion", "Missing binary response; excluded from all panels")
write_data(ex, "excluded_baseline_samples.csv")

# 2. 所有面板都从同一个 d 生成；D/E 的直方图分母为 N。
countries <- c("France", "USA", "China", "United Kingdom", "Germany", "Netherlands", "Spain")
country_cols <- setNames(c("#716A92", "#6E9DB3", "#BC856E", "#A68C66", "#59958B", "#A77485", "#9AA975"), countries)
cancers <- c("Melanoma", "NSCLC", "B cell lymphoma", "RCC", "Gastrointestinal cancer")
cancer_labels <- c("Melanoma", "NSCLC", "B cell lymphoma", "RCC", "Gastrointestinal\ncancer")
cohort_n <- table(d$cohort)
cohorts <- names(cohort_n)[order(-as.integer(cohort_n), names(cohort_n))]
country_n <- table(factor(d$country_plot, levels = countries))
sex_n <- table(factor(d$cohort, levels = cohorts), factor(d$sex, levels = c("Male", "Female", "Missing")))
cancer_n <- table(factor(d$cancer_type_plot, levels = cancers))
response_n <- table(factor(d$cancer_type_plot, levels = cancers), factor(d$response, levels = c("R", "NR")))
ages <- d$age_years[is.finite(d$age_years)]
bmi <- d$bmi_num[is.finite(d$bmi_num)]
age_hist <- hist(ages, breaks = seq(15, 95, 5), right = FALSE, include.lowest = TRUE, plot = FALSE)
bmi_hist <- hist(bmi, breaks = seq(14, 60, 2), right = FALSE, include.lowest = TRUE, plot = FALSE)
stopifnot(sum(country_n) == N, sum(sex_n) == N, sum(cancer_n) == N,
          sum(response_n) == N, sum(age_hist$counts) == length(ages), sum(bmi_hist$counts) == length(bmi))
write_data(data.frame(country = countries, samples = as.integer(country_n), percent = as.integer(country_n) / N * 100), "panel_A_country.csv")
write_data(data.frame(cohort = cohorts, as.data.frame.matrix(sex_n), total = rowSums(sex_n)), "panel_B_sex.csv")
write_data(data.frame(cancer_type = cancers, samples = as.integer(cancer_n)), "panel_C_cancer.csv")
write_hist <- function(h, numeric_n, name) {
  result <- data.frame(category = "Numeric interval", bin_left = head(h$breaks, -1),
    bin_right = tail(h$breaks, -1), samples = h$counts, percent_of_all_samples = h$counts / N * 100,
    denominator = N)
  result <- rbind(result, data.frame(category = "No numeric value", bin_left = NA,
    bin_right = NA, samples = N - numeric_n, percent_of_all_samples = (N - numeric_n) / N * 100, denominator = N))
  stopifnot(sum(result$samples) == N, abs(sum(result$percent_of_all_samples) - 100) < 1e-8)
  write_data(result, name)
}
write_hist(age_hist, length(ages), "panel_D_age.csv")
write_hist(bmi_hist, length(bmi), "panel_E_BMI.csv")
write_data(data.frame(cancer_type = cancers, R = response_n[, "R"], NR = response_n[, "NR"],
  total = rowSums(response_n), R_percent = response_n[, "R"] / rowSums(response_n) * 100,
  NR_percent = response_n[, "NR"] / rowSums(response_n) * 100), "panel_F_response.csv")
write_data(data.frame(panel = LETTERS[1:6], included_samples = N), "panel_sample_counts.csv")

# 3. 地图轮廓随包提供，不联网下载；x/y 为参考图的绘图坐标。
world <- read_data("map_outline.csv")
stopifnot(all(c("polygon_id", "point_order", "x", "y", "country") %in% names(world)))
world <- world[order(world$polygon_id, world$point_order), ]
polys <- split(seq_len(nrow(world)), world$polygon_id)
ink <- "#28363E"; secondary <- "#6A7378"; blue <- "#507C93"; plum <- "#A26980"
pale <- "#E7D8DF"; missing_col <- "#D9DEE0"; grid_col <- "#E9ECEE"; land <- "#DEE3E3"
fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)

# 单画布坐标使六个面板的标题、图例和文字位置一致；图形均为矢量对象。
draw_figure <- function() {
  par(mar = rep(0, 4), family = "sans", xpd = NA, lend = "butt")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  label <- function(x, y, value, pt = 8, col = ink, adj = c(0, .5), font = 1, srt = 0) {
    text(x, y, value, cex = pt / 12, col = col, adj = adj, font = font, srt = srt)
  }
  line <- function(x, y, col = grid_col, lwd = .65, lty = 1) lines(x, y, col = col, lwd = lwd, lty = lty)
  box <- function(x, y, w, h, col, border = NA, lwd = .5) rect(x, y, x + w, y + h, col = col, border = border, lwd = lwd)
  heading <- function(letter, title, x, y) {
    label(x, y, letter, 12, font = 2); label(x + .025, y, title, 9.2, font = 2)
  }
  xaxis <- function(b, max_x, ticks, title) {
    line(c(b[1], b[1] + b[3]), rep(b[2], 2), col = "#758086", lwd = .55)
    for (v in ticks) {
      xp <- b[1] + b[3] * v / max_x
      line(rep(xp, 2), c(b[2] - .005, b[2]), col = "#758086", lwd = .5)
      label(xp, b[2] - .016, v, 7.2, adj = c(.5, .5))
    }
    label(b[1] + b[3] / 2, b[2] - .046, title, 8, adj = c(.5, .5))
  }
  label(.038, .967, "Supplementary Figure 1", 13, font = 2)
  label(.962, .967, "WGS baseline samples", 10, adj = c(1, .5))
  label(.962, .939, sprintf("%s samples across %d cohorts and %d countries", fmt(N), length(cohorts), length(countries)), 8, secondary, adj = c(1, .5))
  line(c(.038, .962), c(.920, .920))
  heading("A", "Geographic distribution", .038, .882)
  heading("B", "Sex by cohort", .373, .882)
  heading("C", "Cancer type", .696, .882)
  heading("D", "Age distribution", .038, .472)
  heading("E", "BMI distribution", .373, .472)
  heading("F", "Response by cancer type", .696, .472)

  # A: all countries have equal treatment; only the common legend carries labels.
  for (ix in polys) {
    w <- world[ix, ]; country <- w$country[1]
    fill <- if (!is.na(country) && country %in% countries) country_cols[[country]] else land
    polygon(.035 + (w$x - 14) / 196 * .287, .662 + (143 - w$y) / 113 * .203,
            col = fill, border = "white", lwd = .15)
  }
  for (i in seq_along(countries)) {
    xx <- .041 + ((i - 1) %/% 4) * .155; yy <- .642 - ((i - 1) %% 4) * .023
    box(xx, yy - .004, .009, .011, country_cols[[countries[i]]])
    label(xx + .014, yy, countries[i], 7)
    label(xx + .137, yy, fmt(as.integer(country_n[i])), 7, adj = c(1, .5))
  }
  label(.041, .544, paste0("n = ", fmt(N), "; country annotation complete"), 7, secondary)

  # B: within-cohort percentages, using the common retained population.
  b <- c(.469, .586, .179, .239)
  for (v in c(0, 25, 50, 75, 100)) line(rep(b[1] + v / 100 * b[3], 2), c(b[2], b[2] + b[4]))
  for (i in seq_along(cohorts)) {
    yy <- b[2] + b[4] * (1 - (i - .5) / length(cohorts)); hh <- .70 * b[4] / length(cohorts)
    totals <- sum(sex_n[i, ]); props <- sex_n[i, ] / totals
    start <- 0
    for (j in 1:3) { box(b[1] + start * b[3], yy - hh / 2, props[j] * b[3], hh, c(blue, plum, missing_col)[j]); start <- start + props[j] }
    label(b[1] - .006, yy, cohorts[i], 6.8, adj = c(1, .5))
    label(b[1] + b[3] + .008, yy, totals, 6.7)
  }
  label(b[1] + b[3] + .008, b[2] + b[4] + .010, "n", 7, secondary)
  for (j in 1:3) { xx <- .410 + (j - 1) * .064; box(xx, .847, .011, .011, c(blue, plum, missing_col)[j]); label(xx + .015, .852, c("Male", "Female", "Missing")[j], 7) }
  xaxis(b, 100, c(0, 25, 50, 75, 100), "Composition (%)")

  # C: cancer sample counts.
  b <- c(.814, .586, .145, .239)
  for (v in c(0, 250, 500)) line(rep(b[1] + v / 550 * b[3], 2), c(b[2], b[2] + b[4]))
  for (i in seq_along(cancers)) {
    yy <- b[2] + b[4] * (1 - (i - .5) / 5); hh <- .59 * b[4] / 5
    box(b[1], yy - hh / 2, as.integer(cancer_n[i]) / 550 * b[3], hh, plum)
    label(b[1] - .006, yy, cancer_labels[i], 7, adj = c(1, .5))
    label(b[1] + as.integer(cancer_n[i]) / 550 * b[3] + .006, yy, as.integer(cancer_n[i]), 7.3)
  }
  label(.718, .848, paste0("n = ", fmt(N)), 7.4, secondary)
  xaxis(b, 550, c(0, 250, 500), "Samples")

  # D/E: a numeric-value histogram plus an availability bar, both based on N.
  histogram_panel <- function(b, h, values, color, title, ticks) {
    numeric_n <- length(values); yy <- h$counts / N * 100
    ylim <- max(yy) * 1.12; xmax <- tail(h$breaks, 1); xmin <- h$breaks[1]
    sx <- function(v) b[1] + (v - xmin) / (xmax - xmin) * b[3]
    sy <- function(v) b[2] + v / ylim * b[4]
    for (i in seq_along(h$counts)) box(sx(h$breaks[i]), b[2], sx(h$breaks[i + 1]) - sx(h$breaks[i]), yy[i] / ylim * b[4], color, "white", .4)
    line(rep(sx(median(values)), 2), c(b[2], b[2] + b[4]), ink, .9, 2)
    line(rep(b[1], 2), c(b[2], b[2] + b[4]), "#758086", .55)
    line(c(b[1], b[1] + b[3]), rep(b[2], 2), "#758086", .55)
    for (v in pretty(c(0, ylim), n = 4)) if (v >= 0 && v <= ylim) {
      line(c(b[1] - .003, b[1]), rep(sy(v), 2), "#758086", .5)
      label(b[1] - .007, sy(v), v, 7.2, adj = c(1, .5))
    }
    for (v in ticks) {
      line(rep(sx(v), 2), c(b[2] - .005, b[2]), "#758086", .5)
      label(sx(v), b[2] - .016, v, 7.2, adj = c(.5, .5))
    }
    label(b[1] + b[3] / 2, b[2] - .046, title, 8, adj = c(.5, .5))
    label(b[1] - .040, b[2] + b[4] / 2, "All samples (%)", 8, adj = c(.5, .5), srt = 90)
    label(b[1] - .020, .435, sprintf("n = %s; median = %g", fmt(N), median(values)), 7.6, secondary)
    box(b[1], .412, b[3] * numeric_n / N, .011, color)
    box(b[1] + b[3] * numeric_n / N, .412, b[3] * (N - numeric_n) / N, .011, missing_col)
    label(b[1], .388, sprintf("Numeric: %s (%.1f%%)", fmt(numeric_n), numeric_n / N * 100), 6.7, color)
    label(b[1] + b[3], .388, paste0("No numeric value: ", fmt(N - numeric_n)), 6.7, secondary, adj = c(1, .5))
  }
  histogram_panel(c(.083, .157, .236, .205), age_hist, ages, plum, "Age (years)", c(20, 40, 60, 80))
  histogram_panel(c(.417, .157, .236, .205), bmi_hist, bmi, blue, "BMI (kg/m²)", c(20, 30, 40, 50, 60))

  # F: labels are counts; segment widths are within-cancer response proportions.
  b <- c(.814, .157, .145, .233)
  for (v in c(0, 25, 50, 75, 100)) line(rep(b[1] + v / 100 * b[3], 2), c(b[2], b[2] + b[4]))
  for (i in seq_along(cancers)) {
    yy <- b[2] + b[4] * (1 - (i - .5) / 5); hh <- .59 * b[4] / 5
    prop <- response_n[i, "R"] / sum(response_n[i, ])
    box(b[1], yy - hh / 2, b[3] * prop, hh, plum)
    box(b[1] + b[3] * prop, yy - hh / 2, b[3] * (1 - prop), hh, pale)
    label(b[1] - .006, yy, cancer_labels[i], 7, adj = c(1, .5))
    label(b[1] + b[3] * prop / 2, yy, response_n[i, "R"], 7, "white", adj = c(.5, .5))
    label(b[1] + b[3] * (prop + (1 - prop) / 2), yy, response_n[i, "NR"], 7, adj = c(.5, .5))
  }
  box(.711, .408, .011, .011, plum); label(.727, .413, "R", 7.3)
  box(.754, .408, .011, .011, pale); label(.770, .413, "NR", 7.3)
  label(.962, .435, paste0("n = ", fmt(N)), 7.6, secondary, adj = c(1, .5))
  xaxis(b, 100, c(0, 25, 50, 75, 100), "Composition (%)")
  line(c(.038, .962), c(.098, .098))
  label(.038, .075, sprintf("All panels use the same %s baseline samples with known response. Age and BMI percentages use all %s as the denominator.", fmt(N), fmt(N)), 7, secondary)
  label(.038, .053, sprintf("Age intervals use midpoints (%d samples). Gray availability segments denote no numeric value; dashed lines show medians.", sum(d$age_uses_interval_midpoint)), 7, secondary)
}

# 4. 导出 PDF、600 dpi PNG 和可编辑 SVG；设备支持取决于本机 R 编译配置。
pdf(file.path(OUT_DIR, "sFigure_1_WGS_baseline_R.pdf"), width = 10.5, height = 6.9,
    family = "Helvetica", pointsize = 12, useDingbats = FALSE)
tryCatch(draw_figure(), finally = dev.off())
if (capabilities("png")) {
  png(file.path(OUT_DIR, "sFigure_1_WGS_baseline_R.png"), width = 10.5, height = 6.9,
      units = "in", res = 600, pointsize = 12, bg = "white")
  tryCatch(draw_figure(), finally = dev.off())
} else warning("PNG device unavailable; the vector PDF was exported.")
if (capabilities("cairo")) {
  svg(file.path(OUT_DIR, "sFigure_1_WGS_baseline_R.svg"), width = 10.5, height = 6.9,
      pointsize = 12, family = "sans", bg = "white")
  tryCatch(draw_figure(), finally = dev.off())
} else warning("Cairo device unavailable; SVG was skipped. The vector PDF was exported.")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Finished: ", normalizePath(OUT_DIR), "\nIncluded samples: ", N,
        " | R: ", sum(d$response == "R"), " | NR: ", sum(d$response == "NR"))
