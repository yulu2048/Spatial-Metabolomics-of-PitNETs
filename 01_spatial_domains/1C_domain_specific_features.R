## ============================================================================
## 1C_domain_specific_features.R
## 基于 1A 和 1B 的双轨 domain-specific feature 分析。
##
## Part A：样本内 domain-associated features（主分析）
##   - 使用 1A 的独立 validation features；每个物理样本独立解释。
## Part B：跨样本 recurrent meta-domain-associated features（补充分析）
##   - 使用 1B 中未参与共同 domain 建模的 features；物理样本为重复单位。
## Part C：连接两套结果
##   - 用稳定核心像素的重叠描述每个样本内 D-domain 与共同 M-domain 的关系。
## Part D：严格验证和稳健性分析
##   - 样本内空间吻合指标、空间分块 bootstrap 和核心阈值敏感性。
##   - M1/M2/M3 物理样本层面的两两比较、专属性分类和样本 bootstrap。
##   - leave-one-sample/diagnosis-out 与低像素样本敏感性。
##
## 本脚本不重新聚类，不需要重跑 1A 或 1B。
## ============================================================================

required_packages <- c("data.table", "Matrix", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

## ==== 1) 路径和参数 ====
data_root <- "/home/gaga/huzhan"
input_dir <- file.path(data_root, "MSI_ALL_aligned_10samples_consensus_tuned225")
part1_dir <- file.path(input_dir,
  "intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225")
part2_dir <- file.path(input_dir, "cross_sample_common_domains_10samples_1B_tuned225")
out_dir <- file.path(input_dir, "domain_specific_features_10samples_1C_tuned225")

within_dir <- file.path(out_dir, "partA_within_sample")
cross_dir <- file.path(out_dir, "partB_cross_sample")
integration_dir <- file.path(out_dir, "partC_within_cross_integration")
robustness_dir <- file.path(out_dir, "partD_specificity_spatial_robustness")
within_spatial_dir <- file.path(robustness_dir, "within_sample_spatial_validation")
cross_pairwise_dir <- file.path(robustness_dir, "cross_sample_pairwise_specificity")
relation_dir <- file.path(robustness_dir, "feature_relation_networks")
within_map_dir <- file.path(within_dir, "top_feature_ion_maps_with_domain_boundary")
cross_map_dir <- file.path(cross_dir, "top_feature_ion_maps_with_domain_boundary")
for (d in c(out_dir, within_dir, cross_dir, integration_dir, robustness_dir,
            within_spatial_dir, cross_pairwise_dir, relation_dir,
            within_map_dir, cross_map_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

input_files <- c(
  aligned_logCPM = file.path(input_dir,
    "pixel_by_feature_matrix_aligned_logCPM.rds"),
  feature_mz = file.path(input_dir, "feature_mz_aligned_filtered.rds"),
  within_effect = file.path(part1_dir,
    "independent_validation_feature_effects_sample_summary.csv"),
  within_partition_effect = file.path(part1_dir,
    "export_partition_validation_feature_effects_supplement.csv"),
  within_model_summary = file.path(part1_dir,
    "within_sample_consensus_model_summary.csv"),
  within_pixels = file.path(part1_dir,
    "tumor_pixels_within_sample_consensus_domains.rds"),
  cross_recurrent = file.path(part2_dir,
    "recurrent_common_domain_features_physical_sample_level.csv"),
  cross_sample_effect = file.path(part2_dir,
    "physical_sample_common_domain_feature_effects_stable_core.csv"),
  nfgnh3_partition_consistency = file.path(part2_dir,
    "NF-GnH-3_export_partition_direction_consistency_supplement.csv"),
  cross_pixels = file.path(part2_dir,
    "tumor_pixels_with_common_meta_domains.rds")
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))

random_seed <- 20250907L

## Part D 计算参数。
spatial_neighbor_radius <- 1L
spatial_block_width <- 10L
spatial_bootstrap_repeats <- 200L
spatial_bootstrap_block_fraction <- 0.80
high_intensity_quantile <- 0.90
within_core_probability_thresholds <- c(0.80, 0.90, 0.95)
cross_margin_thresholds <- c(0.05, 0.10, 0.20)
sample_bootstrap_repeats <- 1000L
pairwise_min_median_effect <- 0.15
pairwise_min_direction_consistency <- 0.80
pairwise_min_prevalence_difference <- 0.01
pairwise_strong_effect_without_prevalence <- 0.35
pairwise_min_samples_evaluable <- 8L
pairwise_max_wilcoxon_FDR <- 0.05
pairwise_min_bootstrap_direction_fraction <- 0.95
shared_domain_similarity_effect <- 0.15
low_mz_background_threshold <- 200

## 样本内候选阈值。效应是 1A 的 residualized z mean difference。
within_min_effect <- 0.50
within_strong_effect_without_prevalence <- 0.75
within_min_prevalence_difference <- 0.05
within_min_pixels_each_group <- 50L
within_min_partition_direction_consistency <- 2 / 3
within_tier_A_effect <- 0.75
within_tier_A_prevalence <- 0.10

## 跨样本候选阈值。效应是各物理样本 logCPM difference 的中位数。
cross_min_samples_evaluable <- 8L
cross_min_direction_consistency <- 0.80
cross_tier_A_direction_consistency <- 0.90
cross_min_median_effect <- 0.15
cross_tier_A_median_effect <- 0.25
cross_min_prevalence_difference <- 0.01
## 新版共识峰在肿瘤像素中近乎普遍检出；当检出率差异饱和为 0 时，
## 允许达到预设主效应阈值的 feature 进入独立验证。
cross_strong_effect_without_prevalence <- 0.15
cross_min_nfgnh3_partition_consistency <- 2 / 3

## 质量关系去冗余。
ion_mode <- "negative"
## 共识峰中心比0.1 Da分配窗口更精确；质量关系使用更严格的0.02 Da，
## 避免把刚好落在0.05 Da边界的无关峰误当成加合物。
isotope_mass_tolerance_da <- 0.02
adduct_mass_tolerance_da <- 0.02
neutral_loss_mass_tolerance_da <- 0.02
within_redundancy_correlation <- 0.90
cross_redundancy_correlation <- 0.90
adduct_redundancy_correlation <- 0.95
collapse_isotope_relations <- TRUE
## 已确认负离子模式。只有质量差匹配且相关性 >= 0.95 才合并潜在加合物。
collapse_adduct_relations <- TRUE
## 中性丢失可能代表真实代谢物或源内碎片，只标记，不自动删除。
collapse_neutral_loss_relations <- FALSE

top_n_within_per_domain <- 30L
top_n_cross_per_domain <- 30L
make_ion_maps <- TRUE
within_ion_maps_per_domain <- 5L
cross_ion_maps_per_domain <- 10L
ion_map_upper_quantile <- 0.99
save_png_dpi <- 220

## ==== 2) 工具函数 ====
safe_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

mass_relation <- function(delta) {
  delta <- abs(as.numeric(delta))
  isotope_delta <- 1.003355 * 1:3
  isotope_name <- paste0("isotope_", 1:3, "x1.003355")
  ## 以下质量差均相对于负离子模式常见基准离子 [M-H]-。
  ## 只根据质量差不能证明加合物身份，后续还要求高相关性。
  negative_adduct_delta <- c(
    "[M+Na-2H]-_vs_[M-H]-" = 21.981943,
    "[M+K-2H]-_vs_[M-H]-" = 37.955882,
    "[M+Cl]-_vs_[M-H]-" = 35.976678,
    "[M+HCOO]-_vs_[M-H]-" = 46.005479,
    "[M+CH3COO]-_vs_[M-H]-" = 60.021129
  )
  neutral_loss_delta <- c(
    "H2O_loss" = 18.010565,
    "NH3_loss" = 17.026549,
    "CO2_loss" = 43.989829
  )
  di <- abs(delta - isotope_delta)
  da <- abs(delta - negative_adduct_delta)
  dl <- abs(delta - neutral_loss_delta)
  if (length(di) && min(di) <= isotope_mass_tolerance_da) {
    ii <- which.min(di)
    return(list(type = "isotope", name = isotope_name[ii],
      observed_delta_da = delta, expected_delta_da = isotope_delta[ii],
      mass_error_da = di[ii]))
  }
  if (length(da) && min(da) <= adduct_mass_tolerance_da) {
    ia <- which.min(da)
    return(list(type = "negative_mode_adduct",
      name = names(negative_adduct_delta)[ia],
      observed_delta_da = delta, expected_delta_da = negative_adduct_delta[ia],
      mass_error_da = da[ia]))
  }
  if (length(dl) && min(dl) <= neutral_loss_mass_tolerance_da) {
    il <- which.min(dl)
    return(list(type = "possible_neutral_loss",
      name = names(neutral_loss_delta)[il],
      observed_delta_da = delta, expected_delta_da = neutral_loss_delta[il],
      mass_error_da = dl[il]))
  }
  NULL
}

safe_cor <- function(x, y) {
  if (length(x) != length(y) || length(x) < 3L ||
      stats::sd(x, na.rm = TRUE) == 0 || stats::sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
}

## 在同一个 domain 内按 score 从高到低保留代表峰。
## 只有质量关系和表达/效应高相关同时成立时才视作冗余。
deduplicate_ranked_group <- function(dt, correlation_matrix,
                                     correlation_threshold,
                                     score_col = "candidate_score") {
  z <- copy(dt)
  setorderv(z, score_col, -1L, na.last = TRUE)
  z[, `:=`(
    nonredundant_representative = TRUE,
    representative_feature_index = as.integer(matrix_feature_index),
    representative_mz = as.numeric(mz),
    mass_relation_to_representative = NA_character_,
    mass_relation_type_to_representative = NA_character_,
    observed_delta_to_representative_da = NA_real_,
    expected_delta_to_representative_da = NA_real_,
    mass_relation_error_da = NA_real_,
    correlation_to_representative = NA_real_,
    possible_mass_relation_flag = FALSE
  )]
  kept_rows <- integer()
  for (i in seq_len(nrow(z))) {
    if (!length(kept_rows)) {
      kept_rows <- i
      next
    }
    candidate_relations <- lapply(kept_rows, function(j) mass_relation(z$mz[i] - z$mz[j]))
    related <- which(!vapply(candidate_relations, is.null, logical(1)))
    if (!length(related)) {
      kept_rows <- c(kept_rows, i)
      next
    }
    cors <- vapply(related, function(k) {
      j <- kept_rows[k]
      fi <- as.character(z$matrix_feature_index[i])
      fj <- as.character(z$matrix_feature_index[j])
      if (is.null(correlation_matrix) ||
          !fi %in% rownames(correlation_matrix) ||
          !fj %in% colnames(correlation_matrix)) return(NA_real_)
      correlation_matrix[fi, fj]
    }, numeric(1))
    best_local <- related[which.max(replace(cors, !is.finite(cors), -Inf))]
    j <- kept_rows[best_local]
    rel <- candidate_relations[[best_local]]
    best_cor <- cors[which(related == best_local)]
    z$possible_mass_relation_flag[i] <- TRUE
    z$mass_relation_to_representative[i] <- rel$name
    z$mass_relation_type_to_representative[i] <- rel$type
    z$observed_delta_to_representative_da[i] <- rel$observed_delta_da
    z$expected_delta_to_representative_da[i] <- rel$expected_delta_da
    z$mass_relation_error_da[i] <- rel$mass_error_da
    z$correlation_to_representative[i] <- best_cor
    collapse_relation <-
      (rel$type == "isotope" && collapse_isotope_relations) ||
      (rel$type == "negative_mode_adduct" && collapse_adduct_relations) ||
      (rel$type == "possible_neutral_loss" && collapse_neutral_loss_relations)
    required_correlation <- if (rel$type == "negative_mode_adduct") {
      adduct_redundancy_correlation
    } else {
      correlation_threshold
    }
    if (collapse_relation && is.finite(best_cor) &&
        best_cor >= required_correlation) {
      z$nonredundant_representative[i] <- FALSE
      z$representative_feature_index[i] <- z$matrix_feature_index[j]
      z$representative_mz[i] <- z$mz[j]
    } else {
      kept_rows <- c(kept_rows, i)
    }
  }
  z
}

exact_sign_test_p <- function(x) {
  x <- x[is.finite(x) & x != 0]
  n <- length(x)
  if (!n) return(1)
  k <- min(sum(x > 0), sum(x < 0))
  min(1, 2 * stats::pbinom(k, n, 0.5))
}

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
    error = function(e) NA_real_)
}

minimum_leave_one_out_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) return(NA_real_)
  min(vapply(seq_along(x), function(i) stats::median(x[-i]), numeric(1)))
}

normalize_ion_intensity <- function(x, upper_quantile = ion_map_upper_quantile) {
  x <- as.numeric(x)
  upper <- stats::quantile(x[x > 0], upper_quantile, na.rm = TRUE, names = FALSE)
  if (!is.finite(upper) || upper <= 0) return(rep(0, length(x)))
  pmin(x / upper, 1)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probability, na.rm = TRUE, names = FALSE))
}

auc_rank_statistic <- function(x_in, x_out) {
  x_in <- as.numeric(x_in); x_out <- as.numeric(x_out)
  x_in <- x_in[is.finite(x_in)]; x_out <- x_out[is.finite(x_out)]
  if (!length(x_in) || !length(x_out)) return(NA_real_)
  r <- rank(c(x_in, x_out), ties.method = "average")
  u <- sum(r[seq_along(x_in)]) - length(x_in) * (length(x_in) + 1) / 2
  u / (length(x_in) * length(x_out))
}

make_spatial_edges <- function(x, y, radius = spatial_neighbor_radius) {
  key <- paste(x, y, sep = "::")
  if (anyDuplicated(key)) stop("同一imzML来源内存在重复坐标。")
  offsets <- data.table::CJ(dx = -radius:radius, dy = -radius:radius)[
    !(dx == 0L & dy == 0L) & (dx > 0L | (dx == 0L & dy > 0L))]
  z <- lapply(seq_len(nrow(offsets)), function(i) {
    to <- match(paste(x + offsets$dx[i], y + offsets$dy[i], sep = "::"), key)
    from <- which(!is.na(to))
    if (!length(from)) return(NULL)
    cbind(from = from, to = to[from])
  })
  z <- Filter(Negate(is.null), z)
  if (!length(z)) return(matrix(integer(), ncol = 2L))
  ans <- do.call(rbind, z); storage.mode(ans) <- "integer"; ans
}

continuous_morans_i <- function(value, edges) {
  value <- as.numeric(value)
  if (!nrow(edges) || length(value) < 3L) return(NA_real_)
  z <- value - mean(value, na.rm = TRUE)
  z[!is.finite(z)] <- 0
  den <- sum(z ^ 2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  length(value) / nrow(edges) *
    sum(z[edges[, 1]] * z[edges[, 2]]) / den
}

make_block_ids <- function(dt, width = spatial_block_width) {
  paste(dt$source_id,
    floor((dt$x - ave(dt$x, dt$source_id, FUN = min)) / width),
    floor((dt$y - ave(dt$y, dt$source_id, FUN = min)) / width), sep = "::")
}

mean_and_prevalence_difference <- function(x, group_a, group_b) {
  xa <- x[group_a]; xb <- x[group_b]
  c(mean_difference = mean(xa) - mean(xb),
    prevalence_difference = mean(xa > 0) - mean(xb > 0))
}

summarize_signed_effect <- function(effect, prevalence = NULL) {
  effect <- as.numeric(effect)
  ok <- is.finite(effect)
  effect <- effect[ok]
  if (!length(effect)) return(data.table())
  if (!is.null(prevalence)) prevalence <- as.numeric(prevalence)[ok]
  loo <- if (length(effect) < 3L) NA_real_ else vapply(seq_along(effect),
    function(i) stats::median(effect[-i]), numeric(1))
  full_sign <- sign(stats::median(effect))
  data.table(
    n_physical_samples_evaluable = length(effect),
    median_effect = stats::median(effect),
    median_absolute_effect = stats::median(abs(effect)),
    median_prevalence_difference = if (is.null(prevalence)) NA_real_ else
      stats::median(prevalence, na.rm = TRUE),
    positive_sample_fraction = mean(effect > 0),
    negative_sample_fraction = mean(effect < 0),
    direction_consistency = max(mean(effect > 0), mean(effect < 0)),
    effect_Q1 = safe_quantile(effect, 0.25),
    effect_Q3 = safe_quantile(effect, 0.75),
    minimum_leave_one_sample_out_absolute_median = {
      if (all(!is.finite(loo))) NA_real_ else min(abs(loo), na.rm = TRUE)
    },
    minimum_leave_one_sample_out_signed_margin = {
      if (all(!is.finite(loo))) NA_real_ else min(full_sign * loo, na.rm = TRUE)
    },
    leave_one_sample_out_same_direction_fraction =
      mean(sign(loo) == full_sign, na.rm = TRUE),
    sign_test_p = exact_sign_test_p(effect),
    wilcoxon_signed_rank_p = safe_wilcox_p(effect)
  )
}

bootstrap_sample_median <- function(effect, repeats, seed) {
  effect <- effect[is.finite(effect)]
  if (length(effect) < 3L) return(c(NA_real_, NA_real_, NA_real_))
  set.seed(seed)
  boot <- replicate(repeats, stats::median(sample(effect,
    length(effect), replace = TRUE)))
  c(bootstrap_median_CI_low = safe_quantile(boot, 0.025),
    bootstrap_median_CI_high = safe_quantile(boot, 0.975),
    bootstrap_same_direction_fraction = mean(sign(boot) == sign(stats::median(effect))))
}

domain_boundary_flag <- function(dt, label_col, target_label) {
  ans <- logical(nrow(dt))
  for (sid in unique(as.character(dt$source_id))) {
    idx <- which(as.character(dt$source_id) == sid)
    edge <- make_spatial_edges(dt$x[idx], dt$y[idx])
    if (!nrow(edge)) next
    lab <- as.character(dt[[label_col]][idx])
    target <- lab == target_label
    crossing <- target[edge[, 1]] != target[edge[, 2]]
    if (!any(crossing)) next
    nodes <- unique(c(edge[crossing, 1], edge[crossing, 2]))
    nodes <- nodes[target[nodes]]
    ans[idx[nodes]] <- TRUE
  }
  ans
}

plot_faceted_ion_map <- function(dt, title, subtitle) {
  p <- ggplot(dt, aes(x = x, y = y, fill = normalized_intensity)) +
    geom_tile() + scale_y_reverse() +
    facet_wrap(~source_id, scales = "free") +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = "x (pixel)", y = "y (pixel)",
      fill = "Relative intensity") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), panel.spacing = grid::unit(1, "lines"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"))
  if ("domain_boundary" %in% names(dt) && any(dt$domain_boundary)) {
    p <- p + geom_point(data = dt[domain_boundary == TRUE],
      aes(x = x, y = y), inherit.aes = FALSE,
      shape = 15, size = 0.35, color = "#00BFC4")
  }
  p
}

## ==== 3) 读入结果 ====
message("[读入] 1A、1B结果和对齐logCPM矩阵。")
feature_mz <- as.numeric(readRDS(input_files[["feature_mz"]]))
aligned_logCPM <- readRDS(input_files[["aligned_logCPM"]])
if (!inherits(aligned_logCPM, "dgCMatrix")) {
  aligned_logCPM <- as(aligned_logCPM, "dgCMatrix")
}
if (ncol(aligned_logCPM) != length(feature_mz)) stop("feature m/z 与矩阵列数不一致。")

within_effect <- fread(input_files[["within_effect"]])
within_partition_effect <- fread(input_files[["within_partition_effect"]])
within_model_summary <- fread(input_files[["within_model_summary"]])
within_pixels <- as.data.table(readRDS(input_files[["within_pixels"]]))
cross_recurrent_input <- fread(input_files[["cross_recurrent"]])
cross_sample_effect <- fread(input_files[["cross_sample_effect"]])
nfgnh3_consistency <- fread(input_files[["nfgnh3_partition_consistency"]])
cross_pixels <- as.data.table(readRDS(input_files[["cross_pixels"]]))

required_within <- c("sample_id", "domain", "matrix_feature_index", "mz",
  "combined_residualized_z_difference", "combined_prevalence_difference",
  "n_pixels_in", "n_pixels_out", "positive_export_partition_fraction")
required_cross <- c("common_meta_domain", "matrix_feature_index", "mz",
  "passed_strict_technical_filter", "selected_for_common_domain_model")
required_within_pixels <- c("sample_id", "source_id", "x", "y",
  "original_pixel_row", "consensus_domain", "consensus_probability",
  "stable_core_pixel")
required_cross_pixels <- c("sample_id", "diagnosis", "source_id", "x", "y",
  "original_pixel_row", "common_meta_domain", "domain_assignment_margin",
  "stable_core_pixel")
if (!all(required_within %in% names(within_effect))) stop("1A effect表字段不完整。")
if (!all(required_cross %in% names(cross_recurrent_input))) stop("1B recurrent表字段不完整。")
if (!all(required_within_pixels %in% names(within_pixels))) {
  stop("1A像素表缺少Part D需要的字段：",
    paste(setdiff(required_within_pixels, names(within_pixels)), collapse = ", "))
}
if (!all(required_cross_pixels %in% names(cross_pixels))) {
  stop("1B像素表缺少Part D需要的字段：",
    paste(setdiff(required_cross_pixels, names(cross_pixels)), collapse = ", "))
}

## ==== 4) Part A：样本内独立验证 features ====
message("[Part A] 筛选每个物理样本内部的独立验证feature。")
within_candidates <- copy(within_effect)
within_candidates <- within_candidates[
  feature_role == "independent_validation" &
    combined_residualized_z_difference >= within_min_effect &
    n_pixels_in >= within_min_pixels_each_group &
    n_pixels_out >= within_min_pixels_each_group &
    (combined_prevalence_difference >= within_min_prevalence_difference |
       combined_residualized_z_difference >= within_strong_effect_without_prevalence) &
    positive_export_partition_fraction >= within_min_partition_direction_consistency
]
within_candidates[, evidence_tier := fifelse(
  combined_residualized_z_difference >= within_tier_A_effect &
    (combined_prevalence_difference >= within_tier_A_prevalence |
       combined_residualized_z_difference >= 1),
  "Tier_A", "Tier_B"
)]
within_candidates[, candidate_score :=
  combined_residualized_z_difference *
    (1 + pmax(combined_prevalence_difference, 0)) *
    positive_export_partition_fraction]
within_candidates <- merge(
  within_candidates,
  within_model_summary[, .(sample_id, selected_k, structure_type,
    stable_core_fraction, spatial_organization_type,
    spatial_organization_label)],
  by = "sample_id", all.x = TRUE, sort = FALSE
)

## 每个样本计算候选feature在1A稳定核心像素中的相关矩阵，再进行质量关系去冗余。
within_deduplicated_list <- list(); wi <- 0L
for (sample_i in unique(as.character(within_candidates$sample_id))) {
  candidate_i <- within_candidates[sample_id == sample_i]
  feature_i <- sort(unique(candidate_i$matrix_feature_index))
  pixel_rows_i <- within_pixels[
    sample_id == sample_i & stable_core_pixel == TRUE,
    as.integer(original_pixel_row)
  ]
  x_i <- as.matrix(aligned_logCPM[pixel_rows_i, feature_i, drop = FALSE])
  cor_i <- stats::cor(x_i, use = "pairwise.complete.obs")
  rownames(cor_i) <- colnames(cor_i) <- as.character(feature_i)
  for (domain_i in unique(as.character(candidate_i$domain))) {
    wi <- wi + 1L
    within_deduplicated_list[[wi]] <- deduplicate_ranked_group(
      candidate_i[domain == domain_i], cor_i,
      within_redundancy_correlation, "candidate_score"
    )
  }
  rm(x_i, cor_i); invisible(gc())
}
within_all <- rbindlist(within_deduplicated_list, use.names = TRUE, fill = TRUE)
within_all[, `:=`(
  low_mz_background_risk = mz < low_mz_background_threshold,
  requires_matrix_blank_review = mz < low_mz_background_threshold
)]
setorder(within_all, sample_id, domain, -candidate_score)
within_all[, rank_before_redundancy := seq_len(.N), by = .(sample_id, domain)]
within_all[nonredundant_representative == TRUE,
  rank_nonredundant := seq_len(.N), by = .(sample_id, domain)]

within_top <- within_all[
  nonredundant_representative == TRUE & rank_nonredundant <= top_n_within_per_domain
]
fwrite(within_all,
  file.path(within_dir, "within_sample_validation_candidates_all_with_redundancy_flags.csv"))
fwrite(within_top,
  file.path(within_dir, "within_sample_top_nonredundant_validation_features.csv"))

within_counts <- within_all[, .(
  n_candidates_before_redundancy = .N,
  n_nonredundant_candidates = sum(nonredundant_representative),
  n_Tier_A_nonredundant = sum(nonredundant_representative & evidence_tier == "Tier_A"),
  n_Tier_B_nonredundant = sum(nonredundant_representative & evidence_tier == "Tier_B")
), by = .(sample_id, diagnosis, domain, selected_k,
          spatial_organization_type, spatial_organization_label)]
fwrite(within_counts,
  file.path(within_dir, "within_sample_candidate_counts_by_domain.csv"))

p_within_counts <- ggplot(within_counts,
    aes(x = interaction(sample_id, domain, sep = " / "),
        y = n_nonredundant_candidates, fill = spatial_organization_label)) +
  geom_col() + coord_flip() +
  labs(title = "Within-tumor domain-associated validation features",
    x = NULL, y = "Nonredundant candidate features", fill = "Spatial organization") +
  theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(within_dir, "within_sample_candidate_counts.png"),
  p_within_counts, width = 9, height = 7.5, dpi = save_png_dpi, bg = "white")

## ==== 5) Part B：跨样本独立验证 features ====
message("[Part B] 以10个物理样本为重复单位筛选共同meta-domain feature。")
cross_stats <- cross_sample_effect[, {
  eff <- mean_logCPM_difference
  prev <- prevalence_difference
  list(
    n_physical_samples_evaluable = uniqueN(sample_id),
    median_logCPM_difference = median(eff, na.rm = TRUE),
    median_prevalence_difference = median(prev, na.rm = TRUE),
    positive_sample_fraction = mean(eff > 0, na.rm = TRUE),
    negative_sample_fraction = mean(eff < 0, na.rm = TRUE),
    effect_Q1 = stats::quantile(eff, 0.25, na.rm = TRUE),
    effect_Q3 = stats::quantile(eff, 0.75, na.rm = TRUE),
    minimum_leave_one_sample_out_median = minimum_leave_one_out_median(eff),
    sign_test_p = exact_sign_test_p(eff),
    wilcoxon_signed_rank_p = safe_wilcox_p(eff)
  )
}, by = .(common_meta_domain, matrix_feature_index, mz)]
cross_stats[, direction_consistency := pmax(
  positive_sample_fraction, negative_sample_fraction
)]
cross_stats[, `:=`(
  sign_test_FDR = p.adjust(sign_test_p, method = "BH"),
  wilcoxon_signed_rank_FDR = p.adjust(wilcoxon_signed_rank_p, method = "BH")
), by = common_meta_domain]

cross_stats <- merge(
  cross_stats,
  unique(cross_recurrent_input[, .(
    matrix_feature_index, common_meta_domain,
    passed_strict_technical_filter,
    selected_for_common_domain_model
  )]),
  by = c("matrix_feature_index", "common_meta_domain"), all.x = TRUE, sort = FALSE
)
cross_stats <- merge(
  cross_stats, nfgnh3_consistency,
  by = c("common_meta_domain", "matrix_feature_index", "mz"),
  all.x = TRUE, sort = FALSE
)
cross_stats[is.na(n_partitions_evaluable), `:=`(
  n_partitions_evaluable = 0L,
  partition_direction_consistency = NA_real_
)]

cross_candidates <- cross_stats[
  passed_strict_technical_filter == TRUE &
    selected_for_common_domain_model == FALSE &
    n_physical_samples_evaluable >= cross_min_samples_evaluable &
    direction_consistency >= cross_min_direction_consistency &
    median_logCPM_difference >= cross_min_median_effect &
    minimum_leave_one_sample_out_median > 0 &
    (median_prevalence_difference >= cross_min_prevalence_difference |
       median_logCPM_difference >= cross_strong_effect_without_prevalence) &
    (is.na(partition_direction_consistency) |
       partition_direction_consistency >= cross_min_nfgnh3_partition_consistency)
]
cross_candidates[, evidence_tier := fifelse(
  direction_consistency >= cross_tier_A_direction_consistency &
    median_logCPM_difference >= cross_tier_A_median_effect &
    partition_direction_consistency >= cross_min_nfgnh3_partition_consistency,
  "Tier_A", "Tier_B"
)]
cross_candidates[, candidate_score :=
  median_logCPM_difference * direction_consistency *
    sqrt(n_physical_samples_evaluable) *
    (1 + pmax(median_prevalence_difference, 0))]

## 用10个物理样本中的feature效应相关性辅助质量关系去冗余。
cross_deduplicated_list <- list(); ci <- 0L
for (domain_i in unique(as.character(cross_candidates$common_meta_domain))) {
  candidate_i <- cross_candidates[common_meta_domain == domain_i]
  effect_i <- cross_sample_effect[
    common_meta_domain == domain_i &
      matrix_feature_index %in% candidate_i$matrix_feature_index,
    .(sample_id, matrix_feature_index, mean_logCPM_difference)
  ]
  effect_wide <- dcast(effect_i, sample_id ~ matrix_feature_index,
    value.var = "mean_logCPM_difference")
  feature_cols <- setdiff(names(effect_wide), "sample_id")
  cor_i <- stats::cor(as.matrix(effect_wide[, ..feature_cols]),
    use = "pairwise.complete.obs")
  ci <- ci + 1L
  cross_deduplicated_list[[ci]] <- deduplicate_ranked_group(
    candidate_i, cor_i, cross_redundancy_correlation, "candidate_score"
  )
}
cross_all <- rbindlist(cross_deduplicated_list, use.names = TRUE, fill = TRUE)
cross_all[, `:=`(
  low_mz_background_risk = mz < low_mz_background_threshold,
  requires_matrix_blank_review = mz < low_mz_background_threshold
)]
setorder(cross_all, common_meta_domain, -candidate_score)
cross_all[, rank_before_redundancy := seq_len(.N), by = common_meta_domain]
cross_all[nonredundant_representative == TRUE,
  rank_nonredundant := seq_len(.N), by = common_meta_domain]
cross_top <- cross_all[
  nonredundant_representative == TRUE & rank_nonredundant <= top_n_cross_per_domain
]

fwrite(cross_stats,
  file.path(cross_dir, "cross_sample_all_validation_feature_statistics.csv"))
fwrite(cross_all,
  file.path(cross_dir, "cross_sample_recurrent_validation_candidates_with_redundancy_flags.csv"))
fwrite(cross_top,
  file.path(cross_dir, "cross_sample_top_nonredundant_recurrent_validation_features.csv"))

cross_counts <- cross_all[, .(
  n_candidates_before_redundancy = .N,
  n_nonredundant_candidates = sum(nonredundant_representative),
  n_Tier_A_nonredundant = sum(nonredundant_representative & evidence_tier == "Tier_A"),
  n_Tier_B_nonredundant = sum(nonredundant_representative & evidence_tier == "Tier_B")
), by = common_meta_domain]
fwrite(cross_counts,
  file.path(cross_dir, "cross_sample_candidate_counts_by_meta_domain.csv"))

p_cross_counts <- ggplot(cross_counts,
    aes(x = common_meta_domain, y = n_nonredundant_candidates,
        fill = common_meta_domain)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c(M1 = "#1B9E77", M2 = "#D95F02", M3 = "#7570B3")) +
  labs(title = "Recurrent meta-domain-associated validation features",
    x = "Common meta-domain", y = "Nonredundant candidate features") +
  theme_minimal(base_size = 11) + theme(panel.grid.major.x = element_blank())
ggsave(file.path(cross_dir, "cross_sample_candidate_counts.png"),
  p_cross_counts, width = 6.5, height = 5, dpi = save_png_dpi, bg = "white")

## 物理样本效应热图。
cross_heat_features <- cross_top[rank_nonredundant <= 15L,
  .(common_meta_domain, matrix_feature_index, mz, rank_nonredundant)]
cross_heat <- merge(
  cross_sample_effect,
  cross_heat_features,
  by = c("common_meta_domain", "matrix_feature_index", "mz"),
  all = FALSE
)
cross_heat[, feature_label := paste0(
  common_meta_domain, " | m/z ", sprintf("%.4f", mz)
)]
cross_heat[, feature_label := factor(feature_label,
  levels = rev(unique(feature_label[order(common_meta_domain, rank_nonredundant)])))]
p_cross_heat <- ggplot(cross_heat,
    aes(x = sample_id, y = feature_label, fill = mean_logCPM_difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0) +
  labs(title = "Physical-sample effects of recurrent validation features",
    x = NULL, y = NULL, fill = "logCPM difference") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(cross_dir, "cross_sample_top_feature_effect_heatmap.png"),
  p_cross_heat, width = 10, height = 11, dpi = save_png_dpi, bg = "white")

## ==== 6) Part C：样本内D-domain与共同M-domain的稳定核心重叠 ====
message("[Part C] 计算1A D-domain与1B M-domain的稳定核心重叠。")
within_labels <- within_pixels[, .(
  original_pixel_row = as.integer(original_pixel_row),
  sample_id = as.character(sample_id),
  within_domain = as.character(consensus_domain),
  within_stable_core = as.logical(stable_core_pixel)
)]
cross_labels <- cross_pixels[, .(
  original_pixel_row = as.integer(original_pixel_row),
  common_meta_domain = as.character(common_meta_domain),
  cross_stable_core = as.logical(stable_core_pixel)
)]
overlap_pixels <- merge(within_labels, cross_labels,
  by = "original_pixel_row", all = FALSE, sort = FALSE)
overlap_pixels <- overlap_pixels[within_stable_core & cross_stable_core]

overlap_table <- overlap_pixels[, .(intersection_pixels = .N),
  by = .(sample_id, within_domain, common_meta_domain)]
within_size <- overlap_pixels[, .(within_domain_pixels = .N),
  by = .(sample_id, within_domain)]
common_size <- overlap_pixels[, .(common_domain_pixels = .N),
  by = .(sample_id, common_meta_domain)]
overlap_table <- merge(overlap_table, within_size,
  by = c("sample_id", "within_domain"), all.x = TRUE)
overlap_table <- merge(overlap_table, common_size,
  by = c("sample_id", "common_meta_domain"), all.x = TRUE)
overlap_table[, `:=`(
  fraction_of_within_domain = intersection_pixels / within_domain_pixels,
  fraction_of_common_domain = intersection_pixels / common_domain_pixels,
  jaccard_index = intersection_pixels /
    (within_domain_pixels + common_domain_pixels - intersection_pixels)
)]
setorder(overlap_table, sample_id, within_domain, -fraction_of_within_domain)
overlap_table[, best_matching_common_domain :=
  common_meta_domain[which.max(fraction_of_within_domain)],
  by = .(sample_id, within_domain)]
overlap_table[, is_best_match := common_meta_domain == best_matching_common_domain]
fwrite(overlap_table,
  file.path(integration_dir, "within_D_to_common_M_stable_core_overlap.csv"))

best_mapping <- overlap_table[is_best_match == TRUE]
fwrite(best_mapping,
  file.path(integration_dir, "best_common_meta_domain_for_each_within_domain.csv"))

p_overlap <- ggplot(overlap_table,
    aes(x = common_meta_domain, y = within_domain,
        fill = fraction_of_within_domain)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", fraction_of_within_domain)), size = 2.7) +
  facet_wrap(~sample_id, scales = "free_y") +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Overlap between within-tumor and common meta-domains",
    subtitle = "Fraction of each within-tumor stable-core domain assigned to M1-M3",
    x = "Common meta-domain", y = "Within-tumor domain", fill = "Fraction") +
  theme_minimal(base_size = 9) + theme(panel.grid = element_blank())
ggsave(file.path(integration_dir, "within_D_to_common_M_overlap_heatmap.png"),
  p_overlap, width = 11, height = 8, dpi = save_png_dpi, bg = "white")

## feature层面的交集只说明同一m/z在两套分析中都出现；D标签不能直接跨样本合并。
feature_overlap <- merge(
  within_top[, .(sample_id, within_domain = domain,
    matrix_feature_index, mz, within_effect = combined_residualized_z_difference,
    within_rank = rank_nonredundant)],
  cross_top[, .(common_meta_domain, matrix_feature_index, mz,
    cross_median_effect = median_logCPM_difference,
    cross_direction_consistency = direction_consistency,
    cross_rank = rank_nonredundant)],
  by = c("matrix_feature_index", "mz"), all = FALSE
)
feature_overlap <- merge(
  feature_overlap,
  best_mapping[, .(sample_id, within_domain,
    best_matching_common_domain, fraction_of_within_domain, jaccard_index)],
  by = c("sample_id", "within_domain"), all.x = TRUE
)
feature_overlap[, domain_mapping_concordant :=
  common_meta_domain == best_matching_common_domain]
fwrite(feature_overlap,
  file.path(integration_dir, "top_feature_overlap_between_1A_and_1B.csv"))

## ==== 7) Part D1：样本内候选的空间吻合与分块稳健性 ====
message("[Part D1] 计算样本内候选feature的AUC、空间自相关、空间富集和分块bootstrap。")
within_spatial_list <- list(); within_boot_list <- list(); within_sensitivity_list <- list()
wsi <- wbi <- wti <- 0L
within_nr <- within_all[nonredundant_representative == TRUE]

for (sample_i in unique(as.character(within_nr$sample_id))) {
  meta_i <- copy(within_pixels[sample_id == sample_i & !is.na(consensus_domain)])
  meta_i[, block_id_1C := make_block_ids(.SD, spatial_block_width)]
  candidate_i <- within_nr[sample_id == sample_i]
  feature_i <- sort(unique(as.integer(candidate_i$matrix_feature_index)))
  x_all_i <- aligned_logCPM[as.integer(meta_i$original_pixel_row), feature_i, drop = FALSE]
  core_i <- which(meta_i$stable_core_pixel == TRUE)
  x_core_i <- x_all_i[core_i, , drop = FALSE]
  meta_core_i <- meta_i[core_i]

  ## 每个feature的Moran's I按导出区域独立计算，再按像素数加权汇总。
  moran_source <- list(); msi <- 0L
  for (sid in unique(as.character(meta_i$source_id))) {
    lr <- which(as.character(meta_i$source_id) == sid)
    edge <- make_spatial_edges(meta_i$x[lr], meta_i$y[lr])
    if (!nrow(edge)) next
    x_s <- as.matrix(x_all_i[lr, , drop = FALSE])
    msi <- msi + 1L
    moran_source[[msi]] <- data.table(
      matrix_feature_index = feature_i,
      source_id = sid,
      source_n_pixels = length(lr),
      feature_morans_i = vapply(seq_along(feature_i),
        function(j) continuous_morans_i(x_s[, j], edge), numeric(1))
    )
  }
  if (length(moran_source)) {
    moran_i <- rbindlist(moran_source, use.names = TRUE, fill = TRUE)[,
      .(weighted_feature_morans_i = stats::weighted.mean(
          feature_morans_i, source_n_pixels, na.rm = TRUE),
        minimum_source_feature_morans_i = min(feature_morans_i, na.rm = TRUE)),
      by = matrix_feature_index]
  } else {
    moran_i <- data.table(matrix_feature_index = feature_i,
      weighted_feature_morans_i = NA_real_,
      minimum_source_feature_morans_i = NA_real_)
  }
  moran_i[!is.finite(weighted_feature_morans_i), weighted_feature_morans_i := NA_real_]
  moran_i[!is.finite(minimum_source_feature_morans_i), minimum_source_feature_morans_i := NA_real_]

  for (domain_i in unique(as.character(candidate_i$domain))) {
    cand_d <- candidate_i[domain == domain_i]
    pos <- match(cand_d$matrix_feature_index, feature_i)
    in_core <- as.character(meta_core_i$consensus_domain) == domain_i
    out_core <- !in_core
    if (sum(in_core) < within_min_pixels_each_group ||
        sum(out_core) < within_min_pixels_each_group) next
    x_d <- as.matrix(x_core_i[, pos, drop = FALSE])

    ## 高表达阈值在每个imzML来源内部定义，避免来源间绝对强度差异。
    high <- matrix(FALSE, nrow(x_d), ncol(x_d))
    for (sid in unique(as.character(meta_core_i$source_id))) {
      lr <- which(as.character(meta_core_i$source_id) == sid)
      if (!length(lr)) next
      threshold <- apply(x_d[lr, , drop = FALSE], 2L,
        function(v) safe_quantile(v[v > 0], high_intensity_quantile))
      threshold[!is.finite(threshold)] <- Inf
      high[lr, ] <- x_d[lr, , drop = FALSE] > 0 &
        sweep(x_d[lr, , drop = FALSE], 2L, threshold, ">=")
    }
    high[!is.finite(x_d)] <- FALSE

    spatial_d <- data.table(
      sample_id = sample_i,
      domain = domain_i,
      matrix_feature_index = as.integer(cand_d$matrix_feature_index),
      mz = as.numeric(cand_d$mz),
      stable_core_n_in = sum(in_core),
      stable_core_n_out = sum(out_core),
      stable_core_auc = vapply(seq_len(ncol(x_d)),
        function(j) auc_rank_statistic(x_d[in_core, j], x_d[out_core, j]), numeric(1)),
      stable_core_mean_logCPM_difference = colMeans(x_d[in_core, , drop = FALSE]) -
        colMeans(x_d[out_core, , drop = FALSE]),
      stable_core_median_logCPM_difference = apply(x_d[in_core, , drop = FALSE], 2L, median) -
        apply(x_d[out_core, , drop = FALSE], 2L, median),
      high_intensity_fraction_in_domain = colMeans(high[in_core, , drop = FALSE]),
      high_intensity_fraction_outside_domain = colMeans(high[out_core, , drop = FALSE]),
      high_intensity_domain_enrichment = colMeans(high[in_core, , drop = FALSE]) /
        pmax(colMeans(high), .Machine$double.eps),
      fraction_high_intensity_pixels_in_target_domain = colSums(high[in_core, , drop = FALSE]) /
        pmax(colSums(high), 1)
    )
    spatial_d <- merge(spatial_d, moran_i,
      by = "matrix_feature_index", all.x = TRUE, sort = FALSE)
    wsi <- wsi + 1L; within_spatial_list[[wsi]] <- spatial_d

    ## 空间块子抽样bootstrap：保持局部空间相关，不把单个像素当独立重复。
    block_values <- unique(as.character(meta_core_i$block_id_1C))
    block_index <- match(as.character(meta_core_i$block_id_1C), block_values)
    block_sum_in <- rowsum(x_d * as.numeric(in_core), block_index, reorder = FALSE)
    block_sum_out <- rowsum(x_d * as.numeric(out_core), block_index, reorder = FALSE)
    block_n_in <- as.numeric(rowsum(as.numeric(in_core), block_index, reorder = FALSE))
    block_n_out <- as.numeric(rowsum(as.numeric(out_core), block_index, reorder = FALSE))
    boot_effect <- matrix(NA_real_, spatial_bootstrap_repeats, ncol(x_d))
    set.seed(random_seed + match(sample_i, unique(as.character(within_nr$sample_id))) * 100L +
      match(domain_i, unique(as.character(candidate_i$domain))))
    for (b in seq_len(spatial_bootstrap_repeats)) {
      keep_blocks <- sample(block_values,
        max(1L, ceiling(length(block_values) * spatial_bootstrap_block_fraction)),
        replace = FALSE)
      keep_block_index <- match(keep_blocks, block_values)
      nin <- sum(block_n_in[keep_block_index])
      nout <- sum(block_n_out[keep_block_index])
      if (nin < 10L || nout < 10L) next
      boot_effect[b, ] <- colSums(block_sum_in[keep_block_index, , drop = FALSE]) / nin -
        colSums(block_sum_out[keep_block_index, , drop = FALSE]) / nout
    }
    boot_d <- data.table(
      sample_id = sample_i, domain = domain_i,
      matrix_feature_index = as.integer(cand_d$matrix_feature_index), mz = as.numeric(cand_d$mz),
      spatial_block_bootstrap_valid_repeats = colSums(is.finite(boot_effect)),
      spatial_block_bootstrap_effect_median = apply(boot_effect, 2L, median, na.rm = TRUE),
      spatial_block_bootstrap_CI_low = apply(boot_effect, 2L, safe_quantile, probability = 0.025),
      spatial_block_bootstrap_CI_high = apply(boot_effect, 2L, safe_quantile, probability = 0.975),
      spatial_block_bootstrap_positive_fraction = colMeans(boot_effect > 0, na.rm = TRUE)
    )
    wbi <- wbi + 1L; within_boot_list[[wbi]] <- boot_d

    ## 不同核心概率阈值下重新计算效应，评估候选对核心定义的敏感性。
    for (probability_i in within_core_probability_thresholds) {
      keep <- is.finite(meta_i$consensus_probability) &
        meta_i$consensus_probability >= probability_i
      sin <- keep & as.character(meta_i$consensus_domain) == domain_i
      sout <- keep & as.character(meta_i$consensus_domain) != domain_i
      if (sum(sin) < 10L || sum(sout) < 10L) next
      sens <- data.table(
        sample_id = sample_i, domain = domain_i,
        matrix_feature_index = as.integer(cand_d$matrix_feature_index), mz = as.numeric(cand_d$mz),
        consensus_probability_threshold = probability_i,
        n_pixels_in = sum(sin), n_pixels_out = sum(sout),
        mean_logCPM_difference = Matrix::colMeans(x_all_i[sin, pos, drop = FALSE]) -
          Matrix::colMeans(x_all_i[sout, pos, drop = FALSE])
      )
      wti <- wti + 1L; within_sensitivity_list[[wti]] <- sens
    }
  }
  rm(x_all_i, x_core_i); invisible(gc())
}

within_spatial_validation <- rbindlist(within_spatial_list, use.names = TRUE, fill = TRUE)
within_spatial_bootstrap <- rbindlist(within_boot_list, use.names = TRUE, fill = TRUE)
within_core_sensitivity <- rbindlist(within_sensitivity_list, use.names = TRUE, fill = TRUE)
within_core_sensitivity_summary <- within_core_sensitivity[, .(
  n_thresholds_evaluable = .N,
  minimum_effect_across_core_thresholds = min(mean_logCPM_difference, na.rm = TRUE),
  median_effect_across_core_thresholds = median(mean_logCPM_difference, na.rm = TRUE),
  positive_core_threshold_fraction = mean(mean_logCPM_difference > 0, na.rm = TRUE)
), by = .(sample_id, domain, matrix_feature_index, mz)]

within_validated <- Reduce(function(x, y) merge(x, y,
  by = c("sample_id", "domain", "matrix_feature_index", "mz"),
  all.x = TRUE, sort = FALSE),
  list(within_nr, within_spatial_validation, within_spatial_bootstrap,
    within_core_sensitivity_summary))
within_validated[, strong_spatial_validation :=
  stable_core_auc >= 0.70 &
  spatial_block_bootstrap_positive_fraction >= 0.95 &
  positive_core_threshold_fraction >= 1]
within_validated[, spatial_validation_score := candidate_score *
  pmax(2 * stable_core_auc - 1, 0) *
  spatial_block_bootstrap_positive_fraction *
  pmax(high_intensity_domain_enrichment, 0)]
setorder(within_validated, sample_id, domain, -spatial_validation_score)
within_validated[strong_spatial_validation == TRUE,
  spatial_validation_rank := seq_len(.N), by = .(sample_id, domain)]

fwrite(within_spatial_validation,
  file.path(within_spatial_dir, "within_feature_spatial_agreement_metrics.csv"))
fwrite(within_spatial_bootstrap,
  file.path(within_spatial_dir, "within_feature_spatial_block_bootstrap.csv"))
fwrite(within_core_sensitivity,
  file.path(within_spatial_dir, "within_feature_core_threshold_sensitivity.csv"))
fwrite(within_validated,
  file.path(within_spatial_dir, "within_nonredundant_candidates_with_spatial_validation.csv"))
fwrite(within_validated[strong_spatial_validation == TRUE & spatial_validation_rank <= 30L],
  file.path(within_spatial_dir,
    "top30_within_strong_spatially_validated_features_per_domain.csv"))

within_validation_counts <- within_validated[, .(
  n_nonredundant_candidates = .N,
  n_strong_spatially_validated = sum(strong_spatial_validation, na.rm = TRUE),
  median_auc = median(stable_core_auc, na.rm = TRUE),
  median_feature_morans_i = median(weighted_feature_morans_i, na.rm = TRUE),
  median_block_bootstrap_positive_fraction =
    median(spatial_block_bootstrap_positive_fraction, na.rm = TRUE)
), by = .(sample_id, domain)]
fwrite(within_validation_counts,
  file.path(within_spatial_dir, "within_spatial_validation_summary_by_domain.csv"))

p_within_auc <- ggplot(within_validated,
    aes(x = interaction(sample_id, domain, sep = " / "), y = stable_core_auc,
      color = strong_spatial_validation)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey60") +
  geom_boxplot(outlier.shape = NA, color = "grey35") +
  geom_jitter(width = 0.18, alpha = 0.55, size = 0.9) +
  coord_flip() + scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#999999")) +
  labs(title = "Spatial agreement of within-tumor domain-associated features",
    x = NULL, y = "AUC: target domain versus other domains",
    color = "Strong spatial validation") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(within_spatial_dir, "within_feature_domain_auc.png"),
  p_within_auc, width = 9, height = 8, dpi = save_png_dpi, bg = "white")

## ==== 8) Part D2：共同meta-domain两两比较和专属性分类 ====
message("[Part D2] 计算M1/M2/M3两两比较、样本bootstrap和专属性分类。")
validation_features <- sort(unique(as.integer(cross_recurrent_input[
  passed_strict_technical_filter == TRUE & selected_for_common_domain_model == FALSE,
  matrix_feature_index])))
cross_core <- copy(cross_pixels[stable_core_pixel == TRUE & !is.na(common_meta_domain)])
sample_domain_summary_list <- list(); sdsi <- 0L
for (sample_i in unique(as.character(cross_core$sample_id))) {
  meta_s <- cross_core[sample_id == sample_i]
  for (domain_i in sort(unique(as.character(meta_s$common_meta_domain)))) {
    rows <- as.integer(meta_s[common_meta_domain == domain_i, original_pixel_row])
    if (!length(rows)) next
    x <- aligned_logCPM[rows, validation_features, drop = FALSE]
    sdsi <- sdsi + 1L
    sample_domain_summary_list[[sdsi]] <- data.table(
      sample_id = sample_i,
      diagnosis = as.character(meta_s$diagnosis[1]),
      common_meta_domain = domain_i,
      matrix_feature_index = validation_features,
      mz = feature_mz[validation_features],
      n_pixels = length(rows),
      mean_logCPM = as.numeric(Matrix::colMeans(x)),
      prevalence = as.numeric(Matrix::colMeans(x > 0))
    )
  }
}
sample_domain_summary <- rbindlist(sample_domain_summary_list)
fwrite(sample_domain_summary,
  file.path(cross_pairwise_dir, "physical_sample_domain_feature_means_validation_features.csv"))

domain_pairs <- list(c("M1", "M2"), c("M1", "M3"), c("M2", "M3"))
pairwise_effect_list <- list(); pei <- 0L
for (pair_i in domain_pairs) {
  a <- sample_domain_summary[common_meta_domain == pair_i[1],
    .(sample_id, diagnosis, matrix_feature_index, mz,
      n_pixels_A = n_pixels, mean_A = mean_logCPM, prevalence_A = prevalence)]
  b <- sample_domain_summary[common_meta_domain == pair_i[2],
    .(sample_id, matrix_feature_index,
      n_pixels_B = n_pixels, mean_B = mean_logCPM, prevalence_B = prevalence)]
  z <- merge(a, b, by = c("sample_id", "matrix_feature_index"), all = FALSE)
  z[, `:=`(
    contrast = paste0(pair_i[1], "_vs_", pair_i[2]),
    domain_A = pair_i[1], domain_B = pair_i[2],
    mean_logCPM_difference = mean_A - mean_B,
    prevalence_difference = prevalence_A - prevalence_B,
    total_contrast_pixels = n_pixels_A + n_pixels_B
  )]
  pei <- pei + 1L; pairwise_effect_list[[pei]] <- z
}
pairwise_sample_effect <- rbindlist(pairwise_effect_list, use.names = TRUE, fill = TRUE)
fwrite(pairwise_sample_effect,
  file.path(cross_pairwise_dir, "physical_sample_pairwise_meta_domain_effects.csv"))

pairwise_stats <- pairwise_sample_effect[, {
  s <- summarize_signed_effect(mean_logCPM_difference, prevalence_difference)
  boot <- bootstrap_sample_median(mean_logCPM_difference,
    sample_bootstrap_repeats,
    random_seed + .GRP + matrix_feature_index[1])
  c(as.list(s[1]), as.list(boot))
}, by = .(contrast, domain_A, domain_B, matrix_feature_index, mz)]
pairwise_stats[, `:=`(
  sign_test_FDR = p.adjust(sign_test_p, method = "BH"),
  wilcoxon_signed_rank_FDR = p.adjust(wilcoxon_signed_rank_p, method = "BH")
), by = contrast]

## NF-GnH-3三个导出区域只作为同一物理样本内的技术一致性检查。
nfg_source <- cross_core[sample_id == "NF-GnH-3"]
nfg_source_summary_list <- list(); nssi <- 0L
for (sid in unique(as.character(nfg_source$source_id))) {
  meta_s <- nfg_source[source_id == sid]
  for (domain_i in sort(unique(as.character(meta_s$common_meta_domain)))) {
    rows <- as.integer(meta_s[common_meta_domain == domain_i, original_pixel_row])
    if (!length(rows)) next
    x <- aligned_logCPM[rows, validation_features, drop = FALSE]
    nssi <- nssi + 1L
    nfg_source_summary_list[[nssi]] <- data.table(
      source_id = sid, common_meta_domain = domain_i,
      matrix_feature_index = validation_features,
      mean_logCPM = as.numeric(Matrix::colMeans(x)))
  }
}
nfg_source_summary <- rbindlist(nfg_source_summary_list, use.names = TRUE, fill = TRUE)
nfg_pairwise_list <- list(); npi <- 0L
for (pair_i in domain_pairs) {
  a <- nfg_source_summary[common_meta_domain == pair_i[1],
    .(source_id, matrix_feature_index, mean_A = mean_logCPM)]
  b <- nfg_source_summary[common_meta_domain == pair_i[2],
    .(source_id, matrix_feature_index, mean_B = mean_logCPM)]
  z <- merge(a, b, by = c("source_id", "matrix_feature_index"), all = FALSE)
  z[, `:=`(contrast = paste0(pair_i[1], "_vs_", pair_i[2]),
    partition_effect = mean_A - mean_B)]
  npi <- npi + 1L; nfg_pairwise_list[[npi]] <- z
}
nfg_pairwise <- rbindlist(nfg_pairwise_list, use.names = TRUE, fill = TRUE)
nfg_pairwise_consistency <- nfg_pairwise[, .(
  n_nfgnh3_partitions_evaluable = .N,
  nfgnh3_positive_partition_fraction = mean(partition_effect > 0),
  nfgnh3_negative_partition_fraction = mean(partition_effect < 0),
  nfgnh3_partition_direction_consistency = max(
    mean(partition_effect > 0), mean(partition_effect < 0))
), by = .(contrast, matrix_feature_index)]
pairwise_stats <- merge(pairwise_stats, nfg_pairwise_consistency,
  by = c("contrast", "matrix_feature_index"), all.x = TRUE, sort = FALSE)

## 留一诊断组和低像素样本敏感性。
pairwise_robustness <- pairwise_sample_effect[, {
  effect <- mean_logCPM_difference
  full_sign <- sign(stats::median(effect))
  diagnoses <- unique(as.character(diagnosis))
  lodo <- vapply(diagnoses, function(dx) {
    z <- effect[as.character(diagnosis) != dx]
    if (!length(z)) return(NA_real_)
    stats::median(z)
  }, numeric(1))
  pixel_cut <- safe_quantile(total_contrast_pixels, 0.20)
  high_pixel_effect <- effect[total_contrast_pixels >= pixel_cut]
  list(
    minimum_leave_one_diagnosis_out_absolute_median = min(abs(lodo), na.rm = TRUE),
    leave_one_diagnosis_out_same_direction_fraction =
      mean(sign(lodo) == full_sign, na.rm = TRUE),
    low_pixel_excluded_median_effect = stats::median(high_pixel_effect),
    low_pixel_excluded_direction_consistency = max(
      mean(high_pixel_effect > 0), mean(high_pixel_effect < 0))
  )
}, by = .(contrast, matrix_feature_index)]
pairwise_stats <- merge(pairwise_stats, pairwise_robustness,
  by = c("contrast", "matrix_feature_index"), all.x = TRUE, sort = FALSE)

pairwise_stats[, enriched_domain := fifelse(median_effect > 0, domain_A, domain_B)]
pairwise_stats[, `:=`(
  low_mz_background_risk = mz < low_mz_background_threshold,
  requires_matrix_blank_review = mz < low_mz_background_threshold
)]
pairwise_stats[, recurrent_pairwise_difference :=
  n_physical_samples_evaluable >= pairwise_min_samples_evaluable &
  direction_consistency >= pairwise_min_direction_consistency &
  abs(median_effect) >= pairwise_min_median_effect &
  minimum_leave_one_sample_out_signed_margin > 0 &
  wilcoxon_signed_rank_FDR <= pairwise_max_wilcoxon_FDR &
  bootstrap_same_direction_fraction >= pairwise_min_bootstrap_direction_fraction &
  fifelse(median_effect > 0,
    bootstrap_median_CI_low > 0, bootstrap_median_CI_high < 0) &
  (abs(median_prevalence_difference) >= pairwise_min_prevalence_difference |
     abs(median_effect) >= pairwise_strong_effect_without_prevalence) &
  (is.na(nfgnh3_partition_direction_consistency) |
     (n_nfgnh3_partitions_evaluable >= 2L &
      nfgnh3_partition_direction_consistency >=
        cross_min_nfgnh3_partition_consistency))]
pairwise_stats[, robust_pairwise_difference := recurrent_pairwise_difference &
  leave_one_diagnosis_out_same_direction_fraction >= 1 &
  sign(low_pixel_excluded_median_effect) == sign(median_effect)]
fwrite(pairwise_stats,
  file.path(cross_pairwise_dir, "pairwise_meta_domain_feature_statistics_all.csv"))
fwrite(pairwise_stats[recurrent_pairwise_difference == TRUE],
  file.path(cross_pairwise_dir, "recurrent_pairwise_meta_domain_differences.csv"))
## 记录原1C one-versus-rest非冗余候选；两两比较可能发现被合并对照稀释的额外feature。
original_cross_feature_universe <- unique(cross_all[nonredundant_representative == TRUE,
  .(matrix_feature_index, mz)])
classify_one_feature <- function(feature_index_i) {
  assoc <- sort(unique(as.character(cross_all[
    nonredundant_representative == TRUE & matrix_feature_index == feature_index_i,
    common_meta_domain])))
  ps <- pairwise_stats[matrix_feature_index == feature_index_i]
  signed_test <- function(target, other) {
    z <- ps[(domain_A == target & domain_B == other) |
      (domain_A == other & domain_B == target)]
    if (!nrow(z)) return(list(pass = FALSE, robust = FALSE, effect = NA_real_))
    effect <- if (z$domain_A[1] == target) z$median_effect[1] else -z$median_effect[1]
    robust <- isTRUE(z$robust_pairwise_difference[1])
    list(pass = robust && effect > 0, robust = robust,
      effect = effect)
  }
  domains <- c("M1", "M2", "M3")
  specific <- domains[vapply(domains, function(d) {
    all(vapply(setdiff(domains, d), function(o) signed_test(d, o)$pass, logical(1)))
  }, logical(1))]
  if (length(specific) == 1L) {
    return(data.table(matrix_feature_index = feature_index_i,
      specificity_class = paste0(specific, "_specific"),
      primary_domain = specific, associated_one_vs_rest_domains = paste(assoc, collapse = ",")))
  }
  shared_labels <- character()
  for (low in domains) {
    high <- setdiff(domains, low)
    high_over_low <- all(vapply(high, function(d) signed_test(d, low)$pass, logical(1)))
    between <- signed_test(high[1], high[2])
    similar <- !between$robust || !is.finite(between$effect) ||
      abs(between$effect) < shared_domain_similarity_effect
    if (high_over_low && similar) {
      shared_labels <- c(shared_labels,
        paste0(paste(high, collapse = "_"), "_shared_over_", low))
    }
  }
  if (length(shared_labels) == 1L) {
    cls <- shared_labels; primary <- sub("_shared.*$", "", shared_labels)
  } else if (length(assoc) > 1L) {
    cls <- "multi_domain_associated_unresolved"; primary <- NA_character_
  } else if (!length(assoc)) {
    cls <- "pairwise_pattern_unresolved"; primary <- NA_character_
  } else {
    cls <- "domain_associated_not_exclusive"; primary <- assoc[1]
  }
  data.table(matrix_feature_index = feature_index_i,
    specificity_class = cls, primary_domain = primary,
    associated_one_vs_rest_domains = paste(assoc, collapse = ","))
}
## 共同domain核心margin阈值敏感性，覆盖全部独立validation features。
cross_sensitivity_features <- validation_features
cross_margin_sensitivity_list <- list(); cmsi <- 0L
for (margin_i in cross_margin_thresholds) {
  meta_margin <- cross_pixels[is.finite(domain_assignment_margin) &
    domain_assignment_margin >= margin_i & !is.na(common_meta_domain)]
  margin_sample_domain <- list(); mdi <- 0L
  for (sample_i in unique(as.character(meta_margin$sample_id))) {
    meta_s <- meta_margin[sample_id == sample_i]
    for (domain_i in sort(unique(as.character(meta_s$common_meta_domain)))) {
      rows <- as.integer(meta_s[common_meta_domain == domain_i, original_pixel_row])
      if (length(rows) < 10L) next
      mdi <- mdi + 1L
      margin_sample_domain[[mdi]] <- data.table(
        sample_id = sample_i, common_meta_domain = domain_i,
        matrix_feature_index = cross_sensitivity_features,
        mean_logCPM = as.numeric(Matrix::colMeans(
          aligned_logCPM[rows, cross_sensitivity_features, drop = FALSE])))
    }
  }
  margin_summary <- rbindlist(margin_sample_domain)
  for (pair_i in domain_pairs) {
    a <- margin_summary[common_meta_domain == pair_i[1],
      .(sample_id, matrix_feature_index, mean_A = mean_logCPM)]
    b <- margin_summary[common_meta_domain == pair_i[2],
      .(sample_id, matrix_feature_index, mean_B = mean_logCPM)]
    z <- merge(a, b, by = c("sample_id", "matrix_feature_index"), all = FALSE)
    z[, effect := mean_A - mean_B]
    zz <- z[, .(
      n_samples = .N,
      median_effect = median(effect),
      direction_consistency = max(mean(effect > 0), mean(effect < 0))
    ), by = matrix_feature_index]
    zz[, `:=`(margin_threshold = margin_i,
      contrast = paste0(pair_i[1], "_vs_", pair_i[2]))]
    cmsi <- cmsi + 1L; cross_margin_sensitivity_list[[cmsi]] <- zz
  }
}
cross_margin_sensitivity <- rbindlist(cross_margin_sensitivity_list,
  use.names = TRUE, fill = TRUE)
fwrite(cross_margin_sensitivity,
  file.path(cross_pairwise_dir, "cross_pairwise_core_margin_sensitivity.csv"))
pairwise_margin_summary <- merge(
  cross_margin_sensitivity,
  pairwise_stats[, .(contrast, matrix_feature_index, full_core_effect = median_effect)],
  by = c("contrast", "matrix_feature_index"), all.x = TRUE, sort = FALSE
)[, .(
  n_margin_thresholds_evaluable = .N,
  core_margin_same_direction_fraction = mean(
    sign(median_effect) == sign(full_core_effect), na.rm = TRUE),
  minimum_absolute_effect_across_margin_thresholds =
    min(abs(median_effect), na.rm = TRUE),
  minimum_direction_consistency_across_margin_thresholds =
    min(direction_consistency, na.rm = TRUE)
), by = .(contrast, matrix_feature_index)]
fwrite(pairwise_margin_summary,
  file.path(cross_pairwise_dir, "cross_pairwise_core_margin_sensitivity_summary.csv"))
pairwise_stats <- merge(pairwise_stats, pairwise_margin_summary,
  by = c("contrast", "matrix_feature_index"), all.x = TRUE, sort = FALSE)
pairwise_stats[, robust_pairwise_difference := robust_pairwise_difference &
  core_margin_same_direction_fraction >= 1 &
  minimum_direction_consistency_across_margin_thresholds >=
    pairwise_min_direction_consistency]
## 覆盖前面的统计表，使最终版本同时包含核心margin敏感性字段。
fwrite(pairwise_stats,
  file.path(cross_pairwise_dir, "pairwise_meta_domain_feature_statistics_all.csv"))
fwrite(pairwise_stats[recurrent_pairwise_difference == TRUE],
  file.path(cross_pairwise_dir, "recurrent_pairwise_meta_domain_differences.csv"))

## 完成所有稳健性检查后再给每个不同m/z赋予最终专属性类别。
fwrite(pairwise_stats[robust_pairwise_difference == TRUE],
  file.path(cross_pairwise_dir,
    "robust_pairwise_differences_all_sensitivity_checks.csv"))

classification_feature_indices <- sort(unique(c(
  original_cross_feature_universe$matrix_feature_index,
  pairwise_stats[robust_pairwise_difference == TRUE, matrix_feature_index]
)))
cross_feature_universe <- data.table(
  matrix_feature_index = classification_feature_indices,
  mz = feature_mz[classification_feature_indices]
)
specificity_classification <- rbindlist(lapply(
  cross_feature_universe$matrix_feature_index, classify_one_feature))
specificity_classification <- merge(cross_feature_universe,
  specificity_classification, by = "matrix_feature_index", all.x = TRUE)
cross_feature_priority <- cross_all[nonredundant_representative == TRUE,
  .(one_vs_rest_best_candidate_score = max(candidate_score),
    best_one_vs_rest_evidence_tier = if (any(evidence_tier == "Tier_A"))
      "Tier_A" else "Tier_B"),
  by = matrix_feature_index]
specificity_classification <- merge(specificity_classification,
  cross_feature_priority, by = "matrix_feature_index", all.x = TRUE, sort = FALSE)
pairwise_feature_priority <- pairwise_stats[robust_pairwise_difference == TRUE,
  .(pairwise_best_absolute_median_effect = max(abs(median_effect), na.rm = TRUE),
    pairwise_minimum_direction_consistency = min(direction_consistency, na.rm = TRUE)),
  by = matrix_feature_index]
specificity_classification <- merge(specificity_classification,
  pairwise_feature_priority, by = "matrix_feature_index", all.x = TRUE, sort = FALSE)
specificity_classification[, `:=`(
  originated_from_one_vs_rest_candidate =
    matrix_feature_index %in% original_cross_feature_universe$matrix_feature_index,
  low_mz_background_risk = mz < low_mz_background_threshold,
  requires_matrix_blank_review = mz < low_mz_background_threshold
)]
specificity_classification[, specificity_priority_score :=
  fcoalesce(one_vs_rest_best_candidate_score, 0) +
  fcoalesce(pairwise_best_absolute_median_effect, 0)]
setorder(specificity_classification, specificity_class,
  -specificity_priority_score)
fwrite(specificity_classification,
  file.path(cross_pairwise_dir, "cross_feature_specificity_classification.csv"))

cross_all_classified <- merge(cross_all, specificity_classification[,
  .(matrix_feature_index, specificity_class, primary_domain)],
  by = "matrix_feature_index", all.x = TRUE, sort = FALSE)
fwrite(cross_all_classified,
  file.path(cross_pairwise_dir,
    "cross_candidates_with_pairwise_specificity_classification.csv"))

specificity_counts <- specificity_classification[, .N, by = specificity_class][order(-N)]
fwrite(specificity_counts,
  file.path(cross_pairwise_dir, "cross_feature_specificity_class_counts.csv"))
p_specificity <- ggplot(specificity_counts,
    aes(x = reorder(specificity_class, N), y = N)) +
  geom_col(fill = "#0072B2") + coord_flip() +
  labs(title = "Pairwise specificity of recurrent meta-domain features",
    x = NULL, y = "Unique nonredundant features") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(cross_pairwise_dir, "cross_feature_specificity_class_counts.png"),
  p_specificity, width = 8, height = 5.5, dpi = save_png_dpi, bg = "white")

## 质量关系网络边表，保留被合并和仅标记关系的审计轨迹。
within_relation_network <- within_all[possible_mass_relation_flag == TRUE,
  .(sample_id, domain, feature_index = matrix_feature_index, mz,
    representative_feature_index, representative_mz,
    relation_type = mass_relation_type_to_representative,
    relation_name = mass_relation_to_representative,
    observed_delta_da = observed_delta_to_representative_da,
    expected_delta_da = expected_delta_to_representative_da,
    mass_error_da = mass_relation_error_da,
    correlation = correlation_to_representative,
    auto_collapsed = !nonredundant_representative)]
cross_relation_network <- cross_all[possible_mass_relation_flag == TRUE,
  .(common_meta_domain, feature_index = matrix_feature_index, mz,
    representative_feature_index, representative_mz,
    relation_type = mass_relation_type_to_representative,
    relation_name = mass_relation_to_representative,
    observed_delta_da = observed_delta_to_representative_da,
    expected_delta_da = expected_delta_to_representative_da,
    mass_error_da = mass_relation_error_da,
    correlation = correlation_to_representative,
    auto_collapsed = !nonredundant_representative)]
fwrite(within_relation_network,
  file.path(relation_dir, "within_candidate_mass_relation_network_edges.csv"))
fwrite(cross_relation_network,
  file.path(relation_dir, "cross_candidate_mass_relation_network_edges.csv"))

## ==== 9) 候选feature空间离子图 ====
if (make_ion_maps) {
  message("[离子图] 绘制样本内和跨样本代表feature。")

  within_map_features <- within_validated[
    strong_spatial_validation == TRUE &
      spatial_validation_rank <= within_ion_maps_per_domain]
  for (i in seq_len(nrow(within_map_features))) {
    z <- within_map_features[i]
    dt <- copy(within_pixels[sample_id == z$sample_id])
    dt[, intensity := as.numeric(aligned_logCPM[
      as.integer(original_pixel_row), z$matrix_feature_index, drop = TRUE])]
    dt[, normalized_intensity := normalize_ion_intensity(intensity), by = source_id]
    dt[, domain_boundary := domain_boundary_flag(
      .SD, "consensus_domain", as.character(z$domain))]
    p <- plot_faceted_ion_map(
      dt,
      paste0(z$sample_id, " / ", z$domain, " / m/z ", sprintf("%.4f", z$mz)),
      "Independent validation feature; cyan marks the target-domain boundary"
    )
    fn <- paste0(safe_filename(z$sample_id), "_", z$domain,
      "_spatialRank", z$spatial_validation_rank,
      "_mz", sprintf("%.4f", z$mz), ".png")
    ggsave(file.path(within_map_dir, fn), p,
      width = ifelse(uniqueN(dt$source_id) > 1, 11, 6.5),
      height = 6, dpi = save_png_dpi, bg = "white")
  }

  cross_map_features <- specificity_classification[
    grepl("_specific$", specificity_class)]
  cross_map_features[, common_meta_domain := primary_domain]
  setorder(cross_map_features, common_meta_domain, -specificity_priority_score)
  cross_map_features[, map_rank := seq_len(.N), by = common_meta_domain]
  cross_map_features <- cross_map_features[map_rank <= cross_ion_maps_per_domain]
  for (i in seq_len(nrow(cross_map_features))) {
    z <- cross_map_features[i]
    dt <- copy(cross_pixels)
    dt[, intensity := as.numeric(aligned_logCPM[
      as.integer(original_pixel_row), z$matrix_feature_index, drop = TRUE])]
    dt[, normalized_intensity := normalize_ion_intensity(intensity), by = source_id]
    dt[, domain_boundary := domain_boundary_flag(
      .SD, "common_meta_domain", as.character(z$common_meta_domain))]
    dt[, source_id := factor(source_id, levels = unique(as.character(source_id)))]
    p <- plot_faceted_ion_map(
      dt,
      paste0(z$common_meta_domain, " recurrent feature / m/z ",
        sprintf("%.4f", z$mz)),
      paste0("Pairwise class: ", z$specificity_class,
        "; cyan marks the target meta-domain boundary")
    )
    fn <- paste0(z$common_meta_domain, "_specific_rank", z$map_rank,
      "_mz", sprintf("%.4f", z$mz), ".png")
    ggsave(file.path(cross_map_dir, fn), p,
      width = 13, height = 10, dpi = save_png_dpi, bg = "white")
  }
}

## ==== 10) 汇总和参数 ====
analysis_summary <- rbindlist(list(
  data.table(
    analysis = "within_sample",
    n_analysis_units = uniqueN(within_all$sample_id),
    n_domains = uniqueN(paste(within_all$sample_id, within_all$domain)),
    n_candidate_rows = nrow(within_all),
    n_nonredundant_candidate_rows = sum(within_all$nonredundant_representative),
    n_unique_features = uniqueN(within_all$matrix_feature_index)
  ),
  data.table(
    analysis = "cross_sample",
    n_analysis_units = uniqueN(cross_sample_effect$sample_id),
    n_domains = uniqueN(cross_all$common_meta_domain),
    n_candidate_rows = nrow(cross_all),
    n_nonredundant_candidate_rows = sum(cross_all$nonredundant_representative),
    n_unique_features = uniqueN(cross_all$matrix_feature_index)
  ),
  data.table(
    analysis = "within_strong_spatial_validation",
    n_analysis_units = uniqueN(within_validated$sample_id),
    n_domains = uniqueN(paste(within_validated$sample_id, within_validated$domain)),
    n_candidate_rows = nrow(within_validated),
    n_nonredundant_candidate_rows = sum(
      within_validated$strong_spatial_validation, na.rm = TRUE),
    n_unique_features = uniqueN(
      within_validated[strong_spatial_validation == TRUE, matrix_feature_index])
  ),
  data.table(
    analysis = "cross_pairwise_recurrent",
    n_analysis_units = uniqueN(pairwise_sample_effect$sample_id),
    n_domains = uniqueN(pairwise_stats$contrast),
    n_candidate_rows = nrow(pairwise_stats),
    n_nonredundant_candidate_rows = sum(
      pairwise_stats$recurrent_pairwise_difference, na.rm = TRUE),
    n_unique_features = uniqueN(
      pairwise_stats[recurrent_pairwise_difference == TRUE, matrix_feature_index])
  )
), use.names = TRUE, fill = TRUE)
fwrite(analysis_summary, file.path(out_dir, "domain_feature_analysis_summary.csv"))

analysis_dictionary <- data.table(
  output_concept = c(
    "within_domain_associated",
    "within_strong_spatial_validation",
    "recurrent_one_vs_rest_associated",
    "robust_pairwise_difference",
    "Mx_specific",
    "Mx_My_shared_over_Mz",
    "domain_associated_not_exclusive",
    "low_mz_background_risk"
  ),
  interpretation = c(
    "Within one physical tumor, the feature is enriched in one D-domain versus its other D-domains.",
    "Within-domain association also has AUC >= 0.70, positive spatial-block support >= 0.95, and the same direction at every tested core threshold.",
    "Across physical samples, the feature is enriched in one M-domain versus the pooled other M-domains; this alone does not prove exclusivity.",
    "A direct M-domain pairwise contrast passes effect, FDR, bootstrap, leave-one-sample-out, leave-one-diagnosis-out, low-pixel and NF-GnH-3 partition checks.",
    "One M-domain is robustly higher than each of the other two M-domains.",
    "Two M-domains are both robustly higher than the third, without a robust difference between those two domains.",
    "The one-versus-rest association is reproducible, but direct pairwise evidence is insufficient to call the feature exclusive.",
    "m/z < 200; matrix-only blank and fragmentation review is required before biological interpretation."
  )
)
fwrite(analysis_dictionary, file.path(out_dir, "analysis_output_dictionary.csv"))

saveRDS(list(
  within_thresholds = list(
    min_effect = within_min_effect,
    strong_effect_without_prevalence = within_strong_effect_without_prevalence,
    min_prevalence_difference = within_min_prevalence_difference,
    min_pixels_each_group = within_min_pixels_each_group,
    min_partition_direction_consistency = within_min_partition_direction_consistency
  ),
  cross_thresholds = list(
    min_samples_evaluable = cross_min_samples_evaluable,
    min_direction_consistency = cross_min_direction_consistency,
    min_median_effect = cross_min_median_effect,
    min_prevalence_difference = cross_min_prevalence_difference,
    strong_effect_without_prevalence = cross_strong_effect_without_prevalence,
    min_nfgnh3_partition_consistency = cross_min_nfgnh3_partition_consistency
  ),
  redundancy = list(
    ion_mode = ion_mode,
    isotope_mass_tolerance_da = isotope_mass_tolerance_da,
    adduct_mass_tolerance_da = adduct_mass_tolerance_da,
    neutral_loss_mass_tolerance_da = neutral_loss_mass_tolerance_da,
    within_correlation = within_redundancy_correlation,
    cross_correlation = cross_redundancy_correlation,
    adduct_correlation = adduct_redundancy_correlation,
    collapse_isotope_relations = collapse_isotope_relations,
    collapse_adduct_relations = collapse_adduct_relations,
    collapse_neutral_loss_relations = collapse_neutral_loss_relations,
    negative_mode_reference_ion = "[M-H]-"
  ),
  spatial_validation = list(
    neighbor_radius = spatial_neighbor_radius,
    block_width = spatial_block_width,
    block_bootstrap_repeats = spatial_bootstrap_repeats,
    block_fraction = spatial_bootstrap_block_fraction,
    high_intensity_quantile = high_intensity_quantile,
    within_core_probability_thresholds = within_core_probability_thresholds,
    strong_validation_auc = 0.70,
    strong_validation_bootstrap_positive_fraction = 0.95
  ),
  pairwise_specificity = list(
    sample_bootstrap_repeats = sample_bootstrap_repeats,
    minimum_median_effect = pairwise_min_median_effect,
    minimum_direction_consistency = pairwise_min_direction_consistency,
    minimum_prevalence_difference = pairwise_min_prevalence_difference,
    strong_effect_without_prevalence = pairwise_strong_effect_without_prevalence,
    minimum_samples = pairwise_min_samples_evaluable,
    maximum_wilcoxon_FDR = pairwise_max_wilcoxon_FDR,
    minimum_bootstrap_direction_fraction =
      pairwise_min_bootstrap_direction_fraction,
    shared_domain_similarity_effect = shared_domain_similarity_effect,
    cross_margin_thresholds = cross_margin_thresholds
  ),
  random_seed = random_seed,
  session_info = utils::sessionInfo()
), file.path(out_dir, "run_parameters_and_sessionInfo.rds"), compress = "gzip")

message("完成。")
message("样本内结果：", within_dir)
message("跨样本结果：", cross_dir)
message("两套domain连接结果：", integration_dir)
message("专属性和空间稳健性结果：", robustness_dir)
message("样本内候选：", nrow(within_all),
  "；去冗余后：", sum(within_all$nonredundant_representative))
message("跨样本候选：", nrow(cross_all),
  "；去冗余后：", sum(cross_all$nonredundant_representative))
message("样本内强空间验证：",
  sum(within_validated$strong_spatial_validation, na.rm = TRUE), " / ",
  nrow(within_validated))
message("跨样本显著两两差异：",
  sum(pairwise_stats$recurrent_pairwise_difference, na.rm = TRUE), " / ",
  nrow(pairwise_stats))
message("离子模式：negative；负离子加合物仅在质量差匹配且相关性 >= ",
  adduct_redundancy_correlation, " 时合并。")
message("同位素/加合物质量差容差已收紧为0.02 Da；m/z < ",
  low_mz_background_threshold, " 的候选已标记为需要基质空白复核。")
message("H2O、NH3、CO2 中性丢失关系只标记，不自动删除。")
