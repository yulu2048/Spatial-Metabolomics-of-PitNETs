## ============================================================================
## 0B_tumor_pixel_cluster_map.R
##
## 目的：
##   1) 从 0_readin_align_10samples.R 保存的共识 feature 矩阵继续分析；
##   2) 对 12 个原始 imzML 来源分别根据像素代谢谱聚类；
##   3) 输出 cluster map、TIC/检出峰数辅助图和 cluster 汇总表；
##   4) 本脚本只生成供人工判读的 cluster，不删除或标记肿瘤像素。
##
## 重要：
##   - source_id 是空间分析单位。不同 imzML 的坐标绝不拼接；
##   - cluster 编号只在同一个 source_id 内有效，不能跨来源直接对应；
##   - 完成人工判读后，再把每个 source_id 的肿瘤 cluster 写入下一步筛选脚本。
## ============================================================================

required_packages <- c("data.table", "Matrix", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "。请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

## irlba 不是必需包；安装后大矩阵 PCA 会明显更快。
has_irlba <- requireNamespace("irlba", quietly = TRUE)

## ==== 1) 路径与参数 ====

data_root <- "/home/gaga/huzhan"
input_dir <- file.path(
  data_root, "MSI_ALL_aligned_10samples_consensus_tuned225"
)
out_dir <- file.path(input_dir, "QC", "tumor_pixel_selection_0B_k4")
map_dir <- file.path(out_dir, "cluster_maps")
aux_map_dir <- file.path(out_dir, "auxiliary_maps")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(aux_map_dir, recursive = TRUE, showWarnings = FALSE)

## 默认每个 imzML 来源分成 4 个 cluster。
k_clusters <- 4L

## 每个来源先按 logCPM 方差选择最多 500 个 feature，再做 PCA。
## 这样不会因为 1430 个 feature 中的大量低变异峰拖慢 PCA 或放大噪声。
max_features_for_clustering <- 500L
min_feature_nonzero_prop <- 0.02
min_feature_nonzero_pixels <- 20L

## 最多计算 20 个 PC；用于 k-means 的 PC 数由解释方差自动决定，
## 但不少于 5 个、不多于 20 个。
max_pcs <- 20L
min_pcs_for_clustering <- 5L
pc_variance_target <- 0.80

## 可选：center_only / pareto / unit_variance。
## Pareto scaling 比单位方差缩放温和，能减少高强度 feature 的支配，
## 又不至于让极低丰度 feature 具有过高权重。
pca_scaling <- "pareto"

kmeans_nstart <- 30L
kmeans_iter_max <- 100L
random_seed <- 20250907L

## NULL 表示运行全部 12 个来源；测试时可设为 c("GH", "cushing_0")。
source_ids_to_run <- NULL

## 默认只用代谢谱聚类。不要把 raw_tic、峰数等技术变量混入代谢区域定义。
include_qc_covariates_in_clustering <- FALSE

## 是否另外保存逐来源的 TIC 和检出峰数图，用于辅助判断组织/背景。
save_auxiliary_maps <- TRUE
save_png_dpi <- 240

input_files <- c(
  matrix = file.path(input_dir, "pixel_by_feature_matrix_aligned_logCPM.rds"),
  feature_mz = file.path(input_dir, "feature_mz_aligned_filtered.rds"),
  pixel_meta = file.path(input_dir, "pixel_metadata_sample_source_xy.rds"),
  manifest = file.path(input_dir, "sample_source_manifest.csv")
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))
}

## ==== 2) 读入并核对数据 ====

message("[读入] 共识 feature 矩阵与像素坐标。")
log_mat <- readRDS(input_files[["matrix"]])
feature_mz <- as.numeric(readRDS(input_files[["feature_mz"]]))
pixel_meta <- data.table::as.data.table(readRDS(input_files[["pixel_meta"]]))
manifest <- data.table::fread(input_files[["manifest"]])

if (!inherits(log_mat, "dgCMatrix")) {
  log_mat <- methods::as(log_mat, "dgCMatrix")
}

required_meta_columns <- c(
  "pixel_row", "sample_id", "diagnosis", "source_id", "x", "y",
  "raw_tic", "alignment_nmatch", "alignment_failed", "detected_peak_count"
)
missing_meta_columns <- setdiff(required_meta_columns, names(pixel_meta))
if (length(missing_meta_columns)) {
  stop("像素元数据缺少列：", paste(missing_meta_columns, collapse = ", "))
}

stopifnot(
  nrow(log_mat) == nrow(pixel_meta),
  ncol(log_mat) == length(feature_mz),
  identical(pixel_meta$pixel_row, seq_len(nrow(pixel_meta)))
)

source_levels <- as.character(manifest$source_id)
if (!is.null(source_ids_to_run)) {
  unknown_sources <- setdiff(source_ids_to_run, source_levels)
  if (length(unknown_sources)) {
    stop("source_ids_to_run 中存在未知来源：",
         paste(unknown_sources, collapse = ", "))
  }
  source_levels <- source_levels[source_levels %in% source_ids_to_run]
}

if (!length(source_levels)) stop("没有需要运行的 source_id。")

message(
  "[数据] ", nrow(log_mat), " 个像素 x ", ncol(log_mat),
  " 个 feature；本次处理 ", length(source_levels), " 个来源。"
)

## ==== 3) 工具函数 ====

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_-]+", "_", x)
}

sparse_column_variance <- function(x) {
  n <- nrow(x)
  if (n < 2L) return(rep(NA_real_, ncol(x)))
  mu <- Matrix::colMeans(x)
  mu2 <- Matrix::colMeans(x ^ 2)
  ans <- (mu2 - mu ^ 2) * n / (n - 1)
  ans[ans < 0 & ans > -1e-12] <- 0
  ans
}

prepare_source_matrix <- function(x, source_meta) {
  n <- nrow(x)
  nonzero_count <- Matrix::colSums(x != 0)
  nonzero_prop <- nonzero_count / n
  variance <- sparse_column_variance(x)

  usable <- which(
    nonzero_count >= min_feature_nonzero_pixels &
      nonzero_prop >= min_feature_nonzero_prop &
      is.finite(variance) & variance > 0
  )
  if (length(usable) < 10L) {
    stop("来源 ", unique(source_meta$source_id),
         " 的可用代谢 feature 少于 10 个。")
  }

  ranked <- usable[order(variance[usable], decreasing = TRUE)]
  selected <- ranked[seq_len(min(length(ranked), max_features_for_clustering))]

  feature_qc <- data.table::data.table(
    source_id = as.character(unique(source_meta$source_id)),
    matrix_feature_index = seq_along(feature_mz),
    mz = feature_mz,
    nonzero_pixels = as.integer(nonzero_count),
    nonzero_prop = as.numeric(nonzero_prop),
    logCPM_variance = as.numeric(variance),
    passed_source_filter = seq_along(feature_mz) %in% usable,
    selected_for_clustering = seq_along(feature_mz) %in% selected
  )

  ## 只把选中的数百个 feature 转为普通矩阵，避免整个稀疏矩阵变稠密。
  dense <- as.matrix(x[, selected, drop = FALSE])
  center <- colMeans(dense)
  sds <- apply(dense, 2L, stats::sd)
  keep <- is.finite(sds) & sds > 0
  dense <- dense[, keep, drop = FALSE]
  selected <- selected[keep]
  center <- center[keep]
  sds <- sds[keep]

  dense <- sweep(dense, 2L, center, FUN = "-")
  denominator <- switch(
    pca_scaling,
    center_only = rep(1, length(sds)),
    pareto = sqrt(sds),
    unit_variance = sds,
    stop("未知 pca_scaling：", pca_scaling)
  )
  dense <- sweep(dense, 2L, denominator, FUN = "/")

  if (isTRUE(include_qc_covariates_in_clustering)) {
    qc_x <- cbind(
      log_raw_tic = log1p(source_meta$raw_tic),
      detected_peak_count = source_meta$detected_peak_count,
      alignment_nmatch = source_meta$alignment_nmatch
    )
    qc_x <- scale(qc_x)
    qc_x[!is.finite(qc_x)] <- 0
    dense <- cbind(dense, qc_x)
  }

  list(
    x = dense,
    selected_feature_indices = selected,
    feature_qc = feature_qc
  )
}

run_pca <- function(x) {
  n_components <- min(max_pcs, nrow(x) - 1L, ncol(x) - 1L)
  if (n_components < 2L) stop("PCA 可计算的维度少于 2。")

  if (has_irlba && n_components < min(dim(x))) {
    fit <- irlba::prcomp_irlba(
      x, n = n_components, center = FALSE, scale. = FALSE
    )
    method <- "irlba::prcomp_irlba"
  } else {
    fit <- stats::prcomp(
      x, center = FALSE, scale. = FALSE, rank. = n_components
    )
    method <- "stats::prcomp"
  }

  total_variance <- sum(x ^ 2) / max(1, nrow(x) - 1L)
  variance_fraction <- fit$sdev ^ 2 / total_variance
  cumulative <- cumsum(variance_fraction)
  first_reaching_target <- which(cumulative >= pc_variance_target)[1]
  if (is.na(first_reaching_target)) first_reaching_target <- ncol(fit$x)

  n_pc_cluster <- max(min_pcs_for_clustering, first_reaching_target)
  n_pc_cluster <- min(n_pc_cluster, ncol(fit$x))

  list(
    fit = fit,
    method = method,
    variance_fraction = variance_fraction,
    n_pc_cluster = n_pc_cluster
  )
}

cluster_palette <- c(
  "1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3",
  "4" = "#E7298A", "5" = "#66A61E", "6" = "#E6AB02",
  "7" = "#A6761D", "8" = "#1F78B4", "9" = "#B15928",
  "10" = "#6A3D9A"
)

plot_cluster_map <- function(dt, sample_label, source_label, k_used,
                             n_pc_used) {
  dt_plot <- data.table::copy(dt)
  dt_plot[, cluster_plot := factor(cluster, levels = seq_len(k_used))]

  ggplot(dt_plot, aes(x = x, y = y, fill = cluster_plot)) +
    geom_tile(width = 1, height = 1) +
    scale_y_reverse() +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(
      values = cluster_palette[as.character(seq_len(k_used))],
      drop = FALSE
    ) +
    labs(
      title = paste0(sample_label, " / ", source_label),
      subtitle = paste0(
        "Metabolic-spectrum clustering: k = ", k_used,
        "; PCA dimensions = ", n_pc_used
      ),
      x = "x (pixel)", y = "y (pixel)", fill = "Cluster"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

plot_continuous_map <- function(dt, value_col, title, legend_title,
                                trans = "identity") {
  ggplot(dt, aes(x = x, y = y, fill = .data[[value_col]])) +
    geom_tile(width = 1, height = 1) +
    scale_y_reverse() +
    coord_fixed(expand = FALSE) +
    scale_fill_viridis_c(option = "magma", trans = trans, na.value = "grey90") +
    labs(title = title, x = "x (pixel)", y = "y (pixel)", fill = legend_title) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank()
    )
}

## ==== 4) 按来源聚类 ====

all_pixel_clusters <- vector("list", length(source_levels))
all_cluster_summary <- vector("list", length(source_levels))
all_pca_summary <- vector("list", length(source_levels))
all_feature_qc <- vector("list", length(source_levels))
all_top_loadings <- vector("list", length(source_levels))
all_cluster_plots <- vector("list", length(source_levels))

for (i in seq_along(source_levels)) {
  sid <- source_levels[i]
  row_idx <- which(as.character(pixel_meta$source_id) == sid)
  source_meta <- data.table::copy(pixel_meta[row_idx])
  sample_label <- as.character(unique(source_meta$sample_id))

  if (!length(row_idx)) stop("像素元数据中找不到来源：", sid)
  if (any(!is.finite(source_meta$x)) || any(!is.finite(source_meta$y))) {
    stop("来源 ", sid, " 存在缺失或非有限坐标。")
  }

  message(
    "[聚类] ", i, "/", length(source_levels), ": ",
    sample_label, " / ", sid, "；像素数 = ", length(row_idx)
  )

  source_matrix <- prepare_source_matrix(
    log_mat[row_idx, , drop = FALSE], source_meta
  )
  pca <- run_pca(source_matrix$x)
  pc_scores <- pca$fit$x[, seq_len(pca$n_pc_cluster), drop = FALSE]

  k_used <- min(k_clusters, nrow(pc_scores) - 1L)
  if (k_used < 2L) stop("来源 ", sid, " 的像素太少，不能聚类。")

  set.seed(random_seed + i)
  km <- tryCatch(
    stats::kmeans(
      pc_scores,
      centers = k_used,
      nstart = kmeans_nstart,
      iter.max = kmeans_iter_max,
      algorithm = "Hartigan-Wong"
    ),
    error = function(e) {
      message(
        "[提示] ", sid, " 的 Hartigan-Wong k-means 未成功（",
        conditionMessage(e), "），改用 Lloyd 算法。"
      )
      stats::kmeans(
        pc_scores,
        centers = k_used,
        nstart = kmeans_nstart,
        iter.max = kmeans_iter_max,
        algorithm = "Lloyd"
      )
    }
  )

  source_meta[, cluster := as.integer(km$cluster)]
  source_meta[, cluster_key := paste(source_id, cluster, sep = "::")]

  cluster_summary <- source_meta[, .(
    n_pixels = .N,
    pixel_fraction = .N / nrow(source_meta),
    x_min = min(x), x_max = max(x),
    y_min = min(y), y_max = max(y),
    median_raw_tic = as.numeric(stats::median(raw_tic, na.rm = TRUE)),
    median_detected_peaks = as.numeric(
      stats::median(detected_peak_count, na.rm = TRUE)
    ),
    median_alignment_matches = as.numeric(
      stats::median(alignment_nmatch, na.rm = TRUE)
    ),
    alignment_failure_rate = mean(alignment_failed, na.rm = TRUE)
  ), by = .(sample_id, diagnosis, source_id, cluster, cluster_key)]
  data.table::setorder(cluster_summary, cluster)

  variance_dt <- data.table::data.table(
    sample_id = sample_label,
    source_id = sid,
    PC = paste0("PC", seq_along(pca$variance_fraction)),
    variance_fraction = pca$variance_fraction,
    variance_percent = 100 * pca$variance_fraction,
    cumulative_variance_percent = 100 * cumsum(pca$variance_fraction),
    used_for_clustering = seq_along(pca$variance_fraction) <= pca$n_pc_cluster,
    pca_method = pca$method,
    pca_scaling = pca_scaling,
    n_pixels = nrow(source_meta),
    n_features_used = length(source_matrix$selected_feature_indices)
  )

  n_loading_pc <- min(5L, ncol(pca$fit$rotation))
  loading_dt <- data.table::rbindlist(lapply(seq_len(n_loading_pc), function(j) {
    loading <- pca$fit$rotation[, j]
    ## QC covariates（若开启）位于末尾，不是 m/z feature，故只读取前面的峰。
    n_mz_loading <- min(
      length(source_matrix$selected_feature_indices), length(loading)
    )
    loading <- loading[seq_len(n_mz_loading)]
    top <- order(abs(loading), decreasing = TRUE)[
      seq_len(min(30L, length(loading)))
    ]
    data.table::data.table(
      sample_id = sample_label,
      source_id = sid,
      PC = paste0("PC", j),
      loading_rank = seq_along(top),
      matrix_feature_index = source_matrix$selected_feature_indices[top],
      mz = feature_mz[source_matrix$selected_feature_indices[top]],
      loading = loading[top],
      abs_loading = abs(loading[top])
    )
  }))

  p_cluster <- plot_cluster_map(
    source_meta, sample_label, sid, k_used, pca$n_pc_cluster
  )

  file_stub <- safe_filename(sid)
  ggsave(
    file.path(map_dir, paste0(file_stub, "_cluster_map.png")),
    p_cluster, width = 7.2, height = 7.2, units = "in", dpi = save_png_dpi
  )
  ggsave(
    file.path(map_dir, paste0(file_stub, "_cluster_map.pdf")),
    p_cluster, width = 7.2, height = 7.2, units = "in"
  )

  if (isTRUE(save_auxiliary_maps)) {
    p_tic <- plot_continuous_map(
      source_meta, "raw_tic",
      paste0(sample_label, " / ", sid, " / raw TIC"),
      "Raw TIC", trans = "log10"
    )
    p_peak <- plot_continuous_map(
      source_meta, "detected_peak_count",
      paste0(sample_label, " / ", sid, " / detected peak count"),
      "Detected peaks"
    )
    ggsave(
      file.path(aux_map_dir, paste0(file_stub, "_raw_TIC.png")),
      p_tic, width = 7.2, height = 7.2, units = "in", dpi = save_png_dpi
    )
    ggsave(
      file.path(aux_map_dir, paste0(file_stub, "_detected_peaks.png")),
      p_peak, width = 7.2, height = 7.2, units = "in", dpi = save_png_dpi
    )
  }

  all_pixel_clusters[[i]] <- source_meta
  all_cluster_summary[[i]] <- cluster_summary
  all_pca_summary[[i]] <- variance_dt
  all_feature_qc[[i]] <- source_matrix$feature_qc
  all_top_loadings[[i]] <- loading_dt
  all_cluster_plots[[i]] <- p_cluster
}

## ==== 5) 汇总并保存 ====

pixel_clusters <- data.table::rbindlist(all_pixel_clusters, use.names = TRUE)
cluster_summary <- data.table::rbindlist(all_cluster_summary, use.names = TRUE)
pca_summary <- data.table::rbindlist(all_pca_summary, use.names = TRUE)
feature_qc <- data.table::rbindlist(all_feature_qc, use.names = TRUE)
top_loadings <- data.table::rbindlist(all_top_loadings, use.names = TRUE)

## 恢复原始矩阵行顺序，保证以后可以直接用 pixel_row 筛选矩阵。
data.table::setorder(pixel_clusters, pixel_row)

saveRDS(
  pixel_clusters,
  file.path(out_dir, "pixel_clusters_12sources.rds"),
  compress = "xz"
)
data.table::fwrite(
  pixel_clusters,
  file.path(out_dir, "pixel_clusters_12sources.csv")
)
data.table::fwrite(
  cluster_summary,
  file.path(out_dir, "cluster_summary_12sources.csv")
)
data.table::fwrite(
  pca_summary,
  file.path(out_dir, "source_PCA_variance.csv")
)
data.table::fwrite(
  feature_qc,
  file.path(out_dir, "source_clustering_feature_qc.csv")
)
data.table::fwrite(
  top_loadings,
  file.path(out_dir, "source_PCA_top_loadings.csv")
)

## 一个 14 页 PDF，便于连续检查全部来源。
grDevices::pdf(
  file.path(out_dir, "all_12sources_cluster_maps.pdf"),
  width = 7.2, height = 7.2, onefile = TRUE
)
for (p in all_cluster_plots) print(p)
grDevices::dev.off()

run_parameters <- list(
  input_dir = input_dir,
  out_dir = out_dir,
  source_ids_run = source_levels,
  k_clusters = k_clusters,
  max_features_for_clustering = max_features_for_clustering,
  min_feature_nonzero_prop = min_feature_nonzero_prop,
  min_feature_nonzero_pixels = min_feature_nonzero_pixels,
  max_pcs = max_pcs,
  min_pcs_for_clustering = min_pcs_for_clustering,
  pc_variance_target = pc_variance_target,
  pca_scaling = pca_scaling,
  include_qc_covariates_in_clustering = include_qc_covariates_in_clustering,
  kmeans_nstart = kmeans_nstart,
  kmeans_iter_max = kmeans_iter_max,
  random_seed = random_seed,
  irlba_available = has_irlba,
  session_info = utils::sessionInfo()
)
saveRDS(
  run_parameters,
  file.path(out_dir, "0B_parameters_and_sessionInfo.rds"),
  compress = "xz"
)

message("完成。当前脚本没有删除任何像素，也没有设置 is_tumor。")
message("请查看：", file.path(out_dir, "all_12sources_cluster_maps.pdf"))
message("Cluster 汇总：", file.path(out_dir, "cluster_summary_12sources.csv"))
message("逐像素标签：", file.path(out_dir, "pixel_clusters_12sources.rds"))
message("人工确认各 source_id 的肿瘤 cluster 后，再进入肿瘤像素筛选。")
