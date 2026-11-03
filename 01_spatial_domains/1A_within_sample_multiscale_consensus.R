## ============================================================================
## 1A_within_sample_multiscale_consensus.R
## 改进的 Part 1：每个物理样本独立进行多尺度、空间分块 bootstrap 共识分析。
##
## 主要目标：
## 1) k = 2:4 与多个空间权重同时评估；
## 2) 用空间分块 bootstrap 评估 PCA + 聚类全流程稳定性；
## 3) 输出稳定核心、过渡区、不确定区，而不是只给一次硬标签；
## 4) 输出 PC1/PC2 连续梯度图，避免把连续变化强行解释为离散亚区；
## 5) 聚类 discovery features 与 domain 验证 features 分离，减少循环论证；
## 6) 输入为上游已确定的 10 个样本，Cushing-3 从 imzML 读入阶段即未纳入；
##    NF-GnH-3 三个 imzML 是同一次测定的分块导出：
##    feature QC、标准化、PCA、聚类、bootstrap和指标均按一个物理样本合并计算；
##    仅空间邻接建图和绘图时按三块各自坐标分开。
## ============================================================================

required_packages <- c("data.table", "Matrix", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("缺少 R 包：", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2)
})
has_irlba <- requireNamespace("irlba", quietly = TRUE)

## ==== 1) 路径、开关和参数 ====
data_root <- "/home/gaga/huzhan"
input_dir <- file.path(data_root, "MSI_ALL_aligned_10samples_consensus_tuned225")
cluster_dir <- file.path(input_dir, "QC", "tumor_pixel_selection_0B_k4")
out_dir <- file.path(input_dir, "intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(out_dir, "sample_checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(
  logCPM = file.path(input_dir, "pixel_by_feature_matrix_aligned_logCPM.rds"),
  feature_mz = file.path(input_dir, "feature_mz_aligned_filtered.rds"),
  pixel_meta = file.path(input_dir, "pixel_metadata_sample_source_xy.rds"),
  manifest = file.path(input_dir, "sample_source_manifest.csv"),
  pixel_clusters = file.path(cluster_dir, "pixel_clusters_12sources.rds")
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))

tumor_cluster_map <- list(
  cushing_0 = c(3L, 4L), Cushing = c(1L, 2L), gh_1 = 2L,
  gh_2 = c(1L, 2L, 3L), GH = c(2L, 4L), NF_ACTH = c(1L, 3L),
  nf_cxx_1 = c(2L, 4L), nf_cxx_2 = c(3L, 4L),
  nf_cxx_3 = c(2L, 4L), NF_CXX = 3L,
  nf_cxx = c(2L, 4L), NF_GH = c(1L, 3L)
)
quick_test <- FALSE
parallel_samples <- TRUE
n_sample_workers <- 4L
resume_completed_samples <- TRUE
random_seed <- 20250908L
k_candidates <- 2:4
spatial_alpha_values <- c(0, 0.10, 0.20, 0.35)
n_spatial_bootstrap <- if (quick_test) 5L else 50L
spatial_block_width <- 10L
bootstrap_block_fraction <- 0.80
max_training_pixels <- 1500L
max_pcs <- 20L
n_pcs_for_domains <- 15L
reference_kmeans_nstart <- 30L
bootstrap_kmeans_nstart <- 5L
kmeans_iter_max <- 100L

min_prevalence <- 0.10
min_nonzero_pixels <- 20L
max_detection_qc_correlation <- 0.35
min_source_variance <- 1e-6
minimum_discovery_features <- 100L
maximum_discovery_features <- 150L
minimum_validation_features <- 80L

minimum_reference_domain_fraction_threshold <- 0.08
minimum_bootstrap_median_ari <- 0.60
minimum_mean_consensus_probability <- 0.65
parsimony_score_tolerance <- 0.05
stable_core_probability <- 0.80
transition_probability <- 0.60
stable_core_min_reference_margin <- 0.15

## 空间组织类型在“物理样本”层面定义，而不是在 imzML 来源层面定义。
## 多分块样本（NF-GnH-3）在合并的像素和块内邻接边上直接计算样本指标，
## 然后得到一个共同的样本类型；三个导出块仅保留各自坐标供建图和绘图。
spatial_coherence_moran_threshold <- 0.35
spatial_coherence_mixing_threshold <- 0.60
spatial_coherence_entropy_threshold <- 0.55
spatial_coherence_required_criteria <- 3L

## "automatic"：按稳定性、共识、CH和简约性选择 k。
## "manual"：使用下方 manual_k_override；可用于对特定样本做敏感性分析。
k_selection_mode <- "automatic"
manual_k_override <- c(
  "GH-1" = 2L, "GH-2" = 4L, "GH-3" = 2L, "NF-GH" = 4L,
  "NF-ACTH" = 3L, "Cushing-1" = 4L, "Cushing-2" = 2L,
  "NF-GnH-1" = 2L, "NF-GnH-2" = 3L, "NF-GnH-3" = 3L
)

run_hierarchical_refinement_screen <- TRUE
hierarchical_min_core_pixels <- 250L
hierarchical_min_subdomain_fraction <- 0.15
hierarchical_min_initialization_ari <- 0.80
hierarchical_min_spatial_ari <- 0.60
save_png_dpi <- 240

analysis_signature <- paste(
  "part1_consensus_v2_combined_export_partitions", quick_test, n_spatial_bootstrap,
  paste(k_candidates, collapse = ","), paste(spatial_alpha_values, collapse = ","),
  spatial_block_width, bootstrap_block_fraction, max_training_pixels,
  maximum_discovery_features, k_selection_mode, random_seed, sep = "|")

## ==== 2) 基础函数 ====
safe_filename <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

safe_z <- function(x) {
  x <- as.numeric(x); s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-12) return(rep(0, length(x)))
  z <- (x - mean(x, na.rm = TRUE)) / s; z[!is.finite(z)] <- 0; z
}

sparse_column_variance <- function(x) {
  n <- nrow(x); if (n < 2L) return(rep(NA_real_, ncol(x)))
  mu <- Matrix::colMeans(x)
  ans <- (Matrix::colMeans(x ^ 2) - mu ^ 2) * n / (n - 1)
  ans[ans < 0 & ans > -1e-12] <- 0; as.numeric(ans)
}

qc_depth_score <- function(dt) {
  safe_z((safe_z(log1p(dt$detected_peak_count)) + safe_z(dt$alignment_nmatch)) / sqrt(2))
}

technical_design <- function(dt) {
  zp <- safe_z(log1p(dt$detected_peak_count)); za <- safe_z(dt$alignment_nmatch)
  za_resid <- stats::lm.fit(cbind(1, zp), za)$residuals
  cbind(intercept = 1, detected_peak_count = zp,
    alignment_residual = safe_z(za_resid), raw_tic = safe_z(log1p(dt$raw_tic)))
}

residualize_and_scale <- function(x, dt, cap = 5) {
  x <- as.matrix(x); qfit <- qr(technical_design(dt))
  q <- qr.Q(qfit)[, seq_len(qfit$rank), drop = FALSE]
  x <- x - q %*% crossprod(q, x)
  x <- sweep(x, 2L, colMeans(x), "-")
  s <- sqrt(colSums(x ^ 2) / pmax(1, nrow(x) - 1L))
  s[!is.finite(s) | s < sqrt(min_source_variance)] <- 1
  x <- sweep(x, 2L, s, "/"); x[!is.finite(x)] <- 0
  x[x > cap] <- cap; x[x < -cap] <- -cap; x
}

binary_qc_correlation <- function(x, q) {
  q <- safe_z(q); n <- nrow(x); p <- as.numeric(Matrix::colMeans(x != 0))
  ans <- (as.numeric(Matrix::crossprod(x != 0, q)) / pmax(n, 1L)) /
    sqrt(p * (1 - p))
  ans[!is.finite(ans)] <- 0; pmax(-1, pmin(1, ans))
}

make_spatial_edges <- function(x, y, radius = 1L) {
  key <- paste(x, y, sep = "::")
  if (anyDuplicated(key)) stop("source 内存在重复坐标。")
  offsets <- data.table::CJ(dx = -radius:radius, dy = -radius:radius)[
    !(dx == 0L & dy == 0L) & (dx > 0L | (dx == 0L & dy > 0L))]
  z <- lapply(seq_len(nrow(offsets)), function(i) {
    to <- match(paste(x + offsets$dx[i], y + offsets$dy[i], sep = "::"), key)
    from <- which(!is.na(to)); if (!length(from)) return(NULL)
    cbind(from = from, to = to[from])
  })
  z <- Filter(Negate(is.null), z)
  if (!length(z)) return(matrix(integer(), ncol = 2L))
  ans <- do.call(rbind, z); storage.mode(ans) <- "integer"; ans
}

combine_partition_edges <- function(partition_rows, partition_edges) {
  z <- lapply(names(partition_rows), function(sid) {
    lr <- partition_rows[[sid]]; edge <- partition_edges[[sid]]
    if (is.null(edge) || !nrow(edge)) return(NULL)
    cbind(from = lr[edge[, 1]], to = lr[edge[, 2]])
  })
  z <- Filter(Negate(is.null), z)
  if (!length(z)) return(matrix(integer(), ncol = 2L,
    dimnames = list(NULL, c("from", "to"))))
  ans <- do.call(rbind, z); storage.mode(ans) <- "integer"; ans
}

smooth_scores <- function(scores, edges, alpha) {
  if (alpha <= 0 || !nrow(edges)) return(scores)
  n <- nrow(scores)
  a <- Matrix::sparseMatrix(i = c(edges[, 1], edges[, 2]),
    j = c(edges[, 2], edges[, 1]), x = 1, dims = c(n, n))
  degree <- as.numeric(Matrix::rowSums(a)); nm <- as.matrix(a %*% scores)
  ok <- degree > 0; nm[ok, ] <- nm[ok, , drop = FALSE] / degree[ok]
  nm[!ok, ] <- scores[!ok, , drop = FALSE]
  (1 - alpha) * scores + alpha * nm
}

fit_pca <- function(x) {
  npc <- min(max_pcs, nrow(x) - 1L, ncol(x) - 1L)
  if (npc < 2L) stop("PCA维度不足。")
  if (has_irlba && npc < min(dim(x))) {
    irlba::prcomp_irlba(x, n = npc, center = TRUE, scale. = FALSE)
  } else stats::prcomp(x, center = TRUE, scale. = FALSE, rank. = npc)
}

safe_kmeans <- function(x, k, nstart, seed) {
  set.seed(seed)
  tryCatch(stats::kmeans(x, centers = k, nstart = nstart,
    iter.max = kmeans_iter_max, algorithm = "Hartigan-Wong"),
    error = function(e) stats::kmeans(x, centers = k, nstart = nstart,
      iter.max = kmeans_iter_max, algorithm = "Lloyd"))
}

nearest_centroid <- function(x, centers) {
  d2 <- outer(rowSums(x ^ 2), rowSums(centers ^ 2), "+") - 2 * tcrossprod(x, centers)
  d2[d2 < 0 & d2 > -1e-10] <- 0
  best <- max.col(-d2, ties.method = "first")
  d1 <- d2[cbind(seq_len(nrow(d2)), best)]
  d2[cbind(seq_len(nrow(d2)), best)] <- Inf
  d_second <- apply(d2, 1L, min)
  list(cluster = best, margin = pmax(0, pmin(1,
    (d_second - d1) / pmax(d_second, .Machine$double.eps))))
}

adjusted_rand_index <- function(a, b) {
  tab <- table(a, b); c2 <- function(z) z * (z - 1) / 2; n <- sum(tab)
  if (n < 2L) return(NA_real_)
  obs <- sum(c2(tab)); rs <- sum(c2(rowSums(tab))); cs <- sum(c2(colSums(tab)))
  expected <- rs * cs / c2(n); den <- (rs + cs) / 2 - expected
  if (abs(den) < .Machine$double.eps) return(1)
  (obs - expected) / den
}

all_permutations <- function(x) {
  if (length(x) <= 1L) return(list(x))
  unlist(lapply(seq_along(x), function(i) {
    lapply(all_permutations(x[-i]), function(z) c(x[i], z))
  }), recursive = FALSE)
}

align_numeric_labels <- function(reference, target, k) {
  perms <- all_permutations(seq_len(k))
  agreement <- vapply(perms, function(p) mean(p[target] == reference), numeric(1))
  best <- perms[[which.max(agreement)]]
  list(labels = best[target], agreement = max(agreement), map = best)
}

normalized_entropy_from_votes <- function(votes) {
  p <- votes / pmax(1, rowSums(votes)); z <- p
  z[z <= 0] <- 1
  -rowSums(p * log(z)) / log(ncol(votes))
}

continuous_morans_i <- function(value, edges) {
  if (!nrow(edges)) return(NA_real_)
  z <- value - mean(value); den <- sum(z ^ 2)
  if (den <= 0) return(NA_real_)
  length(value) / nrow(edges) * sum(z[edges[, 1]] * z[edges[, 2]]) / den
}

mean_domain_morans_i <- function(labels, edges, k) {
  if (!nrow(edges)) return(NA_real_)
  mean(vapply(seq_len(k), function(d) {
    z <- as.numeric(labels == d); z <- z - mean(z); den <- sum(z ^ 2)
    if (den <= 0) return(NA_real_)
    length(labels) / nrow(edges) * sum(z[edges[, 1]] * z[edges[, 2]]) / den
  }, numeric(1)), na.rm = TRUE)
}

local_entropy_values <- function(labels, edges, k) {
  n <- length(labels)
  directed <- rbindlist(list(
    data.table(i = seq_len(n), neighbor = labels),
    data.table(i = edges[, 1], neighbor = labels[edges[, 2]]),
    data.table(i = edges[, 2], neighbor = labels[edges[, 1]])))
  counts <- directed[, .N, by = .(i, neighbor)]
  counts[, p := N / sum(N), by = i]
  ent <- counts[, .(entropy = -sum(p * log(p)) / log(k)), by = i]
  ans <- numeric(n); ans[ent$i] <- ent$entropy; ans
}

count_spatial_patches <- function(labels, edges) {
  n <- length(labels); parent <- seq_len(n); rank_ <- integer(n)
  root <- function(a) {
    while (parent[a] != a) {
      parent[a] <<- parent[parent[a]]; a <- parent[a]
    }
    a
  }
  unite <- function(a, b) {
    ra <- root(a); rb <- root(b); if (ra == rb) return(invisible(NULL))
    if (rank_[ra] < rank_[rb]) parent[ra] <<- rb else if (rank_[ra] > rank_[rb]) {
      parent[rb] <<- ra
    } else {
      parent[rb] <<- ra; rank_[ra] <<- rank_[ra] + 1L
    }
  }
  if (nrow(edges)) {
    same <- edges[labels[edges[, 1]] == labels[edges[, 2]], , drop = FALSE]
    if (nrow(same)) for (i in seq_len(nrow(same))) unite(same[i, 1], same[i, 2])
  }
  roots <- vapply(seq_len(n), root, integer(1))
  nrow(unique(data.table(domain = labels, root = roots)))
}

rank_correlation <- function(x, y) suppressWarnings(stats::cor(
  rank(x, ties.method = "average"), rank(y, ties.method = "average"),
  use = "complete.obs"))

eta_squared <- function(labels, x) {
  x <- as.numeric(x); ok <- is.finite(x) & !is.na(labels)
  x <- x[ok]; labels <- labels[ok]
  if (length(x) < 3L || stats::var(x) <= 0) return(NA_real_)
  mu <- mean(x); total <- sum((x - mu) ^ 2)
  sum(vapply(split(x, labels), function(z) length(z) * (mean(z) - mu) ^ 2,
    numeric(1))) / total
}

domain_spatial_metrics <- function(labels, edges, k) {
  p <- as.numeric(prop.table(table(factor(labels, levels = seq_len(k)))))
  simpson <- 1 - sum(p ^ 2)
  hetero <- if (nrow(edges)) mean(labels[edges[, 1]] != labels[edges[, 2]]) else NA_real_
  n_patches <- count_spatial_patches(labels, edges)
  local_entropy <- local_entropy_values(labels, edges, k)
  data.table(normalized_shannon = {
      q <- p[p > 0]; -sum(q * log(q)) / log(k)
    }, maximum_domain_fraction = max(p), heterotypic_edge_fraction = hetero,
    mixing_ratio_to_random = if (simpson > 0) hetero / simpson else NA_real_,
    total_number_of_patches = n_patches,
    patches_per_1000_tumor_pixels = 1000 * n_patches / length(labels),
    mean_domain_morans_i = mean_domain_morans_i(labels, edges, k),
    mean_local_entropy_r1 = mean(local_entropy))
}

make_block_ids <- function(dt, width) {
  paste(dt$source_id,
    floor((dt$x - min(dt$x)) / width),
    floor((dt$y - min(dt$y)) / width), sep = "::")
}

sample_spatial_blocks <- function(dt, fraction, seed) {
  ## block_id 已包含导出分块 source_id，所以不同坐标系不会冲突。
  ## 但所有 block 在同一物理样本内一起抽样，不把导出分块当成批次。
  set.seed(seed)
  blocks <- unique(dt$block_id)
  nb <- max(1L, ceiling(length(blocks) * fraction))
  keep_blocks <- sample(blocks, nb, replace = FALSE)
  which(dt$block_id %chin% keep_blocks)
}

sample_training_indices <- function(candidate, cap_total, seed) {
  set.seed(seed); candidate <- sort(unique(as.integer(candidate)))
  if (length(candidate) > cap_total) candidate <- sample(candidate, cap_total)
  sort(candidate)
}

project_pca <- function(x, pca, npc) {
  sweep(x, 2L, pca$center, "-") %*% pca$rotation[, seq_len(npc), drop = FALSE]
}

plot_tile <- function(dt, fill_col, title, subtitle, discrete = TRUE, limits = NULL) {
  palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
  p <- ggplot(dt, aes(x = x, y = y, fill = .data[[fill_col]])) + geom_tile() +
    scale_y_reverse() + coord_fixed(expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = "x (pixel)", y = "y (pixel)") +
    theme_minimal(base_size = 12) + theme(panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5))
  if (discrete) p + scale_fill_manual(values = palette, drop = FALSE, na.value = "grey85")
  else p + scale_fill_viridis_c(option = "magma", limits = limits)
}

## ==== 3) 读入、肿瘤像素和 feature QC ====
message("[读入] 对齐矩阵、像素元数据和0B肿瘤标注。")
log_mat_all <- readRDS(input_files[["logCPM"]])
if (!inherits(log_mat_all, "dgCMatrix")) log_mat_all <- as(log_mat_all, "dgCMatrix")
feature_mz <- as.numeric(readRDS(input_files[["feature_mz"]]))
pixel_meta_all <- as.data.table(readRDS(input_files[["pixel_meta"]]))
manifest <- fread(input_files[["manifest"]])
cluster_dt <- as.data.table(readRDS(input_files[["pixel_clusters"]]))
stopifnot(nrow(log_mat_all) == nrow(pixel_meta_all),
  ncol(log_mat_all) == length(feature_mz), nrow(cluster_dt) == nrow(pixel_meta_all))

pixel_meta_all <- merge(pixel_meta_all,
  cluster_dt[, .(pixel_row = as.integer(pixel_row), cluster_0B_k4 = as.integer(cluster))],
  by = "pixel_row", all.x = TRUE, sort = FALSE)
setorder(pixel_meta_all, pixel_row)
if (anyNA(pixel_meta_all$cluster_0B_k4)) stop("0B cluster缺失。")
source_levels <- as.character(manifest$source_id)
sample_levels <- unique(as.character(manifest$sample_id))
if (length(source_levels) != 12L || length(sample_levels) != 10L) stop("预期12来源、10样本。")

pixel_meta_all[, is_tumor := FALSE]
for (sid in names(tumor_cluster_map)) pixel_meta_all[
  source_id == sid & cluster_0B_k4 %in% tumor_cluster_map[[sid]], is_tumor := TRUE]
tumor_global_rows <- which(pixel_meta_all$is_tumor)
tumor_meta <- copy(pixel_meta_all[tumor_global_rows])
tumor_meta[, `:=`(original_pixel_row = pixel_row, tumor_pixel_row = seq_len(.N))]
tumor_log <- log_mat_all[tumor_global_rows, , drop = FALSE]
rm(log_mat_all); invisible(gc())

source_rows <- split(seq_len(nrow(tumor_meta)), as.character(tumor_meta$source_id))[source_levels]
sample_rows <- split(seq_len(nrow(tumor_meta)), as.character(tumor_meta$sample_id))[sample_levels]
tumor_meta[, block_id := make_block_ids(.SD, spatial_block_width)]

message("[feature QC] 按物理样本计算出现率、方差和技术相关性。")
sample_prevalence <- sample_nnz <- sample_variance <- sample_qc_corr <- matrix(
  NA_real_, length(sample_levels), length(feature_mz), dimnames = list(sample_levels, NULL))
for (sample_i in sample_levels) {
  rows <- sample_rows[[sample_i]]; x <- tumor_log[rows, , drop = FALSE]
  sample_nnz[sample_i, ] <- Matrix::colSums(x != 0)
  sample_prevalence[sample_i, ] <- sample_nnz[sample_i, ] / length(rows)
  sample_variance[sample_i, ] <- sparse_column_variance(x)
  sample_qc_corr[sample_i, ] <- binary_qc_correlation(
    x, qc_depth_score(tumor_meta[rows]))
}

## ==== 4) 每个样本的多尺度 bootstrap 共识 ====
analyze_one_sample <- function(si) {
  source_metrics_i <- technical_qc_i <- list()
  plots_domain_i <- plots_core_i <- list()
  sample_i <- sample_levels[si]
  sources_i <- as.character(remaining_manifest[sample_id == sample_i, source_id])
  rows_i <- unlist(source_rows[sources_i], use.names = FALSE)
  meta_i <- copy(tumor_meta[rows_i]); n_i <- nrow(meta_i)
  local_source_rows <- split(seq_len(n_i), as.character(meta_i$source_id))[sources_i]
  edges_i <- lapply(sources_i, function(sid) {
    z <- local_source_rows[[sid]]
    make_spatial_edges(meta_i$x[z], meta_i$y[z], 1L)
  }); names(edges_i) <- sources_i
  combined_edges_i <- combine_partition_edges(local_source_rows, edges_i)
  message("[样本] ", sample_i, " / ", paste(sources_i, collapse = ", "),
    "；像素=", n_i)

  ## 一个物理样本只做一次 feature QC。NF-GnH-3 三个导出块在这里已合并。
  sample_prev <- sample_prevalence[sample_i, ]
  sample_nonzero <- sample_nnz[sample_i, ]
  sample_var <- sample_variance[sample_i, ]
  sample_abs_qc_cor <- abs(sample_qc_corr[sample_i, ])
  strict <- sample_prev >= min_prevalence & sample_nonzero >= min_nonzero_pixels &
    is.finite(sample_var) & sample_var > 0 &
    sample_abs_qc_cor <= max_detection_qc_correlation
  eligible <- which(strict)
  if (length(eligible) < minimum_discovery_features + minimum_validation_features) {
    stop(sample_i, "：严格合格feature不足，无法独立划分discovery/validation。")
  }
  ranking_score <- sample_var[eligible] / (1 + 3 * sample_abs_qc_cor[eligible])
  ranked <- eligible[order(ranking_score, decreasing = TRUE)]
  discovery <- head(ranked[seq(1L, length(ranked), by = 2L)], maximum_discovery_features)
  validation <- setdiff(eligible, discovery)
  if (length(discovery) < minimum_discovery_features ||
      length(validation) < minimum_validation_features) stop(sample_i, " feature拆分后数量不足。")
  message("[", sample_i, "] discovery feature=", length(discovery),
    "；independent validation feature=", length(validation))
  feature_qc_i <- data.table(sample_id = sample_i,
    matrix_feature_index = seq_along(feature_mz), mz = feature_mz,
    physical_sample_prevalence = sample_prev,
    physical_sample_nonzero_pixels = sample_nonzero,
    physical_sample_variance = sample_var,
    physical_sample_abs_detection_qc_correlation = sample_abs_qc_cor,
    feature_qc_unit = "physical_sample_combined",
    passed_strict_filter = strict,
    feature_role = fifelse(seq_along(feature_mz) %in% discovery, "discovery",
      fifelse(seq_along(feature_mz) %in% validation, "validation", "not_used")))

  ## 整个物理样本共同残差化和标准化；导出分块不作为批次。
  x_all <- residualize_and_scale(
    tumor_log[rows_i, discovery, drop = FALSE], meta_i)

  full_train <- sample_training_indices(seq_len(n_i), max_training_pixels,
    random_seed + 10000L * si)
  total_training_variance <- sum(apply(x_all[full_train, , drop = FALSE], 2L, stats::var))
  pca_ref <- fit_pca(x_all[full_train, , drop = FALSE])
  npc <- min(n_pcs_for_domains, ncol(pca_ref$rotation))
  score_ref <- project_pca(x_all, pca_ref, npc)
  message("[", sample_i, "] 参考PCA完成；PC=", npc,
    "；开始空间分块bootstrap。")
  reference <- list(); vote_arrays <- list(); candidate_rows <- list()
  for (k in k_candidates) {
    fit <- safe_kmeans(score_ref[full_train, , drop = FALSE], k,
      reference_kmeans_nstart, random_seed + 20000L * si + k)
    assigned <- nearest_centroid(score_ref, fit$centers)
    ord <- order(fit$centers[, 1]); map <- integer(k); map[ord] <- seq_len(k)
    labels <- map[assigned$cluster]
    centers <- fit$centers[ord, , drop = FALSE]
    ch <- (fit$betweenss / (k - 1)) /
      (fit$tot.withinss / max(1, nrow(score_ref[full_train, , drop = FALSE]) - k))
    reference[[as.character(k)]] <- list(labels = labels,
      margin = assigned$margin, centers = centers, ch = ch,
      minimum_domain_fraction = min(prop.table(table(labels[full_train]))))
    vote_arrays[[as.character(k)]] <- matrix(0L, n_i, k)
  }

  bootstrap_rows <- list(); br <- 0L
  for (b in seq_len(n_spatial_bootstrap)) {
    if (b == 1L || b %% 10L == 0L || b == n_spatial_bootstrap) {
      message("[", sample_i, "] bootstrap ", b, "/", n_spatial_bootstrap)
    }
    block_candidates <- sample_spatial_blocks(meta_i, bootstrap_block_fraction,
      random_seed + 1000000L * si + b)
    boot_train <- sample_training_indices(block_candidates, max_training_pixels,
      random_seed + 1100000L * si + b)
    pca_b <- fit_pca(x_all[boot_train, , drop = FALSE])
    npc_b <- min(npc, ncol(pca_b$rotation))
    score_b0 <- project_pca(x_all, pca_b, npc_b)
    for (alpha in spatial_alpha_values) {
      score_b <- score_b0
      for (sid in sources_i) {
        lr <- local_source_rows[[sid]]
        score_b[lr, ] <- smooth_scores(score_b0[lr, , drop = FALSE],
          edges_i[[sid]], alpha)
      }
      for (k in k_candidates) {
        fit_b <- safe_kmeans(score_b[boot_train, , drop = FALSE], k,
          bootstrap_kmeans_nstart,
          random_seed + 2000000L * si + 10000L * b + 100L * k + round(alpha * 100))
        lab_b <- nearest_centroid(score_b, fit_b$centers)$cluster
        aligned <- align_numeric_labels(reference[[as.character(k)]]$labels, lab_b, k)
        lab_b <- aligned$labels
        votes_k <- vote_arrays[[as.character(k)]]
        vote_pos <- cbind(seq_len(n_i), lab_b)
        votes_k[vote_pos] <- votes_k[vote_pos] + 1L
        vote_arrays[[as.character(k)]] <- votes_k
        br <- br + 1L
        bootstrap_rows[[br]] <- data.table(sample_id = sample_i, bootstrap = b,
          spatial_alpha = alpha, k = k, n_training_pixels = length(boot_train),
          reference_ARI = adjusted_rand_index(reference[[as.character(k)]]$labels, lab_b),
          best_exact_label_agreement = aligned$agreement,
          minimum_domain_fraction = min(prop.table(table(lab_b))))
      }
    }
  }
  boot_dt <- rbindlist(bootstrap_rows)

  for (k in k_candidates) {
    votes <- as.matrix(vote_arrays[[as.character(k)]])
    prob <- apply(votes, 1L, max) / rowSums(votes)
    ent <- normalized_entropy_from_votes(votes)
    lab <- max.col(votes, ties.method = "first")
    k_value <- as.integer(k)
    z <- boot_dt[boot_dt$k == k_value]
    candidate_rows[[as.character(k)]] <- data.table(sample_id = sample_i, k = k,
      calinski_harabasz = reference[[as.character(k)]]$ch,
      minimum_reference_domain_fraction = reference[[as.character(k)]]$minimum_domain_fraction,
      median_bootstrap_ARI_alpha0 = median(z[spatial_alpha == 0, reference_ARI]),
      minimum_bootstrap_ARI_alpha0 = min(z[spatial_alpha == 0, reference_ARI]),
      median_multiscale_bootstrap_ARI = median(z$reference_ARI),
      minimum_multiscale_bootstrap_ARI = min(z$reference_ARI),
      mean_consensus_probability = mean(prob), median_consensus_probability = median(prob),
      fraction_consensus_probability_ge_0.8 = mean(prob >= 0.8),
      median_consensus_entropy = median(ent),
      minimum_consensus_domain_fraction = min(prop.table(table(lab))))
    reference[[as.character(k)]]$consensus_labels <- lab
    reference[[as.character(k)]]$consensus_probability <- prob
    reference[[as.character(k)]]$consensus_entropy <- ent
  }
  candidate_dt <- rbindlist(candidate_rows)
  candidate_dt[, ch_scaled := {
    rg <- range(calinski_harabasz)
    if (diff(rg) <= 0) rep(1, .N) else (calinski_harabasz - rg[1]) / diff(rg)
  }]
  candidate_dt[, selection_score :=
    0.35 * median_bootstrap_ARI_alpha0 +
    0.25 * median_multiscale_bootstrap_ARI +
    0.25 * mean_consensus_probability + 0.15 * ch_scaled]
  candidate_dt[, passes_minimum_requirements :=
    minimum_reference_domain_fraction >= minimum_reference_domain_fraction_threshold &
    median_bootstrap_ARI_alpha0 >= minimum_bootstrap_median_ari &
    mean_consensus_probability >= minimum_mean_consensus_probability]

  if (k_selection_mode == "manual") {
    chosen_k <- unname(manual_k_override[[sample_i]])
  } else {
    pool <- candidate_dt[passes_minimum_requirements == TRUE]
    if (!nrow(pool)) pool <- copy(candidate_dt)
    best_score <- max(pool$selection_score)
    chosen_k <- min(pool[selection_score >= best_score - parsimony_score_tolerance, k])
  }
  candidate_dt[, `:=`(is_selected_k = k == chosen_k,
    k_selection_mode_used = k_selection_mode)]
  message("[", sample_i, "] 自动选择 k=", chosen_k,
    "；开始稳定核心、验证feature和绘图。")

  sel <- reference[[as.character(chosen_k)]]
  final_lab <- sel$consensus_labels
  final_prob <- sel$consensus_probability
  final_entropy <- sel$consensus_entropy
  final_margin <- sel$margin
  reference_label_agreement <- final_lab == sel$labels
  confidence_class <- fifelse(final_prob >= stable_core_probability &
      reference_label_agreement &
      final_margin >= stable_core_min_reference_margin, "stable_core",
    fifelse(final_prob >= transition_probability, "transition", "uncertain"))
  stable_core <- confidence_class == "stable_core"
  final_local_entropy <- numeric(n_i)
  for (sid in sources_i) {
    lr <- local_source_rows[[sid]]
    final_local_entropy[lr] <- local_entropy_values(final_lab[lr], edges_i[[sid]], chosen_k)
  }

  pixels <- copy(meta_i)
  pixels[, `:=`(
    selected_k = chosen_k,
    consensus_domain = factor(paste0("D", final_lab), levels = paste0("D", seq_len(chosen_k))),
    consensus_probability = final_prob,
    consensus_entropy = final_entropy,
    reference_assignment_margin = final_margin,
    consensus_matches_full_reference = reference_label_agreement,
    confidence_class = factor(confidence_class,
      levels = c("stable_core", "transition", "uncertain")),
    stable_core_pixel = stable_core,
    local_domain_entropy_r1 = final_local_entropy,
    reference_PC1 = score_ref[, 1], reference_PC2 = score_ref[, 2]
  )]

  selected_qc <- candidate_dt[k == chosen_k]
  structure_type <- if (mean(stable_core) < 0.50 ||
      selected_qc$median_bootstrap_ARI_alpha0 < 0.60) {
    "continuous_or_unstable"
  } else if (chosen_k == 2L) "two_major_states" else "discrete_domains"
  model_summary_i <- data.table(sample_id = sample_i,
    n_export_partitions = length(sources_i), n_tumor_pixels = n_i,
    n_discovery_features = length(discovery), n_validation_features = length(validation),
    selected_k = chosen_k, structure_type = structure_type,
    median_bootstrap_ARI_alpha0 = selected_qc$median_bootstrap_ARI_alpha0,
    median_multiscale_bootstrap_ARI = selected_qc$median_multiscale_bootstrap_ARI,
    mean_consensus_probability = selected_qc$mean_consensus_probability,
    stable_core_fraction = mean(stable_core), transition_fraction = mean(confidence_class == "transition"),
    uncertain_fraction = mean(confidence_class == "uncertain"),
    PC1_variance_fraction = pca_ref$sdev[1]^2 / total_training_variance,
    PC2_variance_fraction = pca_ref$sdev[2]^2 / total_training_variance)

  ## 主空间指标在物理样本上直接计算。各块内有邻接边，块与块之间不连边。
  sample_metrics_i <- domain_spatial_metrics(final_lab, combined_edges_i, chosen_k)
  sample_metrics_i[, `:=`(
    sample_id = sample_i,
    diagnosis = as.character(unique(meta_i$diagnosis))[1L],
    n_export_partitions = length(sources_i),
    export_partition_ids = paste(sources_i, collapse = ", "),
    selected_k = chosen_k,
    n_tumor_pixels = n_i,
    stable_core_fraction = mean(stable_core),
    PC1_morans_i = continuous_morans_i(score_ref[, 1], combined_edges_i),
    PC2_morans_i = continuous_morans_i(score_ref[, 2], combined_edges_i)
  )]
  sample_technical_qc_i <- data.table(
    sample_id = sample_i,
    diagnosis = as.character(unique(meta_i$diagnosis))[1L],
    n_export_partitions = length(sources_i),
    export_partition_ids = paste(sources_i, collapse = ", "),
    selected_k = chosen_k,
    pc1_rho_detected_peak_count = rank_correlation(score_ref[, 1],
      meta_i$detected_peak_count),
    pc1_rho_alignment_nmatch = rank_correlation(score_ref[, 1],
      meta_i$alignment_nmatch),
    pc1_rho_raw_tic = rank_correlation(score_ref[, 1], meta_i$raw_tic),
    domain_eta2_detected_peak_count = eta_squared(final_lab,
      meta_i$detected_peak_count),
    domain_eta2_alignment_nmatch = eta_squared(final_lab,
      meta_i$alignment_nmatch),
    domain_eta2_log_raw_tic = eta_squared(final_lab, log1p(meta_i$raw_tic))
  )
  sample_technical_qc_i[, technical_confounding_flag :=
    abs(pc1_rho_detected_peak_count) > 0.30 |
    abs(pc1_rho_alignment_nmatch) > 0.30 |
    domain_eta2_detected_peak_count > 0.15 |
    domain_eta2_alignment_nmatch > 0.15]

  ## 分块级指标只作为空间一致性补充，不用于定义样本类型。
  for (sid in sources_i) {
    lr <- local_source_rows[[sid]]; lab_s <- final_lab[lr]; edge_s <- edges_i[[sid]]
    met <- domain_spatial_metrics(lab_s, edge_s, chosen_k)
    met[, `:=`(sample_id = sample_i, diagnosis = as.character(unique(meta_i$diagnosis[lr])),
      source_id = sid, selected_k = chosen_k, n_tumor_pixels = length(lr),
      stable_core_fraction = mean(stable_core[lr]),
      PC1_morans_i = continuous_morans_i(score_ref[lr, 1], edge_s),
      PC2_morans_i = continuous_morans_i(score_ref[lr, 2], edge_s))]
    source_metrics_i[[length(source_metrics_i) + 1L]] <- met
    tech <- data.table(sample_id = sample_i,
      diagnosis = as.character(unique(meta_i$diagnosis[lr])), source_id = sid,
      selected_k = chosen_k,
      pc1_rho_detected_peak_count = rank_correlation(score_ref[lr, 1],
        meta_i$detected_peak_count[lr]),
      pc1_rho_alignment_nmatch = rank_correlation(score_ref[lr, 1],
        meta_i$alignment_nmatch[lr]),
      pc1_rho_raw_tic = rank_correlation(score_ref[lr, 1], meta_i$raw_tic[lr]),
      domain_eta2_detected_peak_count = eta_squared(lab_s,
        meta_i$detected_peak_count[lr]),
      domain_eta2_alignment_nmatch = eta_squared(lab_s,
        meta_i$alignment_nmatch[lr]),
      domain_eta2_log_raw_tic = eta_squared(lab_s, log1p(meta_i$raw_tic[lr])))
    tech[, technical_confounding_flag :=
      abs(pc1_rho_detected_peak_count) > 0.30 |
      abs(pc1_rho_alignment_nmatch) > 0.30 |
      domain_eta2_detected_peak_count > 0.15 |
      domain_eta2_alignment_nmatch > 0.15]
    technical_qc_i[[length(technical_qc_i) + 1L]] <- tech
  }

  ## 独立验证 feature：主效应在整个物理样本合并计算。
  ## 导出分块级效应另存为一致性质控，不代替样本级主效应。
  validation_x_all <- residualize_and_scale(
    tumor_log[rows_i, validation, drop = FALSE], meta_i)
  validation_effects <- discovery_effects <- list()
  validation_partition_effects <- discovery_partition_effects <- list()
  ve <- de <- vpe <- dpe <- 0L
  for (d in seq_len(chosen_k)) {
    inside <- which(stable_core & final_lab == d)
    outside <- which(stable_core & final_lab != d)
    if (length(inside) < 20L || length(outside) < 20L) next
    ve <- ve + 1L
    validation_effects[[ve]] <- data.table(sample_id = sample_i,
      diagnosis = as.character(unique(meta_i$diagnosis))[1L],
      analysis_unit = "physical_sample_combined",
      domain = paste0("D", d), feature_role = "independent_validation",
      matrix_feature_index = validation, mz = feature_mz[validation],
      n_pixels_in = length(inside), n_pixels_out = length(outside),
      mean_residualized_z_difference =
        colMeans(validation_x_all[inside, , drop = FALSE]) -
        colMeans(validation_x_all[outside, , drop = FALSE]),
      prevalence_difference = as.numeric(Matrix::colMeans(
        tumor_log[rows_i[inside], validation, drop = FALSE] != 0) -
        Matrix::colMeans(tumor_log[rows_i[outside], validation, drop = FALSE] != 0)))
    de <- de + 1L
    discovery_effects[[de]] <- data.table(sample_id = sample_i,
      diagnosis = as.character(unique(meta_i$diagnosis))[1L],
      analysis_unit = "physical_sample_combined",
      domain = paste0("D", d), feature_role = "descriptive_discovery",
      matrix_feature_index = discovery, mz = feature_mz[discovery],
      n_pixels_in = length(inside), n_pixels_out = length(outside),
      mean_residualized_z_difference = colMeans(x_all[inside, , drop = FALSE]) -
        colMeans(x_all[outside, , drop = FALSE]),
      prevalence_difference = as.numeric(Matrix::colMeans(
        tumor_log[rows_i[inside], discovery, drop = FALSE] != 0) -
        Matrix::colMeans(tumor_log[rows_i[outside], discovery, drop = FALSE] != 0)))

    for (sid in sources_i) {
      lr <- local_source_rows[[sid]]
      inside_part <- lr[stable_core[lr] & final_lab[lr] == d]
      outside_part <- lr[stable_core[lr] & final_lab[lr] != d]
      if (length(inside_part) < 20L || length(outside_part) < 20L) next
      vpe <- vpe + 1L
      validation_partition_effects[[vpe]] <- data.table(sample_id = sample_i,
        diagnosis = as.character(unique(meta_i$diagnosis))[1L],
        export_partition_id = sid, domain = paste0("D", d),
        feature_role = "independent_validation",
        matrix_feature_index = validation, mz = feature_mz[validation],
        n_pixels_in = length(inside_part), n_pixels_out = length(outside_part),
        mean_residualized_z_difference =
          colMeans(validation_x_all[inside_part, , drop = FALSE]) -
          colMeans(validation_x_all[outside_part, , drop = FALSE]),
        prevalence_difference = as.numeric(Matrix::colMeans(
          tumor_log[rows_i[inside_part], validation, drop = FALSE] != 0) -
          Matrix::colMeans(tumor_log[rows_i[outside_part], validation, drop = FALSE] != 0)))
      dpe <- dpe + 1L
      discovery_partition_effects[[dpe]] <- data.table(sample_id = sample_i,
        diagnosis = as.character(unique(meta_i$diagnosis))[1L],
        export_partition_id = sid, domain = paste0("D", d),
        feature_role = "descriptive_discovery",
        matrix_feature_index = discovery, mz = feature_mz[discovery],
        n_pixels_in = length(inside_part), n_pixels_out = length(outside_part),
        mean_residualized_z_difference =
          colMeans(x_all[inside_part, , drop = FALSE]) -
          colMeans(x_all[outside_part, , drop = FALSE]),
        prevalence_difference = as.numeric(Matrix::colMeans(
          tumor_log[rows_i[inside_part], discovery, drop = FALSE] != 0) -
          Matrix::colMeans(tumor_log[rows_i[outside_part], discovery, drop = FALSE] != 0)))
    }
  }
  validation_effects_i <- rbindlist(validation_effects, use.names = TRUE, fill = TRUE)
  discovery_effects_i <- rbindlist(discovery_effects, use.names = TRUE, fill = TRUE)
  validation_partition_effects_i <- rbindlist(validation_partition_effects,
    use.names = TRUE, fill = TRUE)
  discovery_partition_effects_i <- rbindlist(discovery_partition_effects,
    use.names = TRUE, fill = TRUE)

  ## 层次化细分只作为筛查，不自动替换主domain。
  if (run_hierarchical_refinement_screen) {
    hz <- list(); hi <- 0L
    for (d in seq_len(chosen_k)) {
      idx <- which(stable_core & final_lab == d)
      if (length(idx) < hierarchical_min_core_pixels) next
      fit_h <- safe_kmeans(score_ref[idx, , drop = FALSE], 2L,
        reference_kmeans_nstart, random_seed + 3000000L * si + d)
      ari_init <- vapply(seq_len(10L), function(b) {
        z <- safe_kmeans(score_ref[idx, , drop = FALSE], 2L, 1L,
          random_seed + 3100000L * si + 100L * d + b)
        adjusted_rand_index(fit_h$cluster, z$cluster)
      }, numeric(1))
      score_hs <- score_ref
      for (sid in sources_i) {
        lr <- local_source_rows[[sid]]
        score_hs[lr, ] <- smooth_scores(score_ref[lr, , drop = FALSE], edges_i[[sid]], 0.20)
      }
      fit_hs <- safe_kmeans(score_hs[idx, , drop = FALSE], 2L,
        reference_kmeans_nstart, random_seed + 3200000L * si + d)
      spatial_ari <- adjusted_rand_index(fit_h$cluster, fit_hs$cluster)
      frac <- min(prop.table(table(fit_h$cluster)))
      hi <- hi + 1L
      hz[[hi]] <- data.table(sample_id = sample_i, parent_domain = paste0("D", d),
        n_stable_core_pixels = length(idx), minimum_subdomain_fraction = frac,
        median_initialization_ARI = median(ari_init), minimum_initialization_ARI = min(ari_init),
        alpha0_vs_alpha0.2_ARI = spatial_ari,
        recommended_for_secondary_refinement = frac >= hierarchical_min_subdomain_fraction &
          median(ari_init) >= hierarchical_min_initialization_ari &
          spatial_ari >= hierarchical_min_spatial_ari)
    }
    hierarchy_i <- rbindlist(hz)
  } else {
    hierarchy_i <- data.table()
  }

  ## 每个来源分别画图；稳定核心外显示灰色。
  for (sid in sources_i) {
    dt <- pixels[source_id == sid]
    dt[, stable_domain := as.character(consensus_domain)]
    dt[stable_core_pixel == FALSE, stable_domain := NA_character_]
    dt[, stable_domain := factor(stable_domain, levels = levels(consensus_domain))]
    title <- paste0(sample_i, " / ", sid)
    p_domain <- plot_tile(dt, "consensus_domain", title,
      paste0("Within-sample multiscale consensus; selected k = ", chosen_k), TRUE)
    p_core <- plot_tile(dt, "stable_domain", title,
      "Stable-core domains; grey = transition/uncertain", TRUE)
    p_pc1 <- plot_tile(dt, "reference_PC1", title, "Continuous metabolic gradient: PC1", FALSE)
    p_pc2 <- plot_tile(dt, "reference_PC2", title, "Continuous metabolic gradient: PC2", FALSE)
    p_prob <- plot_tile(dt, "consensus_probability", title,
      "Multiscale bootstrap consensus probability", FALSE, c(0, 1))
    p_entropy <- plot_tile(dt, "local_domain_entropy_r1", title,
      "Local domain entropy; radius = 1 pixel", FALSE, c(0, 1))
    plots_domain_i[[length(plots_domain_i) + 1L]] <- p_domain
    plots_core_i[[length(plots_core_i) + 1L]] <- p_core
    for (z in list(
      list(dir = "consensus_domain_maps", suffix = "consensus_domain", plot = p_domain),
      list(dir = "stable_core_maps", suffix = "stable_core", plot = p_core),
      list(dir = "continuous_PC_maps", suffix = "PC1", plot = p_pc1),
      list(dir = "continuous_PC_maps", suffix = "PC2", plot = p_pc2),
      list(dir = "consensus_diagnostic_maps", suffix = "consensus_probability", plot = p_prob),
      list(dir = "consensus_diagnostic_maps", suffix = "local_entropy", plot = p_entropy))) {
      dd <- file.path(out_dir, z$dir); dir.create(dd, recursive = TRUE, showWarnings = FALSE)
      file_stem <- safe_filename(paste(sample_i, sid, sep = "__"))
      ggsave(file.path(dd, paste0(file_stem, "_", z$suffix, ".png")),
        z$plot, width = 7.2, height = 7.2, dpi = save_png_dpi)
    }
  }
  saved_model_i <- list(pca = pca_ref, discovery_features = discovery,
    validation_features = validation, selected_k = chosen_k,
    analysis_unit = "physical_sample_combined",
    export_partition_ids = sources_i,
    candidate_qc = candidate_dt, reference_model = reference[[as.character(chosen_k)]],
    spatial_alpha_values = spatial_alpha_values,
    stable_core_probability = stable_core_probability,
    stable_core_min_reference_margin = stable_core_min_reference_margin)
  rm(x_all, validation_x_all); invisible(gc())
  result_i <- list(analysis_signature = analysis_signature,
    sample_id = sample_i, candidate_qc = candidate_dt,
    bootstrap_qc = boot_dt, model_summary = model_summary_i,
    pixels = pixels, sample_metrics = sample_metrics_i,
    source_metrics = rbindlist(source_metrics_i),
    sample_technical_qc = sample_technical_qc_i,
    technical_qc = rbindlist(technical_qc_i), feature_qc = feature_qc_i,
    validation_effects = validation_effects_i,
    discovery_effects = discovery_effects_i,
    validation_partition_effects = validation_partition_effects_i,
    discovery_partition_effects = discovery_partition_effects_i,
    hierarchy = hierarchy_i,
    saved_model = saved_model_i,
    domain_plots = plots_domain_i, core_plots = plots_core_i)
  saveRDS(result_i, file.path(checkpoint_dir,
    paste0(safe_filename(sample_i), "_completed.rds")), compress = "gzip")
  message("[", sample_i, "] 完成并写入样本检查点。")
  result_i
}

## ==== 5) 汇总输出 ====
sample_indices <- seq_along(sample_levels)
safe_analyze_one_sample <- function(si) {
  tryCatch(
    {
      sample_i <- sample_levels[si]
      checkpoint_file <- file.path(checkpoint_dir,
        paste0(safe_filename(sample_i), "_completed.rds"))
      cached <- NULL
      if (resume_completed_samples && file.exists(checkpoint_file)) {
        cached <- readRDS(checkpoint_file)
        if (!identical(cached$analysis_signature, analysis_signature)) cached <- NULL
      }
      if (!is.null(cached)) {
        message("[续跑] ", sample_i, " 已有匹配检查点，跳过重新计算。")
        list(ok = TRUE, sample_id = sample_i, result = cached)
      } else {
        list(ok = TRUE, sample_id = sample_i, result = analyze_one_sample(si))
      }
    },
    error = function(e) list(ok = FALSE, sample_id = sample_levels[si],
      error_message = conditionMessage(e), error_call = paste(deparse(conditionCall(e)), collapse = " "))
  )
}
if (parallel_samples && .Platform$OS.type != "windows") {
  workers_used <- min(as.integer(n_sample_workers), length(sample_indices))
  message("[并行] 同时分析最多 ", workers_used, " 个物理样本。")
  sample_jobs <- parallel::mclapply(sample_indices, safe_analyze_one_sample,
    mc.cores = workers_used, mc.preschedule = FALSE, mc.set.seed = FALSE)
} else {
  workers_used <- 1L
  if (parallel_samples && .Platform$OS.type == "windows") {
    warning("Windows不支持mclapply多进程，本次改为顺序运行。")
  }
  sample_jobs <- lapply(sample_indices, safe_analyze_one_sample)
}
failed <- which(!vapply(sample_jobs, `[[`, logical(1), "ok"))
if (length(failed)) {
  details <- vapply(sample_jobs[failed], function(z) paste0(
    z$sample_id, ": ", z$error_message,
    if (nzchar(z$error_call)) paste0(" | call: ", z$error_call) else ""), character(1))
  stop("样本任务失败：\n", paste(details, collapse = "\n"))
}
sample_results <- lapply(sample_jobs, `[[`, "result")

candidate_qc <- rbindlist(lapply(sample_results, `[[`, "candidate_qc"))
bootstrap_qc <- rbindlist(lapply(sample_results, `[[`, "bootstrap_qc"))
model_summary <- rbindlist(lapply(sample_results, `[[`, "model_summary"))
pixels_all <- rbindlist(lapply(sample_results, `[[`, "pixels"), use.names = TRUE, fill = TRUE)
sample_spatial_metrics <- rbindlist(lapply(sample_results, `[[`, "sample_metrics"),
  use.names = TRUE, fill = TRUE)
source_metrics <- rbindlist(lapply(sample_results, `[[`, "source_metrics"))
sample_technical_qc <- rbindlist(lapply(sample_results, `[[`, "sample_technical_qc"),
  use.names = TRUE, fill = TRUE)
partition_technical_qc <- rbindlist(lapply(sample_results, `[[`, "technical_qc"),
  use.names = TRUE, fill = TRUE)
feature_qc <- rbindlist(lapply(sample_results, `[[`, "feature_qc"))
validation_effects <- rbindlist(lapply(sample_results, `[[`, "validation_effects"),
  use.names = TRUE, fill = TRUE)
discovery_effects <- rbindlist(lapply(sample_results, `[[`, "discovery_effects"),
  use.names = TRUE, fill = TRUE)
validation_partition_effects <- rbindlist(lapply(
  sample_results, `[[`, "validation_partition_effects"), use.names = TRUE, fill = TRUE)
discovery_partition_effects <- rbindlist(lapply(
  sample_results, `[[`, "discovery_partition_effects"), use.names = TRUE, fill = TRUE)
hierarchy <- rbindlist(lapply(sample_results, `[[`, "hierarchy"),
  use.names = TRUE, fill = TRUE)
saved_models <- setNames(lapply(sample_results, `[[`, "saved_model"), sample_levels)
plot_domains <- unlist(lapply(sample_results, `[[`, "domain_plots"), recursive = FALSE)
plot_cores <- unlist(lapply(sample_results, `[[`, "core_plots"), recursive = FALSE)
if (!nrow(validation_effects)) stop("没有domain满足独立验证feature所需的核心像素数。")

## ==== 6) 样本级空间组织类型 ====
## 分类单位是物理样本。指标已在合并的样本标签和全部块内邻接边上直接计算。
sample_spatial_metrics[, `:=`(
  passes_moran_criterion = mean_domain_morans_i >=
    spatial_coherence_moran_threshold,
  passes_mixing_criterion = mixing_ratio_to_random <=
    spatial_coherence_mixing_threshold,
  passes_entropy_criterion = mean_local_entropy_r1 <=
    spatial_coherence_entropy_threshold
)]
sample_spatial_metrics[, n_spatial_coherence_criteria_met :=
  as.integer(passes_moran_criterion) + as.integer(passes_mixing_criterion) +
  as.integer(passes_entropy_criterion)]
sample_spatial_metrics[, spatial_organization_type := fifelse(
  n_spatial_coherence_criteria_met >= spatial_coherence_required_criteria,
  "spatially_coherent_domains", "spatially_intermixed_states")]
sample_spatial_metrics[, spatial_organization_label := fifelse(
  spatial_organization_type == "spatially_coherent_domains",
  "Spatially coherent domains", "Spatially intermixed states")]
sample_spatial_metrics[, sample_display_order :=
  match(as.character(sample_id), sample_levels)]
setorder(sample_spatial_metrics, sample_display_order)
sample_spatial_metrics[, sample_display_order := NULL]

## 把物理样本类型回填到每个来源和每个像素。
## NF-GnH-3 的三个来源因此始终拥有相同的 spatial_organization_type。
add_sample_classification <- function(dt) {
  dt[sample_spatial_metrics, on = .(sample_id), `:=`(
    spatial_organization_type = i.spatial_organization_type,
    spatial_organization_label = i.spatial_organization_label,
    n_spatial_coherence_criteria_met = i.n_spatial_coherence_criteria_met,
    passes_moran_criterion = i.passes_moran_criterion,
    passes_mixing_criterion = i.passes_mixing_criterion,
    passes_entropy_criterion = i.passes_entropy_criterion)]
  invisible(dt)
}
add_sample_classification(source_metrics)
add_sample_classification(pixels_all)
add_sample_classification(model_summary)

classification_rules <- data.table(
  metric = c("Physical-sample mean domain Moran's I",
    "Physical-sample mixing ratio to random",
    "Physical-sample mean local entropy (radius 1)"),
  variable = c("mean_domain_morans_i", "mixing_ratio_to_random",
    "mean_local_entropy_r1"),
  coherent_direction = c("greater_than_or_equal", "less_than_or_equal",
    "less_than_or_equal"),
  threshold = c(spatial_coherence_moran_threshold,
    spatial_coherence_mixing_threshold, spatial_coherence_entropy_threshold),
  classification_unit = "physical_sample",
  required_number_of_passed_criteria = spatial_coherence_required_criteria,
  multi_source_aggregation = paste0(
    "direct calculation on combined labels and all within-partition edges; ",
    "no artificial edges between export partitions")
)

fwrite(candidate_qc, file.path(out_dir, "within_sample_k_multiscale_candidate_qc.csv"))
fwrite(bootstrap_qc, file.path(out_dir, "within_sample_spatial_block_bootstrap_qc.csv"))
fwrite(model_summary, file.path(out_dir, "within_sample_consensus_model_summary.csv"))
fwrite(source_metrics, file.path(out_dir, "within_sample_consensus_source_metrics.csv"))
fwrite(source_metrics,
  file.path(out_dir, "export_partition_spatial_metrics_supplement.csv"))
fwrite(sample_spatial_metrics,
  file.path(out_dir, "within_sample_consensus_sample_spatial_metrics.csv"))
fwrite(classification_rules,
  file.path(out_dir, "spatial_organization_classification_rules.csv"))
fwrite(sample_technical_qc,
  file.path(out_dir, "within_sample_consensus_technical_confounding_qc.csv"))
fwrite(partition_technical_qc,
  file.path(out_dir, "export_partition_technical_confounding_qc_supplement.csv"))
fwrite(feature_qc, file.path(out_dir, "within_sample_discovery_validation_feature_qc.csv"))
fwrite(pixels_all, file.path(out_dir, "tumor_pixels_within_sample_consensus_domains.csv"))
saveRDS(pixels_all, file.path(out_dir, "tumor_pixels_within_sample_consensus_domains.rds"),
  compress = "gzip")
fwrite(validation_effects,
  file.path(out_dir, "independent_validation_feature_effects_by_domain.csv"))
fwrite(discovery_effects,
  file.path(out_dir, "descriptive_discovery_feature_effects_by_domain.csv"))
fwrite(validation_partition_effects,
  file.path(out_dir, "export_partition_validation_feature_effects_supplement.csv"))
fwrite(discovery_partition_effects,
  file.path(out_dir, "export_partition_discovery_feature_effects_supplement.csv"))
if (ncol(hierarchy)) {
  fwrite(hierarchy, file.path(out_dir, "hierarchical_secondary_refinement_screen.csv"))
} else {
  fwrite(data.table(note = "No parent domain met the minimum core-pixel requirement."),
    file.path(out_dir, "hierarchical_secondary_refinement_screen.csv"))
}
saveRDS(saved_models, file.path(out_dir, "within_sample_multiscale_consensus_models.rds"),
  compress = "gzip")

core_summary <- pixels_all[, .(n_tumor_pixels = .N,
  n_stable_core = sum(stable_core_pixel), stable_core_fraction = mean(stable_core_pixel),
  transition_fraction = mean(confidence_class == "transition"),
  uncertain_fraction = mean(confidence_class == "uncertain"),
  median_consensus_probability = median(consensus_probability),
  median_consensus_entropy = median(consensus_entropy)),
  by = .(sample_id, diagnosis, source_id, selected_k)]
fwrite(core_summary, file.path(out_dir, "within_sample_consensus_core_summary.csv"))

validation_sample_effects <- copy(validation_effects)
setnames(validation_sample_effects,
  c("mean_residualized_z_difference", "prevalence_difference"),
  c("combined_residualized_z_difference", "combined_prevalence_difference"))
partition_consistency <- validation_partition_effects[, .(
  n_export_partitions_evaluable = uniqueN(export_partition_id),
  minimum_partition_residualized_z_difference =
    min(mean_residualized_z_difference),
  positive_export_partition_fraction =
    mean(mean_residualized_z_difference > 0)),
  by = .(sample_id, diagnosis, domain, matrix_feature_index, mz, feature_role)]
validation_sample_effects <- merge(validation_sample_effects,
  partition_consistency,
  by = c("sample_id", "diagnosis", "domain", "matrix_feature_index", "mz",
    "feature_role"), all.x = TRUE, sort = FALSE)
fwrite(validation_sample_effects,
  file.path(out_dir, "independent_validation_feature_effects_sample_summary.csv"))
top_validation <- validation_sample_effects[
  combined_residualized_z_difference > 0 &
    (n_export_partitions_evaluable <= 1 |
      positive_export_partition_fraction >= 2/3)][
  order(sample_id, domain, -combined_residualized_z_difference), head(.SD, 30L),
  by = .(sample_id, domain)]
fwrite(top_validation, file.path(out_dir, "top30_independent_validation_features_per_domain.csv"))

grDevices::pdf(file.path(out_dir, "all_12sources_consensus_domain_maps.pdf"),
  7.2, 7.2, onefile = TRUE)
for (p in plot_domains) print(p)
grDevices::dev.off()
grDevices::pdf(file.path(out_dir, "all_12sources_stable_core_maps.pdf"),
  7.2, 7.2, onefile = TRUE)
for (p in plot_cores) print(p)
grDevices::dev.off()

## ==== 7) 空间组织类型可视化 ====
organization_dir <- file.path(out_dir, "spatial_organization_classification")
dir.create(organization_dir, recursive = TRUE, showWarnings = FALSE)
organization_palette <- c(
  spatially_coherent_domains = "#0072B2",
  spatially_intermixed_states = "#D55E00"
)

## 7a. 三个主指标的样本级点图；虚线是分类阈值。
sample_metric_long <- melt(sample_spatial_metrics,
  id.vars = c("sample_id", "diagnosis", "n_export_partitions", "n_tumor_pixels",
    "spatial_organization_type", "spatial_organization_label"),
  measure.vars = c("mean_domain_morans_i",
    "mixing_ratio_to_random", "mean_local_entropy_r1"),
  variable.name = "metric", value.name = "value")
metric_lookup <- data.table(
  metric = c("mean_domain_morans_i", "mixing_ratio_to_random",
    "mean_local_entropy_r1"),
  metric_label = c("Mean domain Moran's I (higher = more coherent)",
    "Mixing ratio to random (lower = more coherent)",
    "Mean local entropy (lower = more coherent)"),
  threshold = c(spatial_coherence_moran_threshold,
    spatial_coherence_mixing_threshold, spatial_coherence_entropy_threshold),
  direction = c("higher", "lower", "lower")
)
sample_metric_long <- merge(sample_metric_long, metric_lookup, by = "metric",
  all.x = TRUE, sort = FALSE)
sample_metric_long[, sample_id := factor(as.character(sample_id), levels = sample_levels)]
sample_metric_long[, metric_label := factor(metric_label,
  levels = metric_lookup$metric_label)]
p_sample_metrics <- ggplot(sample_metric_long,
    aes(x = sample_id, y = value, color = spatial_organization_type)) +
  geom_hline(data = metric_lookup, aes(yintercept = threshold),
    inherit.aes = FALSE, linetype = 2, color = "grey35", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", value)), vjust = -0.75,
    size = 2.8, show.legend = FALSE) +
  facet_wrap(~metric_label, ncol = 1, scales = "free_y") +
  scale_color_manual(values = organization_palette,
    breaks = names(organization_palette),
    labels = c("Spatially coherent domains", "Spatially intermixed states")) +
  labs(title = "Sample-level spatial organization metrics",
    subtitle = "Physical-sample classification; dashed lines are prespecified thresholds",
    x = NULL, y = "Metric value", color = "Spatial organization") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(
    angle = 45, hjust = 1), strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"))
ggsave(file.path(organization_dir, "sample_spatial_organization_metric_panels.png"),
  p_sample_metrics, width = 10.5, height = 10.5, dpi = save_png_dpi)
ggsave(file.path(organization_dir, "sample_spatial_organization_metric_panels.pdf"),
  p_sample_metrics, width = 10.5, height = 10.5)

## 7b. 阈值距离热图：正值（绿）代表达到连贯标准，负值（红）代表未达到。
sample_metric_long[, coherence_margin := fifelse(direction == "higher",
  value - threshold, threshold - value)]
sample_metric_long[, sample_id_heatmap := factor(as.character(sample_id),
  levels = rev(sample_levels))]
p_metric_heatmap <- ggplot(sample_metric_long,
    aes(x = metric_label, y = sample_id_heatmap, fill = coherence_margin)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.2) +
  scale_fill_gradient2(low = "#D73027", mid = "white", high = "#1A9850",
    midpoint = 0, name = "Margin to\nthreshold") +
  labs(title = "Spatial coherence criteria by physical sample",
    subtitle = "Cell labels are raw metric values; green meets the criterion",
    x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(
    angle = 25, hjust = 1), plot.title = element_text(face = "bold"))
ggsave(file.path(organization_dir, "sample_spatial_organization_heatmap.png"),
  p_metric_heatmap, width = 10.5, height = 6.5, dpi = save_png_dpi)
ggsave(file.path(organization_dir, "sample_spatial_organization_heatmap.pdf"),
  p_metric_heatmap, width = 10.5, height = 6.5)

## 7c. 导出分块级指标作为透明性补充；颜色仍是物理样本的共同类型。
source_metric_long <- melt(source_metrics,
  id.vars = c("sample_id", "diagnosis", "source_id", "n_tumor_pixels",
    "spatial_organization_type", "spatial_organization_label"),
  measure.vars = c("mean_domain_morans_i", "mixing_ratio_to_random",
    "mean_local_entropy_r1"), variable.name = "source_metric", value.name = "value")
source_metric_lookup <- data.table(
  source_metric = c("mean_domain_morans_i", "mixing_ratio_to_random",
    "mean_local_entropy_r1"),
  metric_label = metric_lookup$metric_label,
  threshold = metric_lookup$threshold
)
source_metric_long <- merge(source_metric_long, source_metric_lookup,
  by = "source_metric", all.x = TRUE, sort = FALSE)
source_order <- as.character(remaining_manifest$source_id)
source_metric_long[, source_label := factor(
  paste0(as.character(sample_id), " / ", as.character(source_id)),
  levels = paste0(as.character(remaining_manifest$sample_id), " / ", source_order))]
source_metric_long[, metric_label := factor(metric_label,
  levels = metric_lookup$metric_label)]
p_source_metrics <- ggplot(source_metric_long,
    aes(x = source_label, y = value, color = spatial_organization_type)) +
  geom_hline(data = source_metric_lookup, aes(yintercept = threshold),
    inherit.aes = FALSE, linetype = 2, color = "grey35", linewidth = 0.5) +
  geom_point(size = 2.8) +
  facet_wrap(~metric_label, ncol = 1, scales = "free_y") +
  scale_color_manual(values = organization_palette,
    breaks = names(organization_palette),
    labels = c("Spatially coherent domains", "Spatially intermixed states")) +
  labs(title = "Export-partition spatial metrics with physical-sample classification",
    subtitle = "Export partitions are one joint sample; partitions only preserve coordinates",
    x = NULL, y = "Metric value", color = "Sample-level class") +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(
    angle = 55, hjust = 1), strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"))
ggsave(file.path(organization_dir, "export_partition_spatial_metrics_sample_class.png"),
  p_source_metrics, width = 12.5, height = 10.5, dpi = save_png_dpi)
ggsave(file.path(organization_dir, "export_partition_spatial_metrics_sample_class.pdf"),
  p_source_metrics, width = 12.5, height = 10.5)

## 7d. 按物理样本类型分成两组展示 domain map。
## 同一样本的多个来源会进入同一个 PDF，但仍使用各自坐标绘图。
organization_domain_plots <- setNames(vector("list", 2L), names(organization_palette))
organization_core_plots <- setNames(vector("list", 2L), names(organization_palette))
for (sample_i in sample_levels) {
  class_i <- sample_spatial_metrics[sample_id == sample_i,
    as.character(spatial_organization_type)][1L]
  class_label_i <- sample_spatial_metrics[sample_id == sample_i,
    as.character(spatial_organization_label)][1L]
  sources_i <- as.character(remaining_manifest[sample_id == sample_i, source_id])
  for (sid in sources_i) {
    dt <- copy(pixels_all[sample_id == sample_i & source_id == sid])
    dt[, stable_domain := as.character(consensus_domain)]
    dt[stable_core_pixel == FALSE, stable_domain := NA_character_]
    dt[, stable_domain := factor(stable_domain, levels = levels(consensus_domain))]
    sm <- source_metrics[sample_id == sample_i & source_id == sid][1L]
    metric_text <- paste0("Sample class: ", class_label_i,
      "; source Moran's I = ", sprintf("%.3f", sm$mean_domain_morans_i),
      ", mixing = ", sprintf("%.3f", sm$mixing_ratio_to_random),
      ", entropy = ", sprintf("%.3f", sm$mean_local_entropy_r1))
    title <- paste0(sample_i, " / ", sid)
    pd <- plot_tile(dt, "consensus_domain", title, metric_text, TRUE)
    pc <- plot_tile(dt, "stable_domain", title,
      paste0(class_label_i, "; grey = transition/uncertain"), TRUE)
    organization_domain_plots[[class_i]][[length(
      organization_domain_plots[[class_i]]) + 1L]] <- pd
    organization_core_plots[[class_i]][[length(
      organization_core_plots[[class_i]]) + 1L]] <- pc
    class_dir <- file.path(organization_dir, class_i)
    dir.create(class_dir, recursive = TRUE, showWarnings = FALSE)
    file_stem <- safe_filename(paste(sample_i, sid, sep = "__"))
    ggsave(file.path(class_dir, paste0(file_stem, "_domain_map.png")),
      pd, width = 7.8, height = 7.2, dpi = save_png_dpi)
  }
}

save_multipage_pdf <- function(plots, filename) {
  if (!length(plots)) return(invisible(NULL))
  grDevices::pdf(filename, 7.8, 7.2, onefile = TRUE)
  for (p in plots) print(p)
  grDevices::dev.off()
  invisible(filename)
}
save_multipage_pdf(organization_domain_plots$spatially_coherent_domains,
  file.path(organization_dir, "spatially_coherent_domain_maps.pdf"))
save_multipage_pdf(organization_domain_plots$spatially_intermixed_states,
  file.path(organization_dir, "spatially_intermixed_state_maps.pdf"))
save_multipage_pdf(organization_core_plots$spatially_coherent_domains,
  file.path(organization_dir, "spatially_coherent_stable_core_maps.pdf"))
save_multipage_pdf(organization_core_plots$spatially_intermixed_states,
  file.path(organization_dir, "spatially_intermixed_stable_core_maps.pdf"))

saveRDS(list(
  quick_test = quick_test, parallel_samples = parallel_samples,
  n_sample_workers_requested = n_sample_workers, n_sample_workers_used = workers_used,
  resume_completed_samples = resume_completed_samples,
  analysis_signature = analysis_signature,
  n_spatial_bootstrap = n_spatial_bootstrap,
  k_candidates = k_candidates, spatial_alpha_values = spatial_alpha_values,
  spatial_block_width = spatial_block_width,
  bootstrap_block_fraction = bootstrap_block_fraction,
  k_selection_mode = k_selection_mode, manual_k_override = manual_k_override,
  stable_core_probability = stable_core_probability,
  transition_probability = transition_probability,
  stable_core_min_reference_margin = stable_core_min_reference_margin,
  spatial_coherence_moran_threshold = spatial_coherence_moran_threshold,
  spatial_coherence_mixing_threshold = spatial_coherence_mixing_threshold,
  spatial_coherence_entropy_threshold = spatial_coherence_entropy_threshold,
  spatial_coherence_required_criteria = spatial_coherence_required_criteria,
  export_partition_policy = paste0(
    "combine within physical sample for feature QC, normalization, PCA, clustering, ",
    "bootstrap, feature effects and sample metrics; separate only for adjacency and plotting"),
  input_cohort = "10 physical samples; Cushing-3 absent from upstream imzML manifest",
  tumor_cluster_map = tumor_cluster_map,
  random_seed = random_seed, session_info = utils::sessionInfo()),
  file.path(out_dir, "run_parameters_and_sessionInfo.rds"), compress = "gzip")

message("完成。输出目录：", out_dir)
message("请先检查 within_sample_consensus_model_summary.csv、",
  "within_sample_k_multiscale_candidate_qc.csv 和 within_sample_consensus_core_summary.csv。")
message("独立验证feature结果：independent_validation_feature_effects_by_domain.csv")
message("空间组织类型和可视化：spatial_organization_classification/")
