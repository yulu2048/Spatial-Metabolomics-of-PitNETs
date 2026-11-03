## ============================================================================
## 0_readin_align_10samples.R
##
## 目的：
##   1) 从 12 个 imzML 区块读入 10 个物理样本；
##   2) 始终保留 sample_id（物理样本）与 source_id（原始 imzML）两级信息；
##   3) 预处理后先做质量轴校正/峰对齐，再由 12 个来源的平均谱
##      建立共识峰；共识峰组最大 m/z 跨度不超过 0.1 Da；
##   4) 保存 TIC 矩阵（用于离子图）和 log-CPM 矩阵（用于 PCA/聚类等比较）；
##   5) Cushing-3 及其两个 imzML 来源不在输入清单中，从读入阶段即排除。
##
## 已确认：三个目录只是存储位置，所有数据属于同一次实验批次。因此全部标记为
##         batch_1，不执行强度批次矫正；目录信息仍由 relative_file/source_file 保留。
## ============================================================================

required_packages <- c(
  "MALDIquant", "MALDIquantForeign", "data.table", "ggplot2", "Matrix"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "。请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(MALDIquant)
  library(MALDIquantForeign)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

set.seed(2025)

## ==== 1) 路径、参数和 10 个物理样本映射 ====

data_root <- "/home/gaga/huzhan"
out_dir   <- file.path(data_root, "MSI_ALL_aligned_10samples_consensus_tuned225")
qc_dir    <- file.path(out_dir, "QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,  recursive = TRUE, showWarnings = FALSE)

## 如果内存不足，可设为整数（例如 10000），表示每个 imzML 最多保留多少像素。
## NULL 表示保留全部像素；完整空间图建议使用 NULL。
max_pixels_per_source <- NULL

## 预处理参数
smooth_half_window <- 6L
baseline_method    <- "SNIP"
peak_snr           <- 4
peak_half_window   <- 20L

## 共识 feature 识别参数：只在对齐后的来源平均谱上检峰。
## 允许只在一个物理样本平均谱上出现，以保留潜在样本特异峰；
## 后面再用样本内像素出现率控制假阳性。
consensus_peak_snr       <- 4
consensus_min_samples    <- 1L

## 对齐参数：MALDIquant 的 tolerance 是相对误差，ppm 必须乘 1e-6。
alignment_search_ppm    <- 250
peak_group_max_span_da  <- 0.1
pixel_assignment_half_window_da <- peak_group_max_span_da / 2
reference_min_frequency <- 0.50
warping_method          <- "lowess"

## 先做文件级校正，再做像素级精细校正。像素级步骤较慢，但能处理文件内漂移。
perform_pixel_alignment <- TRUE
alignment_chunk_size    <- 500L

## 特征过滤参数：按物理样本内的出现率过滤，不再按全部像素的1%。
min_within_sample_nonzero_prop <- 0.05
min_nonzero_pixels_in_sample   <- 20L
min_samples_passing_prevalence <- 1L
high_nonzero_prop <- 0.98
max_cv_quasi_const <- 0.10

## FALSE：仅完成 feature 识别、稀疏矩阵建立和基础 QC，然后停止。
## 确认矩阵满意后改为 TRUE，才运行样本级 PCA 并保存 PCA 结果。
run_downstream_pca <- FALSE

## TRUE：跳过 imzML 读入、预处理、对齐、feature 识别和矩阵构建，
## 直接从 out_dir 读取已保存矩阵，续跑第10节样本QC。
## 本发布版本默认从10样本的 imzML 重新开始，因此保持 FALSE。
resume_after_feature_matrix <- FALSE

## 所有文件为同一批次，不能也不需要做 ComBat/removeBatchEffect 类批次矫正。
## 保留下面的开关和函数仅供未来加入真实多批次数据时使用。
perform_batch_correction <- FALSE
batch_feature_chunk_size <- 250L

sample_levels <- c(
  "GH-1", "GH-2", "GH-3",
  "NF-GH", "NF-ACTH",
  "Cushing-1", "Cushing-2",
  "NF-GnH-1", "NF-GnH-2", "NF-GnH-3"
)

sample_manifest <- data.table::data.table(
  sample_id = c(
    "GH-1", "GH-2", "GH-3",
    "NF-GH", "NF-ACTH",
    "Cushing-1", "Cushing-2",
    "NF-GnH-1", "NF-GnH-2",
    "NF-GnH-3", "NF-GnH-3", "NF-GnH-3"
  ),
  diagnosis = c(
    "GH", "GH", "GH",
    "NF-GH", "NF-ACTH",
    "Cushing", "Cushing",
    "NF-GnH", "NF-GnH",
    "NF-GnH", "NF-GnH", "NF-GnH"
  ),
  source_id = c(
    "GH", "gh_1", "gh_2",
    "NF_GH", "NF_ACTH",
    "Cushing", "cushing_0",
    "NF_CXX", "nf_cxx",
    "nf_cxx_1", "nf_cxx_2", "nf_cxx_3"
  ),
  batch_id = rep("batch_1", 12L),
  relative_file = c(
    "MSI/data/GH.imzML",
    "Other_MSI/gh_1 - total ion count.imzML",
    "Other_MSI/gh_2 - total ion count.imzML",
    "MSI/data/NF_GH.imzML",
    "MSI/data/NF_ACTH.imzML",
    "MSI/data/Cushing.imzML",
    "Other_MSI/cushing_0 - total ion count.imzML",
    "MSI/data/NF_CXX.imzML",
    "NF_prediction/nf_cxx - total ion count.imzML",
    "NF_prediction/nf_cxx_1 - total ion count.imzML",
    "NF_prediction/nf_cxx_2 - total ion count.imzML",
    "NF_prediction/nf_cxx_3 - total ion count.imzML"
  )
)
sample_manifest[, sample_id := factor(sample_id, levels = sample_levels)]
sample_manifest[, imzml_file := file.path(data_root, relative_file)]
sample_manifest[, ibd_file := sub("\\.imzML$", ".ibd", imzml_file, ignore.case = TRUE)]

stopifnot(length(unique(sample_manifest$sample_id)) == 10L)
stopifnot(nrow(sample_manifest) == 12L)
if (anyDuplicated(sample_manifest$source_id)) stop("source_id 必须唯一。")

missing_imzml <- sample_manifest[!file.exists(imzml_file), imzml_file]
missing_ibd   <- sample_manifest[!file.exists(ibd_file), ibd_file]
if (length(missing_imzml) || length(missing_ibd)) {
  stop(
    "缺少数据文件。\n",
    if (length(missing_imzml)) paste0("imzML:\n", paste(missing_imzml, collapse = "\n"), "\n") else "",
    if (length(missing_ibd)) paste0("ibd:\n", paste(missing_ibd, collapse = "\n")) else ""
  )
}

data.table::fwrite(sample_manifest, file.path(out_dir, "sample_source_manifest.csv"))
saveRDS(sample_manifest, file.path(out_dir, "sample_source_manifest.rds"))

## ==== 2) 工具函数 ====

get_xy <- function(sp) {
  xy <- tryCatch(MALDIquant::coordinates(sp), error = function(e) NULL)
  if (!is.null(xy) && length(xy) >= 2L) {
    xy <- as.numeric(xy[1:2])
    return(c(x = xy[1], y = xy[2]))
  }

  md <- tryCatch(MALDIquant::metaData(sp), error = function(e) list())
  im <- md$imaging
  if (!is.null(im$pos) && length(im$pos) >= 2L) {
    return(c(x = as.numeric(im$pos[1]), y = as.numeric(im$pos[2])))
  }
  if (!is.null(im$x) && !is.null(im$y)) {
    return(c(x = as.numeric(im$x), y = as.numeric(im$y)))
  }
  c(x = NA_real_, y = NA_real_)
}

safe_tic <- function(sp) {
  z <- MALDIquant::intensity(sp)
  sum(z[is.finite(z) & z > 0], na.rm = TRUE)
}

preprocess_before_alignment <- function(sp) {
  sp <- MALDIquant::transformIntensity(sp, method = "sqrt")
  sp <- MALDIquant::smoothIntensity(
    sp, method = "SavitzkyGolay", halfWindowSize = smooth_half_window
  )
  MALDIquant::removeBaseline(sp, method = baseline_method)
}

detect_one <- function(sp, snr = peak_snr) {
  MALDIquant::detectPeaks(
    sp, method = "MAD", SNR = snr, halfWindowSize = peak_half_window
  )
}

peak_count <- function(p) length(MALDIquant::mass(p))

## 由 12 个来源平均谱的检测峰建立共识 feature：
##   - 不用全部像素的密集峰直接切分质量轴；
##   - 每组候选峰的 max(m/z)-min(m/z) <= max_span_da；
##   - 同一物理样本在一组内只贡献一个峰位，避免多 imzML 样本过度加权；
##   - feature m/z 为各物理样本代表峰位的中位数。
build_consensus_features <- function(source_average_peaks,
                                     source_sample_id,
                                     mz_min,
                                     mz_max,
                                     max_span_da = 0.1,
                                     min_samples = 1L) {
  stopifnot(
    is.list(source_average_peaks), length(source_average_peaks) > 0L,
    length(source_average_peaks) == length(source_sample_id),
    is.finite(mz_min), is.finite(mz_max), mz_min < mz_max,
    is.finite(max_span_da), max_span_da > 0,
    min_samples >= 1L
  )

  candidate_dt <- data.table::rbindlist(lapply(
    seq_along(source_average_peaks), function(i) {
      p <- source_average_peaks[[i]]
      data.table::data.table(
        source_index = i,
        sample_id = as.character(source_sample_id[i]),
        mz = MALDIquant::mass(p),
        intensity = MALDIquant::intensity(p)
      )
    }
  ))
  candidate_dt <- candidate_dt[
    is.finite(mz) & is.finite(intensity) & intensity > 0 &
      mz >= mz_min & mz <= mz_max
  ]
  if (!nrow(candidate_dt)) stop("来源平均谱上没有检测到候选峰。")

  data.table::setorder(candidate_dt, mz)
  mz_sorted <- candidate_dt$mz
  group_sorted <- integer(length(mz_sorted))
  group_start_idx <- integer()
  group_end_idx <- integer()
  numerical_slack <- sqrt(.Machine$double.eps)
  pos <- 1L
  group_id <- 0L
  while (pos <= length(mz_sorted)) {
    group_id <- group_id + 1L
    end <- findInterval(mz_sorted[pos] + max_span_da + numerical_slack,
                        mz_sorted)
    if (end < pos) end <- pos
    group_sorted[pos:end] <- group_id
    group_start_idx[group_id] <- pos
    group_end_idx[group_id] <- end
    pos <- end + 1L
  }
  candidate_dt[, preliminary_group := group_sorted]

  group_source_stats <- candidate_dt[, .(
    consensus_candidate_mz_min = min(mz),
    consensus_candidate_mz_max = max(mz),
    consensus_candidate_span_da = max(mz) - min(mz),
    n_sources_with_average_peak = uniqueN(source_index),
    n_source_average_peaks = .N
  ), by = preliminary_group]

  ## 每个物理样本在同一峰组内只保留平均谱强度最高的候选峰。
  sample_representatives <- candidate_dt[
    candidate_dt[, .I[which.max(intensity)], by = .(preliminary_group, sample_id)]$V1
  ]
  consensus_dt <- sample_representatives[, .(
    mz = stats::median(mz),
    n_samples_with_average_peak = uniqueN(sample_id)
  ), by = preliminary_group]
  consensus_dt <- merge(
    consensus_dt, group_source_stats,
    by = "preliminary_group", all.x = TRUE, sort = FALSE
  )
  consensus_dt <- consensus_dt[n_samples_with_average_peak >= min_samples]
  if (!nrow(consensus_dt)) {
    stop("没有共识 feature 通过 consensus_min_samples 阈值。")
  }
  data.table::setorder(consensus_dt, mz)
  consensus_dt[, feature_id := seq_len(.N)]
  stopifnot(all(
    consensus_dt$consensus_candidate_span_da <= max_span_da + numerical_slack
  ))
  consensus_dt[]
}

## 把每个像素检测峰分配给最近的共识峰中心：
##   - 只接受 abs(pixel peak m/z - feature m/z) <= half_window_da；
##   - half_window_da=0.05 时，同一 feature 中已分配像素峰的理论总跨度 <=0.1 Da；
##   - 同一像素在同一 feature 内出现多个局部峰时取最高强度；
##   - 未检测到峰时保持稀疏矩阵中的 0。
assign_pixel_peaks_to_consensus <- function(peaks_list,
                                            feature_mz,
                                            half_window_da = 0.05) {
  stopifnot(
    is.list(peaks_list), length(peaks_list) > 0L,
    length(feature_mz) > 0L, all(is.finite(feature_mz)),
    is.finite(half_window_da), half_window_da > 0
  )
  feature_mz <- sort(feature_mz)
  peaks_per_pixel <- vapply(peaks_list, peak_count, integer(1))
  peak_pixel <- rep.int(seq_along(peaks_list), peaks_per_pixel)
  peak_mz <- unlist(lapply(peaks_list, MALDIquant::mass), use.names = FALSE)
  peak_intensity <- unlist(lapply(peaks_list, MALDIquant::intensity),
                           use.names = FALSE)
  valid <- is.finite(peak_mz) & is.finite(peak_intensity) & peak_intensity > 0
  peak_pixel <- peak_pixel[valid]
  peak_mz <- peak_mz[valid]
  peak_intensity <- peak_intensity[valid]
  if (!length(peak_mz)) stop("所有像素都没有可分配的检测峰。")

  left <- findInterval(peak_mz, feature_mz)
  left_idx <- pmax(left, 1L)
  right_idx <- pmin(left + 1L, length(feature_mz))
  left_dist <- abs(peak_mz - feature_mz[left_idx])
  right_dist <- abs(peak_mz - feature_mz[right_idx])
  use_right <- right_dist < left_dist
  feature_idx <- ifelse(use_right, right_idx, left_idx)
  distance_da <- pmin(left_dist, right_dist)
  assigned <- distance_da <= half_window_da

  peak_dt <- data.table::data.table(
    pixel = peak_pixel[assigned],
    feature = feature_idx[assigned],
    intensity = peak_intensity[assigned],
    distance_da = distance_da[assigned]
  )
  if (!nrow(peak_dt)) stop("没有像素峰落在共识峰匹配窗口内。")
  assignment_qc <- peak_dt[, .(
    n_assigned_pixel_peaks = .N,
    median_abs_assignment_error_da = stats::median(distance_da),
    p95_abs_assignment_error_da = stats::quantile(
      distance_da, probs = 0.95, names = FALSE
    ),
    max_abs_assignment_error_da = max(distance_da)
  ), by = feature]
  peak_dt <- peak_dt[, .(intensity = max(intensity)), by = .(pixel, feature)]

  sparse_mat <- Matrix::sparseMatrix(
    i = peak_dt$pixel,
    j = peak_dt$feature,
    x = peak_dt$intensity,
    dims = c(length(peaks_list), length(feature_mz)),
    giveCsparse = TRUE
  )
  colnames(sparse_mat) <- sprintf("mz_%.6f", feature_mz)
  list(matrix = sparse_mat, assignment_qc = assignment_qc)
}

## 对 dgCMatrix 的非零元素原位做 log1p-CPM，避免隐式转成稠密矩阵。
sparse_log1p_cpm <- function(mat, scale_factor = 1e6) {
  if (!inherits(mat, "dgCMatrix")) {
    stop("sparse_log1p_cpm 需要 Matrix::dgCMatrix。")
  }
  ans <- mat
  ans@x <- log1p(ans@x * scale_factor)
  ans
}

nearest_reference_residual_ppm <- function(peaks, reference, search_ppm) {
  pm <- MALDIquant::mass(peaks)
  rm <- MALDIquant::mass(reference)
  if (!length(pm) || !length(rm)) return(numeric())
  ans <- vapply(rm, function(m0) {
    j <- which.min(abs(pm - m0))
    d <- (pm[j] - m0) / m0 * 1e6
    if (abs(d) <= search_ppm) d else NA_real_
  }, numeric(1))
  ans[is.finite(ans)]
}

robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - med) / s
}

if (!resume_after_feature_matrix) {

## ==== 3) 逐 imzML 读入；保留原始像素索引、坐标和两级样本身份 ====

source_objects <- vector("list", nrow(sample_manifest))
names(source_objects) <- as.character(sample_manifest$source_id)

for (i in seq_len(nrow(sample_manifest))) {
  row_i <- sample_manifest[i]
  message("[读取] ", row_i$sample_id, " / ", row_i$source_id, ": ", row_i$imzml_file)

  spectra_raw <- MALDIquantForeign::importImzMl(row_i$imzml_file, verbose = FALSE)
  original_index <- seq_along(spectra_raw)
  raw_tic_all <- vapply(spectra_raw, safe_tic, numeric(1))
  nonempty <- !vapply(spectra_raw, MALDIquant::isEmpty, logical(1))

  spectra_raw  <- spectra_raw[nonempty]
  original_index <- original_index[nonempty]
  raw_tic      <- raw_tic_all[nonempty]
  if (!length(spectra_raw)) stop("该文件全部为空谱：", row_i$imzml_file)

  xy <- t(vapply(spectra_raw, get_xy, numeric(2)))
  colnames(xy) <- c("x", "y")

  if (!is.null(max_pixels_per_source) && length(spectra_raw) > max_pixels_per_source) {
    keep <- sort(sample(seq_along(spectra_raw), max_pixels_per_source))
    spectra_raw   <- spectra_raw[keep]
    original_index <- original_index[keep]
    raw_tic       <- raw_tic[keep]
    xy            <- xy[keep, , drop = FALSE]
  }

  pixel_meta <- data.table::data.table(
    sample_id = as.character(row_i$sample_id),
    diagnosis = row_i$diagnosis,
    source_id = row_i$source_id,
    batch_id = row_i$batch_id,
    source_file = row_i$imzml_file,
    pixel_index_in_source = original_index,
    x = xy[, "x"],
    y = xy[, "y"],
    raw_tic = raw_tic
  )

  spectra_pp <- lapply(spectra_raw, preprocess_before_alignment)
  source_objects[[i]] <- list(spectra = spectra_pp, pixel_meta = pixel_meta)
  rm(spectra_raw, spectra_pp)
  invisible(gc())
}

## 所有谱取共同覆盖的 m/z 范围，避免不同文件边缘范围影响平均谱和对齐。
all_pp_for_range <- unlist(lapply(source_objects, `[[`, "spectra"), recursive = FALSE)
common_min <- max(vapply(all_pp_for_range, function(sp) min(MALDIquant::mass(sp)), numeric(1)))
common_max <- min(vapply(all_pp_for_range, function(sp) max(MALDIquant::mass(sp)), numeric(1)))
if (!is.finite(common_min) || !is.finite(common_max) || common_min >= common_max) {
  stop("所有文件之间不存在有效的共同 m/z 范围。")
}
message("共同 m/z 范围：", signif(common_min, 8), " - ", signif(common_max, 8))
for (i in seq_along(source_objects)) {
  source_objects[[i]]$spectra <- MALDIquant::trim(
    source_objects[[i]]$spectra, range = c(common_min, common_max)
  )
}
rm(all_pp_for_range)
invisible(gc())

## ==== 4) 构建参考峰：每个物理样本只贡献一个代表性文件 ====

source_average <- lapply(source_objects, function(z) {
  MALDIquant::averageMassSpectra(z$spectra, method = "mean")
})
source_peaks_pre <- lapply(source_average, detect_one, snr = 3)

## 多 imzML 样本选择峰数最多的文件参与“参考峰”构建，防止一个物理样本重复加权。
representative_idx <- unlist(lapply(sample_levels, function(sid) {
  idx <- which(as.character(sample_manifest$sample_id) == sid)
  idx[which.max(vapply(source_peaks_pre[idx], peak_count, integer(1)))]
}), use.names = FALSE)

reference_peaks <- MALDIquant::referencePeaks(
  source_peaks_pre[representative_idx],
  method = "strict",
  minFrequency = reference_min_frequency,
  tolerance = alignment_search_ppm * 1e-6
)
if (peak_count(reference_peaks) < 5L) {
  stop(
    "共同参考峰少于 5 个，无法可靠对齐。建议检查平均谱，或适当降低 ",
    "reference_min_frequency / 增大 alignment_search_ppm。"
  )
}
saveRDS(reference_peaks, file.path(out_dir, "alignment_reference_peaks.rds"))

## ==== 5) 文件级质量轴对齐，再做像素级精细对齐 ====

source_warp <- MALDIquant::determineWarpingFunctions(
  source_peaks_pre,
  reference = reference_peaks,
  tolerance = alignment_search_ppm * 1e-6,
  method = warping_method,
  allowNoMatches = TRUE
)
source_nmatch <- attr(source_warp, "nmatch")
if (is.null(source_nmatch)) source_nmatch <- rep(NA_integer_, length(source_warp))

source_peaks_post <- MALDIquant::warpMassPeaks(
  source_peaks_pre, source_warp, emptyNoMatches = FALSE
)

source_alignment_qc <- data.table::rbindlist(lapply(seq_along(source_objects), function(i) {
  pre_res  <- nearest_reference_residual_ppm(
    source_peaks_pre[[i]], reference_peaks, alignment_search_ppm
  )
  post_res <- nearest_reference_residual_ppm(
    source_peaks_post[[i]], reference_peaks, alignment_search_ppm
  )
  data.table::data.table(
    sample_id = as.character(sample_manifest$sample_id[i]),
    source_id = sample_manifest$source_id[i],
    batch_id = sample_manifest$batch_id[i],
    n_reference_matches = source_nmatch[i],
    median_abs_ppm_before = if (length(pre_res)) median(abs(pre_res)) else NA_real_,
    median_abs_ppm_after  = if (length(post_res)) median(abs(post_res)) else NA_real_
  )
}))
data.table::fwrite(source_alignment_qc,
                   file.path(qc_dir, "source_alignment_qc.csv"))

for (i in seq_along(source_objects)) {
  ## 先把该 imzML 平均谱得到的文件级漂移函数用于全部像素。
  w_source_i <- rep(list(source_warp[[i]]), length(source_objects[[i]]$spectra))
  spectra_source_aligned <- MALDIquant::warpMassSpectra(
    source_objects[[i]]$spectra, w_source_i, emptyNoMatches = FALSE
  )

  n_pixel <- length(spectra_source_aligned)
  pixel_nmatch <- rep(NA_integer_, n_pixel)
  pixel_alignment_failed <- rep(FALSE, n_pixel)

  if (perform_pixel_alignment) {
    chunk_starts <- seq.int(1L, n_pixel, by = alignment_chunk_size)
    for (start in chunk_starts) {
      end <- min(start + alignment_chunk_size - 1L, n_pixel)
      idx <- start:end
      p_chunk <- lapply(spectra_source_aligned[idx], detect_one, snr = 3)
      w_chunk <- MALDIquant::determineWarpingFunctions(
        p_chunk,
        reference = reference_peaks,
        tolerance = alignment_search_ppm * 1e-6,
        method = warping_method,
        allowNoMatches = TRUE
      )
      nm <- attr(w_chunk, "nmatch")
      if (is.null(nm)) nm <- rep(NA_integer_, length(idx))
      pixel_nmatch[idx] <- nm
      pixel_alignment_failed[idx] <- vapply(
        w_chunk,
        function(w) is.atomic(w) && length(w) == 1L && is.na(w),
        logical(1)
      )
      ## 对齐失败时 emptyNoMatches=FALSE 会保留文件级对齐后的原谱，不会置零。
      spectra_source_aligned[idx] <- MALDIquant::warpMassSpectra(
        spectra_source_aligned[idx], w_chunk, emptyNoMatches = FALSE
      )
      rm(p_chunk, w_chunk)
    }
  }

  ## 对齐完成后再做 TIC 强度归一化。
  spectra_source_aligned <- MALDIquant::calibrateIntensity(
    spectra_source_aligned, method = "TIC"
  )
  source_objects[[i]]$spectra <- spectra_source_aligned
  source_objects[[i]]$pixel_meta[, `:=`(
    alignment_nmatch = pixel_nmatch,
    alignment_failed = pixel_alignment_failed
  )]
  message("[对齐完成] ", sample_manifest$sample_id[i], " / ",
          sample_manifest$source_id[i])
  invisible(gc())
}

## ==== 6) 对齐后建立平均谱共识峰，再构建稀疏像素矩阵 ====

## feature 候选只来自 12 个对齐后来源平均谱，而不是所有像素峰。
## 这样可防止数万像素的密集 m/z 位置把整个质量轴切成连续 0.1 Da 区块。
source_average_aligned <- lapply(source_objects, function(z) {
  MALDIquant::averageMassSpectra(z$spectra, method = "mean")
})
source_average_peaks_aligned <- lapply(
  source_average_aligned, detect_one, snr = consensus_peak_snr
)
names(source_average_aligned) <- as.character(sample_manifest$source_id)
names(source_average_peaks_aligned) <- as.character(sample_manifest$source_id)
source_average_peak_qc <- sample_manifest[, .(
  sample_id = as.character(sample_id), diagnosis, source_id, batch_id,
  source_file = imzml_file
)]
source_average_peak_qc[, n_detected_average_spectrum_peaks := vapply(
  source_average_peaks_aligned, peak_count, integer(1)
)]
data.table::fwrite(
  source_average_peak_qc,
  file.path(qc_dir, "source_average_spectrum_peak_qc.csv")
)
saveRDS(
  source_average_aligned,
  file.path(out_dir, "aligned_source_average_spectra.rds"), compress = "xz"
)
saveRDS(
  source_average_peaks_aligned,
  file.path(out_dir, "aligned_source_average_peaks.rds"), compress = "xz"
)
consensus_features <- build_consensus_features(
  source_average_peaks = source_average_peaks_aligned,
  source_sample_id = sample_manifest$sample_id,
  mz_min = common_min,
  mz_max = common_max,
  max_span_da = peak_group_max_span_da,
  min_samples = consensus_min_samples
)
feature_mz_all <- consensus_features$mz
message("[共识 feature] 平均谱候选分组完成；feature 数：",
        length(feature_mz_all))

spectra_aligned <- unlist(lapply(source_objects, `[[`, "spectra"), recursive = FALSE)
pixel_meta <- data.table::rbindlist(lapply(source_objects, `[[`, "pixel_meta"),
                                    use.names = TRUE, fill = TRUE)
stopifnot(length(spectra_aligned) == nrow(pixel_meta))

peaks_aligned <- lapply(spectra_aligned, detect_one, snr = peak_snr)
pixel_meta[, detected_peak_count := vapply(peaks_aligned, peak_count, integer(1))]

## 后续矩阵只使用检测峰，不再对连续谱插值。尽早释放大型连续谱对象。
rm(
  spectra_aligned, source_objects, source_average, source_average_aligned,
  source_average_peaks_aligned, source_peaks_pre, source_peaks_post, source_warp
)
invisible(gc())

message("[像素峰分配] 开始；匹配到最近共识峰，半窗宽 = ",
        pixel_assignment_half_window_da, " Da。")
pixel_assignment_result <- assign_pixel_peaks_to_consensus(
  peaks_aligned,
  feature_mz = feature_mz_all,
  half_window_da = pixel_assignment_half_window_da
)
message("[像素峰分配] 完成。")
feat_mat_all <- pixel_assignment_result$matrix
rm(peaks_aligned)
invisible(gc())
stopifnot(nrow(feat_mat_all) == nrow(pixel_meta))
stopifnot(ncol(feat_mat_all) == length(feature_mz_all))

## 保存到此为止累积的警告，便于区分普通检峰警告与真正错误。
warning_snapshot <- warnings()
if (!is.null(warning_snapshot)) {
  capture.output(
    print(warning_snapshot),
    file = file.path(qc_dir, "warnings_through_feature_assignment.txt")
  )
}

## ==== 7) 按物理样本内出现率过滤 feature ====

nonzero_prop <- colSums(feat_mat_all != 0) / nrow(feat_mat_all)
n_nonzero_pixels <- colSums(feat_mat_all != 0)
sample_row_index <- split(seq_len(nrow(pixel_meta)), pixel_meta$sample_id)
sample_row_index <- sample_row_index[sample_levels]
sample_nonzero_count <- do.call(rbind, lapply(sample_row_index, function(idx) {
  colSums(feat_mat_all[idx, , drop = FALSE] != 0)
}))
rownames(sample_nonzero_count) <- sample_levels
sample_nonzero_prop <- sample_nonzero_count /
  matrix(
    lengths(sample_row_index),
    nrow = length(sample_levels), ncol = ncol(feat_mat_all)
  )
sample_prevalence_pass <-
  sample_nonzero_count >= min_nonzero_pixels_in_sample &
  sample_nonzero_prop >= min_within_sample_nonzero_prop
n_samples_passing_prevalence <- colSums(sample_prevalence_pass)
max_within_sample_nonzero_prop <- apply(sample_nonzero_prop, 2L, max)
sample_with_max_nonzero_prop <- sample_levels[
  max.col(t(sample_nonzero_prop), ties.method = "first")
]

feature_mean <- colMeans(feat_mat_all)
feature_mean_sq <- colMeans(feat_mat_all ^ 2)
## 与 stats::sd 接近的样本标准差，同时避免把稀疏矩阵转成稠密矩阵。
feature_var <- pmax(
  (feature_mean_sq - feature_mean ^ 2) * nrow(feat_mat_all) /
    max(nrow(feat_mat_all) - 1L, 1L),
  0
)
feature_sd <- sqrt(feature_var)
feature_cv   <- feature_sd / (feature_mean + 1e-12)

remove_sparse <- n_samples_passing_prevalence < min_samples_passing_prevalence
remove_quasi_const <- nonzero_prop > high_nonzero_prop &
  feature_cv < max_cv_quasi_const
keep_feature <- !(remove_sparse | remove_quasi_const)
if (sum(keep_feature) < 50L) {
  stop("按样本内出现率过滤后少于 50 个 feature，请检查平均谱检峰和过滤参数。")
}
if (!any(keep_feature)) stop("特征过滤后没有剩余特征，请检查检峰和过滤参数。")

feat_mat <- feat_mat_all[, keep_feature, drop = FALSE]
feature_mz <- feature_mz_all[keep_feature]
colnames(feat_mat) <- sprintf("mz_%.6f", feature_mz)

feature_qc <- merge(
  consensus_features,
  pixel_assignment_result$assignment_qc,
  by.x = "feature_id", by.y = "feature", all.x = TRUE, sort = FALSE
)
data.table::setorder(feature_qc, feature_id)
feature_qc[, `:=`(
  distance_to_previous_feature_da = c(NA_real_, diff(mz)),
  assignment_half_window_da = pixel_assignment_half_window_da,
  n_nonzero_pixels = n_nonzero_pixels,
  global_nonzero_prop = nonzero_prop,
  max_within_sample_nonzero_prop = max_within_sample_nonzero_prop,
  sample_with_max_nonzero_prop = sample_with_max_nonzero_prop,
  n_samples_passing_prevalence = n_samples_passing_prevalence,
  min_within_sample_nonzero_prop = min_within_sample_nonzero_prop,
  min_nonzero_pixels_in_sample = min_nonzero_pixels_in_sample,
  mean_intensity = feature_mean,
  cv = feature_cv,
  removed_sparse = remove_sparse,
  removed_quasi_constant = remove_quasi_const,
  retained = keep_feature
)]
data.table::fwrite(feature_qc, file.path(qc_dir, "feature_filter_qc.csv"))
feature_build_summary <- data.table::data.table(
  metric = c(
    "n_consensus_candidate_features",
    "n_retained_features",
    "fraction_candidate_groups_span_ge_99pct_limit",
    "median_candidate_group_span_da",
    "fraction_adjacent_feature_centers_lt_0.01Da",
    "median_adjacent_feature_center_distance_da",
    "assignment_half_window_da",
    "min_within_sample_nonzero_prop",
    "min_nonzero_pixels_in_sample"
  ),
  value = c(
    nrow(feature_qc),
    sum(feature_qc$retained),
    mean(feature_qc$consensus_candidate_span_da >=
           0.99 * peak_group_max_span_da),
    stats::median(feature_qc$consensus_candidate_span_da),
    mean(diff(feature_qc$mz) < 0.01),
    stats::median(diff(feature_qc$mz)),
    pixel_assignment_half_window_da,
    min_within_sample_nonzero_prop,
    min_nonzero_pixels_in_sample
  )
)
data.table::fwrite(
  feature_build_summary,
  file.path(qc_dir, "feature_build_summary.csv")
)
data.table::fwrite(
  data.table::data.table(
    sample_id = rep(sample_levels, each = ncol(feat_mat_all)),
    feature_id = rep(seq_len(ncol(feat_mat_all)), times = length(sample_levels)),
    mz = rep(feature_mz_all, times = length(sample_levels)),
    n_nonzero_pixels = as.vector(t(sample_nonzero_count)),
    within_sample_nonzero_prop = as.vector(t(sample_nonzero_prop)),
    passes_prevalence = as.vector(t(sample_prevalence_pass))
  ),
  file.path(qc_dir, "feature_prevalence_by_sample.csv")
)

## ==== 8) 保存对齐后的单批次数据（离子图应使用 TIC 矩阵） ====

pixel_meta[, sample_id := factor(sample_id, levels = sample_levels)]
pixel_meta[, pixel_row := seq_len(.N)]
data.table::setcolorder(pixel_meta, c(
  "pixel_row", "sample_id", "diagnosis", "source_id", "batch_id",
  "source_file", "pixel_index_in_source", "x", "y",
  "raw_tic", "alignment_nmatch", "alignment_failed", "detected_peak_count"
))

saveRDS(feat_mat,
        file.path(out_dir, "pixel_by_feature_matrix_aligned_TIC.rds"), compress = "xz")
saveRDS(feature_mz,
        file.path(out_dir, "feature_mz_aligned_filtered.rds"), compress = "xz")
saveRDS(pixel_meta,
        file.path(out_dir, "pixel_metadata_sample_source_xy.rds"), compress = "xz")
data.table::fwrite(pixel_meta,
                   file.path(out_dir, "pixel_metadata_sample_source_xy.csv"))

## ==== 9) 可选的多批次矫正函数（本次单批次不会执行） ====

batch_correct_log_matrix <- function(mat, meta, feature_chunk_size = 250L) {
  diagnosis <- factor(meta$diagnosis)
  batch_id  <- factor(meta$batch_id)
  design <- stats::model.matrix(~ 0 + diagnosis + batch_id)
  if (qr(design)$rank < ncol(design)) {
    stop(
      "diagnosis 与 batch_id 完全混杂，无法安全估计批次效应。请检查 manifest 中的 ",
      "batch_id；未生成批次矫正矩阵。"
    )
  }
  batch_cols <- grep("^batch_id", colnames(design))
  if (!length(batch_cols)) stop("batch_id 只有一个水平，无需批次矫正。")

  ## 每个物理样本总权重均为 1，避免像素多的样本主导批次估计。
  n_by_sample <- table(meta$sample_id)
  weights <- 1 / as.numeric(n_by_sample[as.character(meta$sample_id)])
  sqrt_w <- sqrt(weights)
  design_w <- design * sqrt_w
  design_qr <- qr(design_w)

  corrected <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  colnames(corrected) <- colnames(mat)
  starts <- seq.int(1L, ncol(mat), by = feature_chunk_size)
  for (start in starts) {
    end <- min(start + feature_chunk_size - 1L, ncol(mat))
    j <- start:end
    y <- as.matrix(sparse_log1p_cpm(mat[, j, drop = FALSE]))
    coef <- qr.coef(design_qr, y * sqrt_w)
    batch_effect <- design[, batch_cols, drop = FALSE] %*%
      coef[batch_cols, , drop = FALSE]
    ## 以等样本权重的平均批次效应为中心，避免整体强度平移到某个参考批次。
    batch_center <- colSums(batch_effect * weights) / sum(weights)
    corrected[, j] <- y - batch_effect +
      matrix(batch_center, nrow = nrow(y), ncol = length(j), byrow = TRUE)
    rm(y, coef, batch_effect)
    invisible(gc())
  }
  corrected
}

if (perform_batch_correction) {
  feat_mat_log_batch_corrected <- batch_correct_log_matrix(
    feat_mat, pixel_meta, feature_chunk_size = batch_feature_chunk_size
  )
  saveRDS(
    feat_mat_log_batch_corrected,
    file.path(out_dir, "pixel_by_feature_matrix_aligned_logCPM_batch_corrected.rds"),
    compress = "xz"
  )
} else {
  feat_mat_log_batch_corrected <- sparse_log1p_cpm(feat_mat)
  saveRDS(
    feat_mat_log_batch_corrected,
    file.path(out_dir, "pixel_by_feature_matrix_aligned_logCPM.rds"),
    compress = "xz"
  )
}

} else {
  message("[续跑] 跳过 imzML 读入、对齐和 feature matrix 构建。")
  resume_required_files <- c(
    file.path(out_dir, "pixel_by_feature_matrix_aligned_TIC.rds"),
    file.path(out_dir, "feature_mz_aligned_filtered.rds"),
    file.path(out_dir, "pixel_metadata_sample_source_xy.rds"),
    file.path(qc_dir, "feature_filter_qc.csv"),
    if (perform_batch_correction) {
      file.path(
        out_dir,
        "pixel_by_feature_matrix_aligned_logCPM_batch_corrected.rds"
      )
    } else {
      file.path(out_dir, "pixel_by_feature_matrix_aligned_logCPM.rds")
    }
  )
  missing_resume_files <- resume_required_files[!file.exists(resume_required_files)]
  if (length(missing_resume_files)) {
    stop(
      "无法从 feature matrix 续跑，缺少文件：\n",
      paste(missing_resume_files, collapse = "\n")
    )
  }

  feat_mat <- readRDS(
    file.path(out_dir, "pixel_by_feature_matrix_aligned_TIC.rds")
  )
  feature_mz <- readRDS(
    file.path(out_dir, "feature_mz_aligned_filtered.rds")
  )
  pixel_meta <- data.table::as.data.table(readRDS(
    file.path(out_dir, "pixel_metadata_sample_source_xy.rds")
  ))
  feature_qc_resume <- data.table::fread(
    file.path(qc_dir, "feature_filter_qc.csv")
  )
  feature_mz_all <- feature_qc_resume$mz
  common_min <- NA_real_
  common_max <- NA_real_
  feat_mat_log_batch_corrected <- readRDS(
    if (perform_batch_correction) {
      file.path(
        out_dir,
        "pixel_by_feature_matrix_aligned_logCPM_batch_corrected.rds"
      )
    } else {
      file.path(out_dir, "pixel_by_feature_matrix_aligned_logCPM.rds")
    }
  )
  stopifnot(
    nrow(feat_mat) == nrow(pixel_meta),
    ncol(feat_mat) == length(feature_mz),
    all(dim(feat_mat_log_batch_corrected) == dim(feat_mat))
  )
  message("[续跑] 已读入 ", nrow(feat_mat), " 个像素 × ",
          ncol(feat_mat), " 个保留 feature。")
}

## ==== 10) 10 个样本的样本级 QC ====

sample_qc <- pixel_meta[, .(
  n_sources = uniqueN(source_id),
  n_pixels = .N,
  median_raw_tic = as.numeric(median(raw_tic, na.rm = TRUE)),
  median_detected_peaks = as.numeric(
    median(detected_peak_count, na.rm = TRUE)
  ),
  median_alignment_matches = as.numeric(
    median(alignment_nmatch, na.rm = TRUE)
  ),
  alignment_failure_rate = mean(alignment_failed, na.rm = TRUE)
), by = .(sample_id, diagnosis)]

source_pixel_qc <- pixel_meta[, .(
  n_pixels = .N,
  median_raw_tic = as.numeric(median(raw_tic, na.rm = TRUE)),
  median_detected_peaks = as.numeric(
    median(detected_peak_count, na.rm = TRUE)
  ),
  median_alignment_matches = as.numeric(
    median(alignment_nmatch, na.rm = TRUE)
  ),
  alignment_failure_rate = mean(alignment_failed, na.rm = TRUE)
), by = .(sample_id, diagnosis, source_id, batch_id, source_file)]
data.table::fwrite(source_pixel_qc,
                   file.path(qc_dir, "source_pixel_qc_12imzML.csv"))

## 输入 imzML 的 TIC 在同一来源内可能已被归一化为常数，因此 TIC 尺度
## 异常只单独记录，不再自动触发样本排除警报。
sample_qc[, z_log_tic := robust_z(log10(median_raw_tic + 1))]
sample_qc[, z_log_peak_count := robust_z(log10(median_detected_peaks + 1))]
sample_qc[, z_alignment_matches := robust_z(median_alignment_matches)]
sample_qc[, tic_scale_outlier := abs(z_log_tic) > 3]
sample_qc[, qc_flag :=
  z_log_peak_count < -3 |
  z_alignment_matches < -3 |
  alignment_failure_rate > 0.20]

data.table::fwrite(sample_qc, file.path(qc_dir, "sample_qc_10samples.csv"))
saveRDS(sample_qc, file.path(qc_dir, "sample_qc_10samples.rds"))

## PCA 默认关闭；先检查 feature_filter_qc 和 feature_prevalence_by_sample。
if (run_downstream_pca) {
  message("[PCA] 开关已开启，开始计算 10 个物理样本的样本均值 PCA。")
  sample_index <- split(
    seq_len(nrow(pixel_meta)), as.character(pixel_meta$sample_id)
  )
  sample_index <- sample_index[sample_levels]
  sample_mean_before <- do.call(rbind, lapply(sample_index, function(idx) {
    colMeans(sparse_log1p_cpm(feat_mat[idx, , drop = FALSE]))
  }))
  sample_mean_after <- do.call(rbind, lapply(sample_index, function(idx) {
    colMeans(feat_mat_log_batch_corrected[idx, , drop = FALSE], na.rm = TRUE)
  }))

  run_sample_pca <- function(x, stage) {
    ok <- apply(x, 2, stats::sd, na.rm = TRUE) > 0
    p <- stats::prcomp(x[, ok, drop = FALSE], center = TRUE, scale. = TRUE)
    ans <- data.table::data.table(
      sample_id = rownames(x),
      PC1 = p$x[, 1],
      PC2 = p$x[, 2],
      stage = stage
    )
    ans[sample_manifest[, .(
      diagnosis = diagnosis[1], batch_id = batch_id[1]
    ), by = sample_id],
    on = "sample_id", `:=`(
      diagnosis = i.diagnosis, batch_id = i.batch_id
    )]
    ans
  }

  pca_dt <- if (perform_batch_correction) {
    data.table::rbindlist(list(
      run_sample_pca(
        sample_mean_before, "Aligned + TIC (before batch correction)"
      ),
      run_sample_pca(
        sample_mean_after, "Aligned + log-CPM batch corrected"
      )
    ))
  } else {
    run_sample_pca(
      sample_mean_before, "Aligned + TIC + log-CPM (single batch)"
    )
  }
  data.table::fwrite(
    pca_dt, file.path(qc_dir, "sample_level_PCA_scores.csv")
  )

  p_pca <- ggplot2::ggplot(
    pca_dt,
    ggplot2::aes(
      PC1, PC2, colour = diagnosis, shape = batch_id, label = sample_id
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(nudge_y = 0.25, check_overlap = TRUE, size = 3) +
    ggplot2::facet_wrap(~stage, scales = "free") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "10 physical samples: aligned-data PCA QC")
  ggplot2::ggsave(
    file.path(qc_dir, "sample_level_PCA_QC.pdf"),
    p_pca, width = 12, height = 6
  )
} else {
  message("[PCA] run_downstream_pca=FALSE；本次不运行 PCA。")
}

## ==== 11) 运行记录 ====

run_parameters <- list(
  created_at = Sys.time(),
  data_root = data_root,
  out_dir = out_dir,
  sample_levels = sample_levels,
  common_mz_range = c(common_min, common_max),
  max_pixels_per_source = max_pixels_per_source,
  smooth_half_window = smooth_half_window,
  baseline_method = baseline_method,
  peak_snr = peak_snr,
  peak_half_window = peak_half_window,
  consensus_peak_snr = consensus_peak_snr,
  consensus_min_samples = consensus_min_samples,
  alignment_search_ppm = alignment_search_ppm,
  peak_group_max_span_da = peak_group_max_span_da,
  pixel_assignment_half_window_da = pixel_assignment_half_window_da,
  min_within_sample_nonzero_prop = min_within_sample_nonzero_prop,
  min_nonzero_pixels_in_sample = min_nonzero_pixels_in_sample,
  min_samples_passing_prevalence = min_samples_passing_prevalence,
  final_bin_rule = paste(
    "consensus features from aligned source-average spectra;",
    "each physical sample contributes at most one representative peak per group;",
    "feature m/z = median of physical-sample representative m/z;",
    "pixel peaks assigned to nearest center within +/- half window;",
    "max intensity per pixel/feature"
  ),
  matrix_storage = "Matrix::dgCMatrix sparse matrix; undetected peaks stored as zero",
  reference_min_frequency = reference_min_frequency,
  warping_method = warping_method,
  perform_pixel_alignment = perform_pixel_alignment,
  perform_batch_correction = perform_batch_correction,
  resume_after_feature_matrix = resume_after_feature_matrix,
  run_downstream_pca = run_downstream_pca,
  session_info = utils::sessionInfo()
)
saveRDS(run_parameters, file.path(out_dir, "run_parameters_and_sessionInfo.rds"))

message("完成。输入清单仅包含 10 个物理样本、12 个 imzML 来源，不读取 Cushing-3。")
message("主要输出目录：", out_dir)
message("保留 feature 数：", length(feature_mz),
        " / 共识候选 feature 数：", length(feature_mz_all))
print(sample_qc)
