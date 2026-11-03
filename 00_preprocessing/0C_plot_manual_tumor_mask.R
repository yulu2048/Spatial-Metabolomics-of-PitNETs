## ============================================================================
## 0C_plot_manual_tumor_mask.R
##
## 目的：
##   根据 0B 的 k=4 cluster 结果和人工标注，绘制二分类肿瘤像素图：
##     - 肿瘤像素：红色
##     - 非肿瘤像素：灰色
##   Cushing-3 未在上游输入清单中；本脚本处理全部 10 个样本、12 个来源。
##
## 本脚本不重新聚类、不修改原始矩阵或 0B 结果；
## 同时保存人工复核后的肿瘤 mask，供 PCA 和后续分析使用。
## ============================================================================

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("缺少 R 包：", paste(missing_packages, collapse = ", "),
       "。请先安装后再运行。")
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ==== 1) 路径 ====

data_root <- "/home/gaga/huzhan"
input_dir <- file.path(
  data_root, "MSI_ALL_aligned_10samples_consensus_tuned225"
)
k4_dir <- file.path(input_dir, "QC", "tumor_pixel_selection_0B_k4")
map_dir <- file.path(k4_dir, "manual_tumor_mask_maps")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

cluster_file <- file.path(k4_dir, "pixel_clusters_12sources.rds")
manifest_file <- file.path(input_dir, "sample_source_manifest.csv")

missing_files <- c(cluster_file, manifest_file)[
  !file.exists(c(cluster_file, manifest_file))
]
if (length(missing_files)) {
  stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))
}

## ==== 2) 人工肿瘤 cluster 标注 ====

## 名字是 source_id；cluster 编号来自 0B 的 k=4 结果。
tumor_cluster_map <- list(
  cushing_0 = c(3L, 4L),
  Cushing = c(1L, 2L),
  gh_1 = 2L,
  gh_2 = c(1L, 2L, 3L),
  GH = c(2L, 4L),
  NF_ACTH = c(1L, 3L),
  nf_cxx_1 = c(2L, 4L),
  nf_cxx_2 = c(3L, 4L),
  nf_cxx_3 = c(2L, 4L),
  NF_CXX = 3L,
  nf_cxx = c(2L, 4L),
  NF_GH = c(1L, 3L)
)

save_png_dpi <- 240

## ==== 3) 读入与核对 ====

pixel_dt <- data.table::as.data.table(readRDS(cluster_file))
manifest <- data.table::fread(manifest_file)

required_columns <- c(
  "pixel_row", "sample_id", "diagnosis", "source_id", "x", "y", "cluster"
)
missing_columns <- setdiff(required_columns, names(pixel_dt))
if (length(missing_columns)) {
  stop("0B 像素结果缺少列：", paste(missing_columns, collapse = ", "))
}

source_levels <- as.character(manifest$source_id)
if (length(source_levels) != 12L || uniqueN(manifest$sample_id) != 10L) {
  stop("预期输入为12个来源、10个物理样本。")
}

missing_annotations <- setdiff(source_levels, names(tumor_cluster_map))
extra_annotations <- setdiff(names(tumor_cluster_map), source_levels)
if (length(missing_annotations) || length(extra_annotations)) {
  stop(
    if (length(missing_annotations)) {
      paste0("缺少人工标注：", paste(missing_annotations, collapse = ", "), "\n")
    } else "",
    if (length(extra_annotations)) {
      paste0("存在多余人工标注：", paste(extra_annotations, collapse = ", "))
    } else ""
  )
}

pixel_dt <- pixel_dt[source_id %chin% source_levels]

pixel_dt[, is_tumor := FALSE]
for (sid in names(tumor_cluster_map)) {
  selected_clusters <- tumor_cluster_map[[sid]]
  pixel_dt[
    source_id == sid & cluster %in% selected_clusters,
    is_tumor := TRUE
  ]
}

mask_file <- file.path(k4_dir, "manual_tumor_mask_12sources.rds")
saveRDS(pixel_dt, mask_file, compress = "xz")
data.table::fwrite(
  pixel_dt[, .(pixel_row, sample_id, diagnosis, source_id, x, y, cluster,
    is_tumor)],
  file.path(k4_dir, "manual_tumor_mask_12sources.csv")
)
saveRDS(
  as.integer(pixel_dt[is_tumor == TRUE, pixel_row]),
  file.path(k4_dir, "tumor_pixel_rows_10samples.rds"),
  compress = "xz"
)

## 因子顺序决定图例顺序：Tumor 在上，Non-tumor 在下。
pixel_dt[, pixel_type := factor(
  ifelse(is_tumor, "Tumor", "Non-tumor"),
  levels = c("Tumor", "Non-tumor")
)]

## ==== 4) 绘图函数 ====

safe_filename <- function(x) gsub("[^A-Za-z0-9_-]+", "_", x)

plot_one_source <- function(dt) {
  sample_label <- as.character(unique(dt$sample_id))
  source_label <- as.character(unique(dt$source_id))

  ggplot(dt, aes(x = x, y = y, fill = pixel_type)) +
    geom_tile(width = 1, height = 1) +
    scale_y_reverse() +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(
      values = c("Tumor" = "#D73027", "Non-tumor" = "#BDBDBD"),
      breaks = c("Tumor", "Non-tumor"),
      drop = FALSE
    ) +
    labs(
      title = paste0(sample_label, " / ", source_label),
      x = "x (pixel)",
      y = "y (pixel)",
      fill = "Pixel type"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

## ==== 5) 保存每个来源和总览 PDF ====

plot_list <- vector("list", length(source_levels))
names(plot_list) <- source_levels

for (i in seq_along(source_levels)) {
  sid <- source_levels[i]
  dt <- pixel_dt[source_id == sid]
  if (!nrow(dt)) stop("没有找到来源 ", sid, " 的像素。")

  message(
    "[绘图] ", i, "/", length(source_levels), ": ",
    unique(dt$sample_id), " / ", sid,
    "；肿瘤像素 = ", sum(dt$is_tumor),
    "；非肿瘤像素 = ", sum(!dt$is_tumor)
  )

  p <- plot_one_source(dt)
  stub <- safe_filename(sid)

  ggsave(
    file.path(map_dir, paste0(stub, "_manual_tumor_mask.png")),
    p, width = 7.2, height = 7.2, units = "in", dpi = save_png_dpi
  )

  ## 优先使用 Cairo PDF，以提高中文图例在 Linux 服务器上的兼容性。
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggsave(
    file.path(map_dir, paste0(stub, "_manual_tumor_mask.pdf")),
    p, width = 7.2, height = 7.2, units = "in", device = pdf_device
  )

  plot_list[[i]] <- p
}

overview_pdf <- file.path(k4_dir, "all_12sources_manual_tumor_mask_maps.pdf")
if (capabilities("cairo")) {
  grDevices::cairo_pdf(overview_pdf, width = 7.2, height = 7.2, onefile = TRUE)
} else {
  grDevices::pdf(overview_pdf, width = 7.2, height = 7.2, onefile = TRUE)
}
for (p in plot_list) print(p)
grDevices::dev.off()

message("完成。输入清单仅包含10个样本、12个来源，不读取Cushing-3。")
message("肿瘤 mask：", mask_file)
message("12 页总览：", overview_pdf)
message("单独图像目录：", map_dir)
