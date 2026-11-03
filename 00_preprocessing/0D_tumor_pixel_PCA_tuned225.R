## Tumor-only PCA after manual k=4 tumor-pixel annotation
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

set.seed(2025)
base_dir <- "/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225"
k4_dir <- file.path(base_dir, "QC", "tumor_pixel_selection_0B_k4")
out_dir <- file.path(base_dir, "QC", "tumor_only_PCA_10samples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_file <- file.path(base_dir, "pixel_by_feature_matrix_aligned_logCPM.rds")
metadata_file <- file.path(base_dir, "pixel_metadata_sample_source_xy.rds")
mask_file <- file.path(k4_dir, "manual_tumor_mask_12sources.rds")
mz_file <- file.path(base_dir, "feature_mz_aligned_filtered.rds")
required_files <- c(matrix_file, metadata_file, mask_file, mz_file)
if (any(!file.exists(required_files))) {
  stop("Missing input files:\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

X <- readRDS(matrix_file)
pixel_meta <- as.data.table(readRDS(metadata_file))
mask <- as.data.table(readRDS(mask_file))
feature_mz <- readRDS(mz_file)
stopifnot(nrow(X) == nrow(pixel_meta), ncol(X) == length(feature_mz))
stopifnot(all(mask$pixel_row >= 1L & mask$pixel_row <= nrow(X)))

analysis_mask <- mask[is_tumor == TRUE]
if (anyDuplicated(analysis_mask$pixel_row)) stop("Duplicated pixel_row in tumor mask.")
analysis_mask[, sample_id := as.character(sample_id)]
analysis_mask[, diagnosis := as.character(diagnosis)]
analysis_mask[, source_id := as.character(source_id)]

sample_order <- c(
  "GH-1", "GH-2", "GH-3", "NF-GH", "NF-ACTH",
  "Cushing-1", "Cushing-2", "NF-GnH-1", "NF-GnH-2", "NF-GnH-3"
)
sample_order <- sample_order[sample_order %in% unique(analysis_mask$sample_id)]
source_order <- unique(as.character(analysis_mask$source_id))

mean_matrix <- function(group_vector, group_order, rows) {
  ans <- lapply(group_order, function(g) {
    ii <- rows[group_vector == g]
    if (!length(ii)) stop("No pixels for group: ", g)
    Matrix::colMeans(X[ii, , drop = FALSE])
  })
  ans <- do.call(rbind, ans)
  rownames(ans) <- group_order
  ans
}

tumor_sample_mean <- mean_matrix(
  analysis_mask$sample_id, sample_order, analysis_mask$pixel_row
)
tumor_source_mean <- mean_matrix(
  analysis_mask$source_id, source_order, analysis_mask$pixel_row
)
nonexcluded_rows <- seq_len(nrow(pixel_meta))
all_sample_mean <- mean_matrix(
  as.character(pixel_meta$sample_id)[nonexcluded_rows],
  sample_order, nonexcluded_rows
)

feature_sd <- apply(tumor_sample_mean, 2L, sd, na.rm = TRUE)
keep_feature <- is.finite(feature_sd) & feature_sd > 1e-8
if (sum(keep_feature) < 3L) stop("Too few variable features for PCA.")

fit_one_pca <- function(mat, scaling) {
  mat <- mat[, keep_feature, drop = FALSE]
  s <- apply(mat, 2L, sd)
  scale_arg <- switch(
    scaling,
    center = FALSE,
    pareto = sqrt(s),
    autoscale = s,
    stop("Unknown scaling: ", scaling)
  )
  fit <- prcomp(mat, center = TRUE, scale. = scale_arg)
  variance <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  list(fit = fit, variance = variance, scaling = scaling)
}

diagnosis_map <- unique(analysis_mask[, .(sample_id, diagnosis)])
source_map <- unique(analysis_mask[, .(sample_id, diagnosis, source_id)])

score_table <- function(obj, map, id_name) {
  z <- as.data.table(
    obj$fit$x[, seq_len(min(5L, ncol(obj$fit$x))), drop = FALSE],
    keep.rownames = id_name
  )
  z <- merge(z, map, by = id_name, all.x = TRUE, sort = FALSE)
  z[, scaling := obj$scaling]
  z
}

pca_objects <- lapply(
  c("center", "pareto", "autoscale"),
  function(s) fit_one_pca(tumor_sample_mean, s)
)
names(pca_objects) <- c("center", "pareto", "autoscale")
sample_scores <- rbindlist(
  lapply(pca_objects, score_table, map = diagnosis_map, id_name = "sample_id"),
  fill = TRUE
)

source_scores <- rbindlist(lapply(pca_objects, function(obj) {
  z <- predict(obj$fit, newdata = tumor_source_mean[, keep_feature, drop = FALSE])
  z <- as.data.table(
    z[, seq_len(min(5L, ncol(z))), drop = FALSE],
    keep.rownames = "source_id"
  )
  z <- merge(z, source_map, by = "source_id", all.x = TRUE, sort = FALSE)
  z[, scaling := obj$scaling]
  z
}), fill = TRUE)

all_pca_objects <- lapply(
  c("center", "pareto", "autoscale"),
  function(s) fit_one_pca(all_sample_mean, s)
)
names(all_pca_objects) <- names(pca_objects)
all_sample_scores <- rbindlist(
  lapply(all_pca_objects, score_table, map = diagnosis_map, id_name = "sample_id"),
  fill = TRUE
)
all_sample_scores[, pixel_selection := "All nonempty pixels"]
sample_scores[, pixel_selection := "Tumor pixels"]
comparison_scores <- rbindlist(list(sample_scores, all_sample_scores), fill = TRUE)

## Analytical standard errors of sample centroids on the primary Pareto axes.
primary <- pca_objects$pareto$fit
pc_rotation <- primary$rotation[, 1:2, drop = FALSE]
center_vec <- primary$center
scale_vec <- primary$scale
centroid_uncertainty <- rbindlist(lapply(sample_order, function(sid) {
  rr <- analysis_mask[sample_id == sid, pixel_row]
  xx <- as.matrix(X[rr, keep_feature, drop = FALSE])
  xx <- sweep(xx, 2L, center_vec, "-")
  if (!identical(scale_vec, FALSE)) xx <- sweep(xx, 2L, scale_vec, "/")
  sc <- xx %*% pc_rotation
  data.table(
    sample_id = sid,
    n_tumor_pixels = length(rr),
    PC1_pixel_sd = sd(sc[, 1]),
    PC2_pixel_sd = sd(sc[, 2]),
    PC1_centroid_se = sd(sc[, 1]) / sqrt(length(rr)),
    PC2_centroid_se = sd(sc[, 2]) / sqrt(length(rr))
  )
}))
centroid_uncertainty <- merge(
  centroid_uncertainty, diagnosis_map, by = "sample_id", all.x = TRUE
)

primary_scores <- sample_scores[scaling == "pareto"]
primary_scores <- merge(
  primary_scores, centroid_uncertainty,
  by = c("sample_id", "diagnosis"), all.x = TRUE
)
primary_scores[, distance_PC1_PC2 := sqrt(PC1^2 + PC2^2)]

qc_columns <- intersect(
  c("raw_tic", "detected_peaks", "alignment_matches"), names(pixel_meta)
)
meta_rows <- pixel_meta[analysis_mask$pixel_row]
technical_qc <- meta_rows[, c(
  list(n_tumor_pixels = .N, n_sources = uniqueN(source_id)),
  lapply(.SD, median, na.rm = TRUE)
), by = .(
  sample_id = as.character(sample_id),
  diagnosis = as.character(diagnosis)
), .SDcols = qc_columns]
technical_qc <- merge(
  technical_qc,
  primary_scores[, .(sample_id, PC1, PC2, distance_PC1_PC2)],
  by = "sample_id", all.x = TRUE
)

qc_predictors <- intersect(
  c("n_tumor_pixels", "raw_tic", "detected_peaks", "alignment_matches"),
  names(technical_qc)
)
pc_qc_correlations <- rbindlist(lapply(qc_predictors, function(v) {
  rbindlist(lapply(c("PC1", "PC2"), function(pc) {
    ok <- is.finite(technical_qc[[v]]) & is.finite(technical_qc[[pc]])
    ct <- suppressWarnings(cor.test(
      technical_qc[[v]][ok], technical_qc[[pc]][ok],
      method = "spearman", exact = FALSE
    ))
    data.table(
      predictor = v, PC = pc, n = sum(ok),
      spearman_rho = unname(ct$estimate), p_value = ct$p.value
    )
  }))
}))

variance_table <- rbindlist(lapply(names(pca_objects), function(s) {
  data.table(
    scaling = s,
    PC = paste0("PC", seq_along(pca_objects[[s]]$variance)),
    variance_percent = pca_objects[[s]]$variance
  )
}))

saveRDS(
  tumor_sample_mean,
  file.path(out_dir, "tumor_sample_mean_logCPM_matrix.rds"),
  compress = "xz"
)
saveRDS(
  tumor_source_mean,
  file.path(out_dir, "tumor_source_mean_logCPM_matrix.rds"),
  compress = "xz"
)
saveRDS(
  pca_objects,
  file.path(out_dir, "tumor_sample_PCA_objects.rds"),
  compress = "xz"
)
fwrite(sample_scores, file.path(out_dir, "tumor_sample_PCA_scores_all_scalings.csv"))
fwrite(source_scores, file.path(out_dir, "tumor_source_projected_PCA_scores.csv"))
fwrite(comparison_scores, file.path(out_dir, "tumor_vs_allpixel_PCA_scores.csv"))
fwrite(centroid_uncertainty, file.path(out_dir, "tumor_PCA_centroid_uncertainty.csv"))
fwrite(technical_qc, file.path(out_dir, "tumor_sample_technical_QC_with_PCA.csv"))
fwrite(pc_qc_correlations, file.path(out_dir, "PCA_technical_metric_correlations.csv"))
fwrite(variance_table, file.path(out_dir, "PCA_explained_variance.csv"))

palette <- c(
  "GH" = "#2C7BB6", "Cushing" = "#D7191C", "NF-GH" = "#1A9641",
  "NF-ACTH" = "#984EA3", "NF-GnH" = "#FDAE61"
)

plot_pca <- function(dt, variance, title, label_col, errorbars = FALSE) {
  p <- ggplot(dt, aes(x = PC1, y = PC2, color = diagnosis))
  if (errorbars) {
    p <- p +
      geom_errorbar(
        aes(ymin = PC2 - 1.96 * PC2_centroid_se,
            ymax = PC2 + 1.96 * PC2_centroid_se),
        width = 0
      ) +
      geom_errorbarh(
        aes(xmin = PC1 - 1.96 * PC1_centroid_se,
            xmax = PC1 + 1.96 * PC1_centroid_se),
        height = 0
      )
  }
  p +
    geom_point(size = 3.2) +
    geom_text(
      aes(label = .data[[label_col]]),
      hjust = -0.08, vjust = -0.55, size = 3.2, show.legend = FALSE
    ) +
    scale_color_manual(values = palette, drop = FALSE) +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", variance[1]),
      y = sprintf("PC2 (%.1f%%)", variance[2]),
      color = "Diagnosis"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf
for (s in names(pca_objects)) {
  dt <- sample_scores[scaling == s]
  if (s == "pareto") dt <- primary_scores
  p <- plot_pca(
    dt, pca_objects[[s]]$variance,
    paste0("Tumor-only sample PCA (", s, ")"),
    "sample_id", s == "pareto"
  )
  ggsave(
    file.path(out_dir, paste0("tumor_sample_PCA_", s, ".png")),
    p, width = 8, height = 6, dpi = 300
  )
  ggsave(
    file.path(out_dir, paste0("tumor_sample_PCA_", s, ".pdf")),
    p, width = 8, height = 6, device = pdf_device
  )
}

source_primary <- source_scores[scaling == "pareto"]
source_primary[, source_label := paste0(sample_id, " / ", source_id)]
p_source <- plot_pca(
  source_primary, pca_objects$pareto$variance,
  "Tumor-only source projection (Pareto)", "source_label"
)
ggsave(
  file.path(out_dir, "tumor_source_projection_PCA_pareto.png"),
  p_source, width = 9, height = 7, dpi = 300
)
ggsave(
  file.path(out_dir, "tumor_source_projection_PCA_pareto.pdf"),
  p_source, width = 9, height = 7, device = pdf_device
)

for (s in names(all_pca_objects)) {
  dt <- comparison_scores[scaling == s]
  p_t <- plot_pca(
    dt[pixel_selection == "Tumor pixels"], pca_objects[[s]]$variance,
    paste0("Tumor pixels (", s, ")"), "sample_id"
  )
  p_a <- plot_pca(
    dt[pixel_selection == "All nonempty pixels"], all_pca_objects[[s]]$variance,
    paste0("All nonempty pixels (", s, ")"), "sample_id"
  )
  ggsave(
    file.path(out_dir, paste0("comparison_tumor_", s, ".png")),
    p_t, width = 8, height = 6, dpi = 300
  )
  ggsave(
    file.path(out_dir, paste0("comparison_allpixels_", s, ".png")),
    p_a, width = 8, height = 6, dpi = 300
  )
}

cat("Tumor-only PCA completed.\n")
cat(
  "Tumor pixels:", nrow(analysis_mask), "from", length(sample_order),
  "physical samples and", length(source_order), "sources.\n"
)
cat("Variable features:", sum(keep_feature), "of", length(feature_mz), "\n")
cat(
  "Primary Pareto variance: PC1=", round(pca_objects$pareto$variance[1], 2),
  "%, PC2=", round(pca_objects$pareto$variance[2], 2), "%\n", sep = ""
)
print(primary_scores[order(-distance_PC1_PC2), .(
  sample_id, diagnosis, PC1, PC2, n_tumor_pixels,
  PC1_centroid_se, PC2_centroid_se, distance_PC1_PC2
)])
cat("Output:", out_dir, "\n")
