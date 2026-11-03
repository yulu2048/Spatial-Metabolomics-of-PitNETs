## Redraw the retained metabolic-module panels with fewer heatmaps.
## This script is intended to run in the RStudio session that contains
## pB_new, rob, profiles, relations and branch_summary.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

module_root <- "/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225/spatial_metabolite_modules_10samples_tuned225"
out_dir <- file.path(module_root, "publication_modules_reduced_heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_levels <- sprintf("M%02d", 1:15)
module_cols <- setNames(grDevices::hcl.colors(15, "Dynamic"), module_levels)
branch_cols <- c(B1 = "#1B9E77", B2 = "#D95F02", B3 = "#7570B3", B4 = "#E7298A")

theme_pub <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle = element_text(size = 10.5, colour = "#444444"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    strip.background = element_rect(fill = "#F3F3F3", colour = "#BDBDBD"),
    strip.text = element_text(face = "bold", colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

save_pdf <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(out_dir, filename), plot = plot,
    width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
}

## A: integrate module cohesion and leave-one-sample-out reproducibility in one bubble map.
module_summary <- fread(file.path(module_root, "step3_modules/module_summary.csv"))
a_dat <- merge(module_summary, as.data.table(rob), by = "module", all.x = TRUE)
a_dat[, reconstruction_fraction := n_LOSO_ge_0.5 / 10]

pA <- ggplot(
  a_dat,
  aes(median_within_consensus_rho, median_LOSO_Jaccard)
) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "#BDBDBD") +
  geom_point(
    aes(size = n_representative_features, fill = reconstruction_fraction),
    shape = 21, stroke = 0.7, colour = "white", alpha = 0.95
  ) +
  geom_text(aes(label = module), nudge_y = 0.035, size = 24 / .pt,
            fontface = "plain", check_overlap = FALSE) +
  scale_size_continuous(range = c(4.0, 11), breaks = c(10, 20, 40),
                        name = "Features") +
  scale_fill_gradient(low = "#C7E9C0", high = "#006D2C",
                      limits = c(0, 1), breaks = c(0, 0.5, 1),
                      labels = scales::percent,
                      name = "LOSO reconstructions\nwith Jaccard >= 0.5") +
  coord_cartesian(xlim = c(0.4, 0.72), ylim = c(0.35, 1.08),
                  expand = FALSE, clip = "off") +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Median within-module consensus Spearman rho",
    y = "Median leave-one-sample-out Jaccard"
  ) +
  theme_pub +
  guides(
    size = guide_legend(
      order = 1, title.position = "top", title.hjust = 0.5,
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = 21, fill = "#A9D39E",
                          colour = "#4D4D4D", alpha = 1, stroke = 0.8)
    ),
    fill = guide_colourbar(
      order = 2, title.position = "top", title.hjust = 0.5,
      barwidth = unit(5.4, "cm"), barheight = unit(0.8, "cm")
    )
  ) +
  theme(
    text = element_text(size = 24, face = "plain"),
    axis.title = element_text(size = 24, face = "plain"),
    axis.text = element_text(size = 24, face = "plain", colour = "black"),
    legend.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24, face = "plain"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.direction = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.spacing.y = unit(0.08, "cm"),
    panel.grid.major = element_line(colour = "#EFEFEF", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )
save_pdf(pA, "Panel_A_module_cohesion_and_robustness.pdf", 11, 9)

## B: retain the spatial maps that demonstrate distinct/complementary territories.
pB <- pB_new +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10.5),
        strip.text = element_text(face = "bold"))
save_pdf(pB, "Panel_B_representative_spatial_module_maps.pdf", 11.5, 7.5)

## C: hierarchically clustered sample-by-module heatmap of Moran's I.
c_dat <- as.data.table(pC2$data)
c_long <- copy(c_dat)[, .(
  sample = as.character(sample),
  module = as.character(module),
  moran_I
)]
c_wide <- dcast(c_long, sample ~ module, value.var = "moran_I")
c_mat <- as.matrix(c_wide[, -1])
rownames(c_mat) <- c_wide$sample
c_mat <- c_mat[, sprintf("M%02d", 1:15), drop = FALSE]

c_row_hc <- hclust(dist(c_mat), method = "complete")
c_col_hc <- hclust(dist(t(c_mat)), method = "complete")
c_col_fun <- colorRamp2(
  c(0, 0.2, 0.4, 0.6),
  c("#0D0887", "#9C179E", "#ED7953", "#F0F921")
)

htC <- Heatmap(
  c_mat,
  name = "Moran's I",
  col = c_col_fun,
  cluster_rows = c_row_hc,
  cluster_columns = c_col_hc,
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE,
  row_dend_side = "left",
  row_names_side = "right",
  column_dend_side = "top",
  column_names_side = "bottom",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 24, fontface = "plain", col = "black"),
  column_names_gp = gpar(fontsize = 24, fontface = "plain", col = "black"),
  column_names_rot = 90,
  row_dend_width = unit(18, "mm"),
  column_dend_height = unit(16, "mm"),
  width = unit(7.5, "in"),
  height = unit(4.2, "in"),
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "Moran's I",
    at = c(0, 0.2, 0.4, 0.6),
    labels = c("0.0", "0.2", "0.4", "0.6"),
    title_position = "topcenter",
    title_gap = unit(4, "mm"),
    title_gp = gpar(fontsize = 24, fontface = "plain"),
    labels_gp = gpar(fontsize = 24, fontface = "plain"),
    legend_height = unit(1.35, "in")
  )
)
panel_c_pdf <- file.path(out_dir, "Panel_C_module_MoransI_distribution.pdf")
cairo_pdf(panel_c_pdf, width = 12.5, height = 7, family = "sans", onefile = TRUE)
draw(htC, heatmap_legend_side = "right", padding = unit(c(6, 6, 6, 6), "mm"))
dev.off()

## D1: retain one heatmap because the complete D-domain x module composition is the key structure.
d1 <- as.data.table(profiles)
d1[, branch := factor(as.character(branch), levels = names(branch_cols))]
d1[, module := factor(module, levels = module_levels)]
d1[, domain_id := factor(as.character(domain_id), levels = levels(profiles$domain_id))]
d1[, activity_plot := pmax(-2, pmin(2, centered_stable_mean))]

pD1 <- ggplot(d1, aes(domain_id, module, fill = activity_plot)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ branch, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-2, 2), name = "Centered activity"
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Sample-specific D-domain", y = "Metabolic module"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    text = element_text(size = 24, face = "plain"),
    axis.title = element_text(size = 24, face = "plain"),
    axis.text.x = element_text(angle = 55, hjust = 1, colour = "black", size = 24, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 24, face = "plain"),
    panel.grid = element_blank(),
    panel.spacing.x = unit(1.5, "mm"),
    strip.background = element_rect(fill = "#F3F3F3", colour = "#BDBDBD"),
    strip.text = element_text(size = 24, face = "plain", colour = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24, face = "plain"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
save_pdf(pD1, "Panel_D1_Ddomain_module_profile_heatmap.pdf", 12, 11)

## D2: show tumor-level eta-squared values as a dot cloud without interval bars.
d2 <- as.data.table(relations)
d2[, module := factor(module, levels = module_levels)]
supported_modules <- c("M01", "M03", "M04", "M06", "M07", "M08", "M10", "M11")
d2_sum <- d2[, .(
  median_eta2 = median(eta2_validation_stable, na.rm = TRUE),
  q25 = quantile(eta2_validation_stable, 0.25, na.rm = TRUE),
  q75 = quantile(eta2_validation_stable, 0.75, na.rm = TRUE),
  supported = as.character(module[1]) %in% supported_modules
), by = module]

pD2 <- ggplot(d2, aes(module, eta2_validation_stable)) +
  geom_hline(yintercept = 0.25, linetype = 2, colour = "#737373") +
  geom_jitter(width = 0.16, height = 0, shape = 16, size = 2.4,
              colour = "#8C8C8C", alpha = 1) +
  geom_point(data = d2_sum, aes(x = module, y = median_eta2, fill = supported),
             shape = 23, size = 9.6, stroke = 0.55, colour = "white",
             inherit.aes = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "#969696"),
                    labels = c(`TRUE` = "Reproducible D-domain dependence", `FALSE` = "Other module"),
                    name = NULL) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Metabolic module", y = expression("Validation-feature " * eta^2)
  ) +
  theme_classic(base_size = 24) +
  theme(
    text = element_text(size = 24, face = "plain"),
    axis.title = element_text(size = 24, face = "plain"),
    axis.text = element_text(size = 24, face = "plain", colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24, face = "plain")
  )
save_pdf(pD2, "Panel_D2_Ddomain_dependence_distribution.pdf", 12, 6.5)

## D3: replace the branch heatmap with branch-specific effect estimates and confidence intervals.
d3 <- as.data.table(branch_summary)
d3[, branch := factor(as.character(branch), levels = names(branch_cols))]
d3[, module := factor(module, levels = rev(module_levels))]
d3[, supported := as.logical(exploratory_branch_support)]

pD3 <- ggplot(d3, aes(y = module, colour = branch)) +
  geom_vline(xintercept = 0, colour = "#969696", linewidth = 0.5) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = module),
               linewidth = 1.05, alpha = 0.75) +
  geom_point(aes(x = mean_centered_activity, shape = supported),
             size = 3.6, stroke = 0.8) +
  facet_wrap(~ branch, nrow = 1) +
  scale_colour_manual(values = branch_cols, guide = "none") +
  scale_shape_manual(
    values = c(`FALSE` = 1, `TRUE` = 16),
    breaks = c(FALSE, TRUE),
    labels = c("Exploratory only", "Supported branch-module association"),
    name = NULL
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Centered module activity", y = "Metabolic module"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(size = 22, face = "plain"),
    axis.title = element_text(size = 22, face = "plain"),
    axis.text = element_text(size = 22, face = "plain", colour = "black"),
    axis.text.y = element_text(size = 20, face = "plain", colour = "black"),
    strip.text = element_text(size = 22, face = "plain", colour = "black"),
    legend.title = element_text(size = 22, face = "plain"),
    legend.text = element_text(size = 22, face = "plain"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    panel.grid.major.y = element_line(colour = "#EFEFEF", linewidth = 0.35)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE))
save_pdf(pD3, "Panel_D3_branch_module_effect_forest.pdf", 12, 7)

## E: MSI-only m/z distributions for the most reproducible modules that
## contain at least 20 representative features. Features are ranked by
## intramodular TOM connectivity; LC-MS/MS data and annotations are not used.
module_object <- readRDS(file.path(module_root, "step3_modules/MSI_consensus_modules.rds"))
representatives <- as.data.table(module_object$representatives)

eligible_modules <- a_dat[
  n_representative_features >= 20,
  .(module, median_LOSO_Jaccard, median_within_consensus_rho)
][order(-median_LOSO_Jaccard, -median_within_consensus_rho)]
selected_modules <- as.character(head(eligible_modules$module, 4))

top_features <- rbindlist(lapply(selected_modules, function(selected_module) {
  module_rows <- which(as.character(representatives$module) == selected_module)
  module_tom <- module_object$TOM[module_rows, module_rows, drop = FALSE]
  diag(module_tom) <- NA_real_
  feature_connectivity <- rowMeans(module_tom, na.rm = TRUE)
  selected <- order(feature_connectivity, decreasing = TRUE)[seq_len(20)]
  data.table(
    module = selected_module,
    mz = representatives$mz[module_rows[selected]],
    connectivity = feature_connectivity[selected],
    rank = seq_along(selected)
  )
}))
top_features[, normalized_connectivity := {
  value_range <- range(connectivity, na.rm = TRUE)
  if (diff(value_range) == 0) rep(1, .N) else
    (connectivity - value_range[1]) / diff(value_range)
}, by = module]
top_features[, module := factor(module, levels = selected_modules)]
top_features[, mz_label := sprintf("%.4f", mz)]
module_positions <- setNames(rev(seq(1, 7, by = 2)), selected_modules)
top_features[, module_y := unname(module_positions[as.character(module)])]

pE <- ggplot(top_features, aes(mz, module_y, colour = module)) +
  geom_hline(yintercept = unname(module_positions), colour = "#E6E6E6", linewidth = 0.35) +
  geom_point(size = 3.0, alpha = 0.92) +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_x_continuous(
    breaks = seq(100, 800, by = 100),
    minor_breaks = seq(50, 800, by = 50),
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  scale_y_continuous(
    limits = c(0.05, 7.95),
    breaks = sort(unname(module_positions)),
    labels = rev(selected_modules),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "E | Representative m/z distributions of robust MSI modules",
    subtitle = "Top 20 MSI features per module ranked by within-module TOM connectivity; colours denote modules",
    x = "MSI m/z", y = "Metabolic module"
  ) +
  theme_pub +
  theme(panel.grid.major = element_line(colour = "#EFEFEF", linewidth = 0.35),
        panel.grid.major.y = element_blank())

if (requireNamespace("ggrepel", quietly = TRUE)) {
  for (current_module in selected_modules) {
    centre <- unname(module_positions[current_module])
    pE <- pE + ggrepel::geom_text_repel(
      data = top_features[module == current_module],
      aes(mz, module_y, label = mz_label), inherit.aes = FALSE,
      size = 2.35, colour = "black", seed = 225,
      direction = "y", ylim = c(centre - 0.85, centre + 0.85),
      max.overlaps = Inf, min.segment.length = 0,
      box.padding = 0.16, point.padding = 0.10, segment.alpha = 0.38
    )
  }
} else {
  pE <- pE + geom_text(aes(label = mz_label), size = 2.35,
                       angle = 45, hjust = -0.05, vjust = -0.25,
                       check_overlap = TRUE, colour = "black")
}

save_pdf(pE, "Panel_E_MSI_top_feature_mz_distribution.pdf", 14.5, 8.5)
fwrite(top_features, file.path(out_dir, "Panel_E_MSI_top_feature_selection.csv"))

pdf_files <- file.path(out_dir, c(
  "Panel_A_module_cohesion_and_robustness.pdf",
  "Panel_B_representative_spatial_module_maps.pdf",
  "Panel_C_module_MoransI_distribution.pdf",
  "Panel_D1_Ddomain_module_profile_heatmap.pdf",
  "Panel_D2_Ddomain_dependence_distribution.pdf",
  "Panel_D3_branch_module_effect_forest.pdf",
  "Panel_E_MSI_top_feature_mz_distribution.pdf"
))
stopifnot(all(file.exists(pdf_files)), length(pdf_files) == 7L)

zip_file <- file.path(out_dir, "Retained_panels_A_B_C_D1_D2_D3_E.zip")
utils::zip(zip_file, pdf_files, flags = "-j")
writeLines(
  base64enc::base64encode(zip_file, linewidth = 0),
  file.path(out_dir, "Retained_panels_bundle_singleline.b64"),
  useBytes = TRUE
)

cat("Created 7 PDFs in:", out_dir, "\n")
cat(paste(basename(pdf_files), collapse = "\n"), "\n")
