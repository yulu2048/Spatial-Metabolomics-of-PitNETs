suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

module_root <- "/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225/spatial_metabolite_modules_10samples_tuned225"
out_dir <- file.path(module_root, "publication_modules_reduced_heatmaps")
module_object <- readRDS(file.path(module_root, "step3_modules/MSI_consensus_modules.rds"))
representatives <- as.data.table(module_object$representatives)
module_summary <- fread(file.path(module_root, "step3_modules/module_summary.csv"))
robustness <- fread(file.path(
  module_root,
  "step7_robustness_spatial_scale/leave_one_sample_out_module_summary.csv"
))
module_metrics <- merge(module_summary, robustness, by = "module", all.x = TRUE)

domain_summary <- fread(file.path(
  module_root,
  "step5B_sample_specific_domain_relation/module_sample_specific_domain_summary.csv"
))
branch_summary <- fread(file.path(
  module_root,
  "domain_validation_sample_specific_4branches/exploratory_branch_module_effects_physical_sample_unit.csv"
))
branch_evidence <- branch_summary[, .(
  n_D3_supported = sum(exploratory_branch_support, na.rm = TRUE),
  supported_branches = paste(branch[exploratory_branch_support], collapse = ",")
), by = module]

selection_table <- Reduce(
  function(x, y) merge(x, y, by = "module", all.x = TRUE),
  list(module_metrics, domain_summary, branch_evidence)
)
selection_table[, D2_supported := module %in% c(
  "M01", "M03", "M04", "M06", "M07", "M08", "M10", "M11"
)]

# Select modules hierarchically: sufficient feature count, reproducible D-domain
# dependence (D2), number of supported branch associations (D3), domain effect,
# and finally leave-one-sample-out reconstruction stability.
selected_modules <- selection_table[
  n_representative_features >= 10 & D2_supported
][order(
  -n_D3_supported,
  -median_sample_specific_eta2_validation,
  -median_LOSO_Jaccard,
  -median_within_consensus_rho
), head(as.character(module), 6)]

top_features <- rbindlist(lapply(selected_modules, function(selected_module) {
  rows <- which(as.character(representatives$module) == selected_module)
  module_tom <- module_object$TOM[rows, rows, drop = FALSE]
  diag(module_tom) <- NA_real_
  connectivity <- rowMeans(module_tom, na.rm = TRUE)
  keep <- order(connectivity, decreasing = TRUE)[seq_len(10)]
  data.table(
    module = selected_module,
    mz = representatives$mz[rows[keep]],
    connectivity = connectivity[keep],
    rank = seq_along(keep)
  )
}))
top_features[, normalized_connectivity := {
  limits <- range(connectivity, na.rm = TRUE)
  if (diff(limits) == 0) rep(1, .N) else (connectivity - limits[1]) / diff(limits)
}, by = module]
top_features[, module := factor(module, levels = selected_modules)]
module_positions <- setNames(
  rev(seq(1, by = 2, length.out = length(selected_modules))),
  selected_modules
)
top_features[, module_y := unname(module_positions[as.character(module)])]

module_cols <- setNames(grDevices::hcl.colors(15, "Dynamic"), sprintf("M%02d", 1:15))
theme_pub <- theme_classic(base_size = 24) +
  theme(
    text = element_text(size = 24, face = "plain"),
    axis.title = element_text(size = 24, face = "plain"),
    axis.text = element_text(size = 24, face = "plain", colour = "black"),
    strip.background = element_rect(fill = "#F3F3F3", colour = "#BDBDBD"),
    strip.text = element_text(size = 24, face = "plain", colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 24, face = "plain"),
    legend.text = element_text(size = 24, face = "plain"),
    panel.grid.major = element_line(colour = "#EFEFEF", linewidth = 0.35)
  )

pE <- ggplot(top_features, aes(mz, module_y, colour = module,
                               size = normalized_connectivity)) +
  geom_hline(yintercept = unname(module_positions), colour = "#E6E6E6", linewidth = 0.35) +
  geom_point(alpha = 0.88) +
  scale_colour_manual(values = module_cols, guide = "none") +
  scale_size_continuous(
    name = "Within-module connectivity",
    range = c(5, 16), limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = c("Low", "Medium", "High")
  ) +
  scale_x_continuous(
    breaks = seq(200, 800, by = 100),
    minor_breaks = seq(200, 800, by = 50),
    limits = c(200, 820),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.05, 2 * length(selected_modules) - 0.05),
    breaks = sort(unname(module_positions)),
    labels = rev(selected_modules),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "MSI m/z", y = "Metabolic module"
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank())

panel_e_pdf <- file.path(out_dir, "Panel_E_MSI_top_feature_mz_distribution.pdf")
ggsave(panel_e_pdf, pE, width = 22, height = 6, units = "in",
       device = grDevices::cairo_pdf, bg = "white")
fwrite(top_features, file.path(out_dir, "Panel_E_MSI_top_feature_selection.csv"))

cat("Selected modules:", paste(selected_modules, collapse = ", "), "\n")
cat("Panel E:", panel_e_pdf, "\n")
