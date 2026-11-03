#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(20260914)

root <- "/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225"
module_root <- file.path(root, "spatial_metabolite_modules_10samples_tuned225")
s5 <- file.path(module_root, "step5B_sample_specific_domain_relation")
s6 <- file.path(module_root, "step6_posthoc_LCMS_annotation")
s7 <- file.path(module_root, "step7_robustness_spatial_scale")
relation_root <- file.path(root, "domain_relationships_21domains_1D_tuned225")
out <- file.path(module_root, "domain_validation_sample_specific_4branches")
fig_out <- file.path(module_root, "publication_figure_metabolic_modules_sample_domains_4branches")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

profiles <- fread(file.path(s5, "sample_domain_module_profiles_validation_only.csv"))
relations <- fread(file.path(s5, "sample_module_specific_domain_relations.csv"))
domain_summary <- fread(file.path(s5, "module_sample_specific_domain_summary.csv"))
cuts <- fread(file.path(relation_root, "exploratory_dendrogram_cuts_not_common_domain_labels.csv"))
rel_obj <- readRDS(file.path(relation_root, "domain_relationships.rds"))

stopifnot(nrow(cuts) == 21L, uniqueN(cuts$domain_id) == 21L)
cuts[, branch := paste0("B", exploratory_cut_k4)]
profiles[, domain_id := paste(sample, domain, sep = " / ")]
profiles <- merge(profiles, cuts[, .(domain_id, sample_id, diagnosis,
                                    spatial_organization_label, branch)],
                  by = "domain_id", all.x = TRUE)
stopifnot(!anyNA(profiles$branch), all(profiles$sample == profiles$sample_id))

# Stable-core, validation-only scores are the inferential values. Center within
# each physical sample so branch comparisons cannot be driven by sample baselines.
profiles[, stable_weight := fifelse(is.finite(mean_z_stable) & n_stable > 0, n_stable, 0)]
profiles[, sample_stable_mean := {
  ok <- stable_weight > 0
  if (any(ok)) weighted.mean(mean_z_stable[ok], stable_weight[ok]) else NA_real_
}, by = .(sample, module)]
profiles[, centered_stable_mean := mean_z_stable - sample_stable_mean]
profiles[, branch := factor(branch, levels = paste0("B", 1:4))]

# Order columns by the fixed domain dendrogram; do not refit domains or branches.
hc <- rel_obj$hclust
domain_order <- if (!is.null(hc$labels)) hc$labels[hc$order] else cuts$domain_id
if (!setequal(domain_order, cuts$domain_id)) domain_order <- cuts$domain_id
profiles[, domain_id := factor(domain_id, levels = domain_order)]
fwrite(profiles, file.path(out, "module_activity_by_21_sample_specific_domains.csv"))

# Physical-sample unit: first aggregate any same-branch domains within a sample.
branch_sample <- profiles[is.finite(centered_stable_mean), .(
  centered_activity = weighted.mean(centered_stable_mean, pmax(n_stable, 1)),
  n_domains = .N,
  n_stable_pixels = sum(n_stable)
), by = .(sample, module, branch)]
fwrite(branch_sample, file.path(out, "module_activity_by_sample_and_exploratory_branch.csv"))

boot_median_ci <- function(x, B = 2000L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(NA_real_, NA_real_))
  b <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(b, c(.025, .975), na.rm = TRUE))
}

branch_summary <- branch_sample[, {
  ci <- boot_median_ci(centered_activity)
  s <- sign(median(centered_activity))
  .(n_samples = .N,
    median_centered_activity = median(centered_activity),
    mean_centered_activity = mean(centered_activity),
    ci_low = ci[1], ci_high = ci[2],
    direction_fraction = if (s == 0) NA_real_ else mean(sign(centered_activity) == s),
    loso_direction_consistent = if (.N < 3L || s == 0) NA else
      all(vapply(seq_len(.N), function(i) sign(median(centered_activity[-i])) == s,
                 logical(1)))
  )
}, by = .(module, branch)]
branch_summary[, exploratory_branch_support :=
                 n_samples >= 3 & abs(median_centered_activity) >= .35 &
                 direction_fraction >= .8 & is.finite(ci_low) &
                 sign(ci_low) == sign(ci_high) & loso_direction_consistent]
fwrite(branch_summary, file.path(out, "exploratory_branch_module_effects_physical_sample_unit.csv"))

paired_one <- function(b1, b2, label) {
  wide <- dcast(branch_sample[branch %in% c(b1, b2)], sample + module ~ branch,
                value.var = "centered_activity")
  wide <- wide[is.finite(get(b1)) & is.finite(get(b2))]
  wide[, difference := get(b1) - get(b2)]
  fwrite(wide, file.path(out, paste0("paired_sample_differences_", label, ".csv")))
  wide[, {
    x <- difference
    ci <- boot_median_ci(x)
    s <- sign(median(x))
    .(branch_1 = b1, branch_2 = b2, n_paired_samples = .N,
      median_difference = median(x), mean_difference = mean(x),
      ci_low = ci[1], ci_high = ci[2],
      direction_fraction = if (s == 0) NA_real_ else mean(sign(x) == s),
      loso_direction_consistent = if (.N < 3L || s == 0) NA else
        all(vapply(seq_len(.N), function(i) sign(median(x[-i])) == s, logical(1))))
  }, by = module]
}
paired <- rbindlist(list(
  paired_one("B1", "B2", "B1_vs_B2"),
  paired_one("B3", "B4", "B3_vs_B4")
))
paired[, exploratory_pair_support :=
         n_paired_samples >= 3 & abs(median_difference) >= .35 &
         direction_fraction >= .8 & is.finite(ci_low) &
         sign(ci_low) == sign(ci_high) & loso_direction_consistent]
fwrite(paired, file.path(out, "paired_branch_module_comparisons_physical_sample_unit.csv"))

# Where does the locally highest domain fall in the four-branch tree?
high_domain <- copy(relations[, .(sample, module, eta2_validation_stable,
                                  domain_relation_class, highest_activity_domain)])
high_domain[, domain_id := paste(sample, highest_activity_domain, sep = " / ")]
high_domain <- merge(high_domain, cuts[, .(domain_id, high_activity_branch = branch)],
                     by = "domain_id", all.x = TRUE)
high_consistency <- high_domain[is.finite(eta2_validation_stable) & eta2_validation_stable >= .25,
  {tb <- sort(table(high_activity_branch), decreasing = TRUE)
   .(n_strong_samples = .N,
     dominant_high_activity_branch = if (length(tb)) names(tb)[1] else NA_character_,
     dominant_branch_fraction = if (length(tb)) as.numeric(tb[1])/.N else NA_real_,
     n_high_activity_branches = uniqueN(high_activity_branch))}, by = module]
fwrite(high_domain, file.path(out, "highest_activity_domain_with_4branch_annotation.csv"))

pair_best <- paired[order(module, -exploratory_pair_support, -abs(median_difference)), .SD[1], by = module]
setnames(pair_best,
         c("branch_1", "branch_2", "n_paired_samples", "median_difference", "ci_low", "ci_high",
           "direction_fraction", "loso_direction_consistent", "exploratory_pair_support"),
         c("best_pair_branch_1", "best_pair_branch_2", "best_pair_n_samples",
           "best_pair_median_difference", "best_pair_ci_low", "best_pair_ci_high",
           "best_pair_direction_fraction", "best_pair_loso_consistent",
           "four_branch_exploratory_support"))

master <- merge(domain_summary, high_consistency, by = "module", all.x = TRUE)
master <- merge(master, pair_best[, .(module, best_pair_branch_1, best_pair_branch_2,
                                     best_pair_n_samples, best_pair_median_difference,
                                     best_pair_ci_low, best_pair_ci_high,
                                     best_pair_direction_fraction,
                                     best_pair_loso_consistent,
                                     four_branch_exploratory_support)],
                by = "module", all.x = TRUE)
master[is.na(four_branch_exploratory_support), four_branch_exploratory_support := FALSE]
master[, domain_program_class := fcase(
  median_sample_specific_eta2_validation < .10 | n_samples_trans_domain >= 5,
    "trans-domain / continuous",
  four_branch_exploratory_support,
    "four-branch structured (exploratory)",
  n_samples_strong >= 5 & dominant_branch_fraction < .70,
    "locally reconfigured D-domain program",
  n_samples_strong >= 5,
    "recurrent D-domain-aligned program",
  default = "partial / heterogeneous D-domain relation"
)]

lcms <- fread(file.path(s6, "module_LCMS_annotation_summary.csv"))
rob <- fread(file.path(s7, "leave_one_sample_out_module_summary.csv"))
scale <- fread(file.path(s7, "module_spatial_scale_summary.csv"))
master <- Reduce(function(x, y) merge(x, y, by = "module", all.x = TRUE),
                 list(master, lcms, rob, scale))
setorder(master, -median_sample_specific_eta2_validation)
fwrite(master, file.path(out, "MASTER_module_evidence_table_sample_domains_4branches.csv"))

# Figure E1: module activity across the 21 sample-specific D-domains.
ann <- unique(profiles[, .(domain_id, sample, branch)])
ann[, x := as.numeric(domain_id)]
ann[, y := 16.15]
branch_cols <- c(B1 = "#1B9E77", B2 = "#D95F02", B3 = "#7570B3", B4 = "#E7298A")
pE1 <- ggplot(profiles, aes(domain_id, factor(module, levels = rev(sprintf("M%02d", 1:15))),
                            fill = centered_stable_mean)) +
  geom_tile(color = "white", linewidth = .15) +
  geom_tile(data = ann, aes(domain_id, y, fill = branch), inherit.aes = FALSE,
            height = .45, color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       na.value = "grey90", name = "Centered\nactivity") +
  labs(title = "E1  Validation-only module activity in 21 sample-specific D-domains",
       subtitle = "Columns follow the fixed domain dendrogram; top strip denotes exploratory B1–B4",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
        panel.grid = element_blank(), plot.title = element_text(face = "bold"))
# Overlay branch strip with its own colours (new scale if ggnewscale is available).
if (requireNamespace("ggnewscale", quietly = TRUE)) {
  pE1 <- ggplot(profiles, aes(domain_id, factor(module, levels = rev(sprintf("M%02d", 1:15))),
                              fill = centered_stable_mean)) +
    geom_tile(color = "white", linewidth = .15) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         na.value = "grey90", name = "Centered\nactivity") +
    ggnewscale::new_scale_fill() +
    geom_tile(data = ann, aes(domain_id, y, fill = branch), inherit.aes = FALSE,
              height = .45, color = "white") +
    scale_fill_manual(values = branch_cols, name = "Exploratory\nbranch") +
    labs(title = "E1  Validation-only module activity in 21 sample-specific D-domains",
         subtitle = "Columns follow the fixed domain dendrogram; top strip denotes exploratory B1–B4",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
          panel.grid = element_blank(), plot.title = element_text(face = "bold"))
}

# Figure E2: per-sample D-domain association (eta squared).
eta <- relations[, .(sample, module, eta2 = eta2_validation_stable)]
pE2 <- ggplot(eta, aes(sample, factor(module, levels = rev(sprintf("M%02d", 1:15))), fill = eta2)) +
  geom_tile(color = "white", linewidth = .2) +
  scale_fill_gradient(low = "white", high = "#542788", limits = c(0, 1), oob = scales::squish) +
  labs(title = "E2  Within-sample D-domain dependence", x = NULL, y = NULL,
       fill = expression(eta^2)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        panel.grid = element_blank(), plot.title = element_text(face = "bold"))

# Figure E3: branch-centred effects, keeping sample as the replication unit.
pE3 <- ggplot(branch_summary,
              aes(branch, factor(module, levels = rev(sprintf("M%02d", 1:15))),
                  fill = median_centered_activity)) +
  geom_tile(color = "white") +
  geom_point(data = branch_summary[exploratory_branch_support == TRUE], shape = 8, size = 1.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "E3  Four-branch organization (exploratory)",
       subtitle = "Median within-sample-centred effect; * passes bootstrap/LOSO rule",
       x = NULL, y = NULL, fill = "Median effect") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

# Figure E4: only adequately paired branch contrasts.
paired_plot <- copy(paired)
paired_plot[, contrast := paste0(branch_1, "−", branch_2)]
pE4 <- ggplot(paired_plot,
              aes(median_difference, factor(module, levels = rev(sprintf("M%02d", 1:15))),
                  xmin = ci_low, xmax = ci_high, color = contrast)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_errorbarh(height = .22, position = position_dodge(width = .45)) +
  geom_point(aes(shape = exploratory_pair_support), position = position_dodge(width = .45), size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
  labs(title = "E4  Paired branch contrasts", subtitle = "Only B1–B2 and B3–B4 have 5 shared samples",
       x = "Median paired difference (bootstrap 95% CI)", y = NULL, color = NULL, shape = "Support") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

save_plot <- function(p, stem, w, h) {
  ggsave(file.path(fig_out, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_out, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
}
save_plot(pE1, "Panel_E1_21_sample_specific_domains", 12, 6.5)
save_plot(pE2, "Panel_E2_sample_domain_eta2", 8, 6)
save_plot(pE3, "Panel_E3_exploratory_4branch_effects", 5.5, 6)
save_plot(pE4, "Panel_E4_paired_branch_contrasts", 7, 6)

if (requireNamespace("patchwork", quietly = TRUE)) {
  pE <- (pE1 / (pE2 | pE3 | pE4)) + patchwork::plot_layout(heights = c(1.2, 1))
  save_plot(pE, "Panel_E_domain_validation_sample_specific_4branches", 18, 12)
}

# Integrated evidence: no common-domain fields.
evidence <- rbindlist(list(
  master[, .(module, evidence = "D-domain dependence",
             value = pmin(median_sample_specific_eta2_validation/.50, 1),
             label = sprintf("eta2=%.2f", median_sample_specific_eta2_validation))],
  master[, .(module, evidence = "Four-branch structure\n(exploratory)",
             value = pmin(abs(best_pair_median_difference)/.75, 1),
             label = fifelse(four_branch_exploratory_support, "supported", "not robust"))],
  master[, .(module, evidence = "Module stability",
             value = pmin(median_LOSO_Jaccard, 1), label = sprintf("J=%.2f", median_LOSO_Jaccard))],
  master[, .(module, evidence = "LC-MS/MS support\n(post hoc)",
             value = fraction_families_annotated,
             label = sprintf("%.0f%%", 100*fraction_families_annotated))]
), use.names = TRUE)
pG <- ggplot(evidence, aes(evidence, factor(module, levels = rev(sprintf("M%02d", 1:15))),
                           fill = value)) +
  geom_tile(color = "white") + geom_text(aes(label = label), size = 2.25) +
  scale_fill_gradient(low = "white", high = "#006D2C", limits = c(0, 1), na.value = "grey90") +
  labs(title = "G  Integrated evidence for frozen MSI metabolic programs",
       subtitle = "Domain evidence uses sample-specific D-domains; B1–B4 remain exploratory annotations",
       x = NULL, y = NULL, fill = "Scaled evidence") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold"))
save_plot(pG, "Panel_G_integrated_evidence_no_common_domain", 8, 6.5)

# Reassemble the full publication figure from valid frozen-module panels plus new E/G.
old_fig <- file.path(module_root, "publication_figure_metabolic_modules")
if (requireNamespace("magick", quietly = TRUE)) try({
  read_img <- function(x) magick::image_read(file.path(old_fig, x))
  b <- read_img("Panel_B_module_structure_stability.png")
  c <- read_img("Panel_C_representative_spatial_programs.png")
  d <- read_img("Panel_D_cross_sample_spatial_organization.png")
  f <- read_img("Panel_F_LCMS_support.png")
  e <- magick::image_read(file.path(fig_out, "Panel_E_domain_validation_sample_specific_4branches.png"))
  g <- magick::image_read(file.path(fig_out, "Panel_G_integrated_evidence_no_common_domain.png"))
  target_w <- 3600
  norm <- function(im) magick::image_resize(im, paste0(target_w, "x"))
  top <- magick::image_append(c(norm(b), norm(c)), stack = FALSE)
  mid <- magick::image_append(c(norm(d), norm(f)), stack = FALSE)
  full <- magick::image_append(c(top, norm(e), mid, norm(g)), stack = TRUE)
  magick::image_write(full, file.path(fig_out, "Main_Figure_modules_sample_domains_4branches_LCMS.png"))
  magick::image_write(full, file.path(fig_out, "Main_Figure_modules_sample_domains_4branches_LCMS.pdf"), format = "pdf")
}, silent = TRUE)

writeLines(c(
  "REANALYSIS COMPLETE: sample-specific D-domains plus exploratory four-branch annotation",
  "",
  "Frozen inputs:",
  "- 15 MSI modules and their pixel scores were not changed.",
  "- LC-MS/MS remained post-hoc annotation/support and did not alter MSI modules.",
  "",
  "Domain framework:",
  "- 21 within-sample D-domains across 10 physical samples.",
  "- Fixed B1-B4 cut from the domain dendrogram; no common-domain model was used.",
  "- Stable-core validation-only module scores were used for domain inference.",
  "- Physical sample, not domain, was the replicate for branch statistics.",
  "- Only B1-B2 and B3-B4 were tested as paired contrasts (5 shared samples each).",
  "",
  "Thresholds:",
  "- strong D-domain relation: eta2 >= 0.25; trans-domain: eta2 < 0.10.",
  "- exploratory branch support: n>=3, |median effect|>=0.35, direction fraction>=0.80,",
  "  sample-bootstrap 95% CI excludes zero, and leave-one-sample-out direction is unchanged.",
  "",
  "Interpretation limit: B1-B4 are exploratory higher-level annotations of sample-specific domains,",
  "not recurrent common-domain labels and not four independently established biological classes."
), file.path(out, "README.txt"))

summary_lines <- c(
  sprintf("n modules: %d", nrow(master)),
  sprintf("n domains: %d", uniqueN(profiles$domain_id)),
  sprintf("n physical samples: %d", uniqueN(profiles$sample)),
  sprintf("paired B1-B2 samples: %d", max(paired[branch_1 == "B1", n_paired_samples])),
  sprintf("paired B3-B4 samples: %d", max(paired[branch_1 == "B3", n_paired_samples])),
  "module classes:",
  paste(capture.output(print(master[, .N, by = domain_program_class][order(-N)])), collapse = "\n"),
  "robust four-branch-structured modules (after prioritizing D-domain eta2):",
  paste(master[domain_program_class == "four-branch structured (exploratory)", module], collapse = ", "),
  "all modules with at least one supported paired contrast (secondary evidence only):",
  paste(master[four_branch_exploratory_support == TRUE, module], collapse = ", ")
)
writeLines(summary_lines, file.path(out, "analysis_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("OUT=", out, "\nFIG=", fig_out, "\n")
