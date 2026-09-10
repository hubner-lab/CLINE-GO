#!/usr/bin/env Rscript
# =============================================================================
# mvp_manuscript_figures.R -- manuscript figures and tables for the MVP
# simulated benchmark (Lotterhos 2023 corpus; offset arm vs Lind & Lotterhos 2025).
#
# Reads ONLY tables already on disk under benchmarks/mvp_eval/. Fits nothing,
# runs no pipeline mode. Every panel writes the exact values it plots next to
# the image, so no number in the manuscript is untraceable.
#
# DELIBERATELY REDUNDANT. Each claim is rendered under several metrics and
# several framings so one can be chosen on the strength of the full analysis
# rather than on what happened to be plotted first. Outputs are grouped:
#
#   A  benchmark design
#   B  detection, per method          (AUC-PR / best-F1 / precision-recall / MAF arm / redundancy)
#   C  structure-correction ladder    (relative / absolute / best-rung / gain)
#   D  combining methods              (gain / window / transfer mode / selection bias)
#   E  window-based scanning (WZA)
#   F  genetic offset                 (panel / method / stability / decoupling / predictors)
#   G  comparison with published baselines
#
# figure_manifest.tsv maps every output to its source tables and its claim.
#
# Generation discipline (see docs/03-results-simulations.md):
#   journal 07 (report07/, 14 replicates)   -> primary detection results
#   journal 06 (report06/, 1 replicate)     -> single-replicate mechanism only
#   journal 05 (report/mvp_ladders.tsv)     -> the ladder AUC-PR grid; its combining
#       conclusion is superseded, its rung numbers are not
#   phase 1 (figures/, offset09/phase1_*)   -> primary offset results, INTERIM
#
# Usage:  Rscript /pipeline/benchmarks/mvp_manuscript_figures.R
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
# Which scored offset root to read, and where the figures land. Both default to the values
# that produced the 32-replicate figure set, so an unparameterised run reproduces it; a new
# cohort passes OFFSET_DIR=offset10 FIG_OUT=... and cannot overwrite the frozen outputs.
OFF  <- Sys.getenv("OFFSET_DIR", "offset09")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ms"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

# Constants from the seed-selection screen. Sourced from
# docs/gea_simulation_benchmarks.md sections 3.1 and 8.2 -- not re-derived here.
N_CORPUS  <- 2250L               # replicates in the deposit
N_QUALIFY <- 55L                 # meeting all three selection criteria, all SS-Mtn
R2_BAND   <- c(0.30, 0.60)       # target R^2(PC1 ~ temperature) band

# Publication labels. The corpus ships "oliogenic" (sic) as a level name.
ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
relabel_arch <- function(x) factor(x, levels = ARCH_LEVELS, labels = ARCH_LABELS)

MANIFEST <- list()
rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) { message("MISSING: ", f); return(NULL) }
    fread(f)
}
emit <- function(dt, stem) {                      # plotted values, beside the image
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    invisible(dt)
}
reg <- function(id, kind, shows, sources, n_note) # figure_manifest row
    MANIFEST[[length(MANIFEST) + 1L]] <<- data.table(
        id = id, kind = kind, shows = shows, sources = sources, coverage = n_note)
save_fig <- function(p, stem, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h)
    # Guarantee a SAME-STEM table beside every image, so figure and values pair
    # mechanically rather than by knowing which emit() call fed which plot.
    # ggplot objects carry their data; cowplot composites do not, and those
    # already have per-panel tables emitted under their own stems.
    tsv <- file.path(OUT, paste0(stem, ".tsv"))
    if (!file.exists(tsv) && !is.null(p$data) && is.data.frame(p$data) && nrow(p$data))
        fwrite(as.data.table(p$data), tsv, sep = "\t")
    message(sprintf("  OK %s", stem))
}
sp <- function(x, y) suppressWarnings(stats::cor(x, y, method = "spearman",
                                                 use = "complete.obs"))

# =============================================================================
# A -- benchmark design
# =============================================================================
message("\n=== A: benchmark design ===")
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
stopifnot(!is.null(seeds))
seeds[, is_control := arm == "control_degenerate"]
seeds[, arch_lab := relabel_arch(arch_level)]
stopifnot(!any(is.na(seeds$arch_lab)))
PRIM <- seeds[is_control == FALSE]

funnel <- data.table(
    stage = factor(c("deposited", "meet criteria", "converted"),
                   levels = c("deposited", "meet criteria", "converted")),
    n     = c(N_CORPUS, N_QUALIFY, nrow(seeds)),
    note  = c("Lotterhos 2023 corpus",
              "Fst >= 0.05, R2(PC1~temp) in [0.30,0.60], R2(PC1~Env2) <= 0.20",
              sprintf("%d primary + %d degenerate controls",
                      nrow(PRIM), sum(seeds$is_control))))
emit(funnel, "A1a_selection_funnel")

a1 <- ggplot(funnel, aes(stage, n)) +
    geom_col(fill = CLINEGO_RETAINED, width = 0.62) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.4, colour = CLINEGO_COL$fg) +
    scale_y_log10(expand = expansion(mult = c(0, 0.18))) +
    labs(title = "A  Replicate selection", x = NULL, y = "replicates (log scale)") +
    theme_clinego()

a1b_src <- seeds[, .(seed, arch = arch_lab, is_control, meanFst, r2_pc1_temp,
                     r2_pc1_sal, n_snps, n_causal_maf01, final_LA, k_best,
                     n_traits, pleiotropy = ispleiotropy, demog = demog_level_sub)]
emit(a1b_src, "A1b_seed_properties")

a2 <- ggplot(a1b_src, aes(meanFst, r2_pc1_temp, colour = arch, shape = is_control)) +
    geom_hline(yintercept = R2_BAND, linetype = "dashed",
               colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
    geom_point(size = 2.6, stroke = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c("primary", "degenerate control"), name = NULL) +
    scale_color_clinego(name = NULL) +
    labs(title = "B  Structure and its alignment with the environment",
         x = expression(mean~italic(F)[ST]),
         y = expression(italic(R)^2*"(PC1 ~ temperature)")) +
    theme_clinego() +
    theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))

a3 <- ggplot(PRIM, aes(arch_lab, n_causal_maf01, colour = arch_lab)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, colour = CLINEGO_COL$secondary) +
    geom_jitter(width = 0.14, height = 0, size = 2, alpha = 0.85) +
    scale_y_log10() + scale_color_clinego(guide = "none") +
    labs(title = "C  Testable causal loci", x = NULL,
         y = "causal loci at MAF > 0.01 (log scale)") +
    theme_clinego() + theme(axis.text.x = element_text(angle = 12, hjust = 1))

save_fig(plot_grid(a1, a2, a3, nrow = 1, rel_widths = c(0.9, 1.25, 0.95)),
         "A1_benchmark_design", 14, 4.8)
reg("A1", "figure", "Selection funnel; structure vs environment; causal-locus range",
    "mvp_seeds.tsv", sprintf("%d replicates", nrow(seeds)))

# A2 -- the covariates that the offset arm later shows matter. Separate figure so
# the design section can be shown with or without it.
a2src <- melt(PRIM[, .(seed, arch = arch_lab, n_snps, final_LA, meanFst, r2_pc1_temp)],
              id.vars = c("seed", "arch"), variable.name = "property", value.name = "value")
a2src[, property := factor(property, levels = c("n_snps", "final_LA", "meanFst", "r2_pc1_temp"),
        labels = c("SNPs in replicate", "local adaptation (final_LA)",
                   "mean Fst", "R2(PC1 ~ temp)"))]
emit(a2src, "A2_replicate_covariates")
save_fig(ggplot(a2src, aes(arch, value, colour = arch)) +
             geom_boxplot(outlier.shape = NA, width = 0.55, colour = CLINEGO_COL$secondary) +
             geom_jitter(width = 0.14, height = 0, size = 1.8, alpha = 0.85) +
             facet_wrap(~ property, scales = "free_y", nrow = 1) +
             scale_color_clinego(guide = "none") +
             labs(x = NULL, y = NULL) + theme_clinego_grid() +
             theme(axis.text.x = element_text(angle = 20, hjust = 1)),
         "A2_replicate_covariates", 13, 4.2)
reg("A2", "figure", "Replicate covariates by architecture (size, local adaptation, Fst, confounding)",
    "mvp_seeds.tsv", sprintf("%d primary replicates", nrow(PRIM)))

# =============================================================================
# B -- detection, per method
#
# Provenance that must reach the caption: method_rank.tsv / aucpr_by_arm.tsv are
# built (mvp_crossseed_report.R:88, 256-270) from redundancy07*/MVP*_c3_method_
# clusters.tsv -- cell c3, the pipeline's DEFAULT structure-correction rung,
# scored on truth_any. They are not a best-rung figure; group C is.
# =============================================================================
message("\n=== B: detection by method ===")
auc_arm <- rd(EVAL, "report07/aucpr_by_arm.tsv")
mrank   <- rd(EVAL, "report07/method_rank.tsv")
lad     <- rd(EVAL, "report/mvp_ladders.tsv")
stopifnot(!is.null(auc_arm), !is.null(mrank))

b1_src <- melt(auc_arm, id.vars = c("seed", "method", "arch", "control"),
               measure.vars = c("aucpr_maf01", "aucpr_maf05"),
               variable.name = "arm", value.name = "auc_pr")
b1_src[, arm := factor(arm, levels = c("aucpr_maf01", "aucpr_maf05"),
                       labels = c("MAF > 0.01 (primary)", "MAF > 0.05 (sensitivity)"))]
b1_src[, arch := relabel_arch(arch)]
ord <- mrank[arm == "base"][order(-med_aucpr), method]
b1_src[, method := factor(method, levels = ord)]
emit(b1_src, "B1_aucpr_by_method")

n_lab <- mrank[arm == "base", .(method, lab = paste0("n=", n_seeds))]
n_lab[, method := factor(method, levels = ord)]

save_fig(ggplot(b1_src, aes(method, auc_pr)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, colour = CLINEGO_COL$secondary) +
    geom_jitter(aes(colour = arch, shape = control), width = 0.16, height = 0,
                size = 1.9, alpha = 0.9) +
    geom_text(data = n_lab, aes(x = method, y = -0.04, label = lab), inherit.aes = FALSE,
              size = 2.6, colour = CLINEGO_COL$muted) +
    facet_wrap(~ arm, ncol = 1) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                       labels = c("primary", "degenerate control"), name = NULL) +
    scale_color_clinego(name = NULL) +
    labs(x = NULL, y = "AUC-PR (default rung c3)") +
    theme_clinego_grid() + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
    "B1_detection_aucpr", 9.5, 7.5)
reg("B1", "figure", "Threshold-free ranking quality (AUC-PR) per method, both MAF arms",
    "report07/aucpr_by_arm.tsv; method_rank.tsv",
    "14 replicates (MLMM 11, SUPER 12); default rung c3")

# B2 -- the SAME comparison under a threshold-dependent metric. Reported because
# AUC-PR and best-F1 disagree for methods whose ranking is informative but whose
# p-value scale is not calibrated (RDA above all).
if (!is.null(lad)) {
    b2 <- lad[axis == "any" & !is.na(best_f1)]
    b2[, arch := relabel_arch(arch_level)]
    b2_src <- b2[, .(seed, method, cell, arch, auc_pr, best_f1)]
    emit(b2_src, "B2_aucpr_vs_bestf1")
    b2m <- melt(b2_src[cell == "c3"], id.vars = c("seed", "method", "arch"),
                measure.vars = c("auc_pr", "best_f1"),
                variable.name = "metric", value.name = "value")
    b2m[, metric := factor(metric, levels = c("auc_pr", "best_f1"),
                           labels = c("AUC-PR (threshold-free)", "best achievable F1"))]
    save_fig(ggplot(b2m, aes(reorder(method, value, median), value)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, colour = CLINEGO_COL$secondary) +
        geom_jitter(aes(colour = arch), width = 0.15, height = 0, size = 1.9, alpha = 0.9) +
        facet_wrap(~ metric, scales = "free_x") + coord_flip() +
        scale_color_clinego(name = NULL) +
        labs(x = NULL, y = NULL) + theme_clinego_grid(),
        "B2_aucpr_vs_bestf1", 10, 5)
    reg("B2", "figure", "Same methods under a threshold-free vs a threshold-dependent metric",
        "report/mvp_ladders.tsv", "14 replicates, 4 methods, default rung c3")
}

# B3 -- MAF arm, paired within replicate. Ratio > 1 means dropping low-frequency
# markers helped that method on that replicate.
b3_src <- auc_arm[!is.na(ratio) & is.finite(ratio)]
b3_src[, arch := relabel_arch(arch)]
emit(b3_src, "B3_maf_arm_ratio")
save_fig(ggplot(b3_src, aes(reorder(method, ratio, median), ratio)) +
    geom_hline(yintercept = 1, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
    geom_boxplot(outlier.shape = NA, width = 0.6, colour = CLINEGO_COL$secondary) +
    geom_jitter(aes(colour = arch), width = 0.15, height = 0, size = 1.9, alpha = 0.9) +
    coord_flip() + scale_y_log10() + scale_color_clinego(name = NULL) +
    labs(x = NULL, y = "AUC-PR ratio, MAF 0.05 arm / MAF 0.01 arm (log scale)") +
    theme_clinego(),
    "B3_maf_arm_effect", 8.5, 5.5)
reg("B3", "figure", "Effect of the MAF filter on ranking quality, paired within replicate",
    "report07/aucpr_by_arm.tsv", "14 replicates; truth denominator refiltered to match")

# B4 -- the low-MAF detection ceiling (single replicate; mechanism, not a rate).
mafr <- rd(EVAL, "report06/maf_recall.tsv"); mafc <- rd(EVAL, "report06/maf_cut.tsv")
if (!is.null(mafr)) {
    emit(mafr, "B4_recall_by_maf_bin")
    save_fig(ggplot(mafr, aes(maf_bin, recall, group = method, colour = method)) +
        geom_line(linewidth = 0.7) + geom_point(size = 1.9) +
        scale_color_manual(values = clinego_cluster_palette(uniqueN(mafr$method)), name = NULL) +
        labs(x = "minor-allele-frequency bin", y = "recall of causal loci") +
        theme_clinego() + theme(axis.text.x = element_text(angle = 20, hjust = 1)),
        "B4_low_maf_ceiling", 9, 5)
    reg("B4", "figure", "Recall as a function of causal-locus allele frequency",
        "report06/maf_recall.tsv", "1 replicate (1232548), 11 methods")
}
if (!is.null(mafc)) emit(mafc, "B4b_recall_split_at_maf005")

# B5 -- method redundancy. Eleven methods are not eleven independent tests.
mfam <- rd(EVAL, "report06/method_families.tsv")
if (!is.null(mfam)) {
    emit(mfam, "B5_method_clusters")
    save_fig(ggplot(mfam, aes(reorder(method, auc_pr), auc_pr,
                              fill = factor(cluster), shape = representative)) +
        geom_col(width = 0.68) + coord_flip() +
        scale_fill_manual(values = clinego_cluster_palette(uniqueN(mfam$cluster)),
                          name = "redundancy cluster") +
        labs(x = NULL, y = "AUC-PR") + theme_clinego(),
        "B5_method_redundancy", 8, 5)
    reg("B5", "figure", "Methods grouped by how similarly they rank SNPs",
        "report06/method_families.tsv", "1 replicate (1232548), 11 methods")
}

# B6 -- WHAT A SCAN ACTUALLY RETURNS, in absolute loci.
# Rates hide the thing a reader most wants to know: of the markers a method
# hands you, how many are causal, how many merely sit in a causal locus's LD
# block, and how many are pure noise. n_called == tp + expected_linked +
# fp_background exactly, so the three stack to the total with nothing missing.
# Scored at the pipeline's own default calling rule (Bonferroni 0.05) at the
# default rung c3, so this is what a user running defaults receives.
message("\n=== B6: composition of called markers ===")
CELL_DEF <- "c3"; RULE_ADJ <- "bonf"; RULE_VAL <- 0.05
comp <- list()
for (a in list(list(tag = "", arm = "base", sfx = ""),
               list(tag = "_m05", arm = "m05", sfx = "M05"))) {
    dirs <- list.dirs(file.path(EVAL, paste0("sweep07", a$tag)), recursive = FALSE)
    for (d in dirs) {
        proj <- basename(d)
        f <- file.path(d, sprintf("%s_%s_any_threshold_sweep.tsv", proj, CELL_DEF))
        if (!file.exists(f)) next
        x <- fread(f)
        x <- x[adjust == RULE_ADJ & abs(value - RULE_VAL) < 1e-9 & combine_window_kb == 0]
        if (!nrow(x)) next
        sd_ <- as.integer(sub("^MVP", "", sub("M05$", "", proj)))
        comp[[length(comp) + 1L]] <- cbind(seed = sd_, arm = a$arm, x)
    }
}
if (length(comp)) {
    CP <- rbindlist(comp, fill = TRUE)
    CP[, family := fifelse(grepl("^combine_", kind), "combination", "single method")]
    CP[, label := sub("^combine_", "", kind)]
    CP <- merge(CP, seeds[, .(seed, arch = arch_lab, is_control, n_causal_maf01)],
                by = "seed", all.x = TRUE)
    # The three components must reconstruct the total, or the table is not a
    # partition and the stacked bar would silently lie.
    bad <- CP[abs(n_called - (tp + expected_linked + fp_background)) > 0]
    if (nrow(bad)) message(sprintf("WARN: %d rows where causal+linked+background != n_called",
                                   nrow(bad)))
    emit(CP[, .(seed, arm, arch, is_control, family, label, n_called,
                causal = tp, linked_neutral = expected_linked,
                background_neutral = fp_background, missed_causal = fn_testable,
                n_testable, precision_strict = precision, precision_all,
                recall_testable, f1, calls_per_tp)],
         "T5_called_marker_composition")

    med <- function(x) as.numeric(stats::median(as.numeric(x), na.rm = TRUE))
    cm <- CP[is_control == FALSE & family == "single method",
             .(causal = med(tp), linked_neutral = med(expected_linked),
               background_neutral = med(fp_background),
               n_called = med(n_called), n_seeds = uniqueN(seed)),
             by = .(arm, arch, label)]
    emit(cm, "B6_composition_median")
    cml <- melt(cm, id.vars = c("arm", "arch", "label", "n_called", "n_seeds"),
                measure.vars = c("causal", "linked_neutral", "background_neutral"),
                variable.name = "class", value.name = "loci")
    cml[, class := factor(class, levels = c("causal", "linked_neutral", "background_neutral"),
                          labels = c("causal", "linked to a causal locus", "background (noise)"))]
    cml[, label := factor(label, levels = cm[arm == "base",
                          .(t = sum(n_called)), by = label][order(t), label])]

    save_fig(ggplot(cml[arm == "base"], aes(label, loci, fill = class)) +
        geom_col(width = 0.72) + coord_flip() + facet_wrap(~ arch, nrow = 1, scales = "free_x") +
        scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_THRESHOLD, CLINEGO_REMOVED),
                          name = NULL) +
        labs(x = NULL, y = "markers called (median across replicates)") +
        theme_clinego_grid(),
        "B6_call_composition_absolute", 13, 5.4)
    reg("B6", "figure",
        "Absolute causal / linked / background composition of what each method calls",
        "sweep07*/MVP*_c3_any_threshold_sweep.tsv",
        "12 primary replicates, Bonferroni 0.05 at the default rung c3, both MAF arms")

    # Same partition on a log axis, because totals span 1 to >1000 markers and a
    # linear stack renders the small, precise methods invisible.
    save_fig(ggplot(cml[arm == "base" & loci > 0], aes(label, loci, fill = class)) +
        geom_col(position = position_dodge(width = 0.75), width = 0.7) +
        coord_flip() + facet_wrap(~ arch, nrow = 1) + scale_y_log10() +
        scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_THRESHOLD, CLINEGO_REMOVED),
                          name = NULL) +
        labs(x = NULL, y = "markers called, log scale (median across replicates)") +
        theme_clinego_grid(),
        "B6b_call_composition_log", 13, 5.4)
    reg("B6b", "figure", "Same partition, dodged on a log axis so small call sets stay legible",
        "sweep07*/MVP*_c3_any_threshold_sweep.tsv", "12 primary replicates, both MAF arms")

    # D8 -- the same accounting for the combination rules, which is the fair way
    # to compare a union against a single method: same units, same replicates.
    dm <- CP[is_control == FALSE,
             .(causal = med(tp), linked_neutral = med(expected_linked),
               background_neutral = med(fp_background), n_called = med(n_called)),
             by = .(arm, family, label)]
    emit(dm, "D8_composition_single_vs_combined")
    dml <- melt(dm[arm == "base"], id.vars = c("family", "label", "n_called"),
                measure.vars = c("causal", "linked_neutral", "background_neutral"),
                variable.name = "class", value.name = "loci")
    dml[, class := factor(class, levels = c("causal", "linked_neutral", "background_neutral"),
                          labels = c("causal", "linked to a causal locus", "background (noise)"))]
    save_fig(ggplot(dml, aes(reorder(label, loci, sum), loci, fill = class)) +
        geom_col(width = 0.72) + coord_flip() + facet_wrap(~ family, scales = "free") +
        scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_THRESHOLD, CLINEGO_REMOVED),
                          name = NULL) +
        labs(x = NULL, y = "markers called (median across replicates)") +
        theme_clinego_grid(),
        "D8_composition_single_vs_combined", 12, 5.4)
    reg("D8", "figure", "Absolute call composition, single methods against combination rules",
        "sweep07*/MVP*_c3_any_threshold_sweep.tsv",
        "12 primary replicates, Bonferroni 0.05, default rung")
}

# =============================================================================
# C -- structure-correction ladder
# =============================================================================
message("\n=== C: structure-correction ladder ===")
if (!is.null(lad)) {
    L <- lad[axis == "any" & !is.na(auc_pr)]
    L[, arch := relabel_arch(arch_level)]
    L[, def := auc_pr[cell == "c3"][1], by = .(seed, method)]
    L <- L[!is.na(def) & def > 0]
    L[, ratio := auc_pr / def]
    c_src <- L[, .(seed, method, cell, param, rung, arch, auc_pr,
                   auc_pr_default = def, ratio, best_f1)]
    emit(c_src, "C_ladder_long")

    c_med <- c_src[, .(ratio = median(ratio), auc_pr = median(auc_pr),
                       n_seeds = uniqueN(seed)), by = .(method, cell, arch)]
    emit(c_med, "C_ladder_median")

    # C1 -- relative. The "what the default costs" framing.
    save_fig(ggplot(c_src, aes(cell, ratio, group = seed)) +
        geom_hline(yintercept = 1, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_vline(xintercept = 3, linetype = "dashed", colour = CLINEGO_THRESHOLD,
                   linewidth = 0.35) +
        geom_line(colour = CLINEGO_NEUTRAL, alpha = 0.35, linewidth = 0.35) +
        geom_line(data = c_med, aes(cell, ratio, group = 1), colour = CLINEGO_REMOVED,
                  linewidth = 0.9) +
        facet_grid(arch ~ method) + scale_y_continuous(trans = "log2") +
        labs(x = "structure-correction rung (c3 = pipeline default)",
             y = "AUC-PR relative to the default rung (log2)") +
        theme_clinego_grid(),
        "C1_ladder_relative", 11, 7)
    reg("C1", "figure", "Cost of the default structure-correction rung, normalised per replicate",
        "report/mvp_ladders.tsv", "14 replicates x 5 rungs x 4 methods, axis truth_any")

    # C2 -- absolute. Same data, no normalisation; shows the architecture spread
    # that C1 deliberately removes.
    save_fig(ggplot(c_src, aes(cell, auc_pr, group = seed)) +
        geom_vline(xintercept = 3, linetype = "dashed", colour = CLINEGO_THRESHOLD,
                   linewidth = 0.35) +
        geom_line(colour = CLINEGO_NEUTRAL, alpha = 0.35, linewidth = 0.35) +
        geom_line(data = c_med, aes(cell, auc_pr, group = 1), colour = CLINEGO_REMOVED,
                  linewidth = 0.9) +
        facet_grid(arch ~ method, scales = "free_y") +
        labs(x = "structure-correction rung (c3 = pipeline default)", y = "AUC-PR") +
        theme_clinego_grid(),
        "C2_ladder_absolute", 11, 7)
    reg("C2", "figure", "Same ladder without normalisation",
        "report/mvp_ladders.tsv", "14 replicates x 5 rungs x 4 methods")

    # C3 -- where the best rung actually sits.
    best <- c_src[, .SD[which.max(auc_pr)], by = .(seed, method, arch)]
    best_tab <- best[, .N, by = .(method, cell, arch)]
    emit(best[, .(seed, method, arch, best_cell = cell, best_rung = rung,
                  best_auc = auc_pr, default_auc = auc_pr_default,
                  gain = auc_pr / auc_pr_default)], "C3_best_rung_per_replicate")
    save_fig(ggplot(best_tab, aes(cell, N, fill = arch)) +
        geom_col(width = 0.7) + facet_wrap(~ method, nrow = 1) +
        scale_fill_clinego(name = NULL) +
        labs(x = "rung achieving the highest AUC-PR (c3 = default)",
             y = "replicates") + theme_clinego_grid(),
        "C3_best_rung_frequency", 12, 4.2)
    reg("C3", "figure", "How often each rung is the best one",
        "report/mvp_ladders.tsv", "14 replicates x 4 methods")

    # C4 -- the size of the gain a defaults-only user forgoes.
    save_fig(ggplot(best, aes(reorder(method, ratio, median), auc_pr / auc_pr_default)) +
        geom_hline(yintercept = 1, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_boxplot(outlier.shape = NA, width = 0.6, colour = CLINEGO_COL$secondary) +
        geom_jitter(aes(colour = arch), width = 0.14, height = 0, size = 2, alpha = 0.9) +
        coord_flip() + scale_y_log10() + scale_color_clinego(name = NULL) +
        labs(x = NULL, y = "AUC-PR at the best rung / at the default rung (log scale)") +
        theme_clinego(),
        "C4_gain_over_default", 8.5, 5)
    reg("C4", "figure", "Gain available from tuning the rung, per method",
        "report/mvp_ladders.tsv", "14 replicates x 4 methods")
}
lbest <- rd(EVAL, "report06/ladder_best_rung.tsv")
if (!is.null(lbest)) emit(lbest, "C5_ladder_best_rung_11methods")

# =============================================================================
# D -- combining methods
# =============================================================================
message("\n=== D: combining methods ===")
tsum <- rd(EVAL, "report07/transfer_summary.tsv")
cwin <- rd(EVAL, "report06/combine_window_summary.tsv")
tarc <- rd(EVAL, "report07/transfer_by_arch.tsv")
sbia <- rd(EVAL, "report07/selection_bias.tsv")
tseed <- rd(EVAL, "report07/transfer_by_seed.tsv")

if (!is.null(tsum)) {
    d1_src <- melt(tsum[mode == "recal"],
                   id.vars = c("scheme_id", "label", "arm", "n_seeds", "n_beat_best",
                               "n_beat_deflt", "med_f1", "med_precision", "med_called"),
                   measure.vars = c("med_gain_default", "med_gain_best"),
                   variable.name = "baseline", value.name = "median_gain_f1")
    d1_src[, baseline := factor(baseline, levels = c("med_gain_default", "med_gain_best"),
            labels = c("vs pipeline default", "vs tuned single method"))]
    d1_src[, arm := factor(arm, levels = c("base", "m05"),
                           labels = c("MAF > 0.01", "MAF > 0.05"))]
    d1_src[, label := factor(label, levels = tsum[mode == "recal" & arm == "base"][
                                        order(med_gain_default), unique(label)])]
    emit(d1_src, "D1_scheme_gain")
    save_fig(ggplot(d1_src, aes(median_gain_f1, label, colour = baseline, shape = arm)) +
        geom_vline(xintercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
        scale_color_manual(values = c(CLINEGO_RETAINED, CLINEGO_REMOVED), name = NULL) +
        scale_shape_manual(values = c(16, 17), name = NULL) +
        labs(x = "median gain in F1", y = NULL) + theme_clinego(),
        "D1_scheme_gain", 9, 5.4)
    reg("D1", "figure", "Gain of each transferred scheme over the default AND over a tuned single method",
        "report07/transfer_summary.tsv", "8 adequately powered replicates, recalibrated arm")

    # D2 -- the operating point each scheme lands on. A scheme with a high F1 and
    # 2000 calls is not the same offer as one with the same F1 and 80 calls.
    d2_src <- tsum[, .(scheme_id, label, mode, arm, med_precision, med_f1, med_called, n_seeds)]
    emit(d2_src, "D2_scheme_operating_point")
    save_fig(ggplot(d2_src, aes(med_precision, med_f1, colour = label,
                                shape = mode, size = med_called)) +
        geom_point(alpha = 0.85) +
        facet_wrap(~ arm, labeller = as_labeller(c(base = "MAF > 0.01", m05 = "MAF > 0.05"))) +
        scale_color_manual(values = clinego_cluster_palette(uniqueN(d2_src$label)), name = NULL) +
        scale_size_continuous(trans = "log10", name = "SNPs called") +
        scale_shape_manual(values = c(verbatim = 1, recal = 16), name = NULL) +
        labs(x = "median precision", y = "median F1") +
        theme_clinego_grid() + theme(legend.box = "vertical"),
        "D2_scheme_operating_point", 12, 6)
    reg("D2", "figure", "Precision / F1 / number of calls for every scheme and transfer mode",
        "report07/transfer_summary.tsv", "8 replicates x 2 arms x 2 modes")

    # D3 -- verbatim vs recalibrated. Isolates the transferable part (structure)
    # from the untransferable part (the operating point).
    d3 <- dcast(tsum, scheme_id + label + arm ~ mode, value.var = "med_gain_default")
    if (all(c("verbatim", "recal") %in% names(d3))) {
        emit(d3, "D3_verbatim_vs_recal")
        save_fig(ggplot(d3, aes(verbatim, recal, colour = label, shape = arm)) +
            geom_abline(slope = 1, intercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
            geom_point(size = 3) +
            scale_color_manual(values = clinego_cluster_palette(uniqueN(d3$label)), name = NULL) +
            scale_shape_manual(values = c(base = 16, m05 = 17), name = NULL) +
            labs(x = "gain over default, operating point imported verbatim",
                 y = "gain over default, operating point recalibrated") +
            theme_clinego(),
            "D3_verbatim_vs_recal", 9, 6)
        reg("D3", "figure", "How much of a transferred scheme survives without recalibration",
            "report07/transfer_summary.tsv", "8 replicates x 2 arms")
    }
}
if (!is.null(cwin)) {
    emit(cwin, "D4_window_effect")
    save_fig(ggplot(cwin, aes(combine_window_kb, median_f1_delta, colour = rule, group = rule)) +
        geom_hline(yintercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
        scale_color_clinego(name = NULL) +
        labs(x = "agreement window (kb)", y = "median change in F1") +
        theme_clinego(),
        "D4_agreement_window", 8, 5)
    reg("D4", "figure", "Effect of scoring method agreement on LD blocks rather than exact SNPs",
        "report06/combine_window_summary.tsv",
        "1 replicate; union must be exactly 0 at every window (correctness check)")
}
if (!is.null(sbia)) {
    s <- sbia[!is.na(selection_bias)]
    s[, arch := relabel_arch(arch)]
    emit(s, "D5_selection_bias")
    save_fig(ggplot(s, aes(selection_bias, arch, colour = arm)) +
        geom_vline(xintercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_jitter(height = 0.16, size = 2, alpha = 0.85) +
        scale_color_manual(values = c(base = CLINEGO_RETAINED, m05 = CLINEGO_REMOVED), name = NULL) +
        labs(x = "F1 advantage of searching on the same replicate (selection bias)",
             y = NULL) + theme_clinego(),
        "D5_selection_bias", 9, 4.6)
    reg("D5", "figure", "How much of an in-sample combining gain is search, not signal",
        "report07/selection_bias.tsv", "14 replicates x 5 rungs x 2 arms")
}
if (!is.null(tarc)) emit(tarc, "D6_recommended_scheme_by_architecture")
if (!is.null(tseed)) emit(tseed, "D7_transfer_by_replicate")

# =============================================================================
# E -- window-based scanning (WZA)
#
# bonf / qval rows ONLY. A `top` cutoff is a rank rule, which defeats the point
# of windowing and is what inflated the journal-06 headline (its own section 13
# walks it back). Never restore `top` here.
# =============================================================================
message("\n=== E: window-based scanning ===")
wza_s <- rd(EVAL, "report06/wza_scores.tsv"); wza_g <- rd(EVAL, "report06/wza_gate.tsv")
if (!is.null(wza_s)) {
    w <- wza_s[adjust %in% c("bonf", "qval")]
    message(sprintf("  WZA rows kept %d of %d (`top` dropped by design)",
                    nrow(w), nrow(wza_s)))
    emit(w, "E1_wza_window_detection")
    save_fig(ggplot(w, aes(recall, precision, colour = method, shape = adjust)) +
        geom_point(size = 2.4, alpha = 0.9) +
        scale_color_manual(values = clinego_cluster_palette(uniqueN(w$method)), name = NULL) +
        scale_shape_manual(values = c(bonf = 16, qval = 17), name = NULL) +
        labs(x = "recall of causal loci", y = "window-level precision") +
        theme_clinego(),
        "E1_wza_precision_recall", 9, 5.6)
    reg("E1", "figure", "Window-level detection at a 1 kb window, Bonferroni and q-value rules only",
        "report06/wza_scores.tsv", "1 replicate (1232548); `top` rows excluded by design")
}
if (!is.null(wza_g)) emit(wza_g, "E2_wza_viability_gate")

# =============================================================================
# F -- genetic offset   (INTERIM: the `all` / `neutral_all` panels are still
# computing, so no "curated vs everything" contrast appears here.)
#
# Accuracy = -Kendall tau, higher is better, per CLAIMS.md. offset09 stores raw
# tau (negative), so the sign is flipped once, here.
# =============================================================================
message("\n=== F: genetic offset ===")
PANEL_MAP <- data.table(
    marker_set = c("adaptive", "gea_best", "gea_strict", "gea_union",
                   "gea_lfmm_only", "gea_rda_only", "gea_emmax_only"),
    set        = c("truth", "best", "intersect3", "union",
                   "solo_lfmm", "solo_rda", "solo_emmax"),
    # Label strings copied verbatim from figures/panel_summary.tsv so the
    # mapping check below is a real join, not a lookalike comparison.
    panel      = c("true QTNs (oracle)", "≥2 methods agree", "all 3 methods agree",
                   "any method (union)", "LFMM only", "RDA only", "EMMAX only"))

tau  <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
pr   <- rd(EVAL, OFF, "panel_pr_recomputed.tsv")
psum <- rd(EVAL, "figures/panel_summary.tsv")
ptst <- rd(EVAL, "figures/panel_paired_tests.tsv")
abs_counts <- rd(EVAL, "figures/table_absolute_counts.tsv")

if (!is.null(tau)) {
    tau <- merge(tau, PANEL_MAP, by = "marker_set")
    tau[, accuracy := -tau]
    tau <- merge(tau, seeds[, .(seed, arch = arch_lab, is_control, meanFst, final_LA,
                                n_causal_maf01, n_snps_rep = n_snps)],
                 by = "seed", all.x = TRUE)
    tau <- tau[is_control == FALSE]        # the two degenerate controls never aggregate
    N_OFF <- uniqueN(tau$seed)
    message(sprintf("  offset arm: %d replicates, %d panels, %d methods",
                    N_OFF, uniqueN(tau$panel), uniqueN(tau$method_label)))
    emit(tau[, .(seed, panel, method_label, arch, accuracy, final_LA, meanFst,
                 n_causal_maf01)], "F_offset_long")

    lev <- tau[, .(m = median(accuracy)), by = panel][order(m), panel]
    tau[, panel := factor(panel, levels = lev)]
    oracle_arch <- tau[panel == "true QTNs (oracle)", .(y = median(accuracy)), by = arch]
    oracle_all  <- tau[panel == "true QTNs (oracle)", median(accuracy)]

    # F1 -- panel x architecture x method. The most granular view.
    f1_src <- tau[, .(accuracy = median(accuracy), n_models = .N),
                  by = .(panel, arch, method_label)]
    emit(f1_src, "F1_accuracy_panel_arch_method")
    save_fig(ggplot(f1_src, aes(panel, accuracy, colour = method_label)) +
        geom_hline(data = oracle_arch, aes(yintercept = y), inherit.aes = FALSE,
                   linetype = "dashed", colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_point(size = 2.4, position = position_dodge(width = 0.5)) +
        facet_wrap(~ arch, nrow = 1) + coord_flip() +
        scale_color_clinego(name = NULL) +
        labs(x = NULL, y = expression("offset accuracy ("*-tau*", higher is better)")) +
        theme_clinego_grid(),
        "F1_offset_panel_by_architecture", 12, 5.6)
    reg("F1", "figure", "Offset accuracy by marker panel, architecture and offset method",
        paste0(OFF, "/phase1_seed_medians_solo.tsv; mvp_seeds.tsv"),
        sprintf("%d replicates x 112 gardens x 4 methods (INTERIM)", N_OFF))

    # F2 -- the same numbers as a distribution over replicates, not medians.
    emit(tau[, .(seed, panel, method_label, arch, accuracy)], "F2_accuracy_distribution")
    save_fig(ggplot(tau, aes(panel, accuracy)) +
        geom_hline(yintercept = oracle_all, linetype = "dashed",
                   colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_boxplot(outlier.shape = NA, width = 0.62, colour = CLINEGO_COL$secondary) +
        geom_jitter(aes(colour = arch), width = 0.16, height = 0, size = 1.4, alpha = 0.55) +
        coord_flip() + scale_color_clinego(name = NULL) +
        labs(x = NULL, y = expression("offset accuracy ("*-tau*")")) +
        theme_clinego(),
        "F2_offset_distribution", 9, 5.4)
    reg("F2", "figure", "Spread of offset accuracy across replicates, per panel",
        paste0(OFF, "/phase1_seed_medians_solo.tsv"),
        sprintf("%d replicates x 4 methods (INTERIM)", N_OFF))

    # F3 -- paired difference from the oracle, within replicate. The comparison
    # the source design cannot make, because its adaptive panel IS the oracle.
    orc <- tau[panel == "true QTNs (oracle)", .(seed, method_label, oracle = accuracy)]
    f3 <- merge(tau[panel != "true QTNs (oracle)"], orc, by = c("seed", "method_label"))
    f3[, delta := accuracy - oracle]
    emit(f3[, .(seed, panel, method_label, arch, accuracy, oracle, delta)],
         "F3_paired_vs_oracle")
    save_fig(ggplot(f3, aes(reorder(panel, delta, median), delta)) +
        geom_hline(yintercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.5) +
        geom_boxplot(outlier.shape = NA, width = 0.6, colour = CLINEGO_COL$secondary) +
        geom_jitter(aes(colour = arch), width = 0.15, height = 0, size = 1.4, alpha = 0.55) +
        coord_flip() + facet_wrap(~ method_label, nrow = 1) +
        scale_color_clinego(name = NULL) +
        labs(x = NULL, y = "accuracy minus the QTN oracle, paired within replicate") +
        theme_clinego_grid(),
        "F3_paired_vs_oracle", 13, 5.2)
    reg("F3", "figure", "Every panel scored against the true-QTN oracle, paired within replicate",
        paste0(OFF, "/phase1_seed_medians_solo.tsv"),
        sprintf("%d replicates x 4 methods (INTERIM)", N_OFF))

    # F4 -- stability. The argument for >=2 agree is the SPREAD across
    # architectures, not the peak, so plot the spread as its own quantity.
    # Two aggregations, because they do not agree and the difference is not
    # negligible. CLAIMS.md Claim 2 quotes a spread of 0.011 for ">=2 agree";
    # taking the median over replicates AND offset methods gives 0.037, which
    # changes that panel from "smallest spread by 5x" to "among the smallest".
    # The direction survives either way; the margin does not, so both are shown
    # and the manuscript must name which one it quotes.
    f4a <- tau[, .(accuracy = median(accuracy)), by = .(panel, arch)
              ][, .(spread = max(accuracy) - min(accuracy), worst = min(accuracy),
                    best = max(accuracy)), by = panel]
    f4a[, aggregation := "median over replicates and offset methods (this script)"]
    f4s <- copy(f4a)
    if (!is.null(abs_counts)) {
        ab2 <- copy(abs_counts)
        ab2[, arch := relabel_arch(architecture)]
        f4b <- ab2[, .(spread = max(accuracy) - min(accuracy), worst = min(accuracy),
                       best = max(accuracy)), by = panel]
        f4b[, aggregation := "as published in figures/table_absolute_counts.tsv"]
        f4s <- rbind(f4a, f4b)
    }
    emit(f4s, "F4_panel_stability")
    save_fig(ggplot(f4s, aes(reorder(panel, -spread), spread, fill = aggregation)) +
        geom_col(width = 0.68, position = position_dodge(width = 0.75)) + coord_flip() +
        scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_THRESHOLD), name = NULL) +
        geom_text(aes(label = sprintf("%.3f", spread)),
                  position = position_dodge(width = 0.75), hjust = -0.15, size = 2.7,
                  colour = CLINEGO_COL$fg) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
        labs(x = NULL, y = "accuracy range across the three architectures") +
        theme_clinego() + theme(legend.position = "bottom", legend.direction = "vertical"),
        "F4_panel_stability", 9.5, 5.2)
    reg("F4", "figure",
        "Architecture-to-architecture spread of each panel, under two aggregations that disagree",
        paste0(OFF, "/phase1_seed_medians_solo.tsv; figures/table_absolute_counts.tsv"),
        sprintf("%d replicates (INTERIM); see the aggregation column", N_OFF))

    # F5 -- detection versus prediction, PAIRED WITHIN REPLICATE.
    # Replicates differ ~70-fold in how many causal loci exist (7 to ~500), so
    # pooling seed x panel pairs into one rank correlation mixes that between-
    # replicate spread into the within-replicate contrast and cancels it --
    # Simpson's paradox. All three estimands are emitted so the difference is
    # visible rather than a choice made silently:
    #   panel_level  -- one point per marker panel. This is what CLAIMS.md and
    #                   DISCUSSION_NOTES.md report ("across panels").
    #   pooled       -- every seed x panel pair, unpaired.
    #   within_repl  -- computed inside each replicate, then the median across
    #                   replicates. The correct paired analogue.
    if (!is.null(pr)) {
        prm <- merge(pr, PANEL_MAP[, .(set, panel)], by = "set")
        prm[, noise_fraction := pmax(0, n - n_causal - n_linked) / n]
        rep_acc <- tau[, .(accuracy = median(accuracy)), by = .(seed, panel, arch)]
        f5 <- merge(rep_acc, prm[, .(seed, panel, n, precision, recall, noise_fraction)],
                    by = c("seed", "panel"))
        f5 <- f5[panel != "true QTNs (oracle)"]      # recall is 1 by construction
        emit(f5, "F5_detection_vs_prediction")

        panel_lvl <- f5[, .(recall = median(recall), noise_fraction = median(noise_fraction),
                            accuracy = median(accuracy)), by = panel]
        within <- f5[, .(rho_recall = sp(recall, accuracy),
                         rho_noise  = sp(noise_fraction, accuracy),
                         rho_precision = sp(precision, accuracy),
                         n_panels = .N), by = seed]
        emit(within, "F5b_within_replicate_rho")

        rho_tab <- data.table(
            estimand = rep(c("panel_level", "pooled_seed_x_panel",
                             "within_replicate_median"), each = 3),
            quantity = rep(c("rho_recall_accuracy", "rho_noise_accuracy",
                             "rho_precision_accuracy"), 3),
            spearman_rho = c(
                sp(panel_lvl$recall, panel_lvl$accuracy),
                sp(panel_lvl$noise_fraction, panel_lvl$accuracy),
                NA_real_,
                sp(f5$recall, f5$accuracy), sp(f5$noise_fraction, f5$accuracy),
                sp(f5$precision, f5$accuracy),
                median(within$rho_recall, na.rm = TRUE),
                median(within$rho_noise,  na.rm = TRUE),
                median(within$rho_precision, na.rm = TRUE)),
            n_units = c(rep(nrow(panel_lvl), 3), rep(nrow(f5), 3), rep(nrow(within), 3)),
            unit = c(rep("panel", 3), rep("seed x panel", 3), rep("replicate", 3)),
            claims_md_reported = rep(c(-0.79, -0.63, NA_real_), 3))
        emit(rho_tab, "F5c_correlation_estimands")
        print(rho_tab)

        f5l <- melt(f5, id.vars = c("seed", "panel", "arch", "accuracy"),
                    measure.vars = c("recall", "noise_fraction", "precision"),
                    variable.name = "x_metric", value.name = "x")
        f5l[, x_metric := factor(x_metric,
                levels = c("recall", "noise_fraction", "precision"),
                labels = c("causal loci recovered (recall)", "background-marker fraction",
                           "precision"))]
        ann <- data.table(
            x_metric = factor(levels(f5l$x_metric), levels = levels(f5l$x_metric)),
            lab = sprintf("within-replicate rho = %.2f\npooled rho = %.2f",
                rho_tab[estimand == "within_replicate_median", spearman_rho],
                rho_tab[estimand == "pooled_seed_x_panel", spearman_rho]))

        save_fig(ggplot(f5l, aes(x, accuracy)) +
            geom_line(aes(group = seed), colour = CLINEGO_NEUTRAL, alpha = 0.3,
                      linewidth = 0.3) +
            geom_point(aes(colour = arch), size = 1.6, alpha = 0.8) +
            geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
                      hjust = 1.05, vjust = 1.25, size = 2.9, colour = CLINEGO_COL$fg) +
            facet_wrap(~ x_metric, scales = "free_x", nrow = 1) +
            scale_color_clinego(name = NULL) +
            labs(x = NULL, y = expression("offset accuracy ("*-tau*")")) +
            theme_clinego_grid(),
            "F5_detection_vs_prediction", 13, 5)
        reg("F5", "figure",
            "Detection metrics against prediction accuracy, one line per replicate",
            paste0(OFF, "/panel_pr_recomputed.tsv; phase1_seed_medians_solo.tsv"),
            sprintf("%d replicates, %d points (INTERIM)", uniqueN(f5$seed), nrow(f5)))
    }

    # F6 -- what actually predicts offset accuracy. Local adaptation does; Fst
    # does not. This is the "is my dataset suitable" screen.
    f6 <- tau[, .(accuracy = median(accuracy)),
              by = .(seed, arch, final_LA, meanFst, n_causal_maf01, panel)]
    emit(f6, "F6_accuracy_vs_replicate_properties")
    f6l <- melt(f6, id.vars = c("seed", "arch", "panel", "accuracy"),
                measure.vars = c("final_LA", "meanFst", "n_causal_maf01"),
                variable.name = "property", value.name = "x")
    f6l[, property := factor(property,
            levels = c("final_LA", "meanFst", "n_causal_maf01"),
            labels = c("local adaptation (final_LA)", "mean Fst", "causal loci in replicate"))]
    rho6 <- f6l[, .(rho = sp(x, accuracy)), by = .(property, panel)]
    emit(rho6, "F6b_property_correlations")
    rho6a <- f6l[, .(rho = sp(x, accuracy)), by = property]
    save_fig(ggplot(f6l, aes(x, accuracy)) +
        geom_point(aes(colour = arch), size = 1.5, alpha = 0.7) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                    colour = CLINEGO_THRESHOLD, linewidth = 0.6) +
        geom_text(data = rho6a, aes(x = Inf, y = Inf, label = sprintf("rho = %.2f", rho)),
                  inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 3,
                  colour = CLINEGO_COL$fg) +
        facet_wrap(~ property, scales = "free_x", nrow = 1) +
        scale_color_clinego(name = NULL) +
        labs(x = NULL, y = expression("offset accuracy ("*-tau*")")) +
        theme_clinego_grid(),
        "F6_what_predicts_accuracy", 13, 4.8)
    reg("F6", "figure", "Replicate properties against offset accuracy (local adaptation vs Fst)",
        paste0(OFF, "/phase1_seed_medians_solo.tsv; mvp_seeds.tsv"),
        sprintf("%d replicates x 7 panels (INTERIM)", N_OFF))

    # Cross-check the panel-name mapping against the independently built summary
    # rather than trusting it: median panel size must agree.
    # Two checks, because they answer different questions.
    #  (a) snp_sets_summary covers the same replicate subset the published
    #      panel_summary was built on, so its medians must match EXACTLY. This
    #      is the real proof that marker_set -> panel is the right mapping.
    #  (b) panel_pr_recomputed covers more replicates, so its medians are
    #      expected to differ. A difference there is coverage, not a mis-mapping.
    ss <- rd(EVAL, OFF, "snp_sets_summary.tsv")
    if (!is.null(ss) && !is.null(psum)) {
        ssm <- merge(ss, PANEL_MAP[, .(set, panel)], by = "set")
        # Restrict to replicates that carry ALL SEVEN panels. That is the subset
        # figures/panel_summary.tsv was built on, so the medians must match to
        # the digit; including partially covered replicates shifts the median
        # for coverage reasons and would make a correct mapping look wrong.
        complete <- ssm[, .N, by = seed][N == uniqueN(PANEL_MAP$panel), seed]
        complete <- setdiff(complete, seeds[is_control == TRUE, seed])
        chk_a <- merge(ssm[seed %in% complete,
                           .(median_markers = as.numeric(median(n_snps))), by = panel],
                       psum[, .(panel, published = n_snps)], by = "panel")
        chk_a[, exact_match := abs(median_markers - published) < 1e-9]
        chk_a[, check := sprintf("a: %d replicates with all 7 panels, must match exactly",
                                 length(complete))]
        if (!all(chk_a$exact_match))
            message("WARN: panel mapping check (a) failed -- marker_set -> panel is wrong")
        else message(sprintf("  panel mapping verified exactly on %d panels over %d replicates",
                             nrow(chk_a), length(complete)))
        # Coverage fact the manuscript must state: offset ACCURACY aggregates
        # over more replicates than panel COMPOSITION does.
        emit(data.table(
            quantity = c("replicates with an accuracy estimate",
                         "replicates with a complete 7-panel composition"),
            n = c(uniqueN(tau$seed), length(complete))), "F_offset_coverage")
    } else chk_a <- NULL
    if (!is.null(pr) && !is.null(psum)) {
        chk_b <- merge(prm[, .(median_markers = as.numeric(median(n))), by = panel],
                       psum[, .(panel, published = n_snps)], by = "panel", all = TRUE)
        chk_b[, exact_match := NA]
        chk_b[, check := "b: wider replicate coverage, difference expected"]
        emit(rbind(chk_a, chk_b, fill = TRUE), "F_panel_mapping_check")
        print(rbind(chk_a, chk_b, fill = TRUE))
    }
}
# F8 -- the offset arm's marker panels in the same absolute units as B6/D8.
#
# Built from the PER-REPLICATE table, not from figures/table_absolute_counts.tsv.
# That published table reports a median for each class independently, and medians
# do not sum: its causal+linked+noise misses its own `markers` column by -12..+8
# on 12 of 21 rows. Small, but it means a stacked bar drawn from it is not a
# partition. Means over the same replicates do sum exactly, so the bar here is a
# true partition by construction. The published medians are still emitted, as
# T4a, for anyone who wants to quote them.
if (!is.null(pr)) {
    P <- merge(pr, PANEL_MAP[, .(set, panel)], by = "set")
    P <- merge(P, seeds[, .(seed, arch = arch_lab, is_control)], by = "seed")
    P <- P[is_control == FALSE]
    P[, noise := pmax(0, n - n_causal - n_linked)]
    f8 <- P[, .(markers = mean(n), causal = mean(n_causal),
                linked = mean(n_linked), noise = mean(noise),
                causal_loci_total = mean(causal_total),
                n_replicates = uniqueN(seed)), by = .(panel, arch)]
    f8[, resid := markers - (causal + linked + noise)]
    stopifnot(max(abs(f8$resid)) < 1e-6)      # a partition, or the bar is a lie
    acc8 <- tau[, .(accuracy = median(accuracy)), by = .(panel, arch)]
    f8 <- merge(f8, acc8, by = c("panel", "arch"), all.x = TRUE)
    emit(f8[, -"resid"], "F8_panel_composition_absolute")

    f8l <- melt(f8, id.vars = c("panel", "arch", "markers", "accuracy"),
                measure.vars = c("causal", "linked", "noise"),
                variable.name = "class", value.name = "loci")
    f8l[, class := factor(class, levels = c("causal", "linked", "noise"),
            labels = c("causal", "linked to a causal locus", "background (noise)"))]
    save_fig(ggplot(f8l, aes(reorder(panel, markers), loci, fill = class)) +
        geom_col(width = 0.72) +
        geom_text(data = unique(f8[, .(panel, arch, markers, accuracy)]),
                  aes(x = panel, y = markers, label = sprintf("acc %.2f", accuracy)),
                  inherit.aes = FALSE, hjust = -0.12, size = 2.8, colour = CLINEGO_COL$fg) +
        coord_flip() + facet_wrap(~ arch, nrow = 1, scales = "free_x") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
        scale_fill_manual(values = c(CLINEGO_RETAINED, CLINEGO_THRESHOLD, CLINEGO_REMOVED),
                          name = NULL) +
        labs(x = NULL, y = "markers in the panel (mean across replicates)") +
        theme_clinego_grid(),
        "F8_panel_composition_absolute", 13.5, 5.4)
    reg("F8", "figure",
        "Absolute causal / linked / background composition of each offset marker panel, with the accuracy it achieves",
        paste0(OFF, "/panel_pr_recomputed.tsv; phase1_seed_medians_solo.tsv"),
        "30 primary replicates, means (exact partition), INTERIM")

    # Per-replicate detail, so the spread behind those means is inspectable.
    emit(P[, .(seed, arch, panel, markers = n, causal = n_causal, linked = n_linked,
               noise, causal_loci_total = causal_total, precision, recall)],
         "T5b_offset_panel_composition_per_replicate")
}

novc <- rd(EVAL, OFF, "novelty_curve.tsv")
if (!is.null(novc)) {
    emit(novc, "F7_climate_novelty_curve")
    nv <- merge(novc, PANEL_MAP[, .(marker_set, panel)], by = "marker_set", all.x = TRUE)
    nv[is.na(panel), panel := marker_set]
    save_fig(ggplot(nv, aes(dev, -tau_med, colour = panel, group = panel)) +
        geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
        facet_wrap(~ method_label, nrow = 1) +
        scale_color_manual(values = clinego_cluster_palette(uniqueN(nv$panel)), name = NULL) +
        labs(x = "climate novelty (deviation from the training range)",
             y = expression("offset accuracy ("*-tau*")")) +
        theme_clinego_grid(),
        "F7_climate_novelty", 13, 4.6)
    reg("F7", "figure", "Offset accuracy as climate moves outside the training range",
        paste0(OFF, "/novelty_curve.tsv"), "12 replicates (earlier generation); open disagreement with the source paper")
}

# =============================================================================
# G -- comparison with published baselines
# =============================================================================
message("\n=== G: published baselines ===")
pubcmp <- rd(EVAL, "report/published_comparison.tsv")
if (!is.null(pubcmp)) {
    emit(pubcmp, "G1_published_vs_ours_detection")
    g <- pubcmp[!is.na(published) & !is.na(ours_default_c3)]
    if (nrow(g)) save_fig(ggplot(g, aes(published, ours_default_c3, colour = method,
                                        shape = axis)) +
        geom_abline(slope = 1, intercept = 0, colour = CLINEGO_THRESHOLD, linewidth = 0.4) +
        geom_point(size = 2.8, alpha = 0.9) +
        scale_color_clinego(name = NULL) +
        labs(x = "AUC-PR reported by Lotterhos 2023",
             y = "AUC-PR, this pipeline at the default rung") +
        theme_clinego(),
        "G1_published_vs_ours", 8, 6)
    reg("G1", "figure", "Our AUC-PR against the values published for the same replicates",
        "report/published_comparison.tsv",
        "per seed and axis; one published cell is flagged as a suspected column swap")
}

# =============================================================================
# Tables -- wide, so a metric can be chosen after the fact rather than before
# =============================================================================
message("\n=== tables ===")
emit(seeds[order(arch_level, -meanFst),
           .(seed, role = fifelse(is_control, "degenerate control", "primary"),
             architecture = arch_lab, n_traits, pleiotropy = ispleiotropy,
             demography = demog_level_sub, meanFst, r2_pc1_temp, r2_pc1_sal,
             k_best, n_snps, n_causal_pre, n_causal_maf01, n_causal_temp,
             n_causal_sal, final_LA)],
     "T1_replicate_manifest")

emit(dcast(mrank, method ~ arm, value.var = c("med_aucpr", "rank", "n_seeds")),
     "T2a_method_rank")
if (!is.null(lad)) {
    T2b <- lad[axis == "any", .(
        aucpr_med = median(auc_pr, na.rm = TRUE),
        aucpr_q25 = quantile(auc_pr, .25, na.rm = TRUE),
        aucpr_q75 = quantile(auc_pr, .75, na.rm = TRUE),
        bestf1_med = median(best_f1, na.rm = TRUE),
        n_seeds = uniqueN(seed)), by = .(method, cell)]
    emit(T2b[order(method, cell)], "T2b_method_by_rung_full")
}
emit(auc_arm, "T2c_aucpr_per_replicate")

anchor <- rd(ROOT, "benchmarks/mvp_anchor_schemes.tsv")
if (!is.null(anchor) && !is.null(tsum))
    emit(merge(anchor[, .(scheme_id, label, subset, rule, combine_window_kb, assignment)],
               tsum[, .(scheme_id, mode, arm, n_seeds, n_beat_best, n_beat_deflt,
                        med_gain_best, med_gain_default, med_f1, med_precision, med_called)],
               by = "scheme_id", all.y = TRUE)[order(arm, mode, -med_gain_default)],
         "T3a_combination_schemes")
rsch <- rd(EVAL, "report07/recommended_schemes.tsv")
if (!is.null(rsch)) emit(rsch, "T3b_recommended_by_architecture")
rsrc <- rd(EVAL, "report07/research.tsv")
if (!is.null(rsrc)) emit(rsrc, "T3c_in_sample_search_full")

if (!is.null(abs_counts)) emit(abs_counts, "T4a_offset_panel_composition")
if (!is.null(psum))       emit(psum,       "T4b_offset_panel_summary")
if (!is.null(ptst))       emit(ptst,       "T4c_offset_paired_tests")
if (!is.null(pr))         emit(pr,         "T4d_offset_panel_precision_recall")
strt <- rd(EVAL, OFF, "strata_tau.tsv")
if (!is.null(strt)) emit(strt, "T4e_offset_by_strata")
mrk <- rd(EVAL, OFF, "marker_set_comparison_by_method.tsv")
if (!is.null(mrk)) emit(mrk, "T4f_offset_marker_set_contrasts")

excl <- fread(file.path(ROOT, "benchmarks/mvp_method_exclusions.tsv"),
              sep = "\t", header = TRUE, fill = TRUE)
emit(excl[!grepl("^#", seed)], "S1_method_exclusions")
pub <- rd(ROOT, "benchmarks/mvp_published_baseline.tsv")
if (!is.null(pub)) emit(pub, "S2a_published_baseline_detection")
hvp <- rd(EVAL, OFF, "headline_vs_published.tsv")
if (!is.null(hvp)) emit(hvp, "S2b_published_baseline_offset")
hl <- rd(EVAL, "report06/headline.tsv")
if (!is.null(hl)) emit(hl, "S4_single_replicate_best_combinations")
svc <- rd(EVAL, "report06/solo_vs_combine.tsv")
if (!is.null(svc)) emit(svc, "S5_solo_vs_combine_single_replicate")

fwrite(rbindlist(MANIFEST), file.path(OUT, "figure_manifest.tsv"), sep = "\t")
message(sprintf("\nDone. %d figures registered. Outputs in %s",
                length(MANIFEST), OUT))
