#!/usr/bin/env Rscript
# =============================================================================
# mvp_main_figure.R -- the single main-text figure for the MVP simulation arm,
# plus the supplementary figure that justifies dropping RDA-corrected.
#
# Reads ONLY scored tables already on disk. Fits nothing, runs no pipeline mode.
# Every panel writes the exact values it plots to a same-stem .tsv.
#
# Cohort: 90 primary replicates (30 per architecture) + 2 degenerate Est-Clines
# controls, excluded from every aggregate. Source dir: mvp_eval/offset11.
#
#   A  what each selection rule returns   -- absolute causal/linked/background
#   B  signal-vs-noise TRAJECTORY         -- usable vs background across architecture
#   C  prediction accuracy per panel      -- the headline
#   D  robustness across architecture     -- with the random panel as null yardstick
#   E  will it work on my data            -- local adaptation vs Fst
#   S-RDA  a/b/c                          -- why RDA-corrected is dropped
#
# Accuracy = -Kendall's tau (offset vs realised common-garden fitness), higher
# better. NEVER abs(tau): verified 0-0.04% positive for GF/LFMM2/RDA-uncorrected
# but 2.65% for RDA-corrected, which is why that method is excluded here.
#
# Usage:
#   Rscript /pipeline/benchmarks/mvp_main_figure.R
#   OFFSET_DIR=offset11 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_main ...
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(cowplot)
    library(ggrepel)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset11")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_main"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

# =============================================================================
# Palette -- ltc "minou" (Theodosiou), hex read verbatim from the package source
# R/ltc_functions.R. Hardcoded rather than depending on the `ltc` package so the
# figure builds in the shipped image without a Dockerfile change.
#
# ONE rule: the same entity is the same colour in every panel. Hue encodes the
# GROUP; lightness separates members inside a group, so the legend is learned once.
#   combined rules -> teal    single methods -> amber    controls -> grey-blue
#   oracle         -> navy (a reference, but the gold standard)
#   red            -> RESERVED for bad/noise/dropped, nothing else
# =============================================================================
MINOU <- c(teal = "#00798c", red = "#d1495b", amber = "#edae49",
           sage = "#66a182", navy = "#2e4057", grey = "#8d96a3")

# Deterministic lighter/darker members of one hue family.
shades <- function(base, n) {
    if (n == 1L) return(base)
    c(base,
      grDevices::colorRampPalette(c(base, "white"))(6)[2:(1 + (n - 1) %/% 2 + (n - 1) %% 2)],
      grDevices::colorRampPalette(c(base, "black"))(6)[2:(1 + (n - 1) %/% 2)])[seq_len(n)]
}
green3 <- shades(MINOU[["sage"]],  3)   # combining rules
amber3 <- shades(MINOU[["amber"]], 3)   # single methods
grey3  <- shades(MINOU[["grey"]],  3)   # uninformative references

# The true-QTN bar was navy at full depth, which swallowed its own error bar.
# Lightened ~45% toward white so the whisker reads against it.
QTN_BAR  <- grDevices::colorRampPalette(c(MINOU[["navy"]], "white"))(10)[5]
# The two reference LINES carry meaning rather than identity:
#   green = the ceiling being aimed at (perfect knowledge of the causal loci)
#   red   = the floor below which curating markers was not worth doing (all SNPs)
REF_GOOD <- MINOU[["sage"]]
REF_BAD  <- MINOU[["red"]]

RULE_COL <- c(
    "2/3 methods"          = green3[1],   # the recommendation
    "3/3 methods"          = green3[2],
    "1/3 methods"          = green3[3],
    "RDA"                  = amber3[1],
    "EMMAX"                = amber3[2],
    "LFMM"                 = amber3[3],
    "causal loci"            = QTN_BAR,
    "all loci"             = grey3[3],
    "neutral"              = MINOU[["grey"]],
    "random, size-matched" = grey3[2])

# Semantic, at the user's direction: green = signal (the same greens the combining
# rules carry in C/D), red = the thing you do not want. Overrides the earlier
# entity-colour scheme, which put causal in navy and neutral in grey.
CLASS_COL <- c("causal"  = shades(MINOU[["sage"]], 3)[3],   # dark green
               "linked"  = MINOU[["sage"]],                 # green
               "neutral" = MINOU[["red"]])                  # red = wrong

METHOD_COL <- c("GFoffset"        = MINOU[["teal"]],
                "LFMM2offset"     = MINOU[["amber"]],
                "RDA-uncorrected" = MINOU[["sage"]],
                "RDA-corrected"   = MINOU[["red"]])   # red = dropped

GROUP_LTY <- c(combined = "solid", single = "22", reference = "11")


# The corpus ships "oliogenic" (sic) as a level name.
ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
relabel_arch <- function(x) factor(x, levels = ARCH_LEVELS, labels = ARCH_LABELS)

# RDA-corrected is pathological (see S-RDA); excluded from every main-figure aggregate.
DROP_METHOD <- "RDA-corrected"

rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f)
}
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t"); invisible(dt)
}
save_fig <- function(p, stem, w, h) {
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h)
    message(sprintf("  OK %s", stem))
}

# ---------------------------------------------------------------- shared data
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
seeds[, arch_lab := relabel_arch(arch_level)]
PRIM  <- mvp_prim(seeds)
stopifnot(nrow(PRIM) == mvp_n_expect(), !any(is.na(PRIM$arch_lab)))
message(sprintf("primary replicates: %d (%s)", nrow(PRIM),
                paste(sprintf("%s=%d", levels(PRIM$arch_lab), table(PRIM$arch_lab)),
                      collapse = ", ")))

# Marker-panel vocabulary. Two naming schemes exist on disk: the offset tables use
# gea_* / adaptive, the precision-recall table uses best / intersect3 / truth.
PANELS <- data.table(
    acc_key = c("gea_best", "gea_strict", "gea_rda_only", "gea_emmax_only",
                "gea_lfmm_only", "gea_union", "adaptive", "all", "neutral",
                "random_matched"),
    pr_key  = c("best", "intersect3", "solo_rda", "solo_emmax",
                "solo_lfmm", "union", "truth", "all", "neutral_all",
                "rand_best1"),
    # One naming style for the combining rules: how many of the three GEA methods
    # had to flag a marker. "union" and "all 3 agree" said the same thing two ways.
    label   = c("2/3 methods", "3/3 methods", "RDA", "EMMAX",
                "LFMM", "1/3 methods", "causal loci",
                "all loci", "neutral", "random, size-matched"),
    role    = c("combined", "combined", "single", "single",
                "single", "combined", "reference", "reference", "reference",
                "reference"))
# Display order: obtainable rules first (combined then single), references last.
PANEL_ORDER <- PANELS$label[c(1, 2, 6, 3, 4, 5, 7, 8, 9, 10)]
OBTAINABLE  <- PANELS$label[PANELS$role != "reference"]

# Per-replicate accuracy, one row per seed x panel x offset method.
acc <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
acc <- acc[method_label != DROP_METHOD]
acc[, accuracy := -tau]                       # NOT abs(): see header
acc <- merge(acc, PRIM[, .(seed, arch_lab)], by = "seed")   # drops the 2 controls
acc <- merge(acc, PANELS[, .(marker_set = acc_key, label, role)], by = "marker_set")
acc[, label := factor(label, levels = PANEL_ORDER)]
stopifnot(all(acc$accuracy > 0))               # no anti-prediction survives here
message(sprintf("accuracy rows: %d (%d seeds, %d offset methods)",
                nrow(acc), uniqueN(acc$seed), uniqueN(acc$method_label)))

# Per-replicate panel composition, one row per seed x panel.
comp <- rd(EVAL, OFF, "panel_pr_recomputed.tsv")
comp <- merge(comp, PRIM[, .(seed, arch_lab)], by = "seed")
comp <- merge(comp, PANELS[, .(set = pr_key, label, role)], by = "set")
comp[, `:=`(usable     = n_causal + n_linked,
            background = n - n_causal - n_linked)]
comp[, bg_frac := fifelse(n > 0, background / n, NA_real_)]
comp[, label := factor(label, levels = PANEL_ORDER)]

# =============================================================================
# A -- what each selection rule actually returns
# =============================================================================
# Obtainable rules ONLY: the oracle and the random floor are references, not
# candidates, and including them compresses the axis (the abs_D convention).
# Bars are ordered by usable markers POOLED over architectures, so the order is
# identical in every facet and rules stay comparable across them.
message("\n=== A: composition of what each rule returns ===")
a_dt <- comp[label %in% OBTAINABLE]
a_ord <- a_dt[, .(u = median(usable)), by = label][order(u)]$label
a_dt[, label := factor(as.character(label), levels = as.character(a_ord))]

a_src <- a_dt[, .(causal = median(n_causal), linked = median(n_linked),
                  background = median(background), total = median(as.numeric(n)),
                  usable = median(usable), n_seeds = .N),
              by = .(arch_lab, label)]
a_long <- melt(a_src, id.vars = c("arch_lab", "label", "total", "usable", "n_seeds"),
               measure.vars = c("causal", "linked", "background"),
               variable.name = "class", value.name = "markers")
a_long[, class := factor(class, levels = c("causal", "linked", "background"),
                         labels = names(CLASS_COL))]   # -> causal / linked / neutral
emit(a_long, "A_rule_composition")

# Facet labels carry the causal-locus count, so "polygenic" is a number not a word.
n_causal_arch <- comp[label == "causal loci",
                      .(nc = round(median(n_causal))), by = arch_lab]
arch_facet <- setNames(sprintf("%s\n(~%d causal loci)", n_causal_arch$arch_lab,
                               n_causal_arch$nc), n_causal_arch$arch_lab)

pA <- ggplot(a_long, aes(label, markers, fill = class)) +
    geom_col(width = 0.7) +
    geom_text(data = a_src, aes(label, total, label = as.integer(usable)),
              inherit.aes = FALSE, hjust = -0.25, size = 2.7, colour = CLINEGO_COL$fg) +
    coord_flip() +
    facet_wrap(~ arch_lab, nrow = 1, labeller = labeller(arch_lab = arch_facet)) +
    scale_fill_manual(values = CLASS_COL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(tag = "A", x = NULL, y = "markers selected (median)") +
    theme_clinego() +
    theme(legend.position = "bottom", strip.text = element_text(size = 8.5),
          axis.text.y = element_text(size = 8),
          plot.tag = element_text(face = "bold", size = 13)) +
    guides(fill = guide_legend(nrow = 1))
save_fig(pA, "A_rule_composition", w = 11, h = 4.6)

# =============================================================================
# C -- HEADLINE: prediction accuracy per marker panel
# =============================================================================
message("\n=== C: prediction accuracy per panel ===")
# Per-replicate median across offset methods first, then across replicates, so a
# replicate carrying more methods does not outweigh one carrying fewer.
c_seed <- acc[, .(accuracy = median(accuracy)), by = .(seed, label, role)]
c_src  <- c_seed[, .(accuracy = median(accuracy),
                     q25 = quantile(accuracy, .25), q75 = quantile(accuracy, .75),
                     n_seeds = .N), by = .(label, role)]

REF <- "2/3 methods"
paired_vs <- function(dt, unit_cols) {
    w <- dcast(dt, as.formula(paste(paste(unit_cols, collapse = " + "), "~ label")),
               value.var = "accuracy")
    rbindlist(lapply(setdiff(PANEL_ORDER, REF), function(other) {
        if (!other %in% names(w)) return(NULL)
        d <- w[[other]]; r <- w[[REF]]; ok <- !is.na(d) & !is.na(r)
        if (!any(ok)) return(NULL)
        data.table(label = other, n_pairs = sum(ok),
                   win_pct = 100 * mean(r[ok] > d[ok]),
                   median_delta = median(r[ok] - d[ok]),
                   p_wilcox = tryCatch(
                       stats::wilcox.test(r[ok], d[ok], paired = TRUE)$p.value,
                       error = function(e) NA_real_))
    }))
}
pairs      <- paired_vs(c_seed, "seed")                     # n = 90 replicates
pairs_fine <- paired_vs(acc, c("seed", "method_label"))     # n = 270 cells
emit(merge(pairs, pairs_fine[, .(label, n_pairs_fine = n_pairs, win_pct_fine = win_pct,
                                 median_delta_fine = median_delta, p_fine = p_wilcox)],
           by = "label", all = TRUE), "B_paired_both_aggregations")
emit(merge(c_src, pairs, by = "label", all.x = TRUE)[order(-accuracy)],
     "B_accuracy_by_panel")
print(merge(c_src, pairs, by = "label", all.x = TRUE)[order(-accuracy)])

# Only two references survive: the ceiling (the causal loci themselves) and the
# floor a practitioner actually faces (all loci, i.e. no curation at all).
# `neutral` and `random, size-matched` were controls for an earlier question --
# whether the GEA signal is real at all -- which is no longer in doubt.
#
# Three grids, as before, so differences are compared WITHIN a group first; the
# two guides run across all three so the groups stay mutually comparable.
GRP_LAB <- c(reference = "reference", single = "single method",
             combined  = "combining rules")
KEEP_B <- c(OBTAINABLE, "causal loci", "all loci")
c_bar  <- c_src[label %in% KEEP_B]
c_bar[, grp := factor(GRP_LAB[role], levels = GRP_LAB)]
setorder(c_bar, grp, accuracy)
c_bar[, label := factor(as.character(label), levels = as.character(label))]
ref_qtn <- c_bar[label == "causal loci", accuracy]
ref_all <- c_bar[label == "all loci",    accuracy]

x_lo <- min(c_bar$q25) - 0.03
pC <- ggplot(c_bar, aes(label, accuracy, fill = label)) +
    geom_hline(yintercept = ref_qtn, linetype = "22", colour = REF_GOOD, linewidth = 0.6) +
    geom_hline(yintercept = ref_all, linetype = "22", colour = REF_BAD,  linewidth = 0.6) +
    geom_col(width = 0.68) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.2,
                  colour = CLINEGO_COL$fg, linewidth = 0.4) +
    geom_text(aes(y = q75, label = sprintf("%.3f", accuracy)),
              hjust = -0.22, size = 2.9, colour = CLINEGO_COL$fg) +
    coord_flip(ylim = c(x_lo, max(c_bar$q75) + 0.055)) +
    facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = RULE_COL, guide = "none") +
    labs(tag = "B", x = NULL, y = "Accuracy") +
    theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          strip.background = element_blank(),
          strip.placement  = "outside",
          strip.text.y.left = element_text(angle = 90, size = 8.5, face = "bold"),
          panel.border = element_blank(),
          panel.spacing.y = unit(0.5, "lines"),
          axis.text.y = element_text(size = 8.5))
save_fig(pC, "B_accuracy_by_panel", w = 7.8, h = 5.2)

# =============================================================================
# D -- robustness across architecture
# =============================================================================
# Same grouping and the same colours as C, with the group carried a second time by
# linetype, so the reader never has to re-learn which series is which.
message("\n=== D: robustness across architecture ===")
d_src <- acc[, .(accuracy = median(accuracy)), by = .(seed, label, role, arch_lab)][
    , .(accuracy = median(accuracy), n_seeds = .N), by = .(label, role, arch_lab)]
d_spread <- d_src[, .(spread = max(accuracy) - min(accuracy),
                      worst = min(accuracy)), by = .(label, role)][order(spread)]
emit(merge(d_src, d_spread, by = c("label", "role")), "C_robustness")
print(d_spread)

# One panel, not three. The six choosable rules are solid coloured lines; the two
# references are dashed and carry the same green/red meaning as panel B, so the
# reader learns the code once. Lines are labelled at their right end -- no legend
# lookup, which is what made the earlier shade-within-family version unreadable.
d_rule <- d_src[label %in% OBTAINABLE]
d_refs <- d_src[label %in% c("causal loci", "all loci")]
d_rule[, label := factor(as.character(label), levels = OBTAINABLE)]
null_spread <- d_spread[label == "random, size-matched", spread]

xn <- nlevels(d_src$arch_lab)
lab_all <- rbind(
    d_rule[arch_lab == levels(arch_lab)[1L]][, face := "plain"],
    d_refs[arch_lab == levels(arch_lab)[1L]][, face := "italic"])
# Left margin sized to the widest label rather than to a fixed fraction of the axis.
.w_in <- max(strwidth(lab_all$label, units = "inches",
                      cex = 2.9 * .pt / par("ps")))          # widest label, inches
.panel_in <- 7.8 - 1.0                                        # plot width minus margins/axis
LAB_PAD <- max(0.14, (.w_in / .panel_in) * (nlevels(d_src$arch_lab) - 1) * 1.25 + 0.05)

pD <- ggplot(mapping = aes(as.integer(arch_lab), accuracy, group = label)) +
    geom_line(data = d_refs, aes(colour = label), linetype = "22", linewidth = 0.75) +
    geom_line(data = d_rule, aes(colour = label), linewidth = 0.95) +
    geom_point(data = d_rule, aes(colour = label), size = 2.2) +
    geom_text(data = lab_all, aes(x = 1 - 0.03, y = accuracy, label = label,
                                  colour = label, fontface = face),
              hjust = 1, size = 2.9, show.legend = FALSE) +
    scale_x_continuous(breaks = seq_len(xn), labels = levels(d_src$arch_lab),
                       expand = expansion(add = c(LAB_PAD, 0.06))) +
    scale_colour_manual(values = c(RULE_COL[OBTAINABLE],
                                   "causal loci" = REF_GOOD, "all loci" = REF_BAD),
                        guide = "none") +
    labs(tag = "C", x = NULL, y = "Accuracy") +
    theme_clinego() +
    theme(plot.tag = element_text(face = "bold", size = 13),
          panel.border = element_blank(),
          axis.text.x = element_text(size = 8.5))
save_fig(pD, "C_robustness", w = 7.8, h = 5.0)

# =============================================================================
# Composite main figure -- A / B / C
# =============================================================================
message("\n=== composing main figure ===")
bottom <- plot_grid(pC, pD, nrow = 1, rel_widths = c(1, 1.02))
main   <- plot_grid(pA, bottom, ncol = 1, rel_heights = c(1, 1.05))
clinego_save_both(file.path(OUT, "MAIN_FIGURE_simulation"), main, w = 13, h = 9.5)
message("  OK MAIN_FIGURE_simulation")

# =============================================================================
# S-RDA -- why RDA-corrected is dropped
# =============================================================================
message("\n=== S-RDA: justifying the RDA-corrected drop ===")
gp <- rd(EVAL, OFF, "garden_performance.tsv")

# (a) fraction of gardens where the model anti-predicts (tau > 0)
s_a <- gp[!is.na(tau), .(n_gardens = .N, n_anti = sum(tau > 0)), by = method_label]
s_a[, pct_anti := 100 * n_anti / n_gardens]
setorder(s_a, -pct_anti)
emit(s_a, "S_RDA_a_anti_prediction")
print(s_a)

pSa <- ggplot(s_a, aes(reorder(method_label, pct_anti), pct_anti)) +
    geom_col(aes(fill = method_label == DROP_METHOD), width = 0.62) +
    geom_text(aes(label = sprintf("%.2f%%", pct_anti)), hjust = -0.15, size = 3.2,
              colour = CLINEGO_COL$fg) +
    coord_flip(ylim = c(0, max(s_a$pct_anti) * 1.35)) +
    scale_fill_manual(values = c(`TRUE` = MINOU[["red"]], `FALSE` = MINOU[["teal"]]),
                      guide = "none") +
    labs(title = "a  Gardens where the model anti-predicts",
         subtitle = "τ > 0: offset rises where fitness rises",
         x = NULL, y = "% of gardens") +
    theme_clinego()

# (b) replicates where the offset is invariant to the marker set
gea_keys <- PANELS[role != "reference" | acc_key == "adaptive", acc_key]
allm <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
allm <- merge(allm, PRIM[, .(seed)], by = "seed")
s_b <- allm[marker_set %in% gea_keys,
            .(n_distinct = uniqueN(round(tau, 12))), by = .(seed, method_label)][
    , .(n_seeds = .N, n_invariant = sum(n_distinct == 1L)), by = method_label]
s_b[, pct_invariant := 100 * n_invariant / n_seeds]
setorder(s_b, -pct_invariant)
emit(s_b, "S_RDA_b_marker_set_invariance")
print(s_b)

pSb <- ggplot(s_b, aes(reorder(method_label, pct_invariant), pct_invariant)) +
    geom_col(aes(fill = method_label == DROP_METHOD), width = 0.62) +
    geom_text(aes(label = sprintf("%d/%d", n_invariant, n_seeds)), hjust = -0.15,
              size = 3.2, colour = CLINEGO_COL$fg) +
    coord_flip(ylim = c(0, max(s_b$pct_invariant) * 1.5 + 1)) +
    scale_fill_manual(values = c(`TRUE` = MINOU[["red"]], `FALSE` = MINOU[["teal"]]),
                      guide = "none") +
    labs(title = "b  Replicates where the offset ignores the marker set",
         subtitle = "identical accuracy across all GEA panels",
         x = NULL, y = "% of replicates") +
    theme_clinego()

# (c) accuracy on the >=2-agree panel, per architecture, all four methods
s_c <- allm[marker_set == "gea_best"]
s_c <- merge(s_c, PRIM[, .(seed, arch_lab)], by = "seed")
s_c <- s_c[, .(accuracy = median(-tau), n_seeds = .N), by = .(method_label, arch_lab)]
emit(s_c, "S_RDA_c_accuracy_by_method")
print(dcast(s_c, method_label ~ arch_lab, value.var = "accuracy"))

pSc <- ggplot(s_c, aes(arch_lab, accuracy, colour = method_label,
                       group = method_label)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
    scale_colour_manual(values = METHOD_COL, name = NULL) +
    labs(title = "c  Accuracy on the ≥2-agree panel",
         x = NULL, y = "accuracy  (−Kendall's τ)") +
    theme_clinego()

srda <- plot_grid(pSa, pSb, pSc, nrow = 1, rel_widths = c(1, 1, 1.15))
clinego_save_both(file.path(OUT, "S_RDA_why_dropped"), srda, w = 15, h = 4.6)
message("  OK S_RDA_why_dropped")

message("\nAll outputs in: ", OUT)
