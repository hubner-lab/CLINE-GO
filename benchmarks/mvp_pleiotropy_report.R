#!/usr/bin/env Rscript
# =============================================================================
# mvp_pleiotropy_report.R -- pleiotropy as its own factor on the SS-Clines corpus.
#
# WHY THIS EXISTS. Every SS-Clines block crosses {pleiotropy, no-pleiotropy} with
# {equal-S, unequal-S} inside each genic level, and the step-12 report shows that
# cross only as 12 heatmap cells (mvp_block_report.R, block_cell_delta_heatmap).
# Nothing isolates pleiotropy: its marginal effect, whether it moves the ORACLE
# (the causal loci themselves) or only the GEA panel, whether it interacts with
# the selection regime, and whether it changes what the panel is made of. This
# script answers those four questions, stratified by genic level throughout --
# a 2 x 2 pooled over genic levels would average a sign change (see the cells).
#
# UNITS, stated once. accuracy = -tau (Kendall, median over the 112 landscape
# gardens per replicate x panel x method -- that is phase1_seed_medians_solo.tsv).
#   PAIR  one row per replicate x working method: delta = accuracy(>=2-of-3 panel)
#         - accuracy(causal loci), paired within the method. This is the unit of
#         mvp_block_report.R's per-cell table and heatmap (n = seeds x 3).
#   SEED  one row per replicate: median over the three working methods of each
#         panel's accuracy, THEN the difference -- the unit of mvp_oracle_stats.R
#         (delta_per_seed.tsv, delta_by_arch.tsv) and of every test here. The
#         three methods of one replicate are not independent observations.
# Points and reach% use PAIR; medians, IQRs and every p-value use SEED.
#
# TESTS. Pleiotropy groups are DIFFERENT replicates, so the contrast is unpaired:
# Mann-Whitney (wilcox.test, two-sided) with the Hodges-Lehmann shift and its
# 95% CI as the effect size, Holm over the three genic levels. The interaction
# (c) is a difference of medians with a percentile-bootstrap CI over replicates
# within cells -- descriptive, not a test.
#
# Reads only tables already on disk. Fits nothing, runs no pipeline mode.
#
# Usage (pooled corpus):
#   OFFSET_DIR=offset12_ssclines_pooled MVP_ARM=primary_ssclines \
#   MVP_ADDED=ssclines_nvar_mvar,ssclines_ncline_ns,ssclines_nequal_mconst,ssclines_ncline_ctredge,ssclines_nequal_mbreaks \
#   MVP_N_EXPECT=600 FIG_OUT=/pipeline/benchmarks/mvp_eval/figures_ssclines_pooled_j15 \
#     Rscript /pipeline/benchmarks/mvp_pleiotropy_report.R
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OFF  <- Sys.getenv("OFFSET_DIR", "offset12_ssclines_pooled")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_pooled"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))

rd <- function(...) {
    f <- file.path(...)
    if (!file.exists(f)) stop("MISSING: ", f)
    fread(f, colClasses = c("seed" = "character"))
}
rd_opt <- function(...) {          # optional cross-check tables; some carry no seed column
    f <- file.path(...)
    if (!file.exists(f)) return(NULL)
    x <- fread(f)
    if ("seed" %in% names(x)) x[, seed := as.character(seed)]
    x
}
emit <- function(dt, stem) {
    fwrite(dt, file.path(OUT, paste0(stem, ".tsv")), sep = "\t")
    message(sprintf("  wrote %s.tsv (%d rows)", stem, nrow(dt)))
    invisible(dt)
}
save_fig <- function(p, stem, w, h) {
    if (!mvp_subtitle_on()) p <- p + labs(subtitle = NULL)   # MVP_SUBTITLE=0: caption carries it
    clinego_save_both(file.path(OUT, stem), p, w = w, h = h); message("  OK ", stem)
}
message("figure subtitles: ", if (mvp_subtitle_on()) "on" else "off (MVP_SUBTITLE=0)")
SUB <- function(extra = NULL) {
    s <- sprintf("%s -- %d replicates", mvp_arm_label(), nrow(PRIM))
    if (!is.null(extra)) s <- paste(s, extra, sep = "; ")
    s
}

ARCH_LEVELS <- c("oliogenic", "mod-polygenic", "highly-polygenic")   # corpus typo
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
WORKING    <- c("GFoffset", "LFMM2offset", "RDA-uncorrected")   # RDA-corrected excluded, as everywhere
ORACLE_SET <- "adaptive"          # the true causal loci
BEST_SET   <- "gea_best"          # the >=2-of-3 recommended panel
PANEL_LAB  <- c("causal loci", ">=2-of-3 panel")
PLEIO_LAB  <- c("no pleiotropy", "pleiotropy")
SYM_SHAPES <- c(`equal-S` = 16, `unequal-S` = 17)
ZERO_LINE  <- geom_hline(yintercept = 0, linetype = "dashed", colour = CLINEGO_THRESHOLD, linewidth = 0.4)

# ------------------------------------------------------------------- manifest
seeds <- rd(ROOT, "benchmarks/mvp_seeds.tsv")
PRIM  <- mvp_prim(seeds)
PRIM[, arch_lab := factor(arch_level, levels = ARCH_LEVELS, labels = ARCH_LABELS)]
PRIM[, sub_level := sub("^[^_]+_", "", architecture)]
# "unequal-S" first: grepl("equal-S") matches both (same trap as mvp_regime_panel.R).
PRIM[, symmetry := factor(fifelse(grepl("unequal-S", architecture), "unequal-S", "equal-S"),
                          levels = names(SYM_SHAPES))]
PRIM[, pleiotropy := factor(ispleiotropy, levels = 0:1, labels = PLEIO_LAB)]
PRIM[, block := sub("^ssclines_", "", group)]
stopifnot(nrow(PRIM) == mvp_n_expect(),
          uniqueN(PRIM$seed) == nrow(PRIM),
          !anyNA(PRIM$arch_lab), !anyNA(PRIM$pleiotropy),
          # the manifest flag and the architecture label must agree on every replicate
          all(PRIM[, (ispleiotropy == 0L) == grepl("no-pleiotropy", sub_level)]))
CELLS <- PRIM[, .N, by = .(arch_lab, pleiotropy, symmetry)]
setorder(CELLS, arch_lab, pleiotropy, symmetry)
stopifnot(nrow(CELLS) == 12L)
if (any(CELLS$N < 50L))
    message("!! per-cell n < 50 (single block?) -- every test below is underpowered")
DESIGN <- PRIM[, .N, by = .(block, added, arch_lab, pleiotropy, symmetry)]
setorder(DESIGN, block, arch_lab, pleiotropy, symmetry)
emit(DESIGN, "pleio_design_counts")
message(sprintf("design: %d blocks x %d genic levels x 2 pleiotropy x 2 symmetry; %d replicates (%d / %d pleiotropy)",
                uniqueN(PRIM$block), nlevels(PRIM$arch_lab), nrow(PRIM),
                sum(PRIM$ispleiotropy == 0L), sum(PRIM$ispleiotropy == 1L)))

# ------------------------------------------------------------------ accuracy
ACC <- rd(EVAL, OFF, "phase1_seed_medians_solo.tsv")
ACC <- ACC[seed %in% PRIM$seed & method_label %in% WORKING &
           marker_set %in% c(ORACLE_SET, BEST_SET)]
mvp_require_rows(ACC, "phase1 medians")
ACC[, accuracy := -tau]             # NEVER abs()
ACC <- merge(ACC, PRIM[, .(seed, arch_lab, sub_level, symmetry, pleiotropy, block, final_LA)], by = "seed")
ACC[, panel := factor(marker_set, levels = c(ORACLE_SET, BEST_SET), labels = PANEL_LAB)]
stopifnot(nrow(ACC) == nrow(PRIM) * length(WORKING) * 2L)   # both panels complete on every seed

# Cross-check 1: phase1_seed_medians_solo IS the landscape-garden median that
# mvp_block_report.R recomputes from garden_performance. If both are in OUT, they
# must agree to floating-point precision, or the two scripts describe different data.
BLA <- rd_opt(OUT, "block_accuracy_vs_LA.tsv")
if (!is.null(BLA)) {
    chk <- merge(ACC[, .(seed, panel = as.character(panel), method_label, accuracy)],
                 BLA[, .(seed, panel, method_label, accuracy_block = accuracy)],
                 by = c("seed", "panel", "method_label"))
    stopifnot(nrow(chk) == nrow(ACC), max(abs(chk$accuracy - chk$accuracy_block)) < 1e-9)
    message(sprintf("cross-check: phase1 medians == block_accuracy_vs_LA (%d rows, max|diff| %.1e)",
                    nrow(chk), max(abs(chk$accuracy - chk$accuracy_block))))
}

# ------------------------------------------------- paired delta, two units
PAIR <- dcast(ACC, seed + arch_lab + sub_level + symmetry + pleiotropy + block + method_label ~ marker_set,
              value.var = "accuracy")
PAIR[, delta := get(BEST_SET) - get(ORACLE_SET)]
stopifnot(all(is.finite(PAIR$delta)), nrow(PAIR) == nrow(PRIM) * length(WORKING))
emit(PAIR[, .(seed, arch_lab, sub_level, symmetry, pleiotropy, block, method_label,
              acc_causal = get(ORACLE_SET), acc_panel = get(BEST_SET), delta)],
     "pleio_delta_per_seed_method")

SEED_ACC <- ACC[, .(accuracy = median(accuracy), n_methods = .N),
                by = .(seed, arch_lab, sub_level, symmetry, pleiotropy, block, final_LA, marker_set)]
stopifnot(all(SEED_ACC$n_methods == length(WORKING)))
SEED <- dcast(SEED_ACC, seed + arch_lab + sub_level + symmetry + pleiotropy + block + final_LA ~ marker_set,
              value.var = "accuracy")
SEED[, delta := get(BEST_SET) - get(ORACLE_SET)]
stopifnot(nrow(SEED) == nrow(PRIM), all(is.finite(SEED$delta)))
emit(SEED[, .(seed, arch_lab, sub_level, symmetry, pleiotropy, block, final_LA,
              acc_causal = get(ORACLE_SET), acc_panel = get(BEST_SET), delta)],
     "pleio_delta_per_seed")

# Cross-check 2: SEED delta must equal mvp_oracle_stats.R's delta_per_seed for the
# >=2-of-3 panel -- same aggregation, so the headline number in delta_by_arch.tsv
# and every median printed here are the same quantity.
DPS <- rd_opt(OUT, "delta_per_seed.tsv")
if (!is.null(DPS)) {
    chk <- merge(SEED[, .(seed, delta)], DPS[label == "2/3 methods", .(seed, delta_os = delta)], by = "seed")
    stopifnot(nrow(chk) == nrow(SEED), max(abs(chk$delta - chk$delta_os)) < 1e-9)
    message(sprintf("cross-check: SEED delta == delta_per_seed.tsv (2/3 methods) (%d seeds, max|diff| %.1e)",
                    nrow(chk), max(abs(chk$delta - chk$delta_os))))
}

# ----------------------------------------------------------------- stats
# Mann-Whitney rank-sum, two-sided; Hodges-Lehmann shift (yes - no) with 95% CI.
mw <- function(yes, no) {
    w <- tryCatch(suppressWarnings(stats::wilcox.test(yes, no, conf.int = TRUE)),
                  error = function(e) NULL)
    if (is.null(w)) return(list(hl = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
    list(hl = unname(w$estimate), lo = w$conf.int[1], hi = w$conf.int[2], p = w$p.value)
}
pleio_test <- function(dt, value, by = NULL) {
    j <- quote({
        yes <- get(value)[pleiotropy == PLEIO_LAB[2]]
        no  <- get(value)[pleiotropy == PLEIO_LAB[1]]
        t <- mw(yes, no)
        .(n_no = length(no), n_yes = length(yes),
          median_no = median(no), median_yes = median(yes),
          diff_medians = median(yes) - median(no),
          hl_shift = t$hl, hl_lo = t$lo, hl_hi = t$hi, p_mw = t$p)
    })
    out <- if (length(by)) dt[, eval(j), by = by] else dt[, eval(j)]
    if ("arch_lab" %in% names(out)) out[, arch_lab := as.character(arch_lab)]
    # Holm over the genic levels within whatever else `by` holds.
    holm_by <- setdiff(by, "arch_lab")
    if (length(holm_by)) out[, p_holm := p.adjust(p_mw, "holm"), by = holm_by]
    else                 out[, p_holm := p.adjust(p_mw, "holm")]
    out
}
summ <- function(dt, value, by, unit) {
    dt[, .(unit = unit, n = .N,
           median = median(get(value)),
           q25 = quantile(get(value), .25), q75 = quantile(get(value), .75),
           reach_pct = 100 * mean(get(value) >= 0)),
       by = by]
}

# =============================================================================
# (a) marginal effect of pleiotropy, by genic level, symmetry pooled
# =============================================================================
A_PAIR <- summ(PAIR, "delta", c("arch_lab", "pleiotropy"), "seed x method")
A_SEED <- summ(SEED, "delta", c("arch_lab", "pleiotropy"), "seed")
A_TAB  <- rbind(A_PAIR, A_SEED); setorder(A_TAB, unit, arch_lab, pleiotropy)
emit(A_TAB, "pleio_delta_marginal")
A_STAT <- rbind(pleio_test(SEED, "delta", "arch_lab"),
                pleio_test(SEED, "delta")[, arch_lab := "pooled"],
                use.names = TRUE)
emit(A_STAT, "pleio_delta_marginal_stats")
message("\n=== (a) >=2-of-3 minus causal loci, pleiotropy vs none, per genic level (SEED unit) ===")
print(A_STAT[, .(arch_lab, n_no, n_yes, median_no = round(median_no, 4), median_yes = round(median_yes, 4),
                 hl_shift = round(hl_shift, 4), hl_lo = round(hl_lo, 4), hl_hi = round(hl_hi, 4),
                 p_mw = signif(p_mw, 3), p_holm = signif(p_holm, 3))])

P_A <- ggplot(PAIR, aes(pleiotropy, delta)) +
    ZERO_LINE +
    geom_jitter(aes(colour = pleiotropy, shape = symmetry), width = 0.22, height = 0,
                alpha = 0.35, size = 1.1, show.legend = c(colour = FALSE, shape = TRUE)) +
    geom_pointrange(data = A_SEED, aes(y = median, ymin = q25, ymax = q75),
                    colour = CLINEGO_COL$fg, size = 0.45, linewidth = 0.7) +
    geom_text(data = A_PAIR, aes(y = Inf, label = sprintf("%.0f%% reach", reach_pct)),
              vjust = 1.4, size = 2.9, colour = CLINEGO_COL$fg) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_color_clinego() + scale_shape_manual(values = SYM_SHAPES) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = NULL, y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         shape = NULL,
         title = "Pleiotropy effect on the recommended panel, by genic level",
         subtitle = SUB("points: replicate x method; median and IQR over replicates"))
save_fig(P_A, "pleio_delta_marginal", w = 9, h = 4.5)

AM_TAB  <- summ(PAIR, "delta", c("arch_lab", "method_label", "pleiotropy"), "seed")
emit(AM_TAB, "pleio_delta_marginal_by_method")
AM_STAT <- pleio_test(PAIR, "delta", c("arch_lab", "method_label"))
emit(AM_STAT, "pleio_delta_marginal_by_method_stats")
P_AM <- ggplot(PAIR, aes(pleiotropy, delta)) +
    ZERO_LINE +
    geom_jitter(aes(colour = pleiotropy, shape = symmetry), width = 0.22, height = 0,
                alpha = 0.35, size = 1.0, show.legend = c(colour = FALSE, shape = TRUE)) +
    geom_pointrange(data = AM_TAB, aes(y = median, ymin = q25, ymax = q75),
                    colour = CLINEGO_COL$fg, size = 0.4, linewidth = 0.6) +
    geom_text(data = AM_TAB, aes(y = Inf, label = sprintf("%.0f%%", reach_pct)),
              vjust = 1.4, size = 2.7, colour = CLINEGO_COL$fg) +
    facet_grid(method_label ~ arch_lab) +
    scale_color_clinego() + scale_shape_manual(values = SYM_SHAPES) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = NULL, y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         shape = NULL,
         title = "Pleiotropy effect by genic level and offset method",
         subtitle = SUB("one point per replicate"))
save_fig(P_AM, "pleio_delta_marginal_by_method", w = 10, h = 7.5)

# =============================================================================
# (b) does pleiotropy move the oracle itself, or only the panel?
# =============================================================================
B_TAB <- SEED_ACC[, .(unit = "seed", n = .N, median = median(accuracy),
                      q25 = quantile(accuracy, .25), q75 = quantile(accuracy, .75)),
                  by = .(arch_lab, panel = factor(marker_set, levels = c(ORACLE_SET, BEST_SET), labels = PANEL_LAB),
                         pleiotropy)]
setorder(B_TAB, arch_lab, panel, pleiotropy)
emit(B_TAB, "pleio_accuracy_absolute")
SEED_ACC[, panel := factor(marker_set, levels = c(ORACLE_SET, BEST_SET), labels = PANEL_LAB)]
B_STAT <- pleio_test(SEED_ACC, "accuracy", c("arch_lab", "panel"))
emit(B_STAT, "pleio_accuracy_absolute_stats")
message("\n=== (b) absolute accuracy, pleiotropy vs none, per genic level x panel (SEED unit) ===")
print(B_STAT[, .(arch_lab, panel, median_no = round(median_no, 4), median_yes = round(median_yes, 4),
                 hl_shift = round(hl_shift, 4), hl_lo = round(hl_lo, 4), hl_hi = round(hl_hi, 4),
                 p_mw = signif(p_mw, 3), p_holm = signif(p_holm, 3))])

P_B <- ggplot(SEED_ACC, aes(pleiotropy, accuracy, fill = panel)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.7,
                 position = position_dodge(width = 0.8), colour = CLINEGO_COL$fg, linewidth = 0.35) +
    geom_point(aes(shape = symmetry, group = panel), alpha = 0.35, size = 0.9,
               colour = CLINEGO_COL$fg,
               position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8)) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_clinego() + scale_shape_manual(values = SYM_SHAPES) +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = NULL, y = expression("accuracy = " * -tau), fill = NULL, shape = NULL,
         title = "Absolute accuracy of causal loci and of the panel, by pleiotropy",
         subtitle = SUB("methods pooled within replicate")) +
    guides(shape = guide_legend(override.aes = list(alpha = 1, size = 2)))
save_fig(P_B, "pleio_accuracy_absolute", w = 9, h = 4.8)

# =============================================================================
# (c) pleiotropy x symmetry interaction, per genic level
# =============================================================================
C_PAIR <- summ(PAIR, "delta", c("arch_lab", "pleiotropy", "symmetry", "sub_level"), "seed x method")
C_SEED <- summ(SEED, "delta", c("arch_lab", "pleiotropy", "symmetry", "sub_level"), "seed")
C_TAB  <- rbind(C_PAIR, C_SEED); setorder(C_TAB, unit, arch_lab, pleiotropy, symmetry)
emit(C_TAB, "pleio_cell_delta_2x2")

# Cross-check 3: the seed x method cells are the 12 cells of block_cell_delta_best_vs_oracle.
BCD <- rd_opt(OUT, "block_cell_delta_best_vs_oracle.tsv")
if (!is.null(BCD)) {
    chk <- merge(C_PAIR[, .(arch_lab = as.character(arch_lab), sub_level, n, median)],
                 BCD[, .(arch_lab, sub_level, n_block = n, median_block = median_delta)],
                 by = c("arch_lab", "sub_level"))
    stopifnot(nrow(chk) == 12L, all(chk$n == chk$n_block),
              max(abs(chk$median - chk$median_block)) < 1e-9)
    message("cross-check: seed x method cells == block_cell_delta_best_vs_oracle (12 cells)")
}

# Difference of differences on the SEED unit:
#   [median(pleio | unequal) - median(no | unequal)] - [median(pleio | equal) - median(no | equal)]
# Percentile bootstrap over replicates within each of the four cells.
set.seed(15)
B_BOOT <- 2000L
did <- function(d) {
    m <- d[, .(m = median(delta)), by = .(pleiotropy, symmetry)]
    g <- function(p, s) m[pleiotropy == p & symmetry == s, m]
    (g(PLEIO_LAB[2], "unequal-S") - g(PLEIO_LAB[1], "unequal-S")) -
    (g(PLEIO_LAB[2], "equal-S")   - g(PLEIO_LAB[1], "equal-S"))
}
C_STAT <- SEED[, {
    obs  <- did(.SD)
    boot <- replicate(B_BOOT, did(.SD[, .SD[sample(.N, replace = TRUE)], by = .(pleiotropy, symmetry)]))
    .(n = .N, pleio_effect_unequal = .SD[symmetry == "unequal-S", median(delta[pleiotropy == PLEIO_LAB[2]]) - median(delta[pleiotropy == PLEIO_LAB[1]])],
      pleio_effect_equal   = .SD[symmetry == "equal-S",   median(delta[pleiotropy == PLEIO_LAB[2]]) - median(delta[pleiotropy == PLEIO_LAB[1]])],
      diff_in_diff = obs,
      boot_lo = unname(quantile(boot, .025)), boot_hi = unname(quantile(boot, .975)),
      n_boot = B_BOOT)
}, by = arch_lab]
emit(C_STAT, "pleio_interaction_stats")
message("\n=== (c) pleiotropy x symmetry: difference of differences (SEED unit, bootstrap 95% CI) ===")
print(C_STAT[, .(arch_lab, n, pleio_effect_unequal = round(pleio_effect_unequal, 4),
                 pleio_effect_equal = round(pleio_effect_equal, 4),
                 diff_in_diff = round(diff_in_diff, 4),
                 boot_lo = round(boot_lo, 4), boot_hi = round(boot_hi, 4))])

P_C <- ggplot(C_SEED, aes(pleiotropy, median, colour = symmetry, group = symmetry)) +
    ZERO_LINE +
    geom_line(position = position_dodge(width = 0.3), linewidth = 0.6) +
    geom_pointrange(aes(ymin = q25, ymax = q75, shape = symmetry),
                    position = position_dodge(width = 0.3), size = 0.5, linewidth = 0.7) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_color_clinego() + scale_shape_manual(values = SYM_SHAPES) +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = NULL, y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         colour = NULL, shape = NULL,
         title = "Pleiotropy x selection symmetry, by genic level",
         subtitle = SUB("median and IQR over replicates; methods pooled within replicate"))
save_fig(P_C, "pleio_interaction_2x2", w = 9, h = 4.2)

CM_TAB <- summ(PAIR, "delta", c("arch_lab", "method_label", "pleiotropy", "symmetry"), "seed")
emit(CM_TAB, "pleio_cell_delta_2x2_by_method")
P_CM <- ggplot(CM_TAB, aes(pleiotropy, median, colour = symmetry, group = symmetry)) +
    ZERO_LINE +
    geom_line(position = position_dodge(width = 0.3), linewidth = 0.6) +
    geom_pointrange(aes(ymin = q25, ymax = q75, shape = symmetry),
                    position = position_dodge(width = 0.3), size = 0.45, linewidth = 0.6) +
    facet_grid(method_label ~ arch_lab) +
    scale_color_clinego() + scale_shape_manual(values = SYM_SHAPES) +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = NULL, y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         colour = NULL, shape = NULL,
         title = "Pleiotropy x selection symmetry, by genic level and offset method",
         subtitle = SUB("median and IQR over replicates"))
save_fig(P_CM, "pleio_interaction_2x2_by_method", w = 10, h = 7.5)

P_CH <- ggplot(C_PAIR, aes(pleiotropy, symmetry, fill = median)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.3f\n%.0f%% reach", median, reach_pct)),
              size = 3, colour = "grey15") +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_gradient2(low = CLINEGO_REMOVED, mid = "grey92", high = CLINEGO_RETAINED,
                         midpoint = 0, name = "median delta") +
    theme_clinego() +
    labs(x = NULL, y = NULL,
         title = "Recommended panel minus causal loci, pleiotropy x symmetry",
         subtitle = SUB("replicate x method pairs per cell"))
save_fig(P_CH, "pleio_interaction_heatmap", w = 9, h = 3.6)

# =============================================================================
# (d) what the >=2-of-3 panel is made of, by pleiotropy
# =============================================================================
# Same classes and the same definition as panel A of the main figure
# (mvp_main_figure_v2.R): usable = causal + linked, background = the rest.
MINOU <- c(red = "#d1495b", sage = "#66a182")
shades <- function(base, n) {
    if (n == 1L) return(base)
    c(base,
      grDevices::colorRampPalette(c(base, "white"))(6)[2:(1 + (n - 1) %/% 2 + (n - 1) %% 2)],
      grDevices::colorRampPalette(c(base, "black"))(6)[2:(1 + (n - 1) %/% 2)])[seq_len(n)]
}
CLASS_COL <- c(causal = shades(MINOU[["sage"]], 3)[3], linked = MINOU[["sage"]],
               neutral = MINOU[["red"]])

comp <- rd(EVAL, OFF, "panel_pr_recomputed.tsv")
comp <- merge(comp, PRIM[, .(seed, arch_lab, pleiotropy, symmetry)], by = "seed")
stopifnot(uniqueN(comp$seed) == nrow(PRIM))
comp[, `:=`(usable = n_causal + n_linked, background = n - n_causal - n_linked)]
truth <- comp[set == "truth", .(seed, n_causal_truth = n_causal)]
best  <- merge(comp[set == "best"], truth, by = "seed")
stopifnot(nrow(best) == nrow(PRIM))

D_SRC <- best[, .(causal = median(n_causal), linked = median(n_linked),
                  background = median(background), total = median(as.numeric(n)),
                  usable = median(usable), precision = median(precision), recall = median(recall),
                  n_causal_truth = median(n_causal_truth),
                  n_seeds = .N, n_unusable = sum(!usable_for_offset)),
              by = .(arch_lab, pleiotropy)]
setorder(D_SRC, arch_lab, pleiotropy)
D_LONG <- melt(D_SRC, id.vars = c("arch_lab", "pleiotropy", "total", "usable", "precision", "recall",
                                  "n_causal_truth", "n_seeds", "n_unusable"),
               measure.vars = c("causal", "linked", "background"),
               variable.name = "class", value.name = "markers")
D_LONG[, class := factor(class, levels = c("causal", "linked", "background"), labels = names(CLASS_COL))]
emit(D_LONG, "pleio_panel_composition")
D_CELLS <- best[, .(causal = median(n_causal), linked = median(n_linked),
                    background = median(background), total = median(as.numeric(n)),
                    precision = median(precision), recall = median(recall), n_seeds = .N),
                by = .(arch_lab, pleiotropy, symmetry)]
setorder(D_CELLS, arch_lab, pleiotropy, symmetry)
emit(D_CELLS, "pleio_panel_composition_cells")
D_STAT <- rbind(pleio_test(best, "precision", "arch_lab")[, quantity := "precision"],
                pleio_test(best, "recall",    "arch_lab")[, quantity := "recall"],
                pleio_test(best, "n_causal",  "arch_lab")[, quantity := "n_causal"],
                pleio_test(best, "n",         "arch_lab")[, quantity := "n_markers"])
emit(D_STAT, "pleio_panel_composition_stats")
message("\n=== (d) >=2-of-3 panel composition, pleiotropy vs none (medians) ===")
print(D_SRC[, .(arch_lab, pleiotropy, causal, linked, background, total, precision = round(precision, 3),
                recall = round(recall, 3), n_causal_truth, n_seeds, n_unusable)])

P_D <- ggplot(D_LONG, aes(pleiotropy, markers, fill = class)) +
    geom_col(width = 0.7) +
    geom_text(data = D_SRC, aes(pleiotropy, total, label = as.integer(usable)),
              inherit.aes = FALSE, hjust = -0.25, size = 2.7, colour = CLINEGO_COL$fg) +
    coord_flip() +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_fill_manual(values = CLASS_COL, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    theme_clinego() + theme(legend.position = "right", legend.justification = "top") +
    labs(x = NULL, y = "markers in the >=2-of-3 panel (median per replicate)",
         title = "Panel composition by pleiotropy",
         subtitle = SUB("label: median usable (causal + linked) markers")) +
    guides(fill = guide_legend(ncol = 1))
save_fig(P_D, "pleio_panel_composition", w = 10.4, h = 3.0)

# =============================================================================
# (e) is any of this one demography's doing? delta by block
# =============================================================================
E_TAB <- summ(SEED, "delta", c("block", "arch_lab"), "seed")
E_TAB <- merge(E_TAB, unique(PRIM[, .(block, added)]), by = "block")
E_POOL <- summ(SEED, "delta", "block", "seed")[, arch_lab := "pooled"]
E_POOL <- merge(E_POOL, unique(PRIM[, .(block, added)]), by = "block")
E_ALL <- rbind(copy(E_TAB)[, arch_lab := as.character(arch_lab)], E_POOL, use.names = TRUE)
setorder(E_ALL, block, arch_lab)
emit(E_ALL, "pleio_delta_by_block")
if (uniqueN(SEED$block) == 1L)
    message("!! one block only -- the by-block figures have a single column")

P_E <- ggplot(E_TAB, aes(block, median)) +
    ZERO_LINE +
    geom_pointrange(aes(ymin = q25, ymax = q75), colour = CLINEGO_COL$fg, size = 0.5, linewidth = 0.7) +
    facet_wrap(~ arch_lab, nrow = 1) +
    theme_clinego() +
    labs(x = "block (demography)", y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         title = "Recommended panel minus causal loci, by block",
         subtitle = SUB("median and IQR over replicates"))
save_fig(P_E, "pleio_delta_by_block", w = 9, h = 4)

EP_TAB <- summ(SEED, "delta", c("block", "arch_lab", "pleiotropy"), "seed")
EP_TAB <- merge(EP_TAB, unique(PRIM[, .(block, added)]), by = "block")
setorder(EP_TAB, block, arch_lab, pleiotropy)
emit(EP_TAB, "pleio_delta_by_block_pleiotropy")
EP_EFF <- dcast(EP_TAB, block + added + arch_lab ~ pleiotropy, value.var = "median")
setnames(EP_EFF, PLEIO_LAB, c("median_no", "median_yes"))
EP_EFF[, pleio_effect := median_yes - median_no]
EP_EFF[, sign := fifelse(pleio_effect < 0, "-", fifelse(pleio_effect > 0, "+", "0"))]
emit(EP_EFF, "pleio_effect_by_block")
message("\n=== (e) pleiotropy effect (median yes - median no) per block x genic level ===")
print(dcast(EP_EFF, block ~ arch_lab, value.var = "pleio_effect"))

P_EP <- ggplot(EP_TAB, aes(block, median, colour = pleiotropy)) +
    ZERO_LINE +
    geom_pointrange(aes(ymin = q25, ymax = q75), position = position_dodge(width = 0.5),
                    size = 0.45, linewidth = 0.6) +
    facet_wrap(~ arch_lab, nrow = 1) +
    scale_color_clinego() +
    theme_clinego() + theme(legend.position = "bottom") +
    labs(x = "block (demography)", y = expression("delta = accuracy(panel) - accuracy(causal loci)"),
         colour = NULL,
         title = "Pleiotropy effect by block and genic level",
         subtitle = SUB("median and IQR over replicates")) +
    guides(colour = guide_legend(override.aes = list(size = 0.6)))
save_fig(P_EP, "pleio_delta_by_block_pleiotropy", w = 9, h = 4.2)

message("\nPleiotropy report written to ", OUT)
