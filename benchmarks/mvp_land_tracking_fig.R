#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_tracking_fig.R -- review variants of one figure: what each method's SNP set holds and follows, beside the
# offset outcome and the structure of each landscape (plan ~/.claude/plans/let-s-run-the-lanpscape-dynamic-pony.md,
# Phase 1). Descriptive only: measured quantities side by side, no explanatory text in the plots.
#   V1 rows       one row per landscape: a composition of the set, b clines of the set and of its test parts,
#                 c % replicates in ranks 4-6 (median of the three offset methods, t_ranks.tsv)
#   V2 plane      |r| with temperature against |r| with salinity, per landscape, sets and test parts as medians with
#                 interquartile ranges; outcome bars below
#   V3 structure  % replicates in ranks 4-6 within each landscape by r2(neutral PC1, temperature) band and by
#                 selection regime
# Inputs: landscape_ms/explain/{set_tracking,set_outcome}.tsv (benchmarks/mvp_land_set_tracking.R),
#         landscape_ms/tables/{t_ranks,set_order}.tsv (benchmarks/mvp_land_tables.R).
# Writes landscape_ms/explain/figures/track_*.{svg,png} (180 mm, PNG 600 dpi), track_*_preview.png (200 dpi),
#        figures.tsv and tracking_fig_numbers.tsv (every plotted value).
# Usage: Rscript benchmarks/mvp_land_tracking_fig.R
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
LMS  <- file.path(ROOT, "benchmarks/mvp_eval/landscape_ms")
EXP  <- file.path(LMS, "explain")
OUT  <- file.path(EXP, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(ROOT, "benchmarks/ms_style.R"))

# ------------------------------------------------------------------ vocabulary (copied from mvp_land_figures.R)
LAND_LEV <- c("SS-Clines", "SS-Mtn", "Est-Clines")
RULE_COL <- c("LFMM" = MINOU[["red"]], "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE)
SET_LAB  <- c("LFMM" = "LFMM2", "RDA" = "RDA", "EMMAX" = "EMMAX")
BANNED <- c(set_name = "/3 methods", panel = "\\bpanels?\\b", rule = "(?<!tie )\\brules?\\b",
            background = "background", linked_neutral = "linked-neutral",
            linkage_group = "linkage group", share = "\\bshares?\\b")
strip_code <- function(x) {
    x <- gsub("`[^`]*`", " ", x)
    x <- gsub("<code>.*?</code>", " ", x, perl = TRUE)
    gsub("[^ ()]*[_/][^ ()]*\\.(R|tsv|md|html|sh|zip|svg|png)\\b", " ", x, perl = TRUE)
}
vocab_check <- function(x, where) {
    x <- strip_code(as.character(unlist(x)))
    hit <- Filter(length, lapply(BANNED, function(b) grep(b, x, ignore.case = TRUE, perl = TRUE, value = TRUE)))
    if (length(hit)) stop(where, ": banned wording for a reader: ",
                          paste(sprintf("%s in \"%s\"", names(hit), vapply(hit, function(h) substr(h[1], 1, 160), "")),
                                collapse = "; "))
    invisible(length(x))
}
grob_text <- function(g) {
    out <- character()
    walk <- function(x) {
        if (!is.null(x$label) && (is.character(x$label) || is.expression(x$label) || is.factor(x$label)))
            out <<- c(out, as.character(x$label))
        if (is.data.frame(x$data) && "label" %in% names(x$data)) out <<- c(out, as.character(x$data$label))
        if (inherits(x, "gtable")) for (y in x$grobs) walk(y)
        if (length(x$children)) for (y in x$children) walk(y)
    }
    walk(g)
    unique(out[!is.na(out) & nzchar(out)])
}
fig_text <- function(p) { grDevices::pdf(NULL); on.exit(grDevices::dev.off()); grob_text(ggplotGrob(p)) }
FIGS <- list()
save_fig <- function(stem, p, w, h, what, texts) {
    n_txt <- vocab_check(c(texts, if (inherits(p, "ggplot")) fig_text(p)), paste0(stem, " (figure text)"))
    ms_save(file.path(OUT, stem), p, w, h)
    ggsave(file.path(OUT, paste0(stem, "_preview.png")), p, width = w, height = h, units = "mm", dpi = 200, bg = "white")
    FIGS[[stem]] <<- data.table(stem = stem, w_mm = w, h_mm = h, what = what)
    message(sprintf("  %-22s %3.0f x %3.0f mm  (%d text strings checked)", stem, w, h, n_txt))
}

# ------------------------------------------------------------------ inputs
TRK <- fread(file.path(EXP, "set_tracking.tsv"), colClasses = c(seed = "character"))
OC  <- fread(file.path(EXP, "set_outcome.tsv"), colClasses = c(seed = "character"))
stopifnot(all(TRK[part == "all", uniqueN(seed), by = landscape]$V1 == 600L), uniqueN(OC$seed) == 1800L)
ORD <- fread(file.path(LMS, "tables/set_order.tsv"))
SET_ORDER <- ORD[key %in% names(SET_LAB)][order(order), key]           # the pooled order, three methods
stopifnot(identical(SET_ORDER, c("RDA", "EMMAX", "LFMM")))
TRK <- TRK[landscape %in% LAND_LEV]
REP <- unique(TRK[, .(landscape, seed, r2_pc1_temp, r2_pc1_sal)])
STR <- REP[, .(t = median(r2_pc1_temp), s = median(r2_pc1_sal)), by = landscape]
STRIP <- setNames(sprintf("%s\nPC1 r² with\ntemperature %.2f\nsalinity %.2f", STR$landscape, STR$t, STR$s), STR$landscape)
STRIP1 <- setNames(sprintf("%s · PC1 r² with temperature %.2f, salinity %.2f", STR$landscape, STR$t, STR$s), STR$landscape)
land_f <- function(x) factor(x, levels = LAND_LEV)
NUM <- list()
num <- function(D, variant, plot) NUM[[length(NUM) + 1L]] <<- copy(D)[, `:=`(variant = variant, plot = plot)]

# predictor colours, one per predictor across every plot
P_COL <- c(temperature = MINOU[["navy"]], salinity = MINOU[["teal"]])
CLS_LEV <- c("Causal, temperature", "Causal, salinity", "Causal, other", "Linked SNPs", "False positives (neutral SNPs)")
CLS_COL <- setNames(c(MINOU[["navy"]], MINOU[["teal"]], MINOU[["sage"]], "#a9c4cb", "#d3d7dc"), CLS_LEV)

# ------------------------------------------------------------------ data of the plots
# a: composition, mean over replicates of each class's % of the set (so each bar adds to 100)
A <- TRK[part == "all", .(landscape, seed, set, n_all = n)]
CMP <- merge(TRK[set %in% SET_ORDER & part %in% c("causal temperature", "causal salinity", "causal both",
                                                  "causal unassigned", "linked", "neutral")],
             A, by = c("landscape", "seed", "set"))
CMP[, cls := fcase(part == "causal temperature", CLS_LEV[1], part == "causal salinity", CLS_LEV[2],
                   part %in% c("causal both", "causal unassigned"), CLS_LEV[3], part == "linked", CLS_LEV[4],
                   part == "neutral", CLS_LEV[5])]
stopifnot(CMP[, all(n_all > 0)])
CM <- CMP[, .(n = sum(n), n_all = n_all[1]), by = .(landscape, seed, set, cls)][
    , .(pct = mean(100 * n / n_all), replicates = .N), by = .(landscape, set, cls)]
stopifnot(CM[, abs(sum(pct) - 100) < 1e-9, by = .(landscape, set)]$V1)
num(CM[, .(landscape, set, part = cls, measure = "mean % of the set", value = pct, n = replicates)], "V1", "a")

# b: clines of the whole set and of the test parts, median over replicates (IQR kept for V2)
ROWS <- data.table(set = c("RDA", "EMMAX", "EMMAX", "EMMAX", "LFMM", "LFMM", "LFMM", "truth"),
                   part = c("all", "all", "temperature test only", "salinity test only", "all",
                            "temperature test only", "salinity test only", "all"),
                   lab = c("RDA", "EMMAX", "temperature test", "salinity test", "LFMM2",
                           "temperature test", "salinity test", "Causal loci"))
ROWS[, key := paste(set, part)][, ord := .I]
CL <- merge(TRK[, .(landscape, seed, set, part, n, r_temp, r_sal)], ROWS, by = c("set", "part"))
CLm <- CL[n > 0, .(temperature = median(r_temp, na.rm = TRUE), salinity = median(r_sal, na.rm = TRUE),
                   t25 = quantile(r_temp, .25, na.rm = TRUE, names = FALSE), t75 = quantile(r_temp, .75, na.rm = TRUE, names = FALSE),
                   s25 = quantile(r_sal, .25, na.rm = TRUE, names = FALSE), s75 = quantile(r_sal, .75, na.rm = TRUE, names = FALSE),
                   replicates = .N, median_snps = median(n)), by = .(landscape, set, part, lab, key, ord)]
num(melt(CLm, id.vars = c("landscape", "set", "part", "lab", "key", "ord", "replicates", "median_snps"),
         variable.name = "measure", value.name = "value")[, .(landscape, set, part, measure = paste("|r|", measure), value, n = replicates)],
    "V1 V2", "b / plane")
BOTH <- TRK[set %in% c("LFMM", "EMMAX") & part == "both tests", .(replicates_with_both = sum(n > 0), max_both = max(n)), by = .(landscape, set)]
print(BOTH)

# c: % replicates in ranks 4-6, the values and CIs of t_ranks.tsv (median of the three offset methods)
TRp <- fread(file.path(LMS, "tables/t_ranks.tsv"))[offset_method == "Median of the three" & statistic == "ranks 4–6" &
                                                   landscape %in% LAND_LEV & set %in% SET_LAB]
TRp[, set := names(SET_LAB)[match(set, SET_LAB)]]
stopifnot(nrow(TRp) == 9L)
chk <- merge(OC[, .(pct = 100 * mean(w456)), by = .(landscape, set)], TRp[, .(landscape, set, value)], by = c("landscape", "set"))
stopifnot(max(abs(chk$pct - chk$value)) < 1e-9)
num(TRp[, .(landscape, set, part = "all", measure = "% replicates in ranks 4–6 [95% CI]", value, n, lo, hi)], "V1 V2", "c")

# ------------------------------------------------------------------ plot builders
set_f <- function(x) factor(x, levels = rev(SET_ORDER))
lab_set <- function(x) unname(SET_LAB[as.character(x)])
p_comp <- function(strips = TRUE) {
    D <- copy(CM)[, `:=`(landscape = land_f(landscape), set = set_f(set), cls = factor(cls, levels = rev(CLS_LEV)))]
    ggplot(D, aes(x = pct, y = set, fill = cls)) +
        geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
        facet_grid(landscape ~ ., switch = "y", labeller = labeller(landscape = if (strips) STRIP else setNames(LAND_LEV, LAND_LEV))) +
        scale_fill_manual(values = CLS_COL, breaks = CLS_LEV, name = NULL) +
        scale_y_discrete(labels = lab_set) +
        scale_x_continuous(limits = c(0, 100.01), breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0.02))) +
        labs(x = "% of the set's SNPs\n(mean over replicates)", y = NULL) +
        guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
        theme_ms() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = MS_TICK),
                           legend.key.size = unit(2.6, "mm"), legend.location = "plot", legend.justification = "left",
                           legend.key.spacing.x = unit(1.5, "mm"))
}
p_clines <- function() {
    D <- copy(CLm)[, `:=`(landscape = land_f(landscape), row = factor(key, levels = rev(ROWS$key)))]
    L <- melt(D, id.vars = c("landscape", "row", "set"), measure.vars = c("temperature", "salinity"),
              variable.name = "predictor", value.name = "r")
    # thin rules between the set groups (RDA | EMMAX and its two tests | LFMM2 and its two tests | causal loci)
    SEP <- data.table(y = length(ROWS$key) + 1 - c(1, 4, 7) - 0.5)
    ggplot(D, aes(y = row)) +
        geom_hline(data = SEP, aes(yintercept = y), colour = "#c9ced4", linewidth = 0.25) +
        geom_segment(aes(x = salinity, xend = temperature, yend = row), colour = "#b7bcc3", linewidth = 0.5) +
        geom_point(data = L, aes(x = r, colour = predictor), size = 1.5) +
        facet_grid(landscape ~ .) +
        scale_colour_manual(values = P_COL, name = "|r| with") +
        scale_y_discrete(labels = setNames(ROWS$lab, ROWS$key)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
        labs(x = "|r| of allele frequency with the\npredictor across demes (median)", y = NULL) +
        theme_ms() + theme(strip.text.y = element_blank(), panel.grid.major.x = element_line(colour = "#eceef1", linewidth = 0.25))
}
p_outcome <- function(facet = c("rows", "cols")) {
    facet <- match.arg(facet)
    D <- copy(TRp)[, `:=`(landscape = land_f(landscape), set = set_f(set))]
    ggplot(D, aes(x = value, y = set, fill = set)) +
        geom_col(width = 0.62) +
        geom_linerange(aes(xmin = lo, xmax = hi), colour = MS_INK, linewidth = 0.3) +
        geom_vline(xintercept = 50, linetype = "22", colour = MS_INK, linewidth = 0.3) +
        (if (facet == "rows") facet_grid(landscape ~ .) else facet_grid(. ~ landscape)) +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        scale_y_discrete(labels = lab_set) +
        scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0.03))) +
        labs(x = if (facet == "rows") "% replicates in ranks 4–6\n(median of 3 offset methods)" else
                     "% replicates in ranks 4–6 (median of the three offset methods)", y = NULL) +
        theme_ms() + theme(strip.text.y = element_blank())
}

# ------------------------------------------------------------------ V1 rows
pA <- p_comp(); pB <- p_clines(); pC <- p_outcome("rows")
v1 <- cowplot::plot_grid(pA, pB, pC, nrow = 1, rel_widths = c(68, 70, 42), align = "h", axis = "tb",
                         labels = c("a", "b", "c"), label_size = 8, label_fontface = "bold")
save_fig("track_v1_rows", v1, 180, 140, "V1: composition | clines | outcome, one row per landscape",
         c(fig_text(pA), fig_text(pB), fig_text(pC)))

# ------------------------------------------------------------------ V2 plane
PART_SHAPE <- c("whole set" = 21, "temperature test" = 24, "salinity test" = 22)
D2 <- copy(CLm)[, `:=`(landscape = land_f(landscape),
                       shape = fifelse(part == "all", "whole set", fifelse(part == "temperature test only", "temperature test", "salinity test")),
                       who = fifelse(set == "truth", "Causal loci", lab_set(set)))]
WHO_COL <- c(setNames(RULE_COL[SET_ORDER], lab_set(SET_ORDER)), "Causal loci" = MS_REF)
p2a <- ggplot(D2, aes(x = temperature, y = salinity)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#9aa1aa", linewidth = 0.3) +
    geom_errorbar(aes(ymin = s25, ymax = s75, colour = who), width = 0, linewidth = 0.3, alpha = 0.8) +
    geom_errorbar(aes(xmin = t25, xmax = t75, colour = who), width = 0, linewidth = 0.3, alpha = 0.8, orientation = "y") +
    geom_point(aes(fill = who, shape = shape), size = 2, colour = MS_INK, stroke = 0.3) +
    facet_grid(. ~ landscape, labeller = labeller(landscape = STRIP1)) +
    scale_fill_manual(values = WHO_COL, name = NULL, breaks = names(WHO_COL)) +
    scale_colour_manual(values = WHO_COL, guide = "none") +
    scale_shape_manual(values = PART_SHAPE, name = NULL, breaks = names(PART_SHAPE)) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    coord_equal() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2)), shape = guide_legend(override.aes = list(fill = "white"))) +
    labs(x = "|r| with temperature (median, interquartile range)", y = "|r| with salinity") +
    theme_ms() + theme(strip.text.x = element_text(size = MS_TICK))
p2b <- p_outcome("cols")
v2 <- cowplot::plot_grid(p2a, p2b, ncol = 1, rel_heights = c(78, 40), align = "v", axis = "lr",
                         labels = c("a", "b"), label_size = 8, label_fontface = "bold")
save_fig("track_v2_plane", v2, 180, 118, "V2: |r| temperature x |r| salinity plane per landscape | outcome",
         c(fig_text(p2a), fig_text(p2b)))
num(D2[, .(landscape, set, part, measure = "plane point (median |r| temperature, salinity)", value = temperature, n = replicates)], "V2", "a")

# ------------------------------------------------------------------ V3 structure within landscapes
BAND_LEV <- c("< 0.2", "0.2–0.6", "≥ 0.6")
OC[, band := factor(fcase(r2_pc1_temp < 0.2, BAND_LEV[1], r2_pc1_temp < 0.6, BAND_LEV[2], default = BAND_LEV[3]), levels = BAND_LEV)]
OC[, sel := factor(fifelse(regime == "equal-S", "equal", "unequal"), levels = c("equal", "unequal"))]
S1 <- OC[set %in% SET_ORDER, .(pct = 100 * mean(w456), n = .N), by = .(landscape, set, x = as.character(band))][
    , kind := "r² of neutral PC1 with temperature"]
S2 <- OC[set %in% SET_ORDER, .(pct = 100 * mean(w456), n = .N), by = .(landscape, set, x = as.character(sel))][
    , kind := "Selection on the two traits"]
S3 <- rbind(S1, S2)
S3[, `:=`(landscape = land_f(landscape), x = factor(x, levels = c(BAND_LEV, "equal", "unequal")),
          kind = factor(kind, levels = c("r² of neutral PC1 with temperature", "Selection on the two traits")),
          set_l = factor(lab_set(set), levels = lab_set(SET_ORDER)))]
NL <- unique(S3[, .(landscape, kind, x, n)])
num(S3[, .(landscape, set, part = paste(kind, x), measure = "% replicates in ranks 4–6", value = pct, n)], "V3", "structure")
p3 <- ggplot(S3, aes(x = x, y = pct, colour = set_l, group = set_l)) +
    geom_hline(yintercept = 50, linetype = "22", colour = MS_INK, linewidth = 0.3) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.4) +
    geom_text(data = NL, aes(x = x, y = -9, label = paste0("n = ", n)), inherit.aes = FALSE, size = MS_LAB * 0.9, colour = "#59606b") +
    facet_grid(landscape ~ kind, scales = "free_x", space = "free_x", switch = "y",
               labeller = labeller(landscape = STRIP)) +
    scale_colour_manual(values = setNames(RULE_COL[SET_ORDER], lab_set(SET_ORDER)), name = NULL) +
    scale_y_continuous(limits = c(-14, 100), breaks = c(0, 50, 100)) +
    labs(x = NULL, y = "% replicates in ranks 4–6\n(median of the three offset methods)") +
    theme_ms() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = MS_TICK))
save_fig("track_v3_structure", p3, 180, 120, "V3: outcome by PC1-temperature band and selection regime, per landscape", fig_text(p3))

# ------------------------------------------------------------------ V4 one scatter (user, 2026-11-07): every set's whole
# set (filled) and its false positives (open) on one plane; colour = set, shape = landscape
SET6 <- ORD[order(order), key]                                        # the pooled order, all six sets
stopifnot(length(SET6) == 6L)
SETCOL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]], "3/3 methods" = MINOU[["navy"]],
            "LFMM" = MINOU[["red"]], "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE, "truth" = MS_REF)
SETLAB <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set", "LFMM" = "LFMM2",
            "RDA" = "RDA", "EMMAX" = "EMMAX", "truth" = "Causal loci")
LAND_SHAPE <- c("SS-Clines" = 21, "SS-Mtn" = 24, "Est-Clines" = 22)
LAND_LAB <- setNames(sprintf("%s (%.2f / %.2f)", STR$landscape, STR$t, STR$s), STR$landscape)[LAND_LEV]
scatter_fig <- function(tcol, scol, stem, variant, measure, xlab, ylab, what_txt) {
    S4 <- TRK[set %in% c(SET6, "truth") & part %in% c("all", "neutral") & n > 0,
              .(temperature = median(get(tcol), na.rm = TRUE), salinity = median(get(scol), na.rm = TRUE), replicates = .N,
                median_snps = as.numeric(median(n))), by = .(landscape, set, part)][!(set == "truth" & part == "neutral")]
    S4[, `:=`(set_f = factor(set, levels = c(SET6, "truth")), land_f = land_f(landscape),
              what = factor(fifelse(part == "all", "whole set", "false positives"), levels = c("whole set", "false positives")))]
    PAIR <- dcast(S4[set != "truth"], landscape + set + set_f + land_f ~ part, value.var = c("temperature", "salinity"))
    num(S4[, .(landscape, set, part = as.character(what), measure = measure,
               value = temperature, value2 = salinity, n = replicates, median_snps)], variant, "scatter")
    print(S4[order(landscape, set_f, what), .(landscape, set, what, temperature, salinity, replicates, median_snps)], digits = 3)
    p4 <- ggplot() +
        geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#9aa1aa", linewidth = 0.3) +
        geom_segment(data = PAIR, aes(x = temperature_all, y = salinity_all, xend = temperature_neutral, yend = salinity_neutral,
                                      colour = set_f), linewidth = 0.35, alpha = 0.55) +
        geom_point(data = S4[what == "false positives"], aes(x = temperature, y = salinity, colour = set_f, shape = land_f, alpha = what),
                   fill = "white", size = 2.1, stroke = 0.7) +
        geom_point(data = S4[what == "whole set"], aes(x = temperature, y = salinity, fill = set_f, shape = land_f, alpha = what),
                   colour = MS_INK, size = 2.1, stroke = 0.3) +
        scale_colour_manual(values = SETCOL, guide = "none") +
        scale_fill_manual(values = SETCOL, labels = SETLAB, name = "Set", breaks = c(SET6, "truth")) +
        scale_shape_manual(values = LAND_SHAPE, labels = LAND_LAB, name = "Landscape\n(PC1 r² with temperature / salinity)") +
        scale_alpha_manual(values = c("whole set" = 1, "false positives" = 1), breaks = c("whole set", "false positives"),
                           name = NULL) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"), expand = expansion(mult = 0.01)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"), expand = expansion(mult = 0.01)) +
        coord_equal() +
        guides(fill  = guide_legend(order = 1, override.aes = list(shape = 21, colour = MS_INK, size = 2)),
               shape = guide_legend(order = 2, override.aes = list(fill = "#9aa1aa", colour = MS_INK, size = 2)),
               alpha = guide_legend(order = 3, override.aes = list(shape = 21, size = 2, fill = c("#5c6470", "white"),
                                                                   colour = "#5c6470", stroke = c(0.3, 0.7)))) +
        labs(x = xlab, y = ylab) +
        theme_ms() + theme(legend.position = "right", legend.box = "vertical", legend.justification = "top",
                           panel.grid.major = element_line(colour = "#eef0f3", linewidth = 0.25))
    save_fig(stem, p4, 180, 112, what_txt, fig_text(p4))
}
if (all(SET6 %in% unique(TRK$set))) {
    scatter_fig("r_temp", "r_sal", "track_v4_scatter", "V4", "median |r| temperature / salinity",
                "|r| of allele frequency with temperature across demes\n(median over SNPs, then over replicates)", "|r| with salinity",
                "V4: median |r| per SNP; whole set (filled) and false positives (open) of the six sets, shape = landscape")
}
if (all(c("pca1_r_temp", "pca1_r_sal") %in% names(TRK))) {
    scatter_fig("pca1_r_temp", "pca1_r_sal", "track_v5_pca", "V5", "|r| of the PC1 of the part's allele frequencies, temperature / salinity",
                "|r| of the first principal component of the SNPs' allele\nfrequencies with temperature (median over replicates)",
                "|r| of that component with salinity",
                "V5: PCA of each set's SNPs; whole set (filled) and false positives (open), shape = landscape")
}

# ------------------------------------------------------------------ V6 (user, 2026-11-07): false positives only, the |r| of
# the PC1 of each set's false positives with the two climate factors; colour = set, shape = landscape; nothing else
if (all(c("pca1_r_temp", "pca1_r_sal") %in% names(TRK)) && all(SET6 %in% unique(TRK$set))) {
    SETV6 <- SET6[SET6 %in% c("RDA", "EMMAX", "LFMM")]               # the three methods only (user, 2026-11-07)
    S6 <- TRK[set %in% SETV6 & part == "neutral" & n > 0,
              .(temperature = median(pca1_r_temp, na.rm = TRUE), salinity = median(pca1_r_sal, na.rm = TRUE),
                replicates = .N, median_snps = as.numeric(median(n))), by = .(landscape, set)]
    S6[, `:=`(set_f = factor(set, levels = SETV6), land_f = land_f(landscape))]
    setorder(S6, set_f)                         # drawn in the pooled order, so LFMM2 (last) sits on top where points coincide
    S6 <- S6[order(set == "LFMM")]
    num(S6[, .(landscape, set, part = "false positives", measure = "|r| of the PC1 of the false positives, temperature / salinity",
               value = temperature, value2 = salinity, n = replicates, median_snps)], "V6", "scatter")
    print(S6[order(landscape, set_f), .(landscape, set, temperature, salinity, replicates, median_snps)], digits = 3)
    LAND_SHAPE6 <- c("SS-Clines" = 21, "SS-Mtn" = 24, "Est-Clines" = 22)
    p6 <- ggplot(S6, aes(x = temperature, y = salinity, fill = set_f, shape = land_f)) +
        geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#9aa1aa", linewidth = 0.3) +
        geom_point(size = 2.6, colour = MS_INK, stroke = 0.35, alpha = 0.92) +
        scale_fill_manual(values = SETCOL, labels = unname(SETLAB[SETV6]), name = "Method", breaks = SETV6) +
        scale_shape_manual(values = LAND_SHAPE6, name = "Landscape") +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1"),
                           expand = expansion(mult = 0.02)) +
        scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.2), labels = c("0", "0.2", "0.4", "0.6"),
                           expand = expansion(mult = 0.02)) +
        guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.4)),
               shape = guide_legend(order = 2, override.aes = list(fill = "#9aa1aa", size = 2.4))) +
        labs(x = "|r| of the false positives' first principal component\nwith temperature (median over replicates)",
             y = "|r| with salinity") +
        theme_ms() + theme(legend.position = "right", legend.box = "vertical", legend.justification = "top",
                           panel.grid.major = element_line(colour = "#eef0f3", linewidth = 0.25))
    save_fig("track_v6_fp_pca", p6, 140, 90, "V6: PC1 of each set's false positives against temperature and salinity, shape = landscape",
             fig_text(p6))

    # ---- the final supplement figure (user, 2026-11-07: V6 plus the diagonal) and its numbers for the registry
    ms_save(file.path(LMS, "figures", "s_fp_pca"), p6, 140, 90)
    FT <- fread(file.path(LMS, "figures", "figures.tsv"))
    FT <- rbind(FT[stem != "s_fp_pca"], data.table(stem = "s_fp_pca", role = "supp", w_mm = 140, h_mm = 90,
                source = "explain/set_tracking.tsv (benchmarks/mvp_land_set_tracking.R), the 1,800 replicates"), fill = TRUE)
    fwrite(FT, file.path(LMS, "figures", "figures.tsv"), sep = "\t")
    Q <- TRK[set %in% SETV6 & part == "neutral" & n > 0, .(
        fp_pc1_r_temp = median(pca1_r_temp, na.rm = TRUE), t25 = quantile(pca1_r_temp, .25, na.rm = TRUE, names = FALSE),
        t75 = quantile(pca1_r_temp, .75, na.rm = TRUE, names = FALSE),
        fp_pc1_r_sal = median(pca1_r_sal, na.rm = TRUE), s25 = quantile(pca1_r_sal, .25, na.rm = TRUE, names = FALSE),
        s75 = quantile(pca1_r_sal, .75, na.rm = TRUE, names = FALSE),
        fp_pc1_pve = median(pca1_pve, na.rm = TRUE), p25 = quantile(pca1_pve, .25, na.rm = TRUE, names = FALSE),
        p75 = quantile(pca1_pve, .75, na.rm = TRUE, names = FALSE),
        fp_n_median = as.numeric(median(n)), n = .N), by = .(landscape, set)]
    W <- dcast(TRK[set %in% SETV6 & part %in% c("causal temperature", "causal salinity", "causal both")],
               landscape + seed + set ~ part, value.var = "n")
    W <- W[, .(temp = `causal temperature` + `causal both`, sal = `causal salinity` + `causal both`), by = .(landscape, seed, set)]
    CT <- W[, .(temp_median = as.numeric(median(temp)), temp_mean = mean(temp), sal_median = as.numeric(median(sal)),
                sal_mean = mean(sal), n = .N), by = .(landscape, set)]
    src <- "landscape_ms/explain/set_tracking.tsv (benchmarks/mvp_land_set_tracking.R)"
    row <- function(D, metric, key, value, lo = NA_real_, hi = NA_real_, figure = "s_fp_pca")
        D[, .(landscape, figure = figure, rule = set, set = unname(SETLAB[set]), stratum = "pooled", engine = NA_character_, metric = metric,
              metric_key = key, value = get(value), lo = if (is.na(lo)) NA_real_ else get(lo), hi = if (is.na(hi)) NA_real_ else get(hi),
              expected = NA_real_, n, p_perm = NA_real_, source = src)]
    NB <- rbind(
        row(Q, "|r| of the false positives' PC1 with temperature, median over replicates [interquartile range]", "fp_pc1_r_temp", "fp_pc1_r_temp", "t25", "t75"),
        row(Q, "|r| of the false positives' PC1 with salinity, median over replicates [interquartile range]", "fp_pc1_r_sal", "fp_pc1_r_sal", "s25", "s75"),
        row(Q, "variance of the false positives' allele frequencies on their PC1, median [interquartile range]", "fp_pc1_pve", "fp_pc1_pve", "p25", "p75"),
        row(Q, "false positives per replicate, median", "fp_n_median", "fp_n_median"),
        row(CT, "causal loci of temperature identified per replicate, median", "causal_temp_median", "temp_median", figure = "text: causal loci by trait"),
        row(CT, "causal loci of temperature identified per replicate, mean", "causal_temp_mean", "temp_mean", figure = "text: causal loci by trait"),
        row(CT, "causal loci of salinity identified per replicate, median", "causal_sal_median", "sal_median", figure = "text: causal loci by trait"),
        row(CT, "causal loci of salinity identified per replicate, mean", "causal_sal_mean", "sal_mean", figure = "text: causal loci by trait"))
    fwrite(NB[order(match(landscape, LAND_LEV), metric_key, match(rule, SETV6))], file.path(LMS, "tables", "fp_pca_numbers.tsv"), sep = "\t")
    message("final: figures/s_fp_pca.{svg,png} (140 x 90 mm), tables/fp_pca_numbers.tsv (", nrow(NB), " rows)")
}

# ------------------------------------------------------------------ write
fwrite(rbindlist(FIGS), file.path(OUT, "figures.tsv"), sep = "\t")
fwrite(rbindlist(NUM, fill = TRUE), file.path(EXP, "tracking_fig_numbers.tsv"), sep = "\t")
message("figures: ", length(FIGS), "; numbers: ", sum(vapply(NUM, nrow, 0L)), " rows")
