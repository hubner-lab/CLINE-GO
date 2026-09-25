#!/usr/bin/env Rscript
# =============================================================================
# mvp_ms_sim_figures.R -- manuscript figures for the simulation section (SS-Clines, 600
# replicates), chosen 2026-09-25 from the journal-16 gallery (16g, mvp_j16_gallery.R).
#
# PLOTTING ONLY. Every number comes from the tables mvp_j16_gallery.R already wrote to
# figures_ssclines_j16_gallery/ (compute once, reuse); medians are re-derived from the
# per-seed tables only as a regression tie against the gallery's stats tables.
#
#   main  A1_detection_plane       TP-vs-FP plane: contours + medians, own x-axis per panel
#         C5_offset_reach          % replicates with delta >= 0, rule x genic level
#         C2_offset_distribution   raincloud of delta (main since 2026-09-25)
#   supp  D5_detection_cells       A1 view, genic level x (pleiotropy x selection regime)
#         D5_offset_cells          violin of delta + median + 95% bootstrap CI, same 12 panels
#         D6_detection_demography  A1 view, the five demography blocks
#
# LOOK: benchmarks/ms_style.R -- the shared manuscript style (ltc "minou" palette, 7 pt, print
# size, editable SVG via ms_save()). Rule colours: combined rules green/blue (2/3 teal, 1/3 sage,
# 3/3 navy), single methods red/orange (LFMM red, RDA amber, EMMAX orange), references grey.
# Every detection plane labels its medians by rule name; every offset distribution is the same
# raincloud (C2 view).
#
# Also writes captions.md, index.html (served via work/journal/sim_figures on port 8099) and
# sim_figures_ms.zip (every SVG + PNG + captions.md) for the page's "Download all" button.
#
#   GALLERY_DIR  default benchmarks/mvp_eval/figures_ssclines_j16_gallery
#   FIG_OUT      default benchmarks/mvp_eval/figures_ssclines_ms
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT via benchmarks/mvp_arm.R (mandatory, 600)
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggdist); library(ggrepel)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
GAL  <- Sys.getenv("GALLERY_DIR", file.path(EVAL, "figures_ssclines_j16_gallery"))
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "figures_ssclines_ms"))
JITTER_SEED <- 16L
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ vocabulary
# Labels and rule order as mvp_j16_gallery.R:57-70.
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c("N variable, m variable", "N cline north-south", "N equal, m constant",
               "N cline centre-edge", "N equal, m breaks")
RULES <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
# Hue = group (combined rules green/blue, single methods red/orange); ms_style.R palette.
RULE_COL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]],
              "3/3 methods" = MINOU[["navy"]], "LFMM" = MINOU[["red"]],
              "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE)
DELTA_LAB <- "Δ accuracy (panel − causal loci)"

TEXT_TONES <- c("#FFFFFF" = "#FFFFFF", setNames(MS_INK, MS_INK))
ZERO <- geom_vline(xintercept = 0, linetype = "dashed", colour = MS_REF, linewidth = 0.35)

FIGS <- list()                                  # stem -> list(role, w, h, source)
save_ms <- function(stem, p, w, h, role, source) {
    ms_save(file.path(OUT, stem), p, w, h)
    FIGS[[stem]] <<- list(role = role, w = w, h = h, source = source)
    message(sprintf("  %-26s %3.0f x %3.0f mm", stem, w, h))
}

# ------------------------------------------------------------------ inputs
man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
rd <- function(stem, ...) fread(file.path(GAL, paste0(stem, ".tsv")), ...)
fac <- function(D) {
    D[, arch := factor(arch, levels = ARCH_LABELS)]
    if ("regime" %in% names(D)) D[, regime := factor(regime, levels = c("equal-S", "unequal-S"))]
    if ("pleio"  %in% names(D)) D[, pleio  := factor(pleio,  levels = c("no pleiotropy", "pleiotropy"))]
    if ("block"  %in% names(D)) D[, block  := factor(block,  levels = BLOCK_LAB)]
    D[, label := factor(label, levels = RULES)]
    stopifnot(!anyNA(D$arch), !anyNA(D$label))
    D[]
}
DET <- fac(rd("detection_per_seed", colClasses = c(seed = "character"))[label %in% RULES])
OFF <- fac(rd("offset_delta_per_seed", colClasses = c(seed = "character"))[label %in% RULES])
DS  <- rd("detection_stats")[label %in% RULES]
CS  <- rd("C_delta_stats")[label %in% RULES]
D5S <- fac(rd("D5_offset_stats")[label %in% RULES])
stopifnot(!anyNA(DET$block), !anyNA(OFF$regime), !anyNA(D5S$pleio))
# One wrapped panel per design cell (gene level first, so ncol = 4 gives rows = gene level and
# columns = pleiotropy x regime). facet_wrap rather than facet_grid because only facet_wrap can
# give every panel its own x-axis without ggh4x, and both D5 figures use it so they pair up.
CELL_LEVELS <- CJ(arch = ARCH_LABELS, pleio = c("no pleiotropy", "pleiotropy"),
                  regime = c("equal-S", "unequal-S"), sorted = FALSE)[
    , paste0(arch, "\n", pleio, ", ", regime)]
for (D in list(DET, OFF, D5S))
    D[, panel := factor(paste0(arch, "\n", pleio, ", ", regime), levels = CELL_LEVELS)]
stopifnot(!anyNA(DET$panel), !anyNA(OFF$panel), !anyNA(D5S$panel))

# Corpus guards: every seed is a primary SS-Clines seed; detection keeps empty panels (all
# 600 per rule); design cells hold 200 / 50 / 120 replicates.
NREP <- nrow(PRIM)
stopifnot(all(DET$seed %in% PRIM$seed), all(OFF$seed %in% PRIM$seed),
          uniqueN(DET$seed) == NREP, uniqueN(OFF$seed) == NREP,
          DET[, .N, by = label][, all(N == nrow(PRIM))],
          DET[, .N, by = .(label, arch)][, all(N == nrow(PRIM) / 3)],
          DET[, .N, by = .(label, arch, pleio, regime)][, all(N == nrow(PRIM) / 12)],
          DET[, .N, by = .(label, block)][, all(N == nrow(PRIM) / 5)])

# Regression tie: medians re-derived from the per-seed tables equal the gallery's stats tables.
tie <- function(a, b, keys, what) {
    x <- merge(a, b, by = keys, all = TRUE)
    stopifnot(nrow(x) == nrow(a), !anyNA(x$m), !anyNA(x$median),
              max(abs(x$m - x$median)) < 1e-9, all(x$n.x == x$n.y))
    message(sprintf("tie %-24s %3d medians, max |diff| %.1e", what, nrow(x), max(abs(x$m - x$median))))
}
tie(OFF[, .(m = median(delta), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    CS[arch != "pooled", .(label, arch, median, n)], c("label", "arch"), "C offset by arch")
tie(OFF[, .(m = median(delta), n = .N), by = .(label = as.character(label))],
    CS[arch == "pooled", .(label, median, n)], "label", "C offset pooled")
tie(OFF[, .(m = median(delta), n = .N), by = .(label, arch, pleio, regime)],
    D5S[, .(label, arch, pleio, regime, median, n)], c("label", "arch", "pleio", "regime"), "D5 offset")
tie(DET[, .(m = median(TP), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    DS[arch != "pooled", .(label, arch, median = TP_median, n = n_seeds)], c("label", "arch"), "detection TP")
tie(DET[, .(m = median(FP), n = .N), by = .(label = as.character(label), arch = as.character(arch))],
    DS[arch != "pooled", .(label, arch, median = FP_median, n = n_seeds)], c("label", "arch"), "detection FP")

# =============================================================================
# Detection plane (A1, D5-det, D6-det)
# =============================================================================
# Explicit sqrt-axis breaks: the defaults (x 0/40/80/120, y 0/100/200) leave the whole
# 0-40 stretch of a sqrt axis -- where most small panels sit -- without a tick.
BR_FP <- c(0, 5, 20, 50, 100)
BR_TP <- c(0, 10, 50, 100, 200)
stopifnot(max(DET$FP) >= max(BR_FP), max(DET$TP) >= max(BR_TP))
# A 2-D density needs spread in both sqrt dimensions (gallery :225); rules failing it in a
# facet are drawn by their median only, and the caption names them.
contour_ok <- function(D, by) D[, .(ok = .N >= 10L && IQR(sqrt(FP)) > 0 && IQR(sqrt(TP)) > 0), by = by]

plane_ms <- function(D, facet) {
    by <- c("label", facet$vars)
    M  <- D[, .(TP = median(TP), FP = median(FP)), by = by]
    DC <- merge(D, contour_ok(D, by)[ok == TRUE], by = by)
    p <- ggplot(D, aes(FP, TP)) +
        geom_density_2d(data = DC, aes(colour = label), contour_var = "ndensity",
                        breaks = c(0.25, 0.6), linewidth = 0.35) +
        geom_point(data = M, aes(fill = label), shape = 21, size = 1.7, colour = "black", stroke = 0.3) +
        # Median labels are the rule names on solid tags in the rule colour; the text is white or
        # dark, whichever contrasts more with that colour (ms_text_on). Tinted tags with dark
        # text (2026-09-25 draft) were barely visible. The text tone rides on the same colour
        # scale as the contours, under its own two keys.
        geom_label_repel(data = M[, tone := ms_text_on(RULE_COL[as.character(label)])],
                         aes(label = label, fill = label, colour = tone), size = MS_LAB,
                         fontface = "bold", label.size = 0, label.padding = 0.12, box.padding = 0.3,
                         min.segment.length = 0, segment.size = 0.25, segment.colour = MS_INK,
                         seed = JITTER_SEED, show.legend = FALSE, max.overlaps = Inf)
    p + facet$layer +
        scale_x_sqrt(breaks = BR_FP) + scale_y_sqrt(breaks = BR_TP) +
        scale_colour_manual(values = c(RULE_COL, TEXT_TONES), guide = "none") +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = "False positives", y = "True positives") +
        theme_ms()
}
# Caption clause naming the rules drawn by their median only; empty when every rule has contours.
median_only <- function(D, vars) {
    x <- contour_ok(D, c("label", vars))[ok == FALSE, .N, by = label][order(match(label, RULES))]
    if (!nrow(x)) return("")
    paste0(" Shown by the median only (too few distinct counts for a density): ",
           paste(sprintf("%s (%d of %d panels)", x$label, x$N,
                         nrow(unique(D[, vars, with = FALSE]))), collapse = ", "), ".")
}

message("=== detection")
# Every detection plane gives each panel its own false-positive axis (scales = "free_x").
pA1 <- plane_ms(DET, list(vars = "arch", layer = facet_grid(cols = vars(arch), scales = "free_x")))
for (pp in ggplot_build(pA1)$layout$panel_params) {
    bx <- pp$x$get_labels(); by <- pp$y$get_labels()
    message("A1 panel axis labels  x: ", paste(bx, collapse = " "), "   y: ", paste(by, collapse = " "))
    stopifnot(!anyNA(bx), "20" %in% bx, length(bx) >= 4L, !anyNA(by), length(by) == length(BR_TP))
}
save_ms("A1_detection_plane", pA1, 180, 62, "main", "detection_per_seed.tsv")

save_ms("D5_detection_cells",
        plane_ms(DET, list(vars = "panel", layer = facet_wrap(vars(panel), ncol = 4, scales = "free_x"))),
        180, 165, "supp", "detection_per_seed.tsv")

save_ms("D6_detection_demography",
        plane_ms(DET, list(vars = "block", layer = facet_wrap(vars(block), nrow = 2, scales = "free_x"))),
        180, 110, "supp", "detection_per_seed.tsv")

# =============================================================================
# Offset (C5, C2, D5-offset)
# =============================================================================
message("=== offset")
ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\npolygenic",
               "highly polygenic" = "highly\npolygenic", "pooled" = "pooled")
CS5 <- copy(CS)[, arch := factor(ARCH_WRAP[arch], levels = ARCH_WRAP)]
stopifnot(!anyNA(CS5$arch), nrow(CS5) == length(RULES) * 4L)
# Quality scale red -> orange -> amber -> green (user, 2026-09-25: higher share = better =
# green), saturating at 50%, parity with the causal loci. Every observed share is 15-44%, so a
# scale centred on 50 (first draft) drew the whole grid red. Tile text is white or dark per tile,
# from the same ramp.
# Full green is pinned at 35% (user, 2026-09-25) and held to the 50% limit; red, orange and amber
# spread evenly below it.
C5_LIM      <- c(0, 50)
C5_GREEN_AT <- 35
C5_COLS <- c(MS_QUALITY, MS_QUALITY[length(MS_QUALITY)])
C5_VALS <- c(seq(0, C5_GREEN_AT, length.out = length(MS_QUALITY)), max(C5_LIM)) / max(C5_LIM)
C5_PAL  <- function(x) scales::gradient_n_pal(C5_COLS, C5_VALS)(scales::rescale(x, from = C5_LIM))
CS5[, tone := ms_text_on(C5_PAL(pmin(reach_pct, max(C5_LIM))))]
C5_OVER <- CS5[reach_pct > max(C5_LIM), .N]
pC5 <- ggplot(CS5, aes(arch, label)) +
    geom_tile(aes(fill = reach_pct), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f", reach_pct), colour = tone), size = MS_LAB) +
    scale_y_discrete(limits = rev(RULES)) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = C5_COLS, values = C5_VALS, limits = C5_LIM,
                         breaks = c(0, 25, 50), labels = c("0", "25", "≥ 50"),
                         oob = scales::squish,
                         name = "Δ ≥ 0 (%)") +
    guides(fill = guide_colourbar(theme = theme(legend.key.width = unit(18, "mm"),
                                                legend.key.height = unit(2, "mm")))) +
    labs(x = NULL, y = NULL) +
    theme_ms() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                       legend.position = "bottom")
save_ms("C5_offset_reach", pC5, 66, 60, "main", "C_delta_stats.tsv")

# One raincloud for every offset distribution (C2 in the main figure, D5 in the supplement):
# density slab above, box (median, IQR) on the line, replicates as dots below.
raincloud_ms <- function(D, facet_layer) {
    RC <- copy(D)[, ypos := as.numeric(factor(label, levels = rev(RULES)))]
    ggplot(RC, aes(delta, ypos)) + ZERO +
        stat_halfeye(aes(fill = label, group = label), orientation = "horizontal", adjust = 0.8,
                     height = 0.55, justification = -0.2, .width = 0, point_colour = NA,
                     slab_alpha = 0.75) +
        geom_boxplot(aes(group = label), orientation = "y", width = 0.14, outlier.shape = NA,
                     fill = "white", colour = MS_INK, linewidth = 0.25) +
        geom_point(aes(y = ypos - 0.25, colour = label),
                   position = position_jitter(height = 0.08, width = 0, seed = JITTER_SEED),
                   size = 0.6, alpha = 0.6, stroke = 0) +
        facet_layer +
        scale_y_continuous(breaks = seq_along(RULES), labels = rev(RULES)) +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        scale_colour_manual(values = RULE_COL, guide = "none") +
        labs(x = DELTA_LAB, y = NULL) + theme_ms()
}
save_ms("C2_offset_distribution",
        raincloud_ms(OFF, facet_grid(cols = vars(arch), labeller = labeller(arch = ARCH_WRAP))),
        114, 60, "main", "offset_delta_per_seed.tsv")

# D5 offset: the C2 raincloud per design cell (bars, then violins, were tried and replaced on
# 2026-09-25). The x-axis stays SHARED: this figure compares cells against one zero line.
# D5_offset_stats.tsv is not plotted; it is only the regression tie above.
save_ms("D5_offset_cells", raincloud_ms(OFF, facet_wrap(vars(panel), ncol = 4)),
        180, 165, "supp", "offset_delta_per_seed.tsv")

# =============================================================================
# Captions (all explanation lives here, not on the figures)
# =============================================================================
nP <- function(rule) CS[label == rule & arch == "pooled", n]
DETECTION_DEF <- paste0(
    "True positives are selected loci that are causal or linked-neutral (linkage groups 1–10); ",
    "false positives are selected background-neutral loci (linkage groups 11–20). ",
    "Both axes use a square-root scale; each panel has its own false-positive axis. ",
    "Contours enclose 25% and 60% of the peak replicate density of each rule; labelled points are the median replicate. ",
    "Operating points are fixed for all replicates: LFMM top 100, RDA top 100, EMMAX p < 10⁻⁴; ",
    "the 1/3, 2/3 and 3/3 rules keep a locus called by at least one, two or all three methods ",
    "within ±5 kb. The causal loci are not shown, because by construction they carry no linked hits.")
OFFSET_DEF <- paste0(
    "Accuracy is −Kendall τ between the predicted genetic offset and the fitness of the 100 ",
    "source populations in a common garden: median over the 100 landscape gardens, then over three ",
    "offset engines (Gradient Forest, LFMM2 geometric offset, RDA without structure correction). ",
    "Δ = accuracy of the panel − accuracy of the causal loci, paired within a replicate. ",
    sprintf("Panels with fewer than 3 SNPs are not scored (3/3 methods: n = %d; EMMAX: n = %d of %d).",
            nP("3/3 methods"), nP("EMMAX"), NREP))
CAP <- list(
    A1_detection_plane = paste0(
        sprintf("Detection of adaptive loci by six marker-selection rules, by genetic architecture (n = %d replicates per level). ", NREP / 3),
        DETECTION_DEF, median_only(DET, "arch")),
    C5_offset_reach = paste0(
        sprintf("Share of replicates in which the genetic offset from the rule's marker panel is at least as accurate as from the causal loci (Δ ≥ 0), by genetic architecture (n = %d per level) and pooled (n = %d). ", NREP / 3, NREP),
        OFFSET_DEF, " At n = 200, shares of 43–57% lie inside the exact binomial 95% interval of 50%. ",
        sprintf("Colour runs from red (0%%) through orange and amber to green, reached at %d%% and held up to 50%% (parity with the causal loci)", C5_GREEN_AT),
        if (C5_OVER > 0) sprintf("; %d cells exceed 50%%.", C5_OVER) else "; no cell reaches 50%."),
    C2_offset_distribution = paste0(
        sprintf("Distribution of Δ accuracy per replicate, by genetic architecture (n = %d per level). ", NREP / 3),
        "Half-violins: density; boxes: median and interquartile range; dots: replicates. ",
        "Dashed line: accuracy equal to that of the causal loci. ", OFFSET_DEF),
    D5_detection_cells = paste0(
        sprintf("Detection plane (as in the main figure) for every design cell: genetic architecture (rows) × pleiotropy × selection regime (columns), n = %d replicates per cell. ", NREP / 12),
        "equal-S: both environmental axes under equal selection; unequal-S: selection on one axis 8-fold weaker. ",
        DETECTION_DEF, median_only(DET, "panel")),
    D5_offset_cells = paste0(
        sprintf("Offset accuracy relative to the causal loci for every design cell, drawn as main panel c (layout as the detection supplement, n = %d replicates per cell). ", NREP / 12),
        "Half-violins: density; boxes: median and interquartile range; dots: replicates; ",
        "dashed line: accuracy equal to that of the causal loci. ", OFFSET_DEF),
    D6_detection_demography = paste0(
        sprintf("Detection plane (as in the main figure) for the five demographic scenarios of the simulation deposit (n = %d replicates each; migration rate 0.03 in all). ", NREP / 5),
        DETECTION_DEF, median_only(DET, "block")))
stopifnot(setequal(names(CAP), names(FIGS)))

ORDER <- c("A1_detection_plane", "C5_offset_reach", "C2_offset_distribution",
           "D5_detection_cells", "D5_offset_cells", "D6_detection_demography")
TAG <- c(A1_detection_plane = "Main a", C5_offset_reach = "Main b",
         C2_offset_distribution = "Main c", D5_detection_cells = "Supplementary 1",
         D5_offset_cells = "Supplementary 2", D6_detection_demography = "Supplementary 3")
TITLE <- c(A1_detection_plane = "Detection: true vs false positives",
           C5_offset_reach = "Offset: replicates reaching the causal loci",
           C2_offset_distribution = "Offset: distribution of Δ accuracy",
           D5_detection_cells = "Detection by design cell",
           D5_offset_cells = "Offset by design cell",
           D6_detection_demography = "Detection by demography")
md <- unlist(lapply(ORDER, function(s) c(
    sprintf("## %s — %s", TAG[[s]], TITLE[[s]]), "",
    sprintf("`%s.svg` / `.png` · %.0f × %.0f mm · source: `%s`", s, FIGS[[s]]$w, FIGS[[s]]$h,
            FIGS[[s]]$source), "", CAP[[s]], "")))
writeLines(c("# Simulation figures — draft captions",
             "", sprintf("SS-Clines corpus, %d replicates. Arial 7 pt, drawn at print size.", NREP), "", md),
           file.path(OUT, "captions.md"))

# =============================================================================
# Download-all zip + index.html
# =============================================================================
files <- c(paste0(ORDER, ".svg"), paste0(ORDER, ".png"), "captions.md")
ZIP <- "sim_figures_ms.zip"
old <- setwd(OUT)
unlink(ZIP)
stopifnot(utils::zip(ZIP, files, flags = "-9Xq") == 0L)
setwd(old)

esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
                     gsub(">", "&gt;", x, fixed = TRUE) }
kb <- function(f) { b <- file.size(file.path(OUT, f)); if (b > 1e6) sprintf("%.1f MB", b / 1e6) else sprintf("%.0f kB", b / 1e3) }
card <- function(s, cls = "") sprintf(paste0(
    '<figure class="card %s" id="%s">',
    '<figcaption class="head"><span class="tag">%s</span><span class="ttl">%s</span>',
    '<span class="dim">%.0f × %.0f mm</span></figcaption>',
    '<div class="img"><img src="%s.svg" alt="%s" style="--w:%.0f" loading="lazy"></div>',
    '<p class="cap">%s</p>',
    '<div class="dl"><a class="btn" href="%s.svg" download>SVG · %s</a>',
    '<a class="btn" href="%s.png" download>PNG 600 dpi · %s</a></div></figure>'),
    cls, s, esc(TAG[[s]]), esc(TITLE[[s]]), FIGS[[s]]$w, FIGS[[s]]$h, s, esc(TITLE[[s]]),
    FIGS[[s]]$w, esc(CAP[[s]]), s, kb(paste0(s, ".svg")), s, kb(paste0(s, ".png")))

html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
'<meta name="viewport" content="width=device-width, initial-scale=1">',
'<title>Simulation figures</title><style>',
':root{--bg:#f5f6f8;--card:#fff;--fg:#1f2933;--muted:#5f6b7a;--line:#dde2e8;--accent:#00798c;--accent-fg:#fff}',
'@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15181d;--card:#1f242b;--fg:#e6e9ee;--muted:#9aa5b3;--line:#323a44;--accent:#3fa58a;--accent-fg:#0d1512}}',
'*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",Arial,sans-serif}',
'header{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--line);padding:14px 24px;display:flex;flex-wrap:wrap;gap:12px;align-items:center;justify-content:space-between}',
'header h1{font-size:18px;margin:0}header p{margin:2px 0 0;color:var(--muted);font-size:13px}',
'main{max-width:1500px;margin:0 auto;padding:20px 24px 60px}',
'h2{font-size:15px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:28px 0 12px}',
'.note{color:var(--muted);font-size:13px;margin:-6px 0 14px}',
'.grid{display:grid;grid-template-columns:66fr 114fr;gap:16px}.full{grid-column:1/-1}',
'.stack{display:grid;gap:16px}',
'.card{margin:0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px;display:flex;flex-direction:column;min-width:0}',
'.head{display:flex;flex-wrap:wrap;gap:8px;align-items:baseline;margin-bottom:10px}',
'.tag{background:var(--accent);color:var(--accent-fg);font-size:12px;font-weight:600;padding:2px 8px;border-radius:999px}',
'.ttl{font-weight:600}.dim{color:var(--muted);font-size:12px;margin-left:auto}',
'.img{background:#fff;border-radius:6px;padding:6px;display:flex;justify-content:center;overflow-x:auto}',
'.img img{width:100%;height:auto;display:block}',
'body.print .img img{width:calc(var(--w) * 1mm);max-width:none}',
'.cap{color:var(--muted);font-size:13px;margin:10px 0 12px;flex:1}',
'.dl{display:flex;flex-wrap:wrap;gap:8px}',
'.btn{display:inline-block;cursor:pointer;font:inherit;text-decoration:none;font-size:13px;padding:6px 12px;border-radius:7px;border:1px solid var(--line);color:var(--fg);background:var(--bg)}',
'.btn:hover{border-color:var(--accent)}.btn.primary{background:var(--accent);color:var(--accent-fg);border-color:var(--accent);font-weight:600}',
'@media (max-width:900px){.grid{grid-template-columns:1fr}main,header{padding-left:16px;padding-right:16px}body:not(.print) .img img{width:calc(100% * var(--w) / 180)}}',
'</style></head><body>',
sprintf('<header><div><h1>Simulation figures — SS-Clines, %d replicates</h1><p>Editable SVG (Arial 7 pt, print size) and PNG 600 dpi; captions are drafts. Built %s.</p></div>', NREP, format(Sys.time(), "%Y-%m-%d %H:%M")),
sprintf('<div class="dl"><button class="btn" id="ps" type="button" aria-pressed="false">Show print size</button><a class="btn" href="captions.md" download>captions.md</a><a class="btn primary" href="%s" download>Download all (.zip, %s)</a></div></header>', ZIP, kb(ZIP)),
'<main><h2>Main figure — proposed layout</h2>',
'<p class="note">Row 1: detection at full width. Row 2: offset reach (b) beside the offset distribution (c), 66 + 114 = 180 mm.</p>',
'<section class="grid">', card("A1_detection_plane", "full"), card("C5_offset_reach"),
card("C2_offset_distribution"), '</section>',
'<h2>Supplementary</h2><section class="stack">', card("D5_detection_cells"), card("D5_offset_cells"),
card("D6_detection_demography"), '</section></main>',
# Every image is shown at one common scale (its width in mm / 180 of the row), so relative panel
# sizes stay true; "print size" shows each at its physical width (CSS mm) to judge 7 pt text.
'<script>document.getElementById("ps").addEventListener("click",function(){var on=document.body.classList.toggle("print");this.setAttribute("aria-pressed",on);this.textContent=on?"Fit to page":"Show print size";});</script>',
'</body></html>')
writeLines(html, file.path(OUT, "index.html"))

message(sprintf("OK: %d figures (SVG + PNG), captions.md, index.html, %s -> %s", length(FIGS), ZIP, OUT))
