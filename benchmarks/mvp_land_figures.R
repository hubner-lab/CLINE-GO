#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_figures.R -- the simulation figures of the manuscript on all three landscapes (Job D of the
# landscapes pack). Plotting only: every number comes from the landscapes' scored tables or from
# mvp_land_tables.R (t_ranks.tsv, t_offset_accuracy.tsv).
#
# Style and code as the published figures: the functions below (plane_ms, bags_of, sqrt_even_breaks,
# check_ticks, rank_dist, heat, vocab_check, save_ms, y_title_fits) are copied from
# work/closeout_lab/scratch/new/sim_figures_ld.R and rda_correction_ld.R, the copies that drew the published
# figures. One set order everywhere: fewest replicates in ranks 4-6 over the 1,800, median of the three
# offset methods (tables/set_order.tsv), also in the per-landscape figures.
#
#   FIG_MODE=ab     SS-Clines alone in the published order, written under the published stems into FIG_OUT:
#                   each SVG must equal work/closeout_lab/results/figures/<stem>.svg byte for byte. The
#                   proof that this code draws what was published; a difference is a bug here.
#   FIG_MODE=main   Fig2_v1, v2, v3, v4, Fig2_per_landscape_{a,b,c}
#   FIG_MODE=supp   s_rda_correction, s_aucpr_methods, s_cutoffs (a row per landscape);
#                   s_detection_cells_{parts_a,b,c, rows}; s_detection_demography_parts_*,
#                   s_rank_architecture_parts_*, s_rank_strata_parts_* (a = SS-Clines, b = SS-Mtn, c = Est-Clines)
#   FIG_MODE=merge  (Data S1 Job F, 2026-11-08; docs/manuscript/lab/data_s1_landscapes_LAB_PROMPT.md) the merged
#                   figures of the story restructure: s_precision_recall (a s_aucpr_methods above b s_rda_correction),
#                   s_detection_design_parts_{a,b,c} (a design cells above b demographic scenarios, landscape title),
#                   s_rank_strata_parts_{a,b,c} (% of replicates in the top three, landscape title). Each member is
#                   also saved alone under its published stem into CHECK_OUT (default <FIG_OUT>/../checks/members), so
#                   its SVG can be compared with the published one; the strata cells go to CHECK_OUT/../strata_top3.tsv.
# Writes PNG (600 dpi) + SVG and appends to figures.tsv (stem, role, w_mm, h_mm, source).
# Usage: FIG_MODE=<mode> FIG_OUT=<dir> TAB_DIR=<tables> Rscript benchmarks/mvp_land_figures.R
# =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggdist); library(ggrepel)
})
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- Sys.getenv("FIG_OUT", file.path(EVAL, "landscape_ms/figures"))
TAB  <- Sys.getenv("TAB_DIR", file.path(EVAL, "landscape_ms/tables"))
MODE <- Sys.getenv("FIG_MODE", "main")
stopifnot(MODE %in% c("ab", "main", "supp", "merge"))
CHECK <- Sys.getenv("CHECK_OUT", file.path(OUT, "..", "checks", "members"))
JITTER_SEED <- 16L
source(file.path(ROOT, "benchmarks/ms_style.R"))
source(file.path(ROOT, "benchmarks/mvp_rank_ties.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------ vocabulary (closeout sim_figures_ld.R)
LANDS <- data.table(
    land  = c("SS-Clines", "SS-Mtn", "Est-Clines"),
    gal   = c("figures_ld_gallery14", "figures_land_ssmtn_gallery", "figures_land_estclines_gallery"),
    rank  = c("figures_ld_offset14/rank", "figures_land_ssmtn_offset15_ssmtn/rank",
              "figures_land_estclines_offset15_estclines/rank"),
    rem   = c("remeasure600", "remeasure600_ssmtn", "remeasure600_estclines"))
if (MODE == "ab") LANDS <- LANDS[land == "SS-Clines"]
LAND_LEV <- c("SS-Clines", "SS-Mtn", "Est-Clines")
PART_OF <- c("SS-Clines" = "a", "SS-Mtn" = "b", "Est-Clines" = "c")
ARCH_LABELS <- c("oligogenic", "moderately polygenic", "highly polygenic")
BLOCK_LAB <- c("N variable, m variable", "N cline north-south", "N equal, m constant",
               "N cline centre-edge", "N equal, m breaks")
BLOCK_SHOW <- setNames(sub("centre-edge", "center-edge", BLOCK_LAB, fixed = TRUE), BLOCK_LAB)
RULES <- c("2/3 methods", "1/3 methods", "3/3 methods", "LFMM", "RDA", "EMMAX")
RULE_COL <- c("2/3 methods" = MINOU[["teal"]], "1/3 methods" = MINOU[["sage"]],
              "3/3 methods" = MINOU[["navy"]], "LFMM" = MINOU[["red"]],
              "RDA" = MINOU[["amber"]], "EMMAX" = MS_ORANGE)
ENG_LAB <- c(gradientForest = "Gradient Forest", LFMM2 = "Geometric offset", RDA = "RDA-based offset",
             `median of 3` = "Median of the three")
N_PANEL <- length(RULES)
SET_LAB <- c("2/3 methods" = "2/3 set", "1/3 methods" = "1/3 set", "3/3 methods" = "3/3 set",
             "LFMM" = "LFMM2", "RDA" = "RDA", "EMMAX" = "EMMAX")
set_lab <- function(x) unname(SET_LAB[as.character(x)])
KEY_OF <- setNames(names(SET_LAB), unname(SET_LAB))
TP_LAB <- "True positives (causal or linked SNPs)"
FP_LAB <- "False positives (neutral SNPs)"
ARCH_WRAP <- c("oligogenic" = "oligogenic", "moderately polygenic" = "moderately\npolygenic",
               "highly polygenic" = "highly\npolygenic")
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
TEXT_TONES <- c("#FFFFFF" = "#FFFFFF", setNames(MS_INK, MS_INK))
FIGS <- list()
fig_text <- function(p) { grDevices::pdf(NULL); on.exit(grDevices::dev.off()); grob_text(ggplotGrob(p)) }
save_ms <- function(stem, p, w, h, role, source) {
    n_txt <- vocab_check(fig_text(p), paste0(stem, " (figure text)"))
    ms_save(file.path(OUT, stem), p, w, h)
    FIGS[[stem]] <<- list(role = role, w = w, h = h, source = source)
    message(sprintf("  %-34s %3.0f x %3.0f mm  (%d text strings checked)", stem, w, h, n_txt))
}

# ------------------------------------------------------------------ inputs
rd <- function(...) fread(file.path(...), colClasses = c(seed = "character"))
fac <- function(D) {
    D[, arch := factor(arch, levels = ARCH_LABELS)]
    if ("regime" %in% names(D)) D[, regime := factor(regime, levels = c("equal-S", "unequal-S"))]
    if ("pleio"  %in% names(D)) D[, pleio  := factor(pleio,  levels = c("no pleiotropy", "pleiotropy"))]
    if ("block"  %in% names(D)) D[, block  := factor(block,  levels = BLOCK_LAB)]
    D[, label := factor(label, levels = RULES)]
    stopifnot(!anyNA(D$arch), !anyNA(D$label))
    D[]
}
CELL_LEVELS <- CJ(arch = ARCH_LABELS, pleio = c("no pleiotropy", "pleiotropy"),
                  regime = c("equal-S", "unequal-S"), sorted = FALSE)[, paste0(arch, "\n", pleio, ", ", regime)]
DET <- rbindlist(lapply(seq_len(nrow(LANDS)), function(i) {
    d <- fac(rd(EVAL, LANDS$gal[i], "detection_per_seed.tsv")[label %in% RULES])
    d[, landscape := factor(LANDS$land[i], levels = LAND_LEV)]
}))
DET[, panel := factor(paste0(arch, "\n", pleio, ", ", regime), levels = CELL_LEVELS)]
stopifnot(!anyNA(DET$panel), !anyNA(DET$block), !anyNA(DET$regime), !anyNA(DET$pleio))
RL <- rbindlist(lapply(seq_len(nrow(LANDS)), function(i) rd(EVAL, LANDS$rank[i], "rank_long.tsv")[, landscape := factor(LANDS$land[i], levels = LAND_LEV)]))
RL[, key := fifelse(empty, -Inf, accuracy)]
RL[, `:=`(arch = factor(arch_lab, levels = ARCH_LABELS), engine_lab = factor(ENG_LAB[engine], levels = ENG_LAB))]
stopifnot(!anyNA(RL$arch), !anyNA(RL$engine_lab))
RX <- rank_split(RL, c("seed", "engine"))
RX <- merge(RX, RL[, .(seed, engine, panel, arch, engine_lab, landscape)], by = c("seed", "engine", "panel"))

if (MODE == "ab") {
    PANEL_ORDER <- c("2/3 methods", "1/3 methods", "LFMM", "RDA", "3/3 methods", "EMMAX")   # the published order
} else {
    ORD <- fread(file.path(TAB, "set_order.tsv"))[order(order)]
    PANEL_ORDER <- ORD$key
    stopifnot(setequal(PANEL_ORDER, RULES))
    message("set order: ", paste(set_lab(PANEL_ORDER), collapse = " < "))
}

# =============================================================================
# Detection plane (closeout sim_figures_ld.R)
# =============================================================================
SQRT_STEPS <- c(1, 2, 5, 10, 20, 25, 50, 100)
sqrt_even_breaks <- function(limits) {
    m <- max(limits, na.rm = TRUE)
    for (c in SQRT_STEPS) {
        n <- floor(sqrt(m / c)) + 1L
        if (n <= 6L) return(c * (seq_len(n) - 1L)^2)
    }
    c(0, m)
}
in_polygon <- function(px, py, x, y, tol = 1e-9) {
    n <- length(px); j <- n; inside <- FALSE
    for (i in seq_len(n)) {
        ex <- px[j] - px[i]; ey <- py[j] - py[i]
        if (abs((x - px[i]) * ey - (y - py[i]) * ex) < tol &&
            x >= min(px[i], px[j]) - tol && x <= max(px[i], px[j]) + tol &&
            y >= min(py[i], py[j]) - tol && y <= max(py[i], py[j]) + tol) return(TRUE)
        if (((py[i] > y) != (py[j] > y)) && (x < px[i] + (y - py[i]) * ex / ey)) inside <- !inside
        j <- i
    }
    inside
}
BAG_SHARE <- 0.5
bag_hull <- function(fp, tp) {
    x <- sqrt(fp); y <- sqrt(tp)
    sx <- IQR(x); sy <- IQR(y)
    if (length(x) < 10L || sx == 0 || sy == 0) return(NULL)
    mx <- sqrt(median(fp)); my <- sqrt(median(tp))
    k <- order(((x - mx) / sx)^2 + ((y - my) / sy)^2)[seq_len(ceiling(BAG_SHARE * length(x)))]
    hx <- c(x[k], mx); hy <- c(y[k], my)
    h <- chull(hx, hy)
    if (length(h) < 3L) return(NULL)
    list(FP = hx[h]^2, TP = hy[h]^2, has_median = in_polygon(hx[h], hy[h], mx, my))
}
bags_of <- function(D, keys) {
    B <- D[, { b <- bag_hull(FP, TP); if (is.null(b)) NULL else b }, by = keys]
    bad <- unique(B[has_median == FALSE, keys, with = FALSE])
    if (nrow(bad)) stop("bag without its median (would be an orphan region): ",
                        paste(do.call(paste, c(bad, sep = " / ")), collapse = "; "))
    B[, grp := do.call(paste, c(.SD, sep = "|")), .SDcols = keys][]
}
plane_ms <- function(D, facet, y_lab = TP_LAB, y0 = TRUE) {
    keys <- c("label", facet$vars)
    M  <- D[, .(TP = median(TP), FP = median(FP)), by = keys]
    B  <- bags_of(D, keys)
    p <- ggplot(D, aes(FP, TP)) +
        geom_polygon(data = B, aes(group = grp, colour = label), fill = NA, linewidth = 0.35) +
        geom_point(data = M, aes(fill = label), shape = 21, size = 1.7, colour = "black", stroke = 0.3) +
        geom_label_repel(data = M[, tone := ms_text_on(RULE_COL[as.character(label)])],
                         aes(label = set_lab(label), fill = label, colour = tone), size = MS_LAB,
                         fontface = "bold", label.size = 0, label.padding = 0.12, box.padding = 0.3,
                         min.segment.length = 0, segment.size = 0.25, segment.colour = MS_INK,
                         seed = JITTER_SEED, show.legend = FALSE, max.overlaps = Inf)
    p + facet$layer + (if (y0) expand_limits(x = 0, y = 0) else expand_limits(x = 0)) +
        scale_x_sqrt(breaks = sqrt_even_breaks) + scale_y_sqrt(breaks = sqrt_even_breaks) +
        scale_colour_manual(values = c(RULE_COL, TEXT_TONES), guide = "none") +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        labs(x = FP_LAB, y = y_lab) +
        theme_ms()
}
check_ticks <- function(p, what, axes = c("x", "y")) {
    for (pp in ggplot_build(p)$layout$panel_params) for (ax in axes) {
        b <- suppressWarnings(as.numeric(pp[[ax]]$get_labels()))
        b <- b[!is.na(b)]
        g <- diff(sqrt(b))
        if (length(b) < 2L || max(abs(g - g[1])) > 1e-9)
            stop(what, ": ", ax, " ticks not evenly spaced on the sqrt scale: ", paste(b, collapse = " "))
    }
    message(sprintf("%-34s every plot: %s ticks evenly spaced on the sqrt scale", what, paste(axes, collapse = " and ")))
}
y_title_fits <- function(p, h_mm, label) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    g <- ggplotGrob(p)
    fixed <- sum(grid::convertHeight(g$heights[grid::unitType(g$heights) != "null"], "mm", valueOnly = TRUE))
    need  <- grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = MS_BASE))),
                                "mm", valueOnly = TRUE)
    need <= h_mm - fixed
}
# Figure 2a: y starts just below the lowest bag, as published; the plots of a landscape keep their own axes
plane_top <- function(D, vars, layer, what) {
    spec <- list(vars = vars, layer = layer)
    p <- plane_ms(D, spec, y0 = FALSE)
    if (!y_title_fits(p, 62, TP_LAB)) p <- plane_ms(D, spec, y_lab = sub(" (", "\n(", TP_LAB, fixed = TRUE), y0 = FALSE)
    low <- min(bags_of(D, c("label", vars))$TP, D[, .(TP = median(TP)), by = c("label", vars)]$TP)
    for (pp in ggplot_build(p)$layout$panel_params) {
        r <- pp$y$continuous_range
        if (r[1] <= 0) message(sprintf("%s: y axis reaches 0 (lowest bag vertex or median %.1f)", what, low))
        else stopifnot(r[1]^2 < low)
    }
    check_ticks(p, what)
    p
}

# =============================================================================
# Within-replicate rank (closeout sim_figures_ld.R)
# =============================================================================
RANK_COL <- setNames(colorRampPalette(c(MINOU[["red"]], "#F2F2F2", MINOU[["teal"]]))(N_PANEL), N_PANEL:1)
rank_dist <- function(D, extra = NULL, facet_layer) {
    R <- D[, .(n = sum(w)), by = c("engine_lab", "panel", "r", extra)]
    R[, pct := 100 * n / sum(n), by = c("engine_lab", "panel", extra)]
    R[, r := factor(r, levels = N_PANEL:1)]
    R[, panel := factor(panel, levels = rev(PANEL_ORDER))]
    ggplot(R, aes(pct, panel, fill = r)) +
        geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
        facet_layer +
        scale_y_discrete(labels = SET_LAB) +
        scale_fill_manual(values = RANK_COL, name = "rank (1 = most accurate)",
                          breaks = as.character(1:N_PANEL)) +
        scale_x_continuous(breaks = c(0, 50, 100), labels = c("0", "50", "100"), expand = c(0, 0)) +
        labs(x = "% replicates at each rank", y = NULL) + theme_ms() +
        theme(panel.spacing.x = unit(4, "mm"))
}
grid2 <- function(a, b, heights, labels = c("a", "b"))
    cowplot::plot_grid(a, b, ncol = 1, rel_heights = heights, labels = labels, label_size = 8, label_fontface = "bold",
                       label_fontfamily = theme_ms()$text$family, align = "v", axis = "lr")
LAND_SHAPE <- c("SS-Clines" = 21, "SS-Mtn" = 24, "Est-Clines" = 22)
LAND_NUDGE <- c("SS-Clines" = 0.22, "SS-Mtn" = 0, "Est-Clines" = -0.22)
# white markers: the boundary between ranks 3 and 4 in each landscape, one row inside the bar per landscape
add_markers <- function(p, MK, name = "Boundary between ranks 3 and 4:") {
    for (L in LAND_LEV) p <- p + geom_point(data = MK[landscape == L], aes(x = x, y = panel, shape = landscape), inherit.aes = FALSE,
                                            fill = "white", colour = MS_INK, size = 1.25, stroke = 0.3,
                                            position = position_nudge(y = LAND_NUDGE[[L]]))
    p + scale_shape_manual(values = LAND_SHAPE, breaks = LAND_LEV, name = name, drop = FALSE) +
        guides(fill = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
        theme(legend.box = "horizontal")
}

# =============================================================================
# FIG_MODE=ab: the published figures, published order, SS-Clines alone
# =============================================================================
pDcells <- function(D) plane_ms(D, list(vars = "panel", layer = facet_wrap(vars(panel), ncol = 4, scales = "free_x")))
pDdemo  <- function(D) plane_ms(D, list(vars = "block", layer = facet_wrap(vars(block), nrow = 2, scales = "free_x",
                                                                            labeller = labeller(block = BLOCK_SHOW))))
HEAT_COL <- c(MINOU[["teal"]], "#FFFFFF", MINOU[["red"]])
make_heat <- function(chance, heat_max, step, name, v, cols = HEAT_COL) {
    at <- scales::rescale(c(0, chance, heat_max))
    heat_fill <- function(x) scales::gradient_n_pal(cols, at)(scales::rescale(x, from = c(0, heat_max)))
    function(D, x, facet_layer = NULL, legend = TRUE) {
        D <- copy(D)[, `:=`(panel = factor(panel, levels = rev(PANEL_ORDER)), least_pct = get(v))]
        D[, tone := ms_text_on(heat_fill(least_pct))]
        ggplot(D, aes(.data[[x]], panel, fill = least_pct)) +
            geom_tile(colour = "white", linewidth = 0.4) +
            geom_text(aes(label = sprintf("%.1f", least_pct), colour = tone), size = MS_LAB) +
            facet_layer +
            scale_fill_gradientn(colours = cols, values = at, limits = c(0, heat_max),
                                 breaks = seq(0, heat_max, step), labels = function(b) paste0(b, "%"), name = name,
                                 guide = if (legend) guide_colourbar(theme = theme(
                                     legend.key.width = grid::unit(45, "mm"), legend.key.height = grid::unit(2.5, "mm"),
                                     legend.title.position = "top")) else "none") +
            scale_colour_identity() +
            scale_y_discrete(labels = SET_LAB) +
            labs(x = NULL, y = NULL) + theme_ms() +
            theme(axis.line = element_blank(), axis.ticks = element_blank())
    }
}
strata_fig <- function(L, heat) {
    M3 <- RL[engine == "median of 3" & landscape == L, .(seed, panel, rank)]
    COV <- unique(DET[landscape == L, .(seed, arch, pleio, regime, block)])
    I3 <- merge(rank_indicators(rank_split(M3, "seed"), "seed", n_set = N_PANEL), COV, by = "seed")
    strat <- function(cols) I3[, .(least_pct = 100 * mean(bottom1), half_pct = 100 * mean(bottom3), top3_pct = 100 * mean(top3), n = .N),
                               by = c("panel", cols)]
    C <- strat(c("arch", "pleio", "regime")); B <- strat("block")
    a <- heat(C, "regime", facet_grid(cols = vars(arch, pleio), labeller = labeller(arch = ARCH_WRAP)), legend = FALSE)
    b <- heat(B, "block") + scale_x_discrete(labels = function(x) unname(BLOCK_SHOW[x]))
    list(fig = grid2(a, b, c(56, 54)), C = C, B = B)
}
single_cloud <- function(D, y_col, labels, lev, colours, facet_layer) {
    ggplot(D, aes(aucpr, ypos)) +
        stat_halfeye(aes(fill = .data[[y_col]], group = .data[[y_col]]), orientation = "horizontal", adjust = 0.8,
                     height = 0.55, justification = -0.2, .width = 0, point_colour = NA, slab_alpha = 0.75) +
        geom_boxplot(aes(group = .data[[y_col]]), orientation = "y", width = 0.14, outlier.shape = NA,
                     fill = "white", colour = MS_INK, linewidth = 0.25) +
        geom_point(aes(y = ypos - 0.25, colour = .data[[y_col]]),
                   position = position_jitter(height = 0.08, width = 0, seed = JITTER_SEED),
                   size = 0.6, alpha = 0.6, stroke = 0) +
        facet_layer +
        scale_y_continuous(breaks = seq_along(lev), labels = labels) +
        scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
        scale_fill_manual(values = colours, guide = "none") + scale_colour_manual(values = colours, guide = "none") +
        labs(x = "AUC-PR (causal loci vs neutral SNPs)", y = NULL) + theme_ms()
}
# S1 input (mvp_ms_rda_correction.R section 1, the _ldb arms hold the same AUC-PR)
S1_SERIES <- c("RDA uncorrected", "RDA corrected", "LFMM")
S1_COL <- c("RDA uncorrected" = MINOU[["amber"]], "RDA corrected" = MINOU[["navy"]], "LFMM" = MS_REF)
s1_data <- function() {
    rbindlist(lapply(seq_len(nrow(LANDS)), function(i) {
        rk_c <- rd(EVAL, LANDS$rem[i], "onefit_rda2x_ldb/rank_metrics.tsv")
        rk_u <- rd(EVAL, LANDS$rem[i], "rdaunc_rda2x_ldb/rank_metrics.tsv")
        A <- rbind(rk_c[method == "RDA", .(seed, series = "RDA corrected", aucpr)], rk_u[method == "RDA", .(seed, series = "RDA uncorrected", aucpr)],
                   rk_c[method == "LFMM", .(seed, series = "LFMM", aucpr)])
        merge(A, unique(DET[landscape == LANDS$land[i], .(seed, arch, landscape)]), by = "seed")
    }))[, series := factor(series, levels = S1_SERIES)][, ypos := as.numeric(factor(series, levels = rev(S1_SERIES)))]
}
s2_data <- function(single) {
    A <- rbindlist(lapply(seq_len(nrow(LANDS)), function(i) rd(EVAL, LANDS$gal[i], "B4_aucpr_per_seed.tsv")[, landscape := factor(LANDS$land[i], levels = LAND_LEV)]))
    A[, `:=`(arch = factor(arch, levels = ARCH_LABELS), label = factor(label, levels = single))]
    stopifnot(!anyNA(A$arch), !anyNA(A$label), !anyNA(A$aucpr))
    A[, ypos := as.numeric(factor(label, levels = rev(single)))][]
}
CUTS <- c(top_0.001 = 0.1, top_0.0025 = 0.25, top_0.005 = 0.5, top_0.01 = 1, top_0.02 = 2)
s3_data <- function() {
    rbindlist(lapply(seq_len(nrow(LANDS)), function(i) {
        arm <- file.path(LANDS$rem[i], "rdaunc_rda2x_ldb")
        cpm <- rd(EVAL, arm, "calls_per_method.tsv"); ccb <- rd(EVAL, arm, "calls_combine.tsv")
        cnt <- function(D, lab) D[, .(seed, rung, label = lab, n = n_called, TP = tp_any_causal + expected_linked, FP = fp_background)]
        CUT <- rbind(cnt(ccb[combine == "union" & window_kb > 0], "1/3 methods"), cnt(ccb[combine == "ge2" & window_kb > 0], "2/3 methods"),
                     cnt(ccb[combine == "all3" & window_kb > 0], "3/3 methods"),
                     cpm[, .(seed, rung, label = method, n = n_called, TP = tp_any_causal + expected_linked, FP = fp_background)])
        CUT[, cutoff := unname(CUTS[rung])]
        CM <- CUT[, .(TP = as.numeric(median(TP)), FP = as.numeric(median(FP)), empty = sum(n == 0L), n_rep = .N), by = .(cutoff, label)]
        CL <- melt(CM, id.vars = c("cutoff", "label"), measure.vars = c("TP", "FP"), variable.name = "what", value.name = "count")
        CL[, `:=`(what = factor(c(TP = TP_LAB, FP = FP_LAB)[as.character(what)], levels = c(TP_LAB, FP_LAB)),
                  label = factor(label, levels = PANEL_ORDER), landscape = factor(LANDS$land[i], levels = LAND_LEV))]
    }))
}
s3_plot <- function(CL, facet_layer) {
    ggplot(CL, aes(cutoff, count, colour = label)) +
        geom_vline(xintercept = 0.25, linetype = "dashed", colour = MS_REF, linewidth = 0.35) +
        geom_line(linewidth = 0.45) + geom_point(size = 1.2) +
        facet_layer +
        expand_limits(y = 0) +
        scale_x_log10(breaks = unname(CUTS), labels = function(x) paste0(x, "%")) +
        scale_y_sqrt(breaks = sqrt_even_breaks) +
        scale_colour_manual(values = RULE_COL, breaks = PANEL_ORDER, labels = SET_LAB, name = NULL) +
        guides(colour = guide_legend(nrow = 1)) +
        labs(x = "Cutoff (% of tested SNPs per method and predictor)", y = "Median number of SNPs") + theme_ms()
}
SRC_DET <- "detection_per_seed.tsv"; SRC_RANK <- "rank_long.tsv"
if (MODE == "ab") {
    pA1 <- plane_top(DET, "arch", facet_grid(cols = vars(arch), scales = "free_x"), "A1")
    pR1 <- rank_dist(RX, NULL, facet_wrap(~ engine_lab, nrow = 1))
    save_ms("Fig2_simulation", grid2(pA1, pR1, c(62, 66)), 180, 128, "ab", SRC_DET)
    save_ms("s_detection_cells", pDcells(DET), 180, 165, "ab", SRC_DET)
    save_ms("s_detection_demography", pDdemo(DET), 180, 110, "ab", SRC_DET)
    save_ms("s_rank_architecture", rank_dist(RX, "arch", facet_grid(rows = vars(arch), cols = vars(engine_lab))), 180, 120, "ab", SRC_RANK)
    S1 <- s1_data()
    save_ms("s_rda_correction", single_cloud(S1, "series", rev(sub("^LFMM$", "LFMM2", S1_SERIES)), S1_SERIES, S1_COL,
                                             facet_grid(cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP))), 180, 60, "ab", "rank_metrics.tsv")
    SINGLE <- intersect(PANEL_ORDER, c("LFMM", "RDA", "EMMAX"))
    A2 <- s2_data(SINGLE)
    save_ms("s_aucpr_methods", single_cloud(A2, "label", rev(set_lab(SINGLE)), SINGLE, RULE_COL,
                                            facet_grid(cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP))), 180, 60, "ab", "B4_aucpr_per_seed.tsv")
    save_ms("s_cutoffs", s3_plot(s3_data(), facet_wrap(vars(what), nrow = 1, scales = "free_y")), 180, 70, "ab", "calls_combine.tsv")
    MAXV <- strata_fig("SS-Clines", make_heat(100 / N_PANEL, 100, 10, "x", "least_pct"))
    hm <- max(20, 10 * ceiling(max(MAXV$B$least_pct, MAXV$C$least_pct) / 10))
    save_ms("s_least_accurate_strata", strata_fig("SS-Clines", make_heat(100 / N_PANEL, hm, 10,
            "Replicates in which the set is least accurate (white: chance, 16.7%)", "least_pct"))$fig, 180, 110, "ab", SRC_RANK)
}

# =============================================================================
# FIG_MODE=main: Figure 2 in four variants, Figure 2 per landscape
# =============================================================================
if (MODE == "main") {
    TR <- fread(file.path(TAB, "t_ranks.tsv"), sep = "\t")
    MK <- TR[statistic == "ranks 1–3" & landscape %in% LAND_LEV, .(landscape = factor(landscape, levels = LAND_LEV), engine_lab = factor(offset_method, levels = ENG_LAB),
                                                                   panel = factor(unname(KEY_OF[set]), levels = rev(PANEL_ORDER)), x = value)]
    stopifnot(nrow(MK) == 3L * 4L * N_PANEL, !anyNA(MK))
    pA1 <- plane_top(DET, "arch", facet_grid(cols = vars(arch), scales = "free_x"), "Figure 2a, all replicates")
    pR1 <- rank_dist(RX, NULL, facet_wrap(~ engine_lab, nrow = 1))
    save_ms("Fig2_v1", grid2(pA1, pR1, c(62, 66)), 180, 128, "main", paste(SRC_DET, "(a),", SRC_RANK, "(b), the 1,800 replicates"))
    save_ms("Fig2_v2", grid2(pA1, add_markers(pR1, MK), c(62, 70)), 180, 132, "main", paste(SRC_DET, "(a),", SRC_RANK, "+ t_ranks.tsv (b)"))
    # v3: a by landscape; b all landscapes, then each, under the median of the three
    pA3 <- plane_top(DET[, scope := landscape], "scope", facet_grid(cols = vars(scope), scales = "free_x"), "Figure 2a, by landscape")
    X3 <- RX[engine == "median of 3"]
    B3 <- rbind(copy(X3)[, engine_lab := "All landscapes"], copy(X3)[, engine_lab := as.character(landscape)])
    B3[, engine_lab := factor(engine_lab, levels = c("All landscapes", LAND_LEV))]
    pB3 <- rank_dist(B3, NULL, facet_wrap(~ engine_lab, nrow = 1))
    save_ms("Fig2_v3", grid2(pA3, pB3, c(62, 66)), 180, 128, "main", paste(SRC_DET, "(a),", SRC_RANK, "(b)"))
    # v4: a as v1; b median of the three with the markers; c median accuracy lost per set
    pB4 <- add_markers(rank_dist(RX[engine == "median of 3"], NULL, facet_wrap(~ engine_lab, nrow = 1)), MK[engine_lab == "Median of the three"]) +
        theme(legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"), legend.box.just = "left")
    LS <- fread(file.path(TAB, "t_offset_accuracy.tsv"), sep = "\t")[offset_method == "Median of the three" & set != "Causal loci"]
    LS[, panel := factor(unname(KEY_OF[set]), levels = rev(PANEL_ORDER))]
    LB <- LS[landscape == "all"]
    LM <- LS[landscape %in% LAND_LEV][, landscape := factor(landscape, levels = LAND_LEV)]
    LR <- LM[, .(lo = min(median_loss), hi = max(median_loss)), by = panel]
    pC4 <- ggplot() +
        geom_col(data = LB, aes(median_loss, panel, fill = panel), width = 0.72) +
        geom_segment(data = LR, aes(x = lo, xend = hi, y = panel, yend = panel), colour = MS_INK, linewidth = 0.3) +
        geom_point(data = LM, aes(median_loss, panel, shape = landscape), fill = "white", colour = MS_INK, size = 1.4, stroke = 0.3) +
        scale_fill_manual(values = RULE_COL, guide = "none") +
        scale_shape_manual(values = LAND_SHAPE, guide = "none") +
        scale_y_discrete(labels = SET_LAB) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
        labs(x = "Median accuracy lost against\nthe replicate's best set", y = NULL) + theme_ms()
    row2 <- cowplot::plot_grid(pB4, pC4, nrow = 1, rel_widths = c(108, 72), labels = c("b", "c"), label_size = 8, label_fontface = "bold",
                               label_fontfamily = theme_ms()$text$family, align = "h", axis = "tb")
    save_ms("Fig2_v4", cowplot::plot_grid(pA1, row2, ncol = 1, rel_heights = c(62, 72), labels = c("a", ""), label_size = 8, label_fontface = "bold",
                                          label_fontfamily = theme_ms()$text$family), 180, 134, "main",
            paste(SRC_DET, "(a),", SRC_RANK, "+ t_ranks.tsv (b), t_offset_accuracy.tsv (c)"))
    for (L in LAND_LEV) {
        D <- DET[landscape == L]; X <- RX[landscape == L]
        pa <- plane_top(D, "arch", facet_grid(cols = vars(arch), scales = "free_x"), paste("Figure 2a,", L))
        pb <- rank_dist(X, NULL, facet_wrap(~ engine_lab, nrow = 1))
        save_ms(paste0("Fig2_per_landscape_", PART_OF[[L]]), grid2(pa, pb, c(62, 66)), 180, 128, "supp part", paste(L, SRC_DET, SRC_RANK))
    }
}

# =============================================================================
# FIG_MODE=supp
# =============================================================================
if (MODE == "supp") {
    # row per landscape
    S1 <- s1_data()[, landscape := factor(landscape, levels = LAND_LEV)]
    save_ms("s_rda_correction", single_cloud(S1, "series", rev(sub("^LFMM$", "LFMM2", S1_SERIES)), S1_SERIES, S1_COL,
                                             facet_grid(rows = vars(landscape), cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP))),
            180, 170, "supp", "{rdaunc,onefit}_rda2x_ldb/rank_metrics.tsv")
    SINGLE <- intersect(PANEL_ORDER, c("LFMM", "RDA", "EMMAX"))
    A2 <- s2_data(SINGLE)
    save_ms("s_aucpr_methods", single_cloud(A2, "label", rev(set_lab(SINGLE)), SINGLE, RULE_COL,
                                            facet_grid(rows = vars(landscape), cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP))),
            180, 170, "supp", "B4_aucpr_per_seed.tsv")
    CL <- s3_data()
    pS3 <- s3_plot(CL, facet_wrap(vars(landscape, what), ncol = 2, scales = "free_y"))
    check_ticks(pS3, "s_cutoffs", axes = "y")
    save_ms("s_cutoffs", pS3, 180, 170, "supp", "rdaunc_rda2x_ldb/{calls_per_method,calls_combine}.tsv")
    # parts, one page per landscape
    HEAT_R <- make_heat(50, 100, 25, "Replicates in ranks 4–6 (white: chance, 50%)", "half_pct")
    PC <- list()
    for (L in LAND_LEV) {
        k <- PART_OF[[L]]; D <- DET[landscape == L]
        p <- pDcells(D); check_ticks(p, paste("s_detection_cells", L)); PC[[L]] <- p
        save_ms(paste0("s_detection_cells_parts_", k), p, 180, 165, "supp part", paste(L, SRC_DET))
        p <- pDdemo(D); check_ticks(p, paste("s_detection_demography", L))
        save_ms(paste0("s_detection_demography_parts_", k), p, 180, 110, "supp part", paste(L, SRC_DET))
        save_ms(paste0("s_rank_architecture_parts_", k), rank_dist(RX[landscape == L], "arch", facet_grid(rows = vars(arch), cols = vars(engine_lab))),
                180, 120, "supp part", paste(L, SRC_RANK))
        save_ms(paste0("s_rank_strata_parts_", k), strata_fig(L, HEAT_R)$fig, 180, 110, "supp part", paste(L, SRC_RANK, "+", SRC_DET))
    }
    # all landscapes in one figure
    save_ms("s_detection_cells_rows", cowplot::plot_grid(plotlist = unname(PC[LAND_LEV]), ncol = 1, labels = paste(c("a", "b", "c"), LAND_LEV), label_size = 8,
                                                         label_fontface = "bold", label_fontfamily = theme_ms()$text$family, hjust = 0, label_x = 0.01,
                                                         vjust = 1), 180, 495, "supp", SRC_DET)
}
# =============================================================================
# FIG_MODE=merge: the merged figures of the story restructure (Data S1 Job F, 2026-11-08)
# =============================================================================
if (MODE == "merge") {
    dir.create(CHECK, recursive = TRUE, showWarnings = FALSE)
    save_member <- function(stem, p, w, h) { ms_save(file.path(CHECK, stem), p, w, h); message("  member ", stem) }
    page_title <- function(p, title, h_body)                 # the landscape name above a page, panel letters unchanged
        cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(title, fontfamily = theme_ms()$text$family, fontface = "bold",
                                                                   size = 9, x = 0.01, hjust = 0),
                           p, ncol = 1, rel_heights = c(7, h_body))
    # F1 s_precision_recall: a = s_aucpr_methods, b = s_rda_correction, each as published (180 x 170 mm)
    SINGLE <- intersect(PANEL_ORDER, c("LFMM", "RDA", "EMMAX"))
    A2 <- s2_data(SINGLE)
    pa <- single_cloud(A2, "label", rev(set_lab(SINGLE)), SINGLE, RULE_COL,
                       facet_grid(rows = vars(landscape), cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP)))
    S1 <- s1_data()[, landscape := factor(landscape, levels = LAND_LEV)]
    pb <- single_cloud(S1, "series", rev(sub("^LFMM$", "LFMM2", S1_SERIES)), S1_SERIES, S1_COL,
                       facet_grid(rows = vars(landscape), cols = vars(arch), scales = "free_x", labeller = labeller(arch = ARCH_WRAP)))
    save_member("s_aucpr_methods", pa, 180, 170); save_member("s_rda_correction", pb, 180, 170)
    save_ms("s_precision_recall", grid2(pa, pb, c(170, 170)), 180, 340, "supp", "B4_aucpr_per_seed.tsv; {rdaunc,onefit}_rda2x_ldb/rank_metrics.tsv")
    # F2 s_detection_design_parts_*: a = design cells (180 x 165), b = demographic scenarios (180 x 110); F5 s_rank_strata_parts_*
    HEAT_R <- make_heat(50, 100, 25, "Replicates in ranks 4–6 (white: chance, 50%)", "half_pct")           # the published key
    HEAT_T <- make_heat(50, 100, 25, "% of replicates in the top three (white: 50%, the average of the six sets)", "top3_pct",
                        cols = rev(HEAT_COL))                                  # teal above 50%, red below
    ST <- list()
    for (L in LAND_LEV) {
        k <- PART_OF[[L]]; D <- DET[landscape == L]
        a <- pDcells(D); check_ticks(a, paste("s_detection_cells", L))
        b <- pDdemo(D); check_ticks(b, paste("s_detection_demography", L))
        save_member(paste0("s_detection_cells_parts_", k), a, 180, 165)
        save_member(paste0("s_detection_demography_parts_", k), b, 180, 110)
        save_ms(paste0("s_detection_design_parts_", k), page_title(grid2(a, b, c(165, 110)), L, 275), 180, 282, "supp part",
                paste(L, SRC_DET))
        save_member(paste0("s_rank_strata_parts_", k), strata_fig(L, HEAT_R)$fig, 180, 110)
        sf <- strata_fig(L, HEAT_T)
        save_ms(paste0("s_rank_strata_parts_", k), page_title(sf$fig, L, 110), 180, 117, "supp part", paste(L, SRC_RANK, "+", SRC_DET))
        ST[[L]] <- rbind(sf$C[, .(landscape = L, stratum_kind = "design cell", stratum = paste(arch, pleio, regime, sep = ", "),
                                  panel, least_pct, half_pct, top3_pct, n)],
                         sf$B[, .(landscape = L, stratum_kind = "demographic scenario", stratum = as.character(BLOCK_SHOW[as.character(block)]),
                                  panel, least_pct, half_pct, top3_pct, n)])
    }
    fwrite(rbindlist(ST), file.path(CHECK, "..", "strata_top3.tsv"), sep = "\t")
}
F <- rbindlist(lapply(names(FIGS), function(k) data.table(stem = k, role = FIGS[[k]]$role, w_mm = FIGS[[k]]$w, h_mm = FIGS[[k]]$h, source = FIGS[[k]]$source)))
fwrite(F, file.path(OUT, sprintf("figures_%s.tsv", MODE)), sep = "\t")
message(sprintf("=== %s: %d figures -> %s", MODE, nrow(F), OUT))
