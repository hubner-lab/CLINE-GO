#!/usr/bin/env Rscript
# Manuscript pack for the genomic-offset result (user 2026-10-15): main figure + supplementary figures + tables.
# Scope fixed by the user: 27-SNP set (EMMAX x RDA agreement), predictors bio_1,bio_3,bio_5, SSP585 2061-2080 mean of
# 3 GCMs, methods GF / geometric (k 5) / RDA (3 axes). No per-GCM or uncertainty decomposition.
# Rank statistic: Kendall tau-b everywhere (user 2026-10-15; 11 sites = 55 site pairs).
# Consensus = family-weighted rank (Gain et al. 2023: RDA and geometric offsets are theoretically equivalent):
#   family ranks: GF rank, linear rank = mean(geometric, RDA rank); consensus score = mean of the two; 1 = highest risk.
#   category: "higher risk" = upper half (rank < 6) under BOTH families, "lower risk" = lower half (rank > 6) under
#   both, else "methods disagree".
# Layout chosen by the user 2026-10-15 (review of the first pack):
#   main fig_main = a consensus map, b climate-only reference map (no test claim), c loss of isothermality vs
#   linear-family rank; supplementary fig_s_random = adaptive vs random 27-SNP sets (Kendall), fig_s_delta = each
#   predictor's projected change vs the GF and linear-family ranks.
# Reads only existing outputs: site_maps/site_offsets.tsv, site_maps/gf_random_null/, site_maps/random_null/,
#   site_maps/kendall_check.tsv (cross-check). Writes paper_pack/: fig_main, fig_s_random, fig_s_delta (png + svg),
#   consensus_sites.tsv, random_loo.tsv, site_climate_shift.tsv, delta_rank_tau.tsv.
# Run inside the image from /pipeline:  Rscript ANALYSIS/Trifolium86MainK5/maladaptationER/build_paper_pack.R
suppressPackageStartupMessages({ library(data.table); library(terra); library(ggplot2); library(patchwork) })
B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R   <- file.path("/pipeline", paste0(Sys.getenv("PROJ", "Trifolium86MainK5ER"), "_results"))
SC  <- "year2061-2080_ssp585"
SET <- "gea.ER_overlap"
P   <- c("bio_1", "bio_3", "bio_5")
OUT <- file.path(B, "paper_pack"); dir.create(OUT, showWarnings = FALSE)
PAL <- c("#095786", "#0090ab", "#e6cea0", "#db4646", "#ae0000")
CAT_COL <- c("higher risk" = "#ae0000", "methods disagree" = "#e6cea0", "lower risk" = "#095786")
HIST_FILL <- "#0090ab"; ADAPT_LINE <- "#ae0000"
# PALETTE=minou (2026-10-16, plan scalable-doodling-allen.md): the whole paper in minou. Categories, histogram and
# adaptive line switch; the gradient PAL of the maps (Fig 1b, S5) stays (user). Default = the v1 colours, unchanged.
if (Sys.getenv("PALETTE") == "minou") {
    MINOU_ENV <- new.env(); sys.source("/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R", envir = MINOU_ENV)
    CAT_COL <- MINOU_ENV$FIG3_RISK; HIST_FILL <- get("MINOU", MINOU_ENV)[["teal"]]; ADAPT_LINE <- get("MINOU", MINOU_ENV)[["red"]]
}
# SITE_NAMES=<site_names.tsv> (2026-10-16, plan fancy-scribbling-garden.md Phase 5): site point labels show the short
# names 1..11 (S->N); tables keep the raw ids. Unset = raw ids, unchanged.
# With it, the map badges (Fig 1a, 1b) carry the site name instead of the risk rank (user 2026-10-16, option B: a number
# on a point is always the site; rank = colour, axis and Table 1), the titles drop "(1 = highest)", and the Fig 1a
# legend keys are filled squares instead of the "a" glyph of geom_label.
site_lab <- identity
MAP_SITE <- nzchar(Sys.getenv("SITE_NAMES")); KEY_GLYPH <- "label"
if (MAP_SITE) {
    SN <- fread(Sys.getenv("SITE_NAMES"), colClasses = "character")
    site_lab <- function(s) { m <- SN$short[match(s, SN$site_id)]; stopifnot(!anyNA(m)); m }
    KEY_GLYPH <- "rect"
}
# Badge text colour (2026-10-17, plan idempotent-conjuring-noodle.md Phase 3; only with SITE_NAMES): white or ink,
# whichever has the higher WCAG contrast on the badge fill (dark digits on navy / dark red were unreadable).
txt_on <- function(bg, ink = "grey10") {
    lum <- function(x) { v <- grDevices::col2rgb(x) / 255; v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
                         colSums(v * c(0.2126, 0.7152, 0.0722)) }
    l <- lum(bg); ifelse(1.05 / (l + 0.05) >= (l + 0.05) / (lum(ink) + 0.05), "white", ink)
}
METH <- c(gf = "Gradient Forest", geometric = "Geometric", rda = "RDA (3 axes)")
tau  <- function(a, b) cor(a, b, method = "kendall")
# change shown so that larger = more change: bio_3 falls at every site, so its loss is -delta
DLAB <- c(bio_1 = "Warming, annual mean temp. bio_1 (°C)", bio_3 = "Loss of isothermality bio_3 (−Δ, %)",
          bio_5 = "Warming, max temp. warmest month bio_5 (°C)")
FAM  <- c(gf_rank = "Gradient Forest rank", linear_rank = "Linear-family rank")
NO_LINEAR <- Sys.getenv("NO_LINEAR") == "1"   # see the consensus block below
# SUPP_MS=1 (2026-10-17, plan idempotent-conjuring-noodle.md Phase 5 follow-up; user: "redraw S15-S19 at final size"):
# fig_s_random, fig_s_delta_heat and fig_s_offset_maps are also saved at final mm size in the manuscript style
# (benchmarks/ms_style.R ms_save, 7/6 pt Arial) as <name>_ms.{svg,png}; the plot builders take ms = FALSE by default,
# so the unchanged files keep their sizes.
SUPP_MS <- Sys.getenv("SUPP_MS") == "1"
if (SUPP_MS) { MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE) }
if (NO_LINEAR) FAM <- c(gf_rank = "Gradient Forest rank", rank_rda = "RDA offset rank", rank_geometric = "Geometric offset rank")

samples <- fread(file.path(R, "_intermediate/samples/metadata_climate_valid.tsv"),
                 colClasses = c(site = "character", sample = "character"))
xy <- samples[, .(x = mean(longitude), y = mean(latitude), n = .N), by = site]

# ── consensus ───────────────────────────────────────────────────────────────
sm <- fread(file.path(B, "site_maps/site_offsets.tsv"), colClasses = c(site = "character"), na.strings = "")
so <- sm[set == SET & method %in% names(METH), .(site, method, offset)]
so[, rank := frank(-offset), by = method]
w <- dcast(so, site ~ method, value.var = "rank")
w[, `:=`(gf_rank = gf, linear_rank = (geometric + rda) / 2)]
w[, consensus_score := (gf_rank + linear_rank) / 2]
w[, consensus_rank := frank(consensus_score, ties.method = "min")]
w[, category := fifelse(gf_rank < 6 & linear_rank < 6, "higher risk",
                 fifelse(gf_rank > 6 & linear_rank > 6, "lower risk", "methods disagree"))]
# NO_LINEAR=1 (2026-10-17, user: "no linear ranking at all"; consensus = "all 3 methods agree"): each method votes with
# its own rank. Consensus score = mean of the three ranks; category = upper half (rank < 6) under GF, RDA AND geometric,
# lower half under all three, else "methods disagree". On this data the categories equal the family rule (asserted);
# consensus_rank changes because the score no longer averages RDA and geometric first. Unset = v1/v2 rule, unchanged.
if (NO_LINEAR) {
    old_cat <- w$category
    w[, consensus_score := (gf + geometric + rda) / 3]
    w[, consensus_rank := frank(consensus_score, ties.method = "min")]
    w[, category := fifelse(gf < 6 & rda < 6 & geometric < 6, "higher risk",
                     fifelse(gf > 6 & rda > 6 & geometric > 6, "lower risk", "methods disagree"))]
    stopifnot(identical(w$category, old_cat))
    w[, linear_rank := NULL]
}
pres_site <- fread(file.path(R, "climate/tables/present/climate_present_site.tsv"), colClasses = c(sample = "character"))
fut_site  <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_site.tsv")),
                   colClasses = c(sample = "character"))
clim <- merge(pres_site[, c("sample", P), with = FALSE], fut_site[, c("sample", P), with = FALSE], by = "sample",
              suffixes = c("_present", "_future"))[samples[, .(sample, site)], on = "sample"]
clim <- clim[, lapply(.SD, mean), by = site, .SDcols = patterns("_present$|_future$")]
for (p in P) clim[, (paste0("d_", p)) := get(paste0(p, "_future")) - get(paste0(p, "_present"))]
fwrite(clim, file.path(OUT, "site_climate_shift.tsv"), sep = "\t")
cons <- merge(merge(w, xy, by = "site"), clim[, c("site", paste0(P, "_present"), paste0("d_", P)), with = FALSE], by = "site")
setnames(cons, c("gf", "geometric", "rda"), c("rank_gf", "rank_geometric", "rank_rda"))
cons[, `:=`(chg_bio_1 = d_bio_1, chg_bio_3 = -d_bio_3, chg_bio_5 = d_bio_5)]
fwrite(cons[order(consensus_score)], file.path(OUT, "consensus_sites.tsv"), sep = "\t")
print(cons[order(consensus_score), intersect(c("site", "rank_gf", "rank_geometric", "rank_rda", "linear_rank", "consensus_score", "category"), names(cons)), with = FALSE])

# ── map panels (a consensus, b climate-only reference) ─────────────────────
r0   <- rast(file.path(R, "climate/rasters/present/climate_present_rasterstack.tif"))[[1]]
land <- as.polygons(!is.na(r0)); land <- land[land[[1]][, 1] == 1]
lg   <- as.data.frame(geom(land)); lg$grp <- paste(lg$geom, lg$part, lg$hole)
ex   <- as.vector(ext(r0))
# The site badge (SITE_NAMES only): same rounded label as the maps, text colour per badge from d$txt
badge_layer <- function(d, size, r = 0.5, pad = 0.25) list(
    ggrepel::geom_label_repel(data = d, aes(x, y, label = lab, fill = fillv, colour = txt), size = size,
                              fontface = "bold", label.r = unit(r, "lines"), label.padding = unit(pad, "lines"),
                              label.size = 0.3, box.padding = 0.2, min.segment.length = 0, segment.colour = "grey30",
                              segment.size = 0.3, seed = 42, max.overlaps = Inf, key_glyph = KEY_GLYPH),
    scale_colour_identity())
map_panel <- function(d, fill_scale, title, badge = 3, r = 0.5, pad = 0.25) {
    ggplot() +
        geom_polygon(data = lg, aes(x, y, group = grp), fill = "grey95", colour = "grey60", linewidth = 0.25) +
        geom_point(data = d, aes(x, y), size = 0.9, colour = "grey10") +
        (if (MAP_SITE) badge_layer(d, size = badge, r = r, pad = pad) else
        ggrepel::geom_label_repel(data = d, aes(x, y, label = lab, fill = fillv), size = 3,
                                  fontface = "bold", colour = "grey10", label.r = unit(0.5, "lines"),
                                  label.padding = unit(0.25, "lines"), label.size = 0.3, box.padding = 0.2,
                                  min.segment.length = 0, segment.colour = "grey30", segment.size = 0.3, seed = 42,
                                  max.overlaps = Inf, key_glyph = KEY_GLYPH)) +
        fill_scale +
        coord_equal(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 10) +
        theme(plot.title = element_text(face = "bold", size = 10),
              panel.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
              legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8))
}
pa <- map_panel(cons[, .(x, y, lab = if (MAP_SITE) site_lab(site) else consensus_rank, fillv = category,
                         txt = txt_on(CAT_COL[category]))],
                scale_fill_manual(values = CAT_COL, name = NULL),
                if (MAP_SITE) "a  Consensus site risk" else "a  Consensus site risk (1 = highest)")
nul <- merge(sm[method == "null", .(site, offset)], xy, by = "site")[, lab := frank(-offset, ties.method = "min")]
if (MAP_SITE) nul[, lab := site_lab(site)]
pb <- map_panel(nul[, .(x, y, lab, fillv = offset, txt = txt_on(scales::gradient_n_pal(PAL)(scales::rescale(offset))))],
                scale_fill_gradientn(colours = PAL, name = "Climate shift (SD units)",
                                     guide = guide_colourbar(barwidth = unit(3.2, "cm"), barheight = unit(0.25, "cm"),
                                                             title.position = "top")),
                if (MAP_SITE) "b  Climate-only reference" else "b  Climate-only reference (1 = largest shift)")

# ── predictor change vs family rank (main c = bio_3 x linear; supplementary grid = all) ──
dr <- rbindlist(lapply(P, function(p) rbindlist(lapply(names(FAM), function(f)
    data.table(predictor = p, family = f, kendall_tau = tau(cons[[paste0("chg_", p)]], -cons[[f]]),
               p_exact = suppressWarnings(cor.test(cons[[paste0("chg_", p)]], -cons[[f]], method = "kendall",
                                                   exact = TRUE)$p.value))))))
fwrite(dr, file.path(OUT, "delta_rank_tau.tsv"), sep = "\t")
print(dr)
delta_panel <- function(p, f, title = NULL, base = 9, ylab = NULL, badge = 2.6, r = 0.5, pad = 0.25) {
    t <- dr[predictor == p & family == f]
    ggplot(cons, aes(.data[[paste0("chg_", p)]], .data[[f]])) +
        (if (MAP_SITE) list(   # SITE_NAMES: the map badge (user review v2: "same shapes for labels as on maps")
            geom_point(size = 0.9, colour = "grey10"),
            badge_layer(cons[, .(x = get(paste0("chg_", p)), y = get(f), lab = site_lab(site), fillv = category,
                                 txt = txt_on(CAT_COL[category]))], size = badge, r = r, pad = pad)) else list(
        geom_point(aes(fill = category), shape = 21, size = 2.6, colour = "grey15"),
        ggrepel::geom_text_repel(aes(label = site_lab(site)), size = 2.4, colour = "grey30", seed = 42, min.segment.length = 0))) +
        scale_fill_manual(values = CAT_COL, guide = "none") + scale_y_reverse(breaks = c(1, 6, 11)) +
        labs(x = DLAB[[p]], y = paste(ylab %||% FAM[[f]], "(1 = highest risk)"),
             title = title %||% sprintf("Kendall tau = %.2f", t$kendall_tau)) +
        theme_classic(base_size = base) + theme(plot.title = element_text(face = "bold", size = base + 1))
}
if (!NO_LINEAR) pc <- delta_panel("bio_3", "linear_rank",
                  sprintf("c  Isothermality loss orders the linear ranking (tau %.2f)",
                          dr[predictor == "bio_3" & family == "linear_rank", kendall_tau]), base = 9)
# NO_LINEAR=1 (2026-10-17, plan idempotent-conjuring-noodle.md; user: "no linear ranking at all, redraw with RDA"):
# panel c plots the RDA offset rank, not the mean of the RDA and geometric ranks. Unset = v1/v2 panel, unchanged.
if (NO_LINEAR) {
    t_rda <- tau(cons$chg_bio_3, -cons$rank_rda)
    pc <- delta_panel("bio_3", "rank_rda", sprintf("c  Isothermality loss orders the RDA ranking (tau %.2f)", t_rda),
                      base = 9, ylab = "RDA offset rank")
}
fm <- pa + pb + pc + plot_layout(widths = c(1, 1, 1.15))
ggsave(file.path(OUT, "fig_main.png"), fm, width = 13.6, height = 5.6, dpi = 250, bg = "white")
ggsave(file.path(OUT, "fig_main.svg"), fm, width = 13.6, height = 5.6, bg = "white")
# MAIN_MS=1 (2026-10-17, plan idempotent-conjuring-noodle.md Phase 5; user: offset Fig 1 becomes its own manuscript
# figure at 180 mm): the same three panels redrawn at final size in the manuscript style (benchmarks/ms_style.R:
# 7/6 pt Arial, ms_save; badge corner radius 0.25 lines, as 0.5 exceeds half a 6 pt badge), no panel titles (the legend carries them, incl. the tau of c), letters a-c as patchwork tags.
# Writes fig_main_ms.{svg,png} beside the unchanged fig_main. Needs SITE_NAMES and NO_LINEAR (the v3 panels).
if (Sys.getenv("MAIN_MS") == "1") {
    stopifnot(MAP_SITE, NO_LINEAR)
    MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE)
    ms_map <- theme(plot.title = element_blank(), legend.title = element_text(size = MSE$MS_TICK),
                    legend.text = element_text(size = MSE$MS_TICK), legend.key.size = unit(3, "mm"),
                    legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 1, 1, 1, "mm"))
    pa_ms <- map_panel(cons[, .(x, y, lab = site_lab(site), fillv = category, txt = txt_on(CAT_COL[category]))],
                       scale_fill_manual(values = CAT_COL, name = NULL), NULL, badge = MSE$MS_LAB, r = 0.25, pad = 0.18) + ms_map
    pb_ms <- map_panel(nul[, .(x, y, lab, fillv = offset, txt = txt_on(scales::gradient_n_pal(PAL)(scales::rescale(offset))))],
                       scale_fill_gradientn(colours = PAL, name = "Climate shift (SD units)",
                                            guide = guide_colourbar(barwidth = unit(2.4, "cm"), barheight = unit(0.2, "cm"),
                                                                    title.position = "top")),
                       NULL, badge = MSE$MS_LAB, r = 0.25, pad = 0.18) + ms_map
    pc_ms <- delta_panel("bio_3", "rank_rda", base = MSE$MS_BASE, ylab = "RDA offset rank", badge = MSE$MS_LAB,
                         r = 0.25, pad = 0.18) +
        labs(title = NULL) + MSE$theme_ms()
    fm_ms <- pa_ms + pb_ms + pc_ms + plot_layout(widths = c(1, 1, 1.15)) +
        plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 8, colour = MSE$MS_INK))
    MSE$ms_save(file.path(OUT, "fig_main_ms"), fm_ms, 180, 80)
}
fsd <- wrap_plots(lapply(names(FAM), function(f) lapply(P, function(p) delta_panel(p, f))) |> unlist(recursive = FALSE),
                  ncol = 3, byrow = TRUE)
ggsave(file.path(OUT, "fig_s_delta.png"), fsd, width = 11, height = 7.2, dpi = 220, bg = "white")
ggsave(file.path(OUT, "fig_s_delta.svg"), fsd, width = 11, height = 7.2, bg = "white")
# DELTA_CORR=1 (2026-10-17, plan idempotent-conjuring-noodle.md Phase 3; user review v2: "we don't join the ranks from
# linear methods, so here RDA and geometric should be separated; better a heatmap"; then "rows as methods, columns as
# bios", no method-method or bio-bio cells): Kendall tau-b over the 11 sites between each method's offset ranking
# (sign: high = high risk) and each predictor's projected change (same sign convention as chg_*, so bio_3 = its loss).
# Fill = tau (blue negative, red positive, the offset-map ends); bold = exact p < 0.05. Writes fig_s_delta_heat.{png,svg}
# + delta_rank_tau_sep.tsv beside the unchanged fig_s_delta.
if (Sys.getenv("DELTA_CORR") == "1") {
    RM <- c(rank_gf = "Gradient Forest", rank_rda = "RDA (3 axes)", rank_geometric = "Geometric")
    CB <- c(bio_1 = "bio_1\nwarming", bio_3 = "bio_3\nisothermality loss", bio_5 = "bio_5\nwarming")
    cr <- rbindlist(lapply(names(RM), function(m) rbindlist(lapply(names(CB), function(p) {
        a <- cons[[paste0("chg_", p)]]; b <- -cons[[m]]
        data.table(method = sub("rank_", "", m), predictor = p, kendall_tau = tau(a, b),
                   p_exact = suppressWarnings(cor.test(a, b, method = "kendall", exact = TRUE)$p.value))
    }))))
    # the GF row must equal delta_rank_tau.tsv's gf_rank rows (same tau, same sign convention)
    chk <- merge(cr[method == "gf", .(predictor, t1 = kendall_tau)], dr[family == "gf_rank", .(predictor, t2 = kendall_tau)],
                 by = "predictor")
    stopifnot(nrow(chk) == 3, chk[, all(abs(t1 - t2) < 1e-12)])
    fwrite(cr, file.path(OUT, "delta_rank_tau_sep.tsv"), sep = "\t")
    print(cr)
    cr[, `:=`(m = factor(RM[paste0("rank_", method)], rev(RM)), b = factor(CB[predictor], CB))]
    heat_plot <- function(ms = FALSE) ggplot(cr, aes(b, m, fill = kendall_tau)) +
        geom_tile(colour = "white", linewidth = 0.8) +
        geom_text(aes(label = sprintf("%.2f", kendall_tau), fontface = ifelse(p_exact < 0.05, "bold", "plain"),
                      colour = txt_on(scales::gradient_n_pal(c(PAL[1], "white", PAL[5]))(scales::rescale(kendall_tau, from = c(-1, 1))))),
                  size = if (ms) MSE$MS_LAB else 3) +
        scale_colour_identity() +
        scale_fill_gradient2(low = PAL[1], mid = "white", high = PAL[5], midpoint = 0, limits = c(-1, 1),
                             name = "Kendall tau") +
        coord_equal() + labs(x = "Projected change at the site", y = "Offset ranking") +
        theme_minimal(base_size = if (ms) 7 else 8) +
        theme(panel.grid = element_blank(), legend.title = element_text(size = if (ms) 6 else 7),
              legend.key.height = unit(4, "mm"))
    pk <- heat_plot()
    if (SUPP_MS) MSE$ms_save(file.path(OUT, "fig_s_delta_heat_ms"), heat_plot(ms = TRUE) +
                             scale_x_discrete(labels = function(x) sub(" loss", "\nloss", x)) +
                             theme(axis.text = element_text(size = 6, colour = MSE$MS_INK)), 100, 50)
    ggsave(file.path(OUT, "fig_s_delta_heat.png"), pk, width = 4.6, height = 2.6, dpi = 300, bg = "white")
    ggsave(file.path(OUT, "fig_s_delta_heat.svg"), pk, width = 4.6, height = 2.6, bg = "white")
}

# NO_LINEAR=1: consensus-category variants without the linear-family mean (user 2026-10-17: "show me several variants
# and order them by how much they change the current story"). Each method votes with its own rank; upper half = rank < 6,
# lower half = rank > 6 (rank 6 or a tie at 6 is neither). Variants are ordered by the number of sites whose category
# differs from the current family rule. Writes consensus_variants.tsv + fig_consensus_variants.{png,svg}; Fig 1a is not
# changed until the user picks one. Picked 2026-10-17: "all 3 methods agree" (the NO_LINEAR consensus above); the sheet
# is kept behind CONSENSUS_VARIANTS=1 as the record of that choice.
if (Sys.getenv("CONSENSUS_VARIANTS") == "1") {
    R3 <- as.matrix(cons[, .(rank_gf, rank_rda, rank_geometric)])
    nu <- rowSums(R3 < 6); nl <- rowSums(R3 > 6)
    vote <- function(hi, lo) fifelse(hi, "higher risk", fifelse(lo, "lower risk", "methods disagree"))
    cv <- data.table(site = cons$site, x = cons$x, y = cons$y,
        current  = cons$category,
        majority = vote(nu >= 2, nl >= 2),
        all3     = vote(nu == 3, nl == 3),
        gf_rda   = vote(cons$rank_gf < 6 & cons$rank_rda < 6, cons$rank_gf > 6 & cons$rank_rda > 6))
    VT <- c(current = "Current: GF and mean(RDA, geometric) agree",
            majority = "Majority: 2 of 3 methods agree",
            all3 = "Unanimous: all 3 methods agree",
            gf_rda = "GF and RDA agree (geometric left out)")
    nchg <- vapply(names(VT), function(v) sum(cv[[v]] != cv$current), 0L)
    ord <- names(VT)[order(-nchg, names(VT) != "current")]   # most change first, current last
    ord <- c(setdiff(ord, "current"), "current")
    fwrite(cbind(cons[, .(site, rank_gf, rank_rda, rank_geometric)], cv[, names(VT), with = FALSE]),
           file.path(OUT, "consensus_variants.tsv"), sep = "\t")
    print(nchg)
    pv <- lapply(ord, function(v) {
        chg <- cv[[v]] != cv$current
        lab <- if (MAP_SITE) site_lab(cv$site) else cv$site
        d <- data.table(x = cv$x, y = cv$y, lab = ifelse(chg, paste0(lab, "*"), lab), fillv = cv[[v]],
                        txt = txt_on(CAT_COL[cv[[v]]]))
        tab <- table(factor(cv[[v]], names(CAT_COL)))
        map_panel(d, scale_fill_manual(values = CAT_COL, name = NULL, drop = FALSE),
                  sprintf("%s\n%d site(s) change*; higher %d / disagree %d / lower %d", VT[[v]], nchg[[v]],
                          tab[["higher risk"]], tab[["methods disagree"]], tab[["lower risk"]])) +
            theme(plot.title = element_text(face = "bold", size = 8))
    })
    fv <- wrap_plots(pv, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    ggsave(file.path(OUT, "fig_consensus_variants.png"), fv, width = 4 * length(pv), height = 5.6, dpi = 220, bg = "white")
    ggsave(file.path(OUT, "fig_consensus_variants.svg"), fv, width = 4 * length(pv), height = 5.6, bg = "white")
}

# ── supplementary: adaptive order vs the order random 27-SNP sets agree on (Kendall) ──
loo_dist <- function(draws, adaptive) {   # draws: draw, site, offset; adaptive: site, offset
    rk <- dcast(draws[, .(draw, site, r = frank(offset)), by = draw][, .(draw, site, r)], draw ~ site, value.var = "r")
    M <- as.matrix(rk[, -1]); tot <- colSums(M); n <- nrow(M)
    loo <- vapply(seq_len(n), function(i) tau(M[i, ], (tot - M[i, ]) / (n - 1)), 0)
    a <- adaptive[match(colnames(M), site), offset]
    list(loo = loo, adaptive = tau(a, tot / n))
}
gfd <- fread(file.path(B, "site_maps/gf_random_null/site_offsets_draws.tsv"), colClasses = c(site = "character"))[
    kind == "random", .(draw = seed, site, offset)]
lrd <- fread(file.path(B, "site_maps/random_null/site_offsets_draws.tsv"), colClasses = c(site = "character"))
L <- rbindlist(lapply(names(METH), function(m) {
    d <- if (m == "gf") gfd else lrd[engine == m, .(draw, site, offset)]
    x <- loo_dist(d, so[method == m, .(site, offset)])
    data.table(method = m, tau = x$loo, adaptive = x$adaptive)
}))
fwrite(L, file.path(OUT, "random_loo.tsv"), sep = "\t")
Ls <- L[, .(adaptive = adaptive[1], med = median(tau), p = (1 + sum(tau <= adaptive[1])) / (1 + .N)), by = method]
# cross-check against kendall_check.py (same draws, scipy tau-b)
# (first run on a new analysis dir: kendall_check.py needs consensus_sites.tsv from this script, so the check is
# skipped once; run order R -> kendall_check.py -> R)
if (file.exists(file.path(B, "site_maps/kendall_check.tsv"))) {
    kc <- fread(file.path(B, "site_maps/kendall_check.tsv"))[test == "random_consensus" & stat == "kendall"]
    chk <- merge(Ls, kc[, .(method, a2 = adaptive, m2 = random_median, p2 = p)], by = "method")
    stopifnot(chk[, all(abs(adaptive - a2) < 1e-3 & abs(med - m2) < 1e-3 & abs(p - p2) < 1e-6)])
} else message("kendall_check.tsv absent: cross-check skipped; rerun after kendall_check.py")
L[, m := factor(METH[method], METH)]; Ls[, m := factor(METH[method], METH)]
random_plot <- function(ms = FALSE) ggplot(L, aes(tau)) +
    geom_histogram(binwidth = 0.04, boundary = 1, fill = HIST_FILL, colour = "white", linewidth = 0.2) +
    geom_vline(data = Ls, aes(xintercept = adaptive), colour = ADAPT_LINE, linewidth = if (ms) 0.6 else 0.9) +
    geom_text(data = Ls, aes(x = -0.95, y = Inf, label = sprintf("adaptive %.2f\nrandom median %.2f\np = %.3f", adaptive, med, p)),
              hjust = 0, vjust = 1.3, size = if (ms) MSE$MS_LAB else 2.7) +
    facet_wrap(~ m, ncol = 1, scales = "free_y") +
    scale_x_continuous(limits = c(-1.04, 1.04), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Kendall tau with the random-SNP consensus order",
         y = paste0("Random ", nrow(fread(file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv"))), "-SNP sets")) +
    theme_classic(base_size = if (ms) 7 else 9) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0))
ps <- random_plot()
if (SUPP_MS) MSE$ms_save(file.path(OUT, "fig_s_random_ms"), random_plot(ms = TRUE) +
                         theme(axis.text = element_text(size = 6, colour = MSE$MS_INK)), 88, 110)
ggsave(file.path(OUT, "fig_s_random.png"), ps, width = 4.2, height = 5.6, dpi = 250, bg = "white")
ggsave(file.path(OUT, "fig_s_random.svg"), ps, width = 4.2, height = 5.6, bg = "white")
print(Ls)

# ── supplementary: why GF orders sites along bio_5 (user 2026-10-15) ────────
# GF offset = sqrt(sum_p (f_p(future) - f_p(present))^2) with one turnover function f_p per predictor, so the offset
# has no interaction terms; a site's bio_5 component is the rise of f_bio5 between its present and future bio_5,
# which depends on where the site sits on the curve. Left: f_bio5 (solid inside the training range, dashed where
# predict(extrap = TRUE) continues it linearly) with present (circle) and future (triangle) site means. Right: the
# rise per site, sites in GF rank order. Title: bio_5-only GF offset ranking vs the full GF ranking (exact, from
# offset_attribution.R).
suppressPackageStartupMessages({ library(gradientForest); library(qs) })
gf  <- qread(file.path(B, "arm3/intermediate/gradient_forest", paste0(SET, "_nospatial"), "adaptive_model.qs"))
trn <- range(gf$X[["bio_5"]])
fx  <- function(v) {   # f_bio5 at bio_5 = v (other predictors fixed; GF transforms each predictor separately)
    nd <- as.data.frame(lapply(setNames(P, P), function(p) rep(median(gf$X[[p]]), length(v))))
    nd$bio_5 <- v
    predict(gf, nd, extrap = TRUE)[, "bio_5"]
}
cv <- data.table(x = seq(min(cons$bio_5_present) - 0.5, max(clim$bio_5_future) + 0.5, length.out = 300))
cv[, `:=`(f = fx(x), part = fifelse(x < trn[1], "below", fifelse(x <= trn[2], "training range", "above")))]
# segments share their end points so the solid and dashed parts join
cv <- rbind(cv, cv[part == "training range"][c(1, .N)][, part := c("below", "above")])[order(part, x)]
cv[, lt := fifelse(part == "training range", "training range", "extrapolated")]
gs <- merge(cons[, .(site, rank_gf, category, bio_5_present)], clim[, .(site, bio_5_future)], by = "site")
gs[, `:=`(f_pres = fx(bio_5_present), f_fut = fx(bio_5_future))]
pa5 <- fread(file.path(B, "site_maps/attribution/predictor_offsets.tsv"), colClasses = c(site = "character"))[
    method == "gf" & predictor == "bio_5", .(site, off5 = offset)]
tau5 <- merge(pa5, so[method == "gf", .(site, offset)], by = "site")[, tau(off5, offset)]
fwrite(gs[, .(site, rank_gf, category, bio_5_present, bio_5_future, f_present = f_pres, f_future = f_fut,
              f_rise = f_fut - f_pres)][order(rank_gf)], file.path(OUT, "gf_bio5_turnover.tsv"), sep = "\t")
g1 <- ggplot() +
    annotate("rect", xmin = trn[1], xmax = trn[2], ymin = -Inf, ymax = Inf, fill = "grey92") +
    geom_line(data = cv, aes(x, f, linetype = lt, group = part), colour = "grey20", linewidth = 0.6) +
    geom_point(data = gs, aes(bio_5_present, f_pres, fill = category), shape = 21, size = 2.4, colour = "grey15") +
    geom_point(data = gs, aes(bio_5_future, f_fut, fill = category), shape = 24, size = 2.4, colour = "grey15") +
    scale_linetype_manual(values = c("training range" = 1, extrapolated = 2), name = NULL) +
    scale_fill_manual(values = CAT_COL, guide = "none") +
    labs(x = "bio_5, max temp. warmest month (°C)", y = "GF cumulative importance f(bio_5)",
         title = "GF turnover along bio_5 (grey = training range)") +
    theme_classic(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
g2 <- ggplot(gs, aes(y = reorder(site, -rank_gf))) +
    geom_segment(aes(x = f_pres, xend = f_fut, yend = reorder(site, -rank_gf), colour = category), linewidth = 1.1,
                 arrow = arrow(length = unit(0.12, "cm"))) +
    scale_colour_manual(values = CAT_COL, name = NULL) +
    labs(x = "f(bio_5): present → future", y = "Site (GF rank 1 at top)",
         title = sprintf("bio_5-only vs full GF ranking: tau %.2f", tau5)) +
    theme_classic(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
gfb <- g1 + g2 + plot_layout(widths = c(1.3, 1))
ggsave(file.path(OUT, "fig_s_gf_bio5.png"), gfb, width = 10, height = 4.8, dpi = 220, bg = "white")
ggsave(file.path(OUT, "fig_s_gf_bio5.svg"), gfb, width = 10, height = 4.8, bg = "white")
message("tau(bio_5-only GF, full GF) = ", round(tau5, 3))
fwrite(data.table(stat = "kendall_tau_bio5_only_gf_vs_full_gf", value = round(tau5, 3)), file.path(OUT, "gf_bio5_only_tau.tsv"), sep = "\t")

# ── supplementary: full offset maps + extrapolation shares (user 2026-10-15/16) ─
# Offset rasters of the 27-SNP set (the snapshots the site values come from), each method on its own scale. The share
# of cells whose FUTURE climate lies outside the sites' present range (the climate every method trains on, plant level
# as in build_site_maps.R's envelope) goes to tables only, per predictor and as a count of predictors.
TIF <- c(gf = file.path(B, "arm3/intermediate/gradient_forest", paste0(SET, "_nospatial"), SC, "offset_raster.tif"),
         geometric = file.path(B, "arm3/intermediate/geometric_offset", paste0(SET, "_nospatial"), SC, "offset_raster.tif"),
         rda = file.path(B, "arm3_ax3/intermediate/rda_offset", paste0(SET, "_nospatial"), SC, "offset_raster.tif"))
stopifnot(all(file.exists(TIF)))
env_pres <- pres_site[sample %in% samples$sample, ..P]
lo <- sapply(env_pres, min); hi <- sapply(env_pres, max)
fut_all <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_all.tsv")))
fm_ <- as.matrix(fut_all[, ..P])
nout <- rowSums(sweep(fm_, 2, lo, "<") | sweep(fm_, 2, hi, ">"))
rdf <- function(r) setnames(as.data.table(as.data.frame(r, xy = TRUE, na.rm = TRUE)), 3, "v")
ras_panel <- function(d, title, fill_scale, ms = FALSE) {
    ggplot() +
        geom_raster(data = d, aes(x, y, fill = v)) +
        geom_polygon(data = lg, aes(x, y, group = grp), fill = NA, colour = "grey40", linewidth = 0.2) +
        geom_point(data = xy, aes(x, y), shape = 21, size = if (ms) 1 else 1.6, fill = "white", colour = "grey10", stroke = 0.4) +
        fill_scale + coord_equal(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) +
        labs(title = title, x = NULL, y = NULL) + theme_void(base_size = if (ms) 6 else 9) +
        theme(plot.title = element_text(face = "bold", size = if (ms) 7 else 9), legend.position = "bottom",
              legend.key.width = unit(if (ms) 0.8 else 1, "cm"), legend.key.height = unit(if (ms) 0.2 else 0.25, "cm"),
              panel.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))
}
mp <- lapply(names(TIF), function(m) ras_panel(rdf(rast(TIF[[m]])), paste0(c(gf = "a  ", geometric = "b  ", rda = "c  ")[[m]],
                                                                          METH[[m]], " offset"),
                                             scale_fill_gradientn(colours = PAL, name = NULL)))
# the extrapolation map (panel d) was dropped (user 2026-10-16): the shares go into the text, from the two tables below
fom <- wrap_plots(mp, nrow = 1)
if (SUPP_MS) {
    mp_ms <- lapply(names(TIF), function(m) ras_panel(rdf(rast(TIF[[m]])), paste0(c(gf = "a  ", geometric = "b  ", rda = "c  ")[[m]],
                                                                                  METH[[m]], " offset"),
                                                     scale_fill_gradientn(colours = PAL, name = NULL), ms = TRUE))
    MSE$ms_save(file.path(OUT, "fig_s_offset_maps_ms"), wrap_plots(mp_ms, nrow = 1), 180, 78)
}
ggsave(file.path(OUT, "fig_s_offset_maps.png"), fom, width = 10.5, height = 5.4, dpi = 220, bg = "white")
ggsave(file.path(OUT, "fig_s_offset_maps.svg"), fom, width = 10.5, height = 5.4, bg = "white")
ext_cells <- data.table(n_predictors_outside = 0:3)[, cells := sapply(n_predictors_outside, function(k) sum(nout == k))][
    , pct := round(100 * cells / length(nout), 1)]
fwrite(ext_cells, file.path(OUT, "extrapolation_cells.tsv"), sep = "\t")
print(ext_cells)
ext_pred <- rbindlist(lapply(P, function(p) data.table(
    predictor = p, site_min = lo[[p]], site_max = hi[[p]],
    pct_below = round(100 * mean(fm_[, p] < lo[[p]]), 1), pct_above = round(100 * mean(fm_[, p] > hi[[p]]), 1))))
ext_pred[, pct_outside := pct_below + pct_above]
fwrite(ext_pred, file.path(OUT, "extrapolation_by_predictor.tsv"), sep = "\t")
print(ext_pred)

# MAIN6_MS=1 (2026-10-18, plan let-s-implement-this-plan-mighty-bubble.md; manuscript Figure 5): the three offset maps
# of fig_s_offset_maps_ms (row 1, a-c) over the three panels of fig_main_ms (row 2, d-f), one figure at 180 mm, letters
# a-f as patchwork tags. Same builders and data; display labels only ("offset", never "risk"; BIO names), set locally so
# METH / CAT_COL / DLAB / FAM and every earlier output stay as they are. Last block of the script on purpose: it cannot
# shift anything drawn or written above. Writes fig5_offset_ms.{svg,png}. Needs SITE_NAMES and NO_LINEAR (the v3 panels).
if (Sys.getenv("MAIN6_MS") == "1") {
    stopifnot(MAP_SITE, NO_LINEAR)
    MS6 <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MS6)
    MAP_T <- c(gf = "Gradient Forest offset", geometric = "Geometric offset", rda = "RDA-based offset")
    CAT_LAB <- c("higher risk" = "higher offset in all three methods", "lower risk" = "lower offset in all three methods",
                 "methods disagree" = "methods disagree")
    stopifnot(setequal(names(CAT_LAB), names(CAT_COL)), setequal(names(MAP_T), names(TIF)))
    ms_map6 <- theme(plot.title = element_blank(), legend.title = element_text(size = MS6$MS_TICK),
                     legend.text = element_text(size = MS6$MS_TICK), legend.key.size = unit(3, "mm"),
                     legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 1, 1, 1, "mm"))
    row1 <- lapply(names(TIF), function(m) ras_panel(rdf(rast(TIF[[m]])), MAP_T[[m]],
                                                     scale_fill_gradientn(colours = PAL, name = NULL), ms = TRUE))
    pd6 <- map_panel(cons[, .(x, y, lab = site_lab(site), fillv = category, txt = txt_on(CAT_COL[category]))],
                     scale_fill_manual(values = CAT_COL, labels = CAT_LAB, name = NULL,
                                       guide = guide_legend(ncol = 1)), NULL, badge = MS6$MS_LAB, r = 0.25, pad = 0.18) + ms_map6
    pe6 <- map_panel(nul[, .(x, y, lab, fillv = offset, txt = txt_on(scales::gradient_n_pal(PAL)(scales::rescale(offset))))],
                     scale_fill_gradientn(colours = PAL, name = "Climate shift (SD units)",
                                          guide = guide_colourbar(barwidth = unit(2.4, "cm"), barheight = unit(0.2, "cm"),
                                                                  title.position = "top")),
                     NULL, badge = MS6$MS_LAB, r = 0.25, pad = 0.18) + ms_map6
    # panel f: the site badge IS the data point (user 2026-10-18: "labels not as annotation of dots but the dots
    # itself"), centred on (BIO3 loss, RDA rank) with no dot or leader; same data and colours as delta_panel()
    f6 <- cons[, .(x = chg_bio_3, y = rank_rda, lab = site_lab(site), fillv = category, txt = txt_on(CAT_COL[category]))]
    stopifnot(nrow(f6) == 11, !anyNA(f6$x), setequal(f6$y, 1:11))
    pf6 <- ggplot(f6, aes(x, y)) +
        # the map badge call (badge_layer) with repulsion off: force 0 and no iterations leave each badge on its point
        ggrepel::geom_label_repel(aes(label = lab, fill = fillv, colour = txt), size = MS6$MS_LAB, fontface = "bold",
                                  label.r = unit(0.25, "lines"), label.padding = unit(0.18, "lines"), label.size = 0.3,
                                  force = 0, force_pull = 0, max.iter = 0, min.segment.length = Inf, seed = 42,
                                  max.overlaps = Inf) +
        scale_colour_identity() + scale_fill_manual(values = CAT_COL, guide = "none") +
        scale_y_reverse(breaks = c(1, 6, 11), expand = expansion(add = 0.9)) +
        scale_x_continuous(expand = expansion(mult = 0.08)) +
        labs(x = "Loss of isothermality BIO3 (−Δ, %)", y = "RDA offset rank") +
        MS6$theme_ms()
    fig5 <- wrap_plots(row1, nrow = 1) / (pd6 + pe6 + pf6 + plot_layout(widths = c(1, 1, 1.15))) +
        plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold", size = 8, colour = MS6$MS_INK))
    MS6$ms_save(file.path(OUT, "fig5_offset_ms"), fig5, 180, 150)
}
# SUPP2_MS=1 (2026-10-25, plan let-s-follow-the-plan-lively-hellman.md, supplement build B6 + B7; every line of this
# block is added): s_predictor_change (fig_s_delta_heat_ms redrawn: tau on the offsets, Holm-bold, manuscript labels)
# and s_random_sets (fig_s_random_ms with manuscript labels and the in-panel text placed clear of the line and the
# bars). Labels, factors and tables are local copies; L, Ls, cr, cons and every earlier output stay as they are. Last
# block of the script on purpose, as MAIN6_MS was: it cannot shift anything drawn or written above. Writes
# paper_pack/s_predictor_change.{svg,png}, s_predictor_change_cells.tsv, s_random_sets.{svg,png}.
if (Sys.getenv("SUPP2_MS") == "1") {
    stopifnot(MAP_SITE, NO_LINEAR, Sys.getenv("DELTA_CORR") == "1")
    MS2 <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MS2)
    minus <- function(s) gsub("-(?=[0-9])", "−", s, perl = TRUE)
    lab_minus <- scales::label_number(accuracy = 0.1, style_negative = "minus")
    M2 <- c(gf = "Gradient Forest", geometric = "Geometric offset", rda = "RDA-based offset")   # rows / facets, top first
    F2 <- c(bio_1 = "Rise in mean\ntemperature\n(BIO1)", bio_3 = "Loss of\nisothermality\n(BIO3)",
            bio_5 = "Rise in maximum\ntemperature\n(BIO5)")
    stopifnot(setequal(names(M2), names(METH)), identical(names(F2), P))
    # ── B6: Kendall's tau-b between each method's site offsets and each projected shift (chg_*: BIO3 = its loss) ──
    cells <- rbindlist(lapply(names(M2), function(m) rbindlist(lapply(P, function(p) {
        x <- merge(so[method == m, .(site, offset)], cons[, .(site, shift = get(paste0("chg_", p)))], by = "site")
        stopifnot(nrow(x) == 11, !anyNA(x$offset), !anyNA(x$shift))
        warn <- character(0)
        ct <- withCallingHandlers(cor.test(x$offset, x$shift, method = "kendall", exact = TRUE),
                                  warning = function(w) { warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning") })
        tied <- anyDuplicated(x$offset) > 0 || anyDuplicated(x$shift) > 0
        stopifnot(identical(tied, length(warn) > 0), all(grepl("ties", warn)))   # exact p unless a tie forced the approximation
        data.table(method = m, predictor = p, kendall_tau = unname(ct$estimate), p = ct$p.value,
                   p_method = if (tied) "normal approximation" else "exact", warning = paste(warn, collapse = "; "))
    }))))
    cells[, p_holm := p.adjust(p, "holm")][, bold := p_holm < 0.05]
    # tau on the offset equals delta_rank_tau_sep.tsv (tau against -rank) in all nine cells
    sep <- fread(file.path(OUT, "delta_rank_tau_sep.tsv"))
    chk6 <- merge(cells[, .(method, predictor, t1 = kendall_tau)], sep[, .(method, predictor, t2 = kendall_tau)],
                  by = c("method", "predictor"))
    stopifnot(nrow(chk6) == 9, chk6[, all(abs(t1 - t2) < 1e-12)],
              cells[predictor == "bio_5", all(p_method == "normal approximation")],
              cells[predictor != "bio_5", all(p_method == "exact")],
              sum(cells$bold) == 1, cells[bold == TRUE, method == "rda" & predictor == "bio_3"],
              cells[bold == TRUE, abs(kendall_tau - 0.927) < 5e-4 & abs(signif(p, 3) - 3.26e-06) < 1e-12 &
                                  abs(p_holm - 2.9e-05) < 0.05e-05])
    for (p in P) {   # ties and float-noise near-ties of the shift (reported, no value changed)
        v <- setNames(cons[[paste0("chg_", p)]], cons$site); pr <- combn(names(v), 2)
        dd <- abs(v[pr[1, ]] - v[pr[2, ]])
        for (k in which(dd < 1e-5)) message(sprintf("s_predictor_change: %s shift of sites %s and %s differs by %.3g%s",
                                                    p, pr[1, k], pr[2, k], dd[k], if (dd[k] == 0) " (exact tie)" else ""))
    }
    fwrite(cells[, .(method = M2[method], factor = sub("bio_", "BIO", predictor), kendall_tau, p, p_method, p_holm, bold)],
           file.path(OUT, "s_predictor_change_cells.tsv"), sep = "\t")
    print(cells)
    hc <- copy(cells)[, `:=`(m = factor(M2[method], rev(M2)), b = factor(F2[predictor], F2))]
    hc[, `:=`(lab = minus(sprintf("%.2f", kendall_tau)),
              txt = txt_on(scales::gradient_n_pal(c(PAL[1], "white", PAL[5]))(scales::rescale(kendall_tau, from = c(-1, 1)))))]
    p6 <- ggplot(hc, aes(b, m, fill = kendall_tau)) +
        geom_tile(colour = "white", linewidth = 0.8) +
        geom_text(aes(label = lab, fontface = ifelse(bold, "bold", "plain"), colour = txt), size = MS2$MS_LAB) +
        scale_colour_identity() +
        scale_fill_gradient2(low = PAL[1], mid = "white", high = PAL[5], midpoint = 0, limits = c(-1, 1),
                             labels = lab_minus, name = "Kendall's τ") +
        coord_equal() + labs(x = "Projected shift at the site", y = NULL) +
        theme_minimal(base_size = 7) +
        theme(panel.grid = element_blank(), legend.title = element_text(size = 6), legend.key.height = unit(4, "mm"),
              axis.text = element_text(size = 6, colour = MS2$MS_INK), text = element_text(colour = MS2$MS_INK))
    MS2$ms_save(file.path(OUT, "s_predictor_change"), p6, 100, 55)
    # ── B7: the associated SNPs' tau against the leave-one-out taus of the random sets ──
    kc7 <- fread(file.path(B, "site_maps/kendall_check.tsv"))[test == "random_consensus" & stat == "kendall"]
    stopifnot(identical(kc7$method, c("gf", "geometric", "rda")),
              all(abs(Ls[match(kc7$method, method), adaptive] - c(0.164, 0.455, -0.491)) < 5e-4),
              all(abs(Ls[match(kc7$method, method), med] - c(0.734, 0.818, 0.927)) < 5e-4),
              all(abs(Ls[match(kc7$method, method), p] - c(0.015, 0.005, 0.005)) < 5e-4))
    n_set <- nrow(fread(file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv")))
    stopifnot(n_set == 30)
    L7 <- copy(L)[, m := factor(M2[method], M2)]
    Ls7 <- copy(Ls)[, m := factor(M2[method], M2)]
    Ls7[, lab := minus(sprintf("associated SNPs %.2f\nrandom median %.2f\np = %.3f", adaptive, med, p))]
    # Text placement: top left at x = -0.95 (as fig_s_random) unless the red line falls inside the text's span
    # (TXT_W, the width of the longest line at 6 pt in tau units of this panel); then just right of the line.
    TXT_W <- 0.62
    Ls7[, `:=`(tx = fifelse(adaptive > -0.95 - 0.04 & adaptive < -0.95 + TXT_W + 0.04, adaptive + 0.05, -0.95))]
    Ls7[, bar_min := L7[, .(v = min(tau)), by = method][match(Ls7$method, method), v]]
    message("s_random_sets text: ", paste(sprintf("%s x = %.2f (line %.2f, lowest random tau %.2f)", Ls7$method, Ls7$tx,
                                                  Ls7$adaptive, Ls7$bar_min), collapse = "; "))
    p7 <- ggplot(L7, aes(tau)) +
        geom_histogram(binwidth = 0.04, boundary = 1, fill = HIST_FILL, colour = "white", linewidth = 0.2) +
        geom_vline(data = Ls7, aes(xintercept = adaptive), colour = ADAPT_LINE, linewidth = 0.6) +
        geom_text(data = Ls7, aes(x = tx, y = Inf, label = lab), hjust = 0, vjust = 1.3, size = MS2$MS_LAB) +
        facet_wrap(~ m, ncol = 1, scales = "free_y") +
        scale_x_continuous(limits = c(-1.04, 1.04), breaks = seq(-1, 1, 0.5), labels = lab_minus) +
        labs(x = "Kendall's τ with the consensus of the random sets", y = paste0("Random ", n_set, "-SNP sets")) +
        theme_classic(base_size = 7) +
        theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
              axis.text = element_text(size = 6, colour = MS2$MS_INK))
    MS2$ms_save(file.path(OUT, "s_random_sets"), p7, 88, 110)
}
message("done")
