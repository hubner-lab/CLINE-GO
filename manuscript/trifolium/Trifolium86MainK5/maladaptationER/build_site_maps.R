#!/usr/bin/env Rscript
# Site-only offset maps for Trifolium86MainK5ER (plan ~/.claude/plans/peaceful-kindling-micali.md, Phase 1).
# The 11 sampling sites as points coloured by their site-mean genetic offset on a white map (land = the
# present climate raster's non-NA mask, no offset raster). Run chosen by the user 2026-10-15:
# predictors bio_1,bio_3,bio_5; GF adaptive + geometric from arm3/, GF random from arm3/gf_random/,
# RDA offset with axes fixed at 3 from arm3_ax3/. SSP585 2061-2080, mean of 3 GCMs (pipeline default).
#   rows = method, cols = SNP set (27 / 86 / 205); fig_site_grid = own colour scale per panel,
#   fig_site_grid_shared = one scale per method row (native units);
#   the climate-only null (scaled Euclidean present -> future distance on the same 3 predictors) has its own.
#   Point label = site risk rank, 1 = highest offset in that panel.
# Outputs (site_maps/): fig_site_grid.{png,svg}, site_offsets.tsv, envelope.tsv
# Run inside the image from /pipeline:
#   Rscript ANALYSIS/Trifolium86MainK5/maladaptationER/build_site_maps.R
suppressPackageStartupMessages({ library(data.table); library(terra); library(ggplot2); library(patchwork) })

B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R   <- file.path("/pipeline", paste0(Sys.getenv("PROJ", "Trifolium86MainK5ER"), "_results"))
SC  <- "year2061-2080_ssp585"
OUT <- file.path(B, "site_maps"); dir.create(OUT, showWarnings = FALSE)
P   <- c("bio_1", "bio_3", "bio_5")
# H. spontaneum paper offset palette (color-hex 111971), as in F13/F16
PAL <- c("#095786", "#0090ab", "#e6cea0", "#db4646", "#ae0000")
POINT_COL <- "#095786"
# PALETTE=minou (2026-10-16): rank scatter points in minou navy; the gradient PAL of the maps stays (user).
if (Sys.getenv("PALETTE") == "minou") POINT_COL <- "#2e4057"
# SITE_NAMES=<site_names.tsv> (2026-10-16, plan fancy-scribbling-garden.md Phase 5): rank-plot site labels show the
# short names 1..11 (S->N); tables keep the raw ids. Unset = raw ids, unchanged.
site_lab <- identity
if (nzchar(Sys.getenv("SITE_NAMES"))) {
    SN <- fread(Sys.getenv("SITE_NAMES"), colClasses = "character")
    site_lab <- function(s) { m <- SN$short[match(s, SN$site_id)]; stopifnot(!anyNA(m)); m }
}
# With SITE_NAMES the map badges carry the site name, not the rank (user 2026-10-16, option B: a number on a point
# is always the site; rank = colour and the rank plots' axes). Unset = rank, unchanged.
map_lab <- if (nzchar(Sys.getenv("SITE_NAMES"))) function(site, rank) site_lab(site) else function(site, rank) rank
# With SITE_NAMES the badge text is white or ink per badge, whichever contrasts more with its fill (2026-10-17, plan
# idempotent-conjuring-noodle.md Phase 3: dark digits on navy / dark red were unreadable). Unset = grey10, unchanged.
txt_on <- function(bg, ink = "grey10") {
    lum <- function(x) { v <- grDevices::col2rgb(x) / 255; v <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
                         colSums(v * c(0.2126, 0.7152, 0.0722)) }
    l <- lum(bg); ifelse(1.05 / (l + 0.05) >= (l + 0.05) / (lum(ink) + 0.05), "white", ink)
}

# labels from the set files (v1: 27 / 86 / 205; the PC5 rebuild changes them)
SETS <- sapply(c("gea.ER_overlap", "gea.ER_rda1e4", "gea.ER_rda1e3"), function(s)
    paste(nrow(fread(file.path(R, "_intermediate/snp_sets", s, "selected_snps.tsv"))), "SNPs"))
METHODS <- c(gf = "GF adaptive", gf_random = "GF random", geometric = "Geometric", rda = "RDA offset (3 axes)")
site_table <- function(m, set) {
    s <- paste0(set, "_nospatial")
    switch(m,
        gf        = file.path(B, "arm3/Maladaptation/tables/gradient_forest", s, SC, "genetic_offset_site.tsv"),
        gf_random = file.path(B, "arm3/gf_random", set, "genetic_offset_site.tsv"),
        geometric = file.path(B, "arm3/Maladaptation/tables/geometric_offset", s, SC, "genetic_offset_site.tsv"),
        rda       = file.path(B, "arm3_ax3/Maladaptation/tables/rda_offset", s, SC, "genetic_offset_site.tsv"))
}

samples <- fread(file.path(R, "_intermediate/samples/metadata_climate_valid.tsv"),
                 colClasses = c(site = "character", sample = "character"))
xy <- samples[, .(x = mean(longitude), y = mean(latitude)), by = site]

# ── offsets: site mean per method x set, plus the climate-only null ─────────
long <- list()
for (m in names(METHODS)) for (set in names(SETS)) {
    f <- site_table(m, set); stopifnot(file.exists(f))
    long[[paste(m, set)]] <- fread(f, colClasses = c(site = "character", sample = "character"))[
        , .(offset = mean(genetic_offset)), by = site][, `:=`(method = m, set = set)]
}
pres_site <- fread(file.path(R, "climate/tables/present/climate_present_site.tsv"), colClasses = c(sample = "character"))
fut_site  <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_site.tsv")),
                   colClasses = c(sample = "character"))
fs  <- fut_site[match(pres_site$sample, fut_site$sample)]
sdv <- sapply(pres_site[, ..P], sd)
dn  <- sqrt(rowSums(((as.matrix(fs[, ..P]) - as.matrix(pres_site[, ..P])) / rep(sdv, each = nrow(pres_site)))^2))
long[["null"]] <- data.table(sample = pres_site$sample, offset = dn)[samples, on = "sample"][
    , .(offset = mean(offset)), by = site][, `:=`(method = "null", set = "climate")]
D <- rbindlist(long, use.names = TRUE)
stopifnot(D[, .N, by = .(method, set)][, all(N == 11)])
D[, rank := frank(-offset, ties.method = "min"), by = .(method, set)]
D <- merge(D, xy, by = "site")

# ── training envelope: the sites' present climate (what every engine trains on) ──
env <- pres_site[sample %in% samples$sample, ..P]
lo <- sapply(env, min); hi <- sapply(env, max)
out_any <- function(tab) {
    m <- as.matrix(tab[, ..P]); ok <- stats::complete.cases(m)
    o <- sweep(m[ok, , drop = FALSE], 2, lo, "<") | sweep(m[ok, , drop = FALSE], 2, hi, ">")
    list(any = 100 * mean(rowSums(o) > 0), by = 100 * colMeans(o))
}
pres_all <- fread(file.path(R, "climate/tables/present/climate_present_all.tsv"))
fut_all  <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_all.tsv")))
ep <- out_any(pres_all); ef <- out_any(fut_all)
fsite <- fs[sample %in% samples$sample]
site_out <- data.table(sample = fsite$sample,
                       n_pred_outside = rowSums(sweep(as.matrix(fsite[, ..P]), 2, lo, "<") |
                                                sweep(as.matrix(fsite[, ..P]), 2, hi, ">")))[
    samples, on = "sample"][, .(future_preds_outside = max(n_pred_outside)), by = site]
envelope <- data.table(scope = c("present cells", "future cells", "future sites"),
                       pct_outside_any = round(c(ep$any, ef$any, 100 * mean(site_out$future_preds_outside > 0)), 1),
                       bio_1 = round(c(ep$by["bio_1"], ef$by["bio_1"], NA), 1),
                       bio_3 = round(c(ep$by["bio_3"], ef$by["bio_3"], NA), 1),
                       bio_5 = round(c(ep$by["bio_5"], ef$by["bio_5"], NA), 1))
fwrite(envelope, file.path(OUT, "envelope.tsv"), sep = "\t")
print(envelope)
D <- merge(D, site_out, by = "site")
fwrite(D[order(method, set, rank)], file.path(OUT, "site_offsets.tsv"), sep = "\t")

# ── land backdrop from the present raster's non-NA mask ─────────────────────
r0   <- rast(file.path(R, "climate/rasters/present/climate_present_rasterstack.tif"))[[1]]
land <- as.polygons(!is.na(r0)); land <- land[land[[1]][, 1] == 1]
lg   <- as.data.frame(geom(land)); lg$grp <- paste(lg$geom, lg$part, lg$hole)
ex   <- as.vector(ext(r0))

panel <- function(d, title, lims, legend = TRUE, ms = FALSE) {
    sz <- if (ms) list(b = MSE$MS_LAB, r = 0.25, pad = 0.18, base = 7, lt = 6, kh = 0.6, kw = 0.2) else
                  list(b = 2.4, r = 0.5, pad = 0.22, base = 9, lt = 7, kh = 1.1, kw = 0.3)
    ggplot() +
        geom_polygon(data = lg, aes(x, y, group = grp), fill = "grey94", colour = "grey60", linewidth = 0.25) +
        # true location = small dot; the coloured badge (rank) is repelled with a leader line so the
        # close pairs (29010/29024 are 0.02 deg apart) never hide each other
        geom_point(data = d, aes(x, y), size = 0.9, colour = "grey10") +
        (if (nzchar(Sys.getenv("SITE_NAMES"))) list(
        ggrepel::geom_label_repel(data = transform(d, txt = txt_on(scales::gradient_n_pal(PAL)(scales::rescale(offset, from = lims)))),
                                  aes(x, y, label = map_lab(site, rank), fill = offset, colour = txt), size = sz$b,
                                  fontface = "bold", label.r = unit(sz$r, "lines"),
                                  label.padding = unit(sz$pad, "lines"), label.size = 0.3,
                                  box.padding = 0.15, point.padding = 0.1, min.segment.length = 0,
                                  segment.colour = "grey30", segment.size = 0.3, seed = 42,
                                  max.overlaps = Inf), scale_colour_identity()) else
        ggrepel::geom_label_repel(data = d, aes(x, y, label = map_lab(site, rank), fill = offset), size = 2.4,
                                  fontface = "bold", colour = "grey10", label.r = unit(0.5, "lines"),
                                  label.padding = unit(0.22, "lines"), label.size = 0.3,
                                  box.padding = 0.15, point.padding = 0.1, min.segment.length = 0,
                                  segment.colour = "grey30", segment.size = 0.3, seed = 42,
                                  max.overlaps = Inf)) +
        scale_fill_gradientn(colours = PAL, limits = lims, name = NULL) +
        coord_equal(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = sz$base) +
        theme(plot.title = element_text(size = sz$base, face = "bold", hjust = 0.5, margin = margin(b = 3)),
              panel.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
              plot.background = element_rect(fill = "white", colour = NA),
              legend.key.height = unit(sz$kh, "cm"), legend.key.width = unit(sz$kw, "cm"),
              legend.text = element_text(size = sz$lt),
              legend.position = if (legend) "right" else "none")
}
# Two versions (user 2026-10-15): fig_site_grid = own scale + legend per panel, to compare the spatial
# PATTERN of each map; fig_site_grid_shared = one scale per method row, to compare magnitudes across sets.
build_grid <- function(shared) {
    rows <- lapply(names(METHODS), function(m) {
        ps <- lapply(names(SETS), function(s) {
            d <- D[method == m & set == s]
            panel(d, paste(METHODS[[m]], SETS[[s]], sep = " · "),
                  if (shared) range(D[method == m, offset]) else range(d$offset),
                  legend = !shared || s == tail(names(SETS), 1))
        })
        ps[[4]] <- if (m == "gf") panel(D[method == "null"], "Climate-only null", range(D[method == "null", offset])) else plot_spacer()
        ps
    })
    wrap_plots(unlist(rows, recursive = FALSE), ncol = 4, byrow = TRUE)
}
for (v in c(own = FALSE, shared = TRUE) |> names()) {
    fig <- build_grid(v == "shared")
    f <- if (v == "own") "fig_site_grid" else "fig_site_grid_shared"
    ggsave(file.path(OUT, paste0(f, ".png")), fig, width = 15, height = 4.6 * length(METHODS), dpi = 200, bg = "white")
    ggsave(file.path(OUT, paste0(f, ".svg")), fig, width = 15, height = 4.6 * length(METHODS), bg = "white")
}
# ── Figure: 27-SNP set, the 3 engines vs the climate-only null (user 2026-10-15) ──
# Top: site maps on own scales. Bottom: site risk rank of each engine vs the null's rank (1 = highest offset),
# Kendall tau-b (user 2026-10-15: Kendall throughout the paper; 55 site pairs) with its exact two-sided p
# (n = 11 sites, no ties). p treats the 11 sites as exchangeable; neighbouring sites share climate, so it is optimistic
# if anything. The figure shows tau only: the climate null is a visual reference, the test is against random SNPs
# (kendall_check.py).
CMP <- c(gf = "GF adaptive", geometric = "Geometric", rda = "RDA offset (3 axes)")
S27 <- "gea.ER_overlap"
nl  <- D[method == "null", .(site, null_rank = rank, null = offset)]
tests <- rbindlist(lapply(c(names(CMP), "gf_random"), function(m) {
    x <- merge(D[method == m & set == S27, .(site, rank, offset)], nl, by = "site")
    ct <- suppressWarnings(cor.test(x$offset, x$null, method = "kendall", exact = TRUE))
    data.table(method = m, n_sites = nrow(x), kendall_tau = round(unname(ct$estimate), 3),
               p_exact = signif(ct$p.value, 3),
               top3_shared_with_null = length(intersect(x[rank <= 3, site], x[null_rank <= 3, site])))
}))
fwrite(tests, file.path(OUT, "null_compare_27.tsv"), sep = "\t")
print(tests)
maps <- c(list(panel(D[method == "null"], "Climate-only null", range(D[method == "null", offset]))),
          lapply(names(CMP), function(m) { d <- D[method == m & set == S27]
              panel(d, paste0(CMP[[m]], " · ", SETS[[S27]]), range(d$offset)) }))
rankplot <- function(m, ms = FALSE) {
    sz <- if (ms) list(pt = 1.2, lab = MSE$MS_LAB, base = 7) else list(pt = 2, lab = 2.3, base = 9)
    x <- merge(D[method == m & set == S27, .(site, rank)], nl, by = "site")
    t <- tests[method == m]
    ggplot(x, aes(null_rank, rank)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
        geom_point(size = sz$pt, colour = POINT_COL) +
        ggrepel::geom_text_repel(aes(label = site_lab(site)), size = sz$lab, colour = "grey30", seed = 42,
                                 min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf) +
        scale_x_continuous(breaks = c(1, 6, 11), limits = c(0.5, 11.5)) +
        scale_y_continuous(breaks = c(1, 6, 11), limits = c(0.5, 11.5)) +
        coord_equal() +
        labs(x = "Climate-only null rank", y = paste(CMP[[m]], "rank"),
             title = sprintf("Kendall tau = %.2f", t$kendall_tau)) +
        theme_classic(base_size = sz$base) +
        theme(plot.title = element_text(size = sz$base, face = "bold", hjust = 0.5))
}
fig2 <- wrap_plots(c(maps, list(plot_spacer()), lapply(names(CMP), rankplot)), ncol = 4, byrow = TRUE,
                   heights = c(1.15, 1))
ggsave(file.path(OUT, "fig_null_compare_27.png"), fig2, width = 15, height = 8.6, dpi = 200, bg = "white")
ggsave(file.path(OUT, "fig_null_compare_27.svg"), fig2, width = 15, height = 8.6, bg = "white")
# SUPP_MS=1 (2026-10-17, plan idempotent-conjuring-noodle.md Phase 5 follow-up; user: "redraw S15-S19 at final size"):
# the same figure at 180 mm in the manuscript style (benchmarks/ms_style.R ms_save: 7/6 pt Arial, badge radius 0.25
# lines). Writes fig_null_compare_27_ms.{svg,png} beside the unchanged figure; panel(ms = FALSE) / rankplot(ms = FALSE)
# keep the default sizes.
if (Sys.getenv("SUPP_MS") == "1") {
    MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE)
    maps_ms <- c(list(panel(D[method == "null"], "Climate-only null", range(D[method == "null", offset]), ms = TRUE)),
                 lapply(names(CMP), function(m) { d <- D[method == m & set == S27]
                     panel(d, CMP[[m]], range(d$offset), ms = TRUE) }))
    fig2_ms <- wrap_plots(c(maps_ms, list(plot_spacer()), lapply(names(CMP), rankplot, ms = TRUE)), ncol = 4,
                          byrow = TRUE, heights = c(1.15, 1)) & theme(plot.margin = margin(1, 1, 1, 1, "mm"))
    MSE$ms_save(file.path(OUT, "fig_null_compare_27_ms"), fig2_ms, 180, 82)
}
# SUPP2_MS=1 (2026-10-25, plan let-s-follow-the-plan-lively-hellman.md, supplement build B5; every line of this block is
# added): s_site_rankings = fig_null_compare_27_ms with the manuscript's display labels only, and tau formatted once from
# the unrounded estimate (null_compare_27.tsv stores round(tau, 3), which "%.2f" turned into 0.34 for 19/55 = 0.3455).
# panel(ms = TRUE) / rankplot(ms = TRUE) read the global MSE, so this block builds it itself and does not depend on
# SUPP_MS. null_compare_27.tsv is unchanged; writes site_maps/s_site_rankings.{svg,png} + s_site_rankings_tau.tsv.
if (Sys.getenv("SUPP2_MS") == "1") {
    stopifnot(nzchar(Sys.getenv("SITE_NAMES")))
    MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE)
    CMP2 <- c(gf = "Gradient Forest", geometric = "Geometric offset", rda = "RDA-based offset")
    YL2  <- c(gf = "Gradient Forest offset rank", geometric = "Geometric offset rank", rda = "RDA-based offset rank")
    stopifnot(identical(names(CMP2), names(CMP)), identical(names(YL2), names(CMP)))
    minus <- function(s) gsub("-(?=[0-9])", "−", s, perl = TRUE)
    tau2 <- rbindlist(lapply(names(CMP2), function(m) {
        x <- merge(D[method == m & set == S27, .(site, offset)], nl, by = "site")
        stopifnot(nrow(x) == 11, !anyDuplicated(x$offset), !anyDuplicated(x$null))   # no ties: the p is exact
        ct <- withCallingHandlers(cor.test(x$offset, x$null, method = "kendall", exact = TRUE),
                                  warning = function(w) stop("cor.test warned: ", conditionMessage(w)))
        data.table(method = m, label = CMP2[[m]], n_sites = nrow(x), n_pairs = choose(nrow(x), 2),
                   kendall_tau = unname(ct$estimate), p_exact = ct$p.value)
    }))
    stopifnot(all(abs(tau2$kendall_tau * 55 - c(19, 1, -17)) < 1e-9),
              all(abs(tau2$kendall_tau - tests[match(tau2$method, method), kendall_tau]) < 5e-4))
    tau2[, title := minus(sprintf("Kendall's τ = %.2f", kendall_tau))]
    stopifnot(identical(tau2$title, c("Kendall's τ = 0.35", "Kendall's τ = 0.02", "Kendall's τ = \u{2212}0.31")))
    fwrite(tau2[, .(method, label, n_sites, n_pairs, kendall_tau, p_exact)], file.path(OUT, "s_site_rankings_tau.tsv"),
           sep = "\t")
    print(tau2)
    maps2 <- c(list(panel(D[method == "null"], "Projected climate shift", range(D[method == "null", offset]), ms = TRUE)),
               lapply(names(CMP2), function(m) { d <- D[method == m & set == S27]
                   panel(d, CMP2[[m]], range(d$offset), ms = TRUE) }))
    ranks2 <- lapply(names(CMP2), function(m) rankplot(m, ms = TRUE) +
        labs(x = "Projected climate shift rank", y = YL2[[m]], title = tau2[method == m, title]))
    fig2_s2 <- wrap_plots(c(maps2, list(plot_spacer()), ranks2), ncol = 4, byrow = TRUE, heights = c(1.15, 1)) &
        theme(plot.margin = margin(1, 1, 1, 1, "mm"))
    MSE$ms_save(file.path(OUT, "s_site_rankings"), fig2_s2, 180, 82)
}
message("done")
