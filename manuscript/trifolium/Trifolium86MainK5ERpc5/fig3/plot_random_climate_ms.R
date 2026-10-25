#!/usr/bin/env Rscript
# s_random_climate (2026-10-25, plan let-s-follow-the-plan-lively-hellman.md, supplement build B8): the random SNP sets
# against the projected climate shift, laid out as s_site_rankings (fig_null_compare_27_ms). Adapted from the reference
# implementation drawn on the other machine on 2026-10-25 (plan Appendix); changes from it:
#   - paths: AN/site_maps (site_offsets.tsv, random_null/ and gf_random_null/ site_offsets_draws.tsv, kendall_check.tsv),
#     fig3/site_names.tsv, the run's present raster; output AN/site_maps;
#   - the three method maps are coloured by the CONSENSUS rank of the random sets: frank(mean_rank, ties = "min") per
#     engine, i.e. each site's place when the sites are ordered by their mean offset rank over the 199 sets (1 = highest);
#   - y title "Median offset rank over the random sets" (two lines); a rank-plot label goes left of its point when
#     another point sits within 1.5 ranks to the right and less than 1 rank above or below; save stem s_random_climate;
#   - random_climate_site_ranks.tsv gains the consensus_rank column; the spec's printed values are asserted too.
# Kendall tau-b per random set between its site offsets and the climate shift; the medians, the associated SNPs' tau
# and the one-sided p are asserted against kendall_check.tsv (climate_null, kendall) before drawing. The rho_* columns
# of the draws and summary tables are Spearman and are never read. The Gradient Forest draws are numbered seed - 100000
# only to align draw ids; each engine is scored on its own, nothing is paired across engines.
# Run inside the image (AN required; never the approved pack folder):
#   AN=ANALYSIS/Trifolium86MainK5ERpc5/fig3/offset_minou Rscript ANALYSIS/Trifolium86MainK5ERpc5/fig3/plot_random_climate_ms.R
# Writes AN/site_maps/s_random_climate.{svg,png}, random_climate_tau.tsv, random_climate_site_ranks.tsv.
suppressPackageStartupMessages({ library(data.table); library(terra); library(ggplot2); library(patchwork) })
stopifnot(nzchar(Sys.getenv("AN")))
W   <- file.path("/pipeline", Sys.getenv("AN"), "site_maps")
OUT <- W
SNF <- "/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/site_names.tsv"
RAS <- "/pipeline/Trifolium86MainK5ERpc5_results/climate/rasters/present/climate_present_rasterstack.tif"
MSE <- new.env(); sys.source("/pipeline/benchmarks/ms_style.R", envir = MSE)
PAL <- c("#095786", "#0090ab", "#e6cea0", "#db4646", "#ae0000")   # the s_site_rankings map gradient
POINT_COL <- MSE$MINOU[["navy"]]
METH <- c(gf = "Gradient Forest", geometric = "Geometric offset", rda = "RDA-based offset")
tau  <- function(a, b) cor(a, b, method = "kendall")
minus <- function(s) gsub("-(?=[0-9])", "−", s, perl = TRUE)

SN <- fread(SNF, colClasses = "character")
site_lab <- function(s) { m <- SN$short[match(s, SN$site_id)]; stopifnot(!anyNA(m)); m }

so   <- fread(file.path(W, "site_offsets.tsv"), colClasses = c(site = "character"))
null <- so[method == "null", .(site, shift = offset, x, y)]
null[, shift_rank := frank(-shift, ties.method = "min")]
stopifnot(nrow(null) == 11, !anyDuplicated(null$shift))
assoc <- so[set == "gea.ER_overlap" & method %in% names(METH), .(site, engine = method, offset)]
stopifnot(assoc[, .N, by = engine][, all(N == 11)], uniqueN(assoc$engine) == 3)

rnd <- fread(file.path(W, "random_null/site_offsets_draws.tsv"), colClasses = c(site = "character"))
gfr <- fread(file.path(W, "gf_random_null/site_offsets_draws.tsv"), colClasses = c(site = "character"))[
    kind == "random", .(site, offset, engine = "gf", draw = seed - 100000L)]
dr <- rbind(rnd[engine %in% c("geometric", "rda"), .(site, offset, engine, draw)], gfr)
stopifnot(dr[, .N, by = .(engine, draw)][, all(N == 11)],
          dr[, uniqueN(draw), by = engine][, all(V1 == 199)])
dr <- merge(dr, null[, .(site, shift)], by = "site")

# ---- tau per random set, the associated SNPs, the one-sided p, checked against the lab's numbers ----
td <- dr[, .(tau = tau(offset, shift)), by = .(engine, draw)]
ta <- merge(assoc, null[, .(site, shift)], by = "site")[, .(tau_assoc = tau(offset, shift)), by = engine]
sm <- td[, .(med = median(tau), q025 = quantile(tau, 0.025), q975 = quantile(tau, 0.975), n = .N), by = engine]
sm <- merge(sm, ta, by = "engine")
sm[, p := vapply(seq_len(.N), function(i) (1 + sum(td[engine == sm$engine[i], tau] <= sm$tau_assoc[i])) / (1 + n[i]), 0)]
kc <- fread(file.path(W, "kendall_check.tsv"))[test == "climate_null" & stat == "kendall",
                                                .(engine = method, a2 = adaptive, m2 = random_median, l2 = q025, h2 = q975, p2 = p)]
chk <- merge(sm, kc, by = "engine")
print(chk)
stopifnot(nrow(chk) == 3, chk[, all(abs(med - m2) < 1e-3 & abs(tau_assoc - a2) < 1e-3 & abs(p - p2) < 1e-6)])
s3 <- sm[match(names(METH), engine)]   # the spec's values, Gradient Forest, geometric, RDA
stopifnot(all(abs(s3$med - c(0.636, 0.527, 0.782)) < 5e-4), all(abs(s3$tau_assoc - c(0.345, 0.018, -0.309)) < 5e-4),
          all(abs(s3$p - c(0.015, 0.005, 0.005)) < 5e-4))
fwrite(chk, file.path(OUT, "random_climate_tau.tsv"), sep = "\t")

# ---- each site's rank (1 = highest offset) in every random set; consensus rank of the 199 sets ----
dr[, rank := frank(-offset, ties.method = "average"), by = .(engine, draw)]
rs <- dr[, .(med_rank = median(rank), lo = quantile(rank, 0.025), hi = quantile(rank, 0.975), mean_rank = mean(rank)),
         by = .(engine, site)]
rs[, consensus_rank := frank(mean_rank, ties.method = "min"), by = engine]
rs <- merge(rs, null[, .(site, shift_rank, x, y)], by = "site")
for (e in names(METH)) {   # report any tie of the consensus rank
    tt <- rs[engine == e][, .N, by = consensus_rank][N > 1]
    message("s_random_climate consensus-rank ties, ", e, ": ",
            if (nrow(tt)) paste(sprintf("rank %d x %d", tt$consensus_rank, tt$N), collapse = ", ") else "none")
}
fwrite(rs[order(engine, shift_rank)], file.path(OUT, "random_climate_site_ranks.tsv"), sep = "\t")

# ---- maps (land = the present raster's non-NA mask, as in s_site_rankings) ----
r0   <- rast(RAS)[[1]]
land <- as.polygons(!is.na(r0)); land <- land[land[[1]][, 1] == 1]
lg   <- as.data.frame(geom(land)); lg$grp <- paste(lg$geom, lg$part, lg$hole)
ex   <- as.vector(ext(r0))
map_panel <- function(d, fill_scale, title, cols) {
    ggplot() +
        geom_polygon(data = lg, aes(x, y, group = grp), fill = "grey94", colour = "grey60", linewidth = 0.25) +
        geom_point(data = d, aes(x, y), size = 0.6, colour = "grey10") +
        ggrepel::geom_label_repel(data = transform(d, txt = MSE$ms_text_on(cols)),
                                  aes(x, y, label = site_lab(site), fill = v, colour = txt), size = MSE$MS_LAB,
                                  fontface = "bold", label.r = unit(0.25, "lines"), label.padding = unit(0.18, "lines"),
                                  label.size = 0.3, box.padding = 0.15, point.padding = 0.1, min.segment.length = 0,
                                  segment.colour = "grey30", segment.size = 0.3, seed = 42, max.overlaps = Inf) +
        scale_colour_identity() + fill_scale +
        coord_equal(xlim = ex[1:2], ylim = ex[3:4], expand = FALSE) +
        labs(title = title, x = NULL, y = NULL) +
        theme_void(base_size = 7) +
        theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 3)),
              panel.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
              legend.key.height = unit(0.6, "cm"), legend.key.width = unit(0.2, "cm"),
              legend.text = element_text(size = 6), legend.title = element_text(size = 6))
}
grad <- scales::gradient_n_pal(PAL)
cs <- null[, .(site, x, y, v = shift)]
p_shift <- map_panel(cs, scale_fill_gradientn(colours = PAL, name = NULL),
                     "Projected climate shift", grad(scales::rescale(cs$v)))
p_maps <- lapply(names(METH), function(m) {
    d <- rs[engine == m, .(site, x, y, v = consensus_rank)]
    map_panel(d, scale_fill_gradientn(colours = rev(PAL), limits = c(1, 11), breaks = c(1, 6, 11),
                                      name = "Consensus\nrank", guide = guide_colourbar(reverse = TRUE)),
              METH[[m]], grad(scales::rescale(12 - d$v, from = c(1, 11))))
})

# ---- each site's rank over the 199 random sets against its climate-shift rank ----
rankplot <- function(m, first) {
    x <- rs[engine == m]; s <- sm[engine == m]
    x[, side := vapply(seq_len(.N), function(i) {
        dx <- shift_rank - shift_rank[i]; dy <- abs(med_rank - med_rank[i])
        if (any(dx > 0 & dx <= 1.5 & dy < 1)) "left" else "right" }, "")]
    stopifnot(x[side == "left", all(shift_rank >= 2)])
    if (any(x$side == "left")) message("s_random_climate ", m, ": label left of site(s) ", paste(site_lab(x[side == "left", site]), collapse = ", "))
    ggplot(x, aes(shift_rank, med_rank)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
        geom_linerange(aes(ymin = lo, ymax = hi), colour = POINT_COL, alpha = 0.45, linewidth = 0.5) +
        geom_point(size = 1.2, colour = POINT_COL) +
        # labels right of the point as in the reference, except where another point sits within 1.5 ranks to the
        # right and less than 1 rank above or below: then left (display fix 2026-10-25, RDA "5" ran into site 4)
        geom_text(data = x[side == "right"], aes(label = site_lab(site)), size = MSE$MS_LAB, colour = "grey30",
                  hjust = 0, nudge_x = 0.28) +
        geom_text(data = x[side == "left"], aes(label = site_lab(site)), size = MSE$MS_LAB, colour = "grey30",
                  hjust = 1, nudge_x = -0.28) +
        scale_x_continuous(breaks = c(1, 6, 11), limits = c(0.5, 11.9)) +
        scale_y_continuous(breaks = c(1, 6, 11), limits = c(0.5, 11.5)) +
        coord_equal() +
        labs(x = "Projected climate shift rank", y = if (first) "Median offset rank\nover the random sets" else NULL,   # two lines: one overran the panel
             title = minus(sprintf("median τ = %.2f", s$med))) +
        MSE$theme_ms() +
        theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5))
}
fig <- wrap_plots(c(list(p_shift), p_maps, list(plot_spacer()),
                    lapply(seq_along(METH), function(i) rankplot(names(METH)[i], i == 1))),
                  ncol = 4, byrow = TRUE, heights = c(1.15, 1)) & theme(plot.margin = margin(1, 1, 1, 1, "mm"))
MSE$ms_save(file.path(OUT, "s_random_climate"), fig, 180, 82)
message("done")
