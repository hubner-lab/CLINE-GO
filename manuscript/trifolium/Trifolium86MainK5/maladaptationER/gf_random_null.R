#!/usr/bin/env Rscript
# Empirical random-SNP null for the 27-SNP GF site ranking (user 2026-10-15; plan
# ~/.claude/plans/peaceful-kindling-micali.md, Phase 1 follow-up).
# Question: GF adaptive (27 SNPs) ranks the 11 sites at Spearman 0.50 with the climate-only null, while the
# pipeline's one size-matched random GF draw gives 0.91. Is the adaptive set's departure from the climate order
# outside what random 27-SNP sets give?
#   - N_DRAWS random 27-SNP sets from the same matrix the pipeline samples (lfmm_imp_full_climate, all SNPs),
#     each fitted exactly as gradient_forest_model.R does (random branch: set.seed(seed); sample(); gradientForest
#     with ntree 1000, corr.threshold 0.5, maxLevel log2(0.368 n / 2), predictors bio_1,bio_3,bio_5, no spatial).
#   - N_ADAPT refits of the adaptive 27 SNPs with different seeds: the forest's own seed noise around 0.50.
#   - Site offset exactly as gradient_forest_offset.R: predict(extrap = TRUE) on the present and future cell of
#     each sample's raster cell, Euclidean distance, then the site mean (as site_maps/site_offsets.tsv).
#   - Statistic: Spearman rho between the site offsets and the climate-only null (site_maps/site_offsets.tsv).
#     Empirical one-sided p = (1 + #random draws with rho <= adaptive rho) / (1 + N), i.e. "random sets agree
#     with climate change at least as weakly as the adaptive set does".
# Checks before the draws: the pipeline's own random draw (seed 1396) and the saved adaptive model are
# re-projected and must reproduce their genetic_offset_site.tsv values.
# Outputs (site_maps/gf_random_null/): draws.tsv, site_offsets_draws.tsv, summary.tsv, check.tsv, fig_gf_random_null.png
# Run inside the image from /pipeline:  Rscript ANALYSIS/Trifolium86MainK5/maladaptationER/gf_random_null.R
suppressPackageStartupMessages({
    library(data.table); library(gradientForest); library(terra); library(qs); library(parallel); library(ggplot2)
    library(patchwork)
})
N_DRAWS <- 199L; N_ADAPT <- 20L; CORES <- as.integer(Sys.getenv("CORES", "40"))   # env CORES (v2 run: 160)
B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R   <- file.path("/pipeline", paste0(Sys.getenv("PROJ", "Trifolium86MainK5ER"), "_results"))
SC  <- "year2061-2080_ssp585"
SET <- "gea.ER_overlap"
# ARM=arm3 (default; bio_1,bio_3,bio_5) or ARM=arm16 (all 16 Climate.predictors, user 2026-10-15). The climate-only
# null is computed here on the arm's own predictors (scaled Euclidean present -> future distance on the site table,
# as build_site_maps.R), plus the bio_1,3,5 null as a second reference (rho_null3).
ARM <- Sys.getenv("ARM", "arm3")
P   <- switch(ARM, arm3 = c("bio_1", "bio_3", "bio_5"),
              arm16 = strsplit("bio_1,bio_2,bio_3,bio_4,bio_5,bio_6,bio_7,bio_8,bio_9,bio_10,bio_11,bio_12,bio_13,bio_15,bio_16,bio_19", ",")[[1]])
OUT <- file.path(B, "site_maps", if (ARM == "arm3") "gf_random_null" else paste0("gf_random_null_", ARM))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
WF  <- file.path(R, "_work/maf0.1_miss0.3_smiss0.8_rel0.999")

lfmm   <- fread(file.path(R, "_intermediate/climate_subset/lfmm_imp_full_climate.lfmm"))
vcfsnp <- fread(file.path(WF, "Trifolium86_chr1-6_K5imp_maf0.1.vcfsnp"), header = FALSE)[, paste0(V1, ":", V2)]
removed <- fread(file.path(WF, "Trifolium86_chr1-6.removed"), header = FALSE)
removed <- if (nrow(removed)) removed[, paste0(V1, ":", V2)] else character()
stopifnot(ncol(lfmm) == length(vcfsnp))
sig <- fread(file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv"))$SNPID
adapt_cols <- names(lfmm)[(!vcfsnp %in% removed) & (vcfsnp %in% sig)]
message("matrix ", nrow(lfmm), " x ", ncol(lfmm), "; adaptive SNPs ", length(adapt_cols))

samples <- fread(file.path(R, "_intermediate/samples/metadata_climate_valid.tsv"),
                 colClasses = c(site = "character", sample = "character"))
env <- fread(file.path(R, "climate/tables/present/climate_present_site.tsv"))[, ..P]
stopifnot(nrow(env) == nrow(lfmm), nrow(samples) == nrow(lfmm))
maxLevel <- log2(0.368 * nrow(lfmm) / 2)

# the present / future cell of every sample (what terra::extract on the offset raster reads)
template <- rast(file.path(R, "climate/rasters/present/climate_present_rasterstack.tif"))[[1]]
cells <- cellFromXY(template, as.matrix(samples[, .(longitude, latitude)]))
pres_all <- fread(file.path(R, "climate/tables/present/climate_present_all.tsv"))
fut_all  <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_all.tsv")))
pres_c <- as.data.frame(pres_all[match(cells, ID), ..P]); fut_c <- as.data.frame(fut_all[match(cells, ID), ..P])
stopifnot(!anyNA(pres_c), !anyNA(fut_c))

site_offset <- function(gf) {
    if (is.null(gf) || inherits(gf, "error") || is.null(gf$result) || !length(gf$result)) return(NULL)
    d <- sqrt(rowSums((predict(gf, fut_c, extrap = TRUE) - predict(gf, pres_c, extrap = TRUE))^2, na.rm = TRUE))
    data.table(site = samples$site, offset = d)[, .(offset = mean(offset)), by = site]
}
fit <- function(cols, seed) {
    set.seed(seed)   # pipeline: set.seed(42 + label) right before sample(); GF consumes the stream after it
    if (is.null(cols)) cols <- sample(names(lfmm), length(adapt_cols))
    gf <- tryCatch(suppressWarnings(gradientForest(cbind(env, lfmm[, ..cols]), predictor.vars = P,
                                                   response.vars = cols, ntree = 1000, maxLevel = maxLevel,
                                                   trace = FALSE, corr.threshold = 0.5)),
                   error = function(e) NULL)
    list(gf = gf, n_used = if (is.null(gf)) 0L else length(gf$result), so = site_offset(gf))
}

# ── checks: reproduce the pipeline's random draw and the saved adaptive model ──
ref <- function(f) fread(f, colClasses = c(site = "character", sample = "character"))[
    , .(ref = mean(genetic_offset)), by = site]
chk_rand  <- merge(fit(NULL, 42L + sum(utf8ToInt(SET)))$so,
                   ref(file.path(B, ARM, "gf_random", SET, "genetic_offset_site.tsv")), by = "site")
chk_adapt <- merge(site_offset(qread(file.path(B, ARM, "intermediate/gradient_forest", paste0(SET, "_nospatial"),
                                               "adaptive_model.qs"))),
                   ref(file.path(B, ARM, "Maladaptation/tables/gradient_forest", paste0(SET, "_nospatial"), SC,
                                 "genetic_offset_site.tsv")), by = "site")
check <- data.table(model = c("pipeline random draw (seed 1396, refit)", "saved adaptive model (re-projected)"),
                     n_sites = c(nrow(chk_rand), nrow(chk_adapt)),
                     max_abs_diff = signif(c(max(abs(chk_rand$offset - chk_rand$ref)),
                                             max(abs(chk_adapt$offset - chk_adapt$ref))), 3))
fwrite(check, file.path(OUT, "check.tsv"), sep = "\t"); print(check)
stopifnot(check$max_abs_diff < 1e-6)

# ── draws ───────────────────────────────────────────────────────────────────
pres_site <- fread(file.path(R, "climate/tables/present/climate_present_site.tsv"), colClasses = c(sample = "character"))
fut_site  <- fread(file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_site.tsv")),
                   colClasses = c(sample = "character"))
clim_null <- function(p) {
    fs <- fut_site[match(pres_site$sample, fut_site$sample)]
    sdv <- sapply(pres_site[, ..p], sd)
    d <- sqrt(rowSums(((as.matrix(fs[, ..p]) - as.matrix(pres_site[, ..p])) / rep(sdv, each = nrow(pres_site)))^2))
    data.table(sample = pres_site$sample, null = d)[samples, on = "sample"][, .(null = mean(null)), by = site]
}
null  <- clim_null(P)
null3 <- clim_null(c("bio_1", "bio_3", "bio_5"))[, .(site, null3 = null)]
fwrite(merge(null, null3, by = "site"), file.path(OUT, "climate_null_sites.tsv"), sep = "\t")
adapt_ref <- chk_adapt[, .(site, adaptive = ref)]
jobs <- c(lapply(seq_len(N_DRAWS), function(i) list(kind = "random", seed = 100000L + i)),
          lapply(seq_len(N_ADAPT), function(i) list(kind = "adaptive refit", seed = 200000L + i)))
res <- mclapply(jobs, function(j) {
    r <- fit(if (j$kind == "random") NULL else adapt_cols, j$seed)
    list(kind = j$kind, seed = j$seed, n_used = r$n_used, so = r$so)
}, mc.cores = CORES, mc.preschedule = FALSE)

draws <- rbindlist(lapply(res, function(r) {
    if (is.null(r$so)) return(data.table(kind = r$kind, seed = r$seed, n_used = r$n_used,
                                         rho_null = NA_real_, rho_null3 = NA_real_, rho_adaptive = NA_real_))
    x <- merge(merge(merge(r$so, null, by = "site"), null3, by = "site"), adapt_ref, by = "site")
    data.table(kind = r$kind, seed = r$seed, n_used = r$n_used,
               rho_null = cor(x$offset, x$null, method = "spearman"),
               rho_null3 = cor(x$offset, x$null3, method = "spearman"),
               rho_adaptive = cor(x$offset, x$adaptive, method = "spearman"))
}))
fwrite(draws, file.path(OUT, "draws.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, function(r) if (!is.null(r$so)) r$so[, `:=`(kind = r$kind, seed = r$seed)])),
       file.path(OUT, "site_offsets_draws.tsv"), sep = "\t")

a0  <- merge(adapt_ref, null, by = "site")[, cor(adaptive, null, method = "spearman")]
a03 <- merge(adapt_ref, null3, by = "site")[, cor(adaptive, null3, method = "spearman")]
rnd <- draws[kind == "random" & !is.na(rho_null)]
adr <- draws[kind == "adaptive refit" & !is.na(rho_null)]
summary <- data.table(
    quantity = c("adaptive rho vs null (pipeline model)", "adaptive refits: median rho vs null [min, max]",
                 "random draws fitted / requested", "random draws empty (no SNP with R2 > 0)",
                 "random: median rho vs null [2.5%, 97.5%]", "random draws with rho <= adaptive rho",
                 "empirical one-sided p", "empirical p using the adaptive refits' median rho",
                 "random: median rho with the adaptive order",
                 "predictors", "adaptive rho vs bio_1,3,5 null", "random: median rho vs bio_1,3,5 null",
                 "empirical one-sided p vs bio_1,3,5 null"),
    value = c(sprintf("%.3f", a0),
              sprintf("%.3f [%.3f, %.3f]", median(adr$rho_null), min(adr$rho_null), max(adr$rho_null)),
              sprintf("%d / %d", nrow(rnd), N_DRAWS), sum(draws[kind == "random", is.na(rho_null)]),
              sprintf("%.3f [%.3f, %.3f]", median(rnd$rho_null), quantile(rnd$rho_null, .025), quantile(rnd$rho_null, .975)),
              sum(rnd$rho_null <= a0),
              sprintf("%.3f", (1 + sum(rnd$rho_null <= a0)) / (1 + nrow(rnd))),
              sprintf("%.3f", (1 + sum(rnd$rho_null <= median(adr$rho_null))) / (1 + nrow(rnd))),
              sprintf("%.3f", median(rnd$rho_adaptive)),
              paste(P, collapse = ","), sprintf("%.3f", a03), sprintf("%.3f", median(rnd$rho_null3)),
              sprintf("%.3f", (1 + sum(rnd$rho_null3 <= a03)) / (1 + nrow(rnd)))))
fwrite(summary, file.path(OUT, "summary.tsv"), sep = "\t"); print(summary)

# ── figure ──────────────────────────────────────────────────────────────────
pa <- ggplot(rnd, aes(rho_null)) +
    geom_histogram(binwidth = 0.05, boundary = 1, fill = "#0090ab", colour = "white", linewidth = 0.2) +
    geom_rug(data = adr, aes(rho_null), colour = "#ae0000", linewidth = 0.6, length = unit(0.06, "npc")) +
    geom_vline(xintercept = a0, colour = "#ae0000", linewidth = 0.8) +
    scale_x_continuous(limits = c(-1, 1.05), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Spearman rho: GF site offsets vs climate-only null (11 sites)", y = "Random 27-SNP draws",
         title = sprintf("a  %d random 27-SNP GF fits; red line = adaptive 27 SNPs (rho %.2f), red ticks = %d adaptive refits",
                         nrow(rnd), a0, nrow(adr))) +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
so <- fread(file.path(OUT, "site_offsets_draws.tsv"), colClasses = c(site = "character"))
so[, rank := frank(-offset), by = .(kind, seed)]
nr <- null[, .(site, null_rank = frank(-null))]
ar <- adapt_ref[, .(site, adaptive_rank = frank(-adaptive))]
ord <- nr[order(null_rank), site]
sr <- merge(so[kind == "random"], nr, by = "site")[, site := factor(site, ord)]
pb <- ggplot(sr, aes(site, rank)) +
    geom_boxplot(fill = "#e6cea0", colour = "grey35", outlier.size = 0.4, linewidth = 0.3, width = 0.6) +
    geom_point(data = nr[, site := factor(site, ord)], aes(site, null_rank), shape = 4, size = 2.6, colour = "black") +
    geom_point(data = ar[, site := factor(site, ord)], aes(site, adaptive_rank), size = 2.4, colour = "#ae0000") +
    scale_y_reverse(breaks = c(1, 6, 11)) +
    labs(x = "Site (ordered by climate-null rank)", y = "Site risk rank (1 = highest)",
         title = "b  Site ranks: boxes = random draws, x = climate-only null, red = adaptive 27 SNPs") +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(OUT, "fig_gf_random_null.png"), pa / pb, width = 10, height = 8, dpi = 200, bg = "white")
message("done")
