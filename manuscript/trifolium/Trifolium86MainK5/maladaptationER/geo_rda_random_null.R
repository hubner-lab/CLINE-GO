#!/usr/bin/env Rscript
# Random-SNP null for the geometric and RDA offsets of the 27-SNP set (user 2026-10-15), the counterpart of
# gf_random_null.R. Each of the 199 draws is the SAME random SNP set as GF draw i (set.seed(100000 + i);
# sample() over the same 177,259 matrix columns), so the three engines are compared on identical draws.
# Each draw runs the pipeline's own scripts with the arguments of their rules (as run_gcm_offsets.sh, but on the
# pipeline's 3-GCM mean future):
#   geometric_offset.R  k = 5, scale TRUE, predictors bio_1,bio_3,bio_5
#   rda_offset.R        axes 3 (fixed), axis_alpha 0.05, condition_pcs 0, seed 42, permutations 9: with fixed
#                       axes the by-axis anova only fills diagnostics, it never touches the offset (rda_offset.R:269-285)
# Check first: the real 27-SNP set through the same calls must reproduce arm3/ (geometric) and arm3_ax3/ (RDA).
# The two engines are deterministic given the SNPs, so there are no "adaptive refits" as for GF.
# Outputs (site_maps/random_null/): draws.tsv (all 3 engines, GF from gf_random_null/), summary.tsv, check.tsv,
#   fig_random_null_3engines.png; per-draw site tables under work/<engine>/<i>/ (rasters deleted)
# Run inside the image from /pipeline:  Rscript ANALYSIS/Trifolium86MainK5/maladaptationER/geo_rda_random_null.R
suppressPackageStartupMessages({ library(data.table); library(parallel); library(ggplot2); library(patchwork) })
N_DRAWS <- 199L; CORES <- as.integer(Sys.getenv("CORES", "40"))   # env CORES (v2 run: 160)
B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R   <- file.path("/pipeline", paste0(Sys.getenv("PROJ", "Trifolium86MainK5ER"), "_results"))
W   <- file.path(R, "_work/maf0.1_miss0.3_smiss0.8_rel0.999")
SC  <- "year2061-2080_ssp585"
SET <- "gea.ER_overlap"
OUT <- file.path(B, "site_maps/random_null"); dir.create(file.path(OUT, "work"), recursive = TRUE, showWarnings = FALSE)
F <- list(lfmm = file.path(R, "_intermediate/climate_subset/lfmm_imp_full_climate.lfmm"),
          vcfsnp = file.path(W, "Trifolium86_chr1-6_K5imp_maf0.1.vcfsnp"),
          removed = file.path(W, "Trifolium86_chr1-6.removed"),
          pres_site = file.path(R, "climate/tables/present/climate_present_site.tsv"),
          fut_site = file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_site.tsv")),
          pres_all = file.path(R, "climate/tables/present/climate_present_all.tsv"),
          fut_all = file.path(R, "climate/tables/future", paste0("climate_future_", SC, "_all.tsv")),
          rast = file.path(R, "climate/rasters/present/climate_present_rasterstack.tif"),
          samp = file.path(R, "_intermediate/samples/metadata_climate_valid.tsv"),
          pca = file.path(W, "ld0.2_win100_step10/Trifolium86_chr1-6.pca/Trifolium86_chr1-6.projections"),
          order = file.path(R, "_intermediate/samples/samples_order.list"),
          cvalid = file.path(R, "_intermediate/samples/climate_valid_samples.list"))
stopifnot(all(file.exists(unlist(F))))
PRED <- "bio_1,bio_3,bio_5"

run_engine <- function(engine, snp_tsv, d) {
    dir.create(file.path(d, "plots"), recursive = TRUE, showWarnings = FALSE)
    site <- file.path(d, "site.tsv")
    if (file.exists(site)) return(site)
    args <- if (engine == "geometric") c(F$lfmm, F$vcfsnp, F$removed, snp_tsv, F$pres_site, F$fut_site,
        F$pres_all, F$fut_all, F$rast, F$samp, PRED, 5, "TRUE", site, file.path(d, "map.tsv"),
        file.path(d, "raster.tif"), file.path(d, "importance.png"), SC, file.path(d, "diagnostics.tsv"))
    else c(F$lfmm, F$vcfsnp, F$removed, snp_tsv, F$pres_site, F$fut_site, F$pres_all, F$fut_all, F$rast, F$samp,
        PRED, 3, 0.05, 9, 0, F$pca, F$order, F$cvalid, 42, 1, site, file.path(d, "map.tsv"),
        file.path(d, "raster.tif"), file.path(d, "importance.png"), file.path(d, "diagnostics.tsv"),
        paste0(file.path(d, "plots"), "/"), SC)
    script <- if (engine == "geometric") "/pipeline/scripts/geometric_offset.R" else "/pipeline/scripts/rda_offset.R"
    st <- system2("Rscript", c(script, shQuote(args)), stdout = file.path(d, "run.log"), stderr = file.path(d, "run.log"))
    unlink(c(file.path(d, "raster.tif"), file.path(d, "map.tsv"), file.path(d, "importance.png"),
             file.path(d, "importance.svg"), file.path(d, "plots")), recursive = TRUE)
    if (st != 0 || !file.exists(site)) return(NA_character_)
    site
}
site_mean <- function(f) fread(f, colClasses = c(site = "character", sample = "character"))[
    , .(offset = mean(genetic_offset)), by = site]

# ── random SNP sets: same draws as gf_random_null.R ─────────────────────────
vcfsnp <- fread(F$vcfsnp, header = FALSE)[, .(SNPID = paste0(V1, ":", V2), chr = V1, pos = V2)]
# size-matched to the adaptive set (v1: 27; was hardcoded, which left the PC5 rebuild's draws at 27 instead of 30)
n_adapt <- nrow(fread(file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv")))
message("random draws of ", n_adapt, " SNPs")
snp_files <- vapply(seq_len(N_DRAWS), function(i) {
    set.seed(100000L + i)
    idx <- sample(seq_len(nrow(vcfsnp)), n_adapt)   # identical indices to sample(names(lfmm), 27)
    f <- file.path(OUT, "work/snps", sprintf("draw_%03d.tsv", i)); dir.create(dirname(f), showWarnings = FALSE)
    fwrite(vcfsnp[idx][, min_pvalue := NA_real_], f, sep = "\t"); f
}, "")

# ── check: the real 27-SNP set reproduces the snapshots ────────────────────
real <- file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv")
ref_geo <- file.path(B, "arm3/Maladaptation/tables/geometric_offset", paste0(SET, "_nospatial"), SC, "genetic_offset_site.tsv")
ref_rda <- file.path(B, "arm3_ax3/Maladaptation/tables/rda_offset", paste0(SET, "_nospatial"), SC, "genetic_offset_site.tsv")
chk <- mclapply(c("geometric", "rda"), function(e) run_engine(e, real, file.path(OUT, "work", e, "real")), mc.cores = 2)
cmp <- function(a, b) { m <- merge(site_mean(a), site_mean(b), by = "site"); signif(max(abs(m$offset.x - m$offset.y)), 3) }
check <- data.table(engine = c("geometric", "rda"),
                    max_abs_diff_vs_snapshot = c(cmp(chk[[1]], ref_geo), cmp(chk[[2]], ref_rda)))
fwrite(check, file.path(OUT, "check.tsv"), sep = "\t"); print(check)
stopifnot(check$max_abs_diff_vs_snapshot < 1e-6)

# ── draws ───────────────────────────────────────────────────────────────────
jobs <- CJ(engine = c("geometric", "rda"), i = seq_len(N_DRAWS))
res <- mclapply(seq_len(nrow(jobs)), function(j) {
    e <- jobs$engine[j]; i <- jobs$i[j]
    run_engine(e, snp_files[i], file.path(OUT, "work", e, sprintf("%03d", i)))
}, mc.cores = CORES, mc.preschedule = FALSE)
jobs[, site_file := unlist(res)]
message("failed runs: ", jobs[is.na(site_file), .N], " (", paste(jobs[is.na(site_file), paste(engine, i)], collapse = ", "), ")")

# ── statistics: all three engines on the same draws ─────────────────────────
null <- fread(file.path(B, "site_maps/site_offsets.tsv"), colClasses = c(site = "character"))[
    method == "null", .(site, null = offset)]
adapt <- rbind(site_mean(ref_geo)[, engine := "geometric"], site_mean(ref_rda)[, engine := "rda"],
    fread(file.path(B, "site_maps/site_offsets.tsv"), colClasses = c(site = "character"))[
        method == "gf" & set == SET, .(site, offset, engine = "gf")])
so <- rbind(rbindlist(lapply(which(!is.na(jobs$site_file)), function(j)
             site_mean(jobs$site_file[j])[, `:=`(engine = jobs$engine[j], draw = jobs$i[j])])),
            fread(file.path(B, "site_maps/gf_random_null/site_offsets_draws.tsv"), colClasses = c(site = "character"))[
                kind == "random", .(site, offset, engine = "gf", draw = seed - 100000L)])
fwrite(so, file.path(OUT, "site_offsets_draws.tsv"), sep = "\t")
rho <- function(x, y) cor(x, y, method = "spearman")
draws <- merge(merge(so, null, by = "site"), adapt[, .(site, engine, adaptive = offset)], by = c("site", "engine"))[
    , .(rho_null = rho(offset, null), rho_adaptive = rho(offset, adaptive)), by = .(engine, draw)]
fwrite(draws, file.path(OUT, "draws.tsv"), sep = "\t")
a0 <- merge(adapt, null, by = "site")[, .(rho_adaptive_vs_null = rho(offset, null)), by = engine]
summary <- merge(a0, draws[, .(n_draws = .N, random_median = median(rho_null),
                               random_q025 = quantile(rho_null, .025), random_q975 = quantile(rho_null, .975),
                               random_median_rho_with_adaptive = median(rho_adaptive)), by = engine], by = "engine")
summary[, n_le := draws[summary, on = "engine", sum(rho_null <= i.rho_adaptive_vs_null), by = .EACHI]$V1]
summary[, p_lower := (1 + n_le) / (1 + n_draws)]                 # adaptive less climate-like than random
summary[, p_two_sided := pmin(1, 2 * pmin(p_lower, (1 + n_draws - n_le) / (1 + n_draws)))]
num <- names(summary)[vapply(summary, is.numeric, TRUE)]
summary[, (num) := lapply(.SD, function(v) round(v, 3)), .SDcols = num]
fwrite(summary, file.path(OUT, "summary.tsv"), sep = "\t"); print(summary)

# ── figure ──────────────────────────────────────────────────────────────────
LAB <- c(gf = "GF", geometric = "Geometric", rda = "RDA offset (3 axes)")
draws[, eng := factor(LAB[engine], LAB)]; a0[, eng := factor(LAB[engine], LAB)]
summary[, eng := factor(LAB[engine], LAB)]
pa <- ggplot(draws, aes(rho_null)) +
    geom_histogram(binwidth = 0.1, boundary = 0, fill = "#0090ab", colour = "white", linewidth = 0.2) +
    geom_vline(data = a0, aes(xintercept = rho_adaptive_vs_null), colour = "#ae0000", linewidth = 0.8) +
    geom_text(data = summary, aes(x = -1, y = Inf, label = sprintf("adaptive rho %.2f\nempirical p (lower) %.3f\ntwo-sided %.3f",
                                                                    rho_adaptive_vs_null, p_lower, p_two_sided)),
              hjust = 0, vjust = 1.2, size = 2.6) +
    facet_wrap(~ eng, nrow = 1) + scale_x_continuous(limits = c(-1.15, 1.15), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Spearman rho: site offsets vs climate-only null (11 sites)", y = "Random 27-SNP draws",
         title = "a  199 random 27-SNP sets (same draws for all engines); red = adaptive 27 SNPs") +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"),
                                         strip.background = element_blank(), strip.text = element_text(face = "bold"))
nr  <- null[, .(site, null_rank = frank(-null))]; ord <- nr[order(null_rank), site]
sr  <- copy(so)[, rank := frank(-offset), by = .(engine, draw)][, `:=`(site = factor(site, ord), eng = factor(LAB[engine], LAB))]
ar  <- copy(adapt)[, rank := frank(-offset), by = engine][, `:=`(site = factor(site, ord), eng = factor(LAB[engine], LAB))]
pb <- ggplot(sr, aes(site, rank)) +
    geom_boxplot(fill = "#e6cea0", colour = "grey35", outlier.size = 0.3, linewidth = 0.3, width = 0.6) +
    geom_point(data = nr[, site := factor(site, ord)], aes(site, null_rank), shape = 4, size = 2.2, colour = "black",
               inherit.aes = FALSE) +
    geom_point(data = ar, aes(site, rank), size = 2, colour = "#ae0000") +
    facet_wrap(~ eng, nrow = 1) + scale_y_reverse(breaks = c(1, 6, 11)) +
    labs(x = "Site (ordered by climate-null rank)", y = "Site risk rank (1 = highest)",
         title = "b  boxes = random draws, x = climate-only null, red = adaptive 27 SNPs") +
    theme_classic(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"),
                                         axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5),
                                         strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUT, "fig_random_null_3engines.png"), pa / pb, width = 13, height = 8, dpi = 200, bg = "white")
message("done")
