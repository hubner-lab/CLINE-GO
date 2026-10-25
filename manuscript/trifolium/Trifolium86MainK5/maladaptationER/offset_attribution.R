#!/usr/bin/env Rscript
# What drives the 27-SNP site risk ranking? (user 2026-10-15, manuscript: link the offset to the GEA biology)
# Predictors bio_1,bio_3,bio_5; SSP585 2061-2080 mean of 3 GCMs (the pipeline's averaged future); methods GF,
# geometric (k 5), RDA (3 axes). Consensus = family-weighted rank: (GF rank + mean(geometric, RDA rank)) / 2,
# because geometric and RDA offsets are theoretically equivalent (Gain et al. 2023).
#
# A. Predictor attribution, one predictor at a time: future climate = present climate with ONLY predictor p moved
#    to its future value. GF: the saved adaptive model re-projected (GF transforms each predictor separately, so
#    offset^2 = sum_p offset_p^2 exactly). Geometric / RDA: their pipeline scripts on the modified future tables
#    (cross terms exist there, so shares are approximate).
# B. Leave-one-region-out: the 27 SNPs minus the SNPs of one of the 13 GEA regions (regions_combined.tsv), all three
#    methods refitted (GF with the pipeline's set.seed(42); the full-set refit must reproduce the saved model).
#    Per region: Spearman of the site offsets with the full set, and the change of the consensus rank per site.
# Outputs (site_maps/attribution/): predictor_offsets.tsv, predictor_shares.tsv, loro_sites.tsv, loro_summary.tsv,
#   check.tsv; work/ holds the per-run site tables.
# Run inside the image from /pipeline:  Rscript ANALYSIS/Trifolium86MainK5/maladaptationER/offset_attribution.R
suppressPackageStartupMessages({ library(data.table); library(gradientForest); library(terra); library(qs); library(parallel) })
CORES <- as.integer(Sys.getenv("CORES", "40"))   # env CORES (v2 run: 160)
B   <- file.path("/pipeline", Sys.getenv("AN", "ANALYSIS/Trifolium86MainK5/maladaptationER"))   # env AN/PROJ: v2 = PC5 rebuild (2026-10-16)
R   <- file.path("/pipeline", paste0(Sys.getenv("PROJ", "Trifolium86MainK5ER"), "_results"))
W   <- file.path(R, "_work/maf0.1_miss0.3_smiss0.8_rel0.999")
SC  <- "year2061-2080_ssp585"
SET <- "gea.ER_overlap"
P   <- c("bio_1", "bio_3", "bio_5"); PRED <- paste(P, collapse = ",")
OUT <- file.path(B, "site_maps/attribution"); dir.create(file.path(OUT, "work"), recursive = TRUE, showWarnings = FALSE)
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
          cvalid = file.path(R, "_intermediate/samples/climate_valid_samples.list"),
          sig = file.path(R, "_intermediate/snp_sets", SET, "selected_snps.tsv"),
          regions = file.path(R, "GEA/tables/regions_combined.tsv"))
stopifnot(all(file.exists(unlist(F))))

samples <- fread(F$samp, colClasses = c(site = "character", sample = "character"))
site_mean <- function(f) fread(f, colClasses = c(site = "character", sample = "character"))[
    , .(offset = mean(genetic_offset)), by = site]

run_engine <- function(engine, snp_tsv, fut_site, fut_all, d) {
    dir.create(file.path(d, "plots"), recursive = TRUE, showWarnings = FALSE)
    site <- file.path(d, "site.tsv")
    if (file.exists(site)) return(site)
    args <- if (engine == "geometric") c(F$lfmm, F$vcfsnp, F$removed, snp_tsv, F$pres_site, fut_site,
        F$pres_all, fut_all, F$rast, F$samp, PRED, 5, "TRUE", site, file.path(d, "map.tsv"),
        file.path(d, "raster.tif"), file.path(d, "importance.png"), SC, file.path(d, "diagnostics.tsv"))
    else c(F$lfmm, F$vcfsnp, F$removed, snp_tsv, F$pres_site, fut_site, F$pres_all, fut_all, F$rast, F$samp,
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

# ── GF helpers (as gf_random_null.R) ────────────────────────────────────────
lfmm   <- fread(F$lfmm)
vcfsnp <- fread(F$vcfsnp, header = FALSE)[, .(SNPID = paste0(V1, ":", V2), chr = as.character(V1), pos = V2)]
env    <- fread(F$pres_site)[, ..P]
maxLevel <- log2(0.368 * nrow(lfmm) / 2)
template <- rast(F$rast)[[1]]
cells    <- cellFromXY(template, as.matrix(samples[, .(longitude, latitude)]))
pres_all <- fread(F$pres_all); fut_all <- fread(F$fut_all)
pres_c <- as.data.frame(pres_all[match(cells, ID), ..P]); fut_c <- as.data.frame(fut_all[match(cells, ID), ..P])
gf_sites <- function(gf, newfut = fut_c) {
    d <- sqrt(rowSums((predict(gf, newfut, extrap = TRUE) - predict(gf, pres_c, extrap = TRUE))^2, na.rm = TRUE))
    data.table(site = samples$site, offset = d)[, .(offset = mean(offset)), by = site]
}
gf_fit <- function(snpids) {
    cols <- names(lfmm)[vcfsnp$SNPID %in% snpids]
    set.seed(42)   # gradient_forest_model.R: the only seed before the adaptive fit
    gradientForest(cbind(env, lfmm[, ..cols]), predictor.vars = P, response.vars = cols, ntree = 1000,
                   maxLevel = maxLevel, trace = FALSE, corr.threshold = 0.5)
}
saved_gf <- qread(file.path(B, "arm3/intermediate/gradient_forest", paste0(SET, "_nospatial"), "adaptive_model.qs"))

# ── A. predictor attribution ────────────────────────────────────────────────
pres_site <- fread(F$pres_site, colClasses = c(sample = "character"))
fut_site  <- fread(F$fut_site, colClasses = c(sample = "character"))
modA <- lapply(P, function(p) {
    d <- file.path(OUT, "work/predictor", p); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    s <- copy(pres_site); s[[p]] <- fut_site[match(s$sample, sample)][[p]]
    a <- copy(pres_all[ID %in% fut_all$ID]); a[[p]] <- fut_all[match(a$ID, ID)][[p]]
    fwrite(s, file.path(d, "fut_site.tsv"), sep = "\t"); fwrite(a, file.path(d, "fut_all.tsv"), sep = "\t")
    nf <- pres_c; nf[[p]] <- fut_c[[p]]
    list(p = p, d = d, gf = gf_sites(saved_gf, nf))
})
resA <- mclapply(CJ(eng = c("geometric", "rda"), i = seq_along(P))[, .(eng, i)] |> split(by = c("eng", "i")), function(j) {
    m <- modA[[j$i]]
    run_engine(j$eng, F$sig, file.path(m$d, "fut_site.tsv"), file.path(m$d, "fut_all.tsv"), file.path(m$d, j$eng))
}, mc.cores = 6)
predA <- rbindlist(c(
    lapply(modA, function(m) m$gf[, `:=`(method = "gf", predictor = m$p)]),
    lapply(modA, function(m) rbindlist(lapply(c("geometric", "rda"), function(e)
        site_mean(file.path(m$d, e, "site.tsv"))[, `:=`(method = e, predictor = m$p)])))))
full <- rbind(gf_sites(saved_gf)[, method := "gf"],
              site_mean(file.path(B, "arm3/Maladaptation/tables/geometric_offset", paste0(SET, "_nospatial"), SC,
                                  "genetic_offset_site.tsv"))[, method := "geometric"],
              site_mean(file.path(B, "arm3_ax3/Maladaptation/tables/rda_offset", paste0(SET, "_nospatial"), SC,
                                  "genetic_offset_site.tsv"))[, method := "rda"])
fwrite(rbind(predA, full[, predictor := "all (full future)"]), file.path(OUT, "predictor_offsets.tsv"), sep = "\t")
shares <- predA[, .(share = offset^2 / sum(offset^2), predictor), by = .(method, site)]
rho_full <- merge(predA, full[, .(site, method, full = offset)], by = c("site", "method"))[
    , .(rho_with_full_order = cor(offset, full, method = "spearman")), by = .(method, predictor)]
psum <- merge(shares[, .(mean_share = mean(share)), by = .(method, predictor)], rho_full, by = c("method", "predictor"))
fwrite(psum[order(method, -mean_share)], file.path(OUT, "predictor_shares.tsv"), sep = "\t")
fwrite(dcast(shares, method + site ~ predictor, value.var = "share"), file.path(OUT, "predictor_shares_by_site.tsv"), sep = "\t")
print(psum[order(method, -mean_share)])

# ── B. leave-one-region-out ─────────────────────────────────────────────────
sig  <- fread(F$sig)
regs <- fread(F$regions, colClasses = c(chr = "character"))
sig[, region := vapply(seq_len(.N), function(i)
    regs[chr == as.character(sig$chr[i]) & start <= sig$pos[i] & end >= sig$pos[i], region_id][1], "")]
stopifnot(!anyNA(sig$region))
fwrite(sig[, .(SNPID, region)], file.path(OUT, "snp_regions.tsv"), sep = "\t")

# check: GF full refit reproduces the saved model's site offsets
gf_full <- gf_sites(gf_fit(sig$SNPID))
chk <- merge(gf_full, gf_sites(saved_gf), by = "site")
check <- data.table(model = "GF full 27-SNP refit vs saved adaptive model",
                    max_abs_diff = signif(max(abs(chk$offset.x - chk$offset.y)), 3))
fwrite(check, file.path(OUT, "check.tsv"), sep = "\t"); print(check)

jobs <- CJ(region = unique(sig$region), eng = c("gf", "geometric", "rda"))
resB <- mclapply(seq_len(nrow(jobs)), function(j) {
    r <- jobs$region[j]; e <- jobs$eng[j]
    keep <- sig[region != r]
    d <- file.path(OUT, "work/loro", r, e); dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (e == "gf") {
        f <- file.path(d, "site_gf.tsv")
        if (!file.exists(f)) fwrite(gf_sites(gf_fit(keep$SNPID)), f, sep = "\t")
        return(fread(f, colClasses = c(site = "character"))[, `:=`(region = r, method = e)])
    }
    snp <- file.path(d, "snps.tsv"); fwrite(keep[, .(SNPID, chr, pos, min_pvalue)], snp, sep = "\t")
    s <- run_engine(e, snp, F$fut_site, F$fut_all, d)
    if (is.na(s)) return(NULL)
    site_mean(s)[, `:=`(region = r, method = e)]
}, mc.cores = CORES, mc.preschedule = FALSE)
loro <- rbindlist(resB, use.names = TRUE)
message("leave-one-region-out runs: ", nrow(unique(loro[, .(region, method)])), " / ", nrow(jobs))

cons <- function(x) {   # x: site, method, offset -> family-weighted consensus rank (1 = highest risk)
    x <- copy(x)[, rk := frank(-offset), by = method]
    w <- dcast(x, site ~ method, value.var = "rk")
    w[, .(site, consensus = (gf + (geometric + rda) / 2) / 2)]
}
full_c <- cons(full[, .(site, method, offset)])[, full_rank := frank(consensus)]
loro_sites <- rbindlist(lapply(unique(loro$region), function(r) {
    x <- loro[region == r]
    cmp <- merge(x, full[, .(site, method, full = offset)], by = c("site", "method"))
    cc  <- cons(x[, .(site, method, offset)])[, drop_rank := frank(consensus)]
    merge(cc, full_c[, .(site, full_rank)], by = "site")[, `:=`(region = r, delta_rank = drop_rank - full_rank)]
}))
fwrite(loro_sites, file.path(OUT, "loro_sites.tsv"), sep = "\t")
loro_sum <- merge(
    sig[, .(n_snps_removed = .N), by = region],
    dcast(merge(loro, full[, .(site, method, full = offset)], by = c("site", "method"))[
        , .(rho = round(cor(offset, full, method = "spearman"), 3)), by = .(region, method)],
        region ~ method, value.var = "rho"), by = "region")
loro_sum <- merge(loro_sum, loro_sites[, .(consensus_rho = round(cor(drop_rank, full_rank, method = "spearman"), 3),
                                           max_abs_rank_shift = max(abs(delta_rank))), by = region], by = "region")
setnames(loro_sum, c("gf", "geometric", "rda"), c("rho_gf", "rho_geometric", "rho_rda"))
fwrite(loro_sum[order(consensus_rho)], file.path(OUT, "loro_summary.tsv"), sep = "\t")
print(loro_sum[order(consensus_rho)])
print(full_c[order(full_rank)])
message("done")
