#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_setup_check.R -- Phase 1 of the landscapes pack: the inputs every later job reads.
#   * 600 distinct seeds per landscape, no seed in two landscapes, 3 x 4 x 5 x 10 balance
#     (architecture x pleiotropy/selection x demography x replicates);
#   * every source table and folder exists and covers exactly the landscape's seeds;
#   * the on-disk SS-Clines tables equal their Data_S2.zip members (sha256).
# Writes <OUT>/setup_check.tsv (check, landscape, ok, detail). Stops if any check fails.
# Usage: OUT_DIR=<dir> Rscript benchmarks/mvp_land_setup_check.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "landscape_ms"))
LANDS <- data.table(
    land  = c("SS-Clines", "SS-Mtn", "Est-Clines"),
    lkey  = c("ssclines", "ssmtn", "estclines"),
    gal   = c("figures_ld_gallery14", "figures_land_ssmtn_gallery", "figures_land_estclines_gallery"),
    rank  = c("figures_ld_offset14/rank", "figures_land_ssmtn_offset15_ssmtn/rank",
              "figures_land_estclines_offset15_estclines/rank"),
    rem   = c("remeasure600", "remeasure600_ssmtn", "remeasure600_estclines"),
    off   = c("offset14", "offset15_ssmtn", "offset15_estclines"),
    det   = c("detection600", "detection600_ssmtn", "detection600_estclines"),
    seeds = c("work/ld_phase3_2026-10-26/scripts/seeds600.txt", "work/mvp_land/seeds_ssmtn.txt",
              "work/mvp_land/seeds_estclines.txt"))
R <- list()
rec <- function(check, land, ok, detail = "") {
    R[[length(R) + 1L]] <<- data.table(check = check, landscape = land, ok = isTRUE(ok), detail = as.character(detail))
    message(sprintf("%-34s %-11s %s %s", check, land, if (isTRUE(ok)) "ok  " else "FAIL", detail))
}
man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
ALL <- list()
for (i in seq_len(nrow(LANDS))) {
    L <- LANDS$land[i]
    s <- trimws(readLines(file.path(ROOT, LANDS$seeds[i]))); s <- s[nzchar(s)]
    ALL[[L]] <- s
    rec("600 distinct seeds", L, length(s) == 600L && !anyDuplicated(s), length(s))
    m <- man[seed %in% s]
    rec("seeds in the manifest", L, nrow(m) == 600L, nrow(m))
    m[, cond := sub("^[^_]*_", "", architecture)]
    bal <- m[, .N, by = .(arch_level, cond, demog_level_sub)]
    rec("3 x 4 x 5 x 10 balance", L, nrow(bal) == 60L && all(bal$N == 10L), sprintf("%d cells, N %d-%d", nrow(bal), min(bal$N), max(bal$N)))
    rec("n_traits 2", L, all(m$n_traits == 2L), paste(unique(m$n_traits), collapse = ","))
    # source tables cover exactly the landscape's seeds
    for (f in c(file.path(LANDS$gal[i], "detection_per_seed.tsv"), file.path(LANDS$gal[i], "offset_delta_per_seed.tsv"),
                file.path(LANDS$gal[i], "B4_aucpr_per_seed.tsv"), file.path(LANDS$rank[i], "rank_long.tsv"),
                file.path(LANDS$rem[i], "rdaunc_rda2x_ldb/rank_metrics.tsv"),
                file.path(LANDS$rem[i], "rdaunc_rda2x_ldb/calls_combine.tsv"),
                file.path(LANDS$rem[i], "rdaunc_rda2x_ldb/calls_per_method.tsv"),
                file.path(LANDS$rem[i], "onefit_rda2x_ldb/rank_metrics.tsv"),
                file.path(LANDS$det[i], "aucpr_per_seed.tsv"))) {
        p <- file.path(EVAL, f)
        if (!file.exists(p)) { rec("source exists", L, FALSE, f); next }
        D <- fread(p, select = "seed", colClasses = c(seed = "character"))
        rec("source covers the seeds", L, setequal(D$seed, s), f)
    }
    rec("offset folder", L, dir.exists(file.path(EVAL, LANDS$off[i], "gardens")),
        sprintf("%d seed folders", length(list.dirs(file.path(EVAL, LANDS$off[i], "gardens"), recursive = FALSE))))
    tree <- if (L == "SS-Clines") sprintf("_archive/mvp_ssclines/results/MVP%shd_results", s) else sprintf("MVP%shd_results", s)
    rec("hd trees", L, all(dir.exists(file.path(ROOT, tree))), sum(dir.exists(file.path(ROOT, tree))))
    if (L != "SS-Clines") {
        rc <- sprintf("MVP%shdrc_results", s)
        rec("hdrc trees", L, all(dir.exists(file.path(ROOT, rc))), sum(dir.exists(file.path(ROOT, rc))))
    }
    rec("data/mvp VCF", L, all(file.exists(file.path(ROOT, sprintf("data/mvp/MVP%s/MVP%s.vcf", s, s)))), "")
}
all_s <- unlist(ALL, use.names = FALSE)
rec("no seed in two landscapes", "all", length(all_s) == 1800L && !anyDuplicated(all_s), length(unique(all_s)))
# LD window table and the VIF source
LDW <- fread(file.path(ROOT, "work/mvp_land/ld_window_by_replicate.tsv"), colClasses = c(seed = "character"))
for (i in seq_len(nrow(LANDS)))
    rec("LD window covers the seeds", LANDS$land[i], setequal(LDW[landscape == LANDS$lkey[i], seed], ALL[[LANDS$land[i]]]), "")
rec("declared 619 bp", "Est-Clines", LDW[distance_source == "declared_median_619", .N] == 39L &&
        LDW[distance_source == "declared_median_619", all(clumping_distance_bp == 619 & landscape == "estclines")],
    LDW[distance_source == "declared_median_619", .N])
VD <- fread(file.path(EVAL, "remeasure600/onefit_diag/rda_diag.tsv"), colClasses = c(seed = "character"))
rec("SS-Clines VIF source", "SS-Clines", setequal(VD[key == "max_vif", seed], ALL[["SS-Clines"]]), "remeasure600/onefit_diag/rda_diag.tsv")
for (s in c("1231114", "1231189", "1231264", "1231669", "1231805")) {
    f <- Sys.glob(file.path(ROOT, sprintf("_archive/mvp_ssclines/results/MVP%shdrc_results/GEA/tables/methods/RDA/RDA_diagnostics_K*.tsv", s)))
    D <- fread(f, sep = "\t", header = TRUE, fill = TRUE)
    a <- D[key == "max_vif", value]; b <- VD[seed == s & key == "max_vif", value]
    rec("VIF source = hdrc tree", s, identical(as.character(a), as.character(b)), paste(a, b))
}
# the on-disk tables are the published Data S2 members
ZIP <- file.path(EVAL, "figures_ld_ms/Data_S2.zip")
sha <- function(cmd) strsplit(system(cmd, intern = TRUE), " +")[[1]][1]
for (m in c("figures_ld_gallery14/detection_per_seed.tsv", "figures_ld_offset14/rank/rank_long.tsv",
            "figures_ld_offset14/rank/rank_summary.tsv", "figures_ld_gallery14/offset_delta_per_seed.tsv",
            "figures_ld_gallery14/C_delta_stats.tsv", "figures_ld_gallery14/B4_aucpr_per_seed.tsv",
            "remeasure600/rdaunc_rda2x_ldb/calls_per_method.tsv", "remeasure600/rdaunc_rda2x_ldb/calls_combine.tsv",
            "remeasure600/rdaunc/rank_metrics.tsv", "remeasure600/onefit/rank_metrics.tsv")) {
    a <- sha(sprintf("unzip -p %s %s | sha256sum", shQuote(ZIP), shQuote(m)))
    b <- sha(sprintf("sha256sum %s", shQuote(file.path(EVAL, m))))
    rec("disk = Data_S2 member", "SS-Clines", identical(a, b), m)
}
R <- rbindlist(R)
fwrite(R, file.path(OUT, "setup_check.tsv"), sep = "\t")
if (!all(R$ok)) stop(sum(!R$ok), " setup check(s) failed: ", paste(unique(R[ok == FALSE, check]), collapse = "; "))
message(sprintf("setup check: %d of %d ok", sum(R$ok), nrow(R)))
