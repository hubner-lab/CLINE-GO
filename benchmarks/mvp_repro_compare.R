#!/usr/bin/env Rscript
# mvp_repro_compare.R -- is a re-run's p-value table the same as the frozen pmax arm's?
#
#   Rscript benchmarks/mvp_repro_compare.R --seeds=1231418[,...] \
#       [--ref=benchmarks/mvp_eval/params] [--new=benchmarks/mvp_eval/params_repro] \
#       [--methods=LFMM,EMMAX,RDA] [--out=benchmarks/mvp_eval/params_repro/repro_compare.tsv]
#
# SS-Clines re-analysis, Phase 0 step 6 (docs/gea-simulation-reanalysis.md). Phase 3 measures
# what one-fit RDA changes by comparing a re-run against the frozen tables, which is only a
# before/after design if an UNCHANGED re-run reproduces them. One row per (seed, method):
#   identical      byte-identical files (cmp) -- the pass condition
#   same_keys      same chr:pos key set, in the same order
#   max_abs_dp     max |p_new - p_ref| over every trait column (NA only if keys differ)
#   max_rel_dp     max |dp| / p_ref over p_ref > 0 -- floating-point noise is ~1e-12 here
#   n_top100_diff  per trait, SNPs in one top-100 set but not the other, summed over traits
#   spearman_min   min over traits of the rank correlation
# Exits 1 if any row is not identical, after writing the table: a non-identical row is a
# result to read, not a crash.

suppressPackageStartupMessages(library(data.table))

A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
ROOT    <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
opt     <- function(k, d) if (is.null(A[[k]])) d else A[[k]]
SEEDS   <- strsplit(opt("seeds", stop("--seeds is required")), ",", fixed = TRUE)[[1]]
REF     <- file.path(ROOT, opt("ref", "benchmarks/mvp_eval/params"))
NEW     <- file.path(ROOT, opt("new", "benchmarks/mvp_eval/params_repro"))
METHODS <- strsplit(opt("methods", "LFMM,EMMAX,RDA"), ",", fixed = TRUE)[[1]]
OUT     <- file.path(ROOT, opt("out", "benchmarks/mvp_eval/params_repro/repro_compare.tsv"))

STRUCT <- c("SNPID", "chr", "pos")

compare_one <- function(seed, m) {
    fr <- file.path(REF, paste0("MVP", seed), "c1", paste0(m, "_pvalues.tsv"))
    fn <- file.path(NEW, paste0("MVP", seed), "c1", paste0(m, "_pvalues.tsv"))
    if (!file.exists(fr) || !file.exists(fn))
        stop("missing table: ", if (!file.exists(fr)) fr else fn)
    identical_bytes <- unname(tools::md5sum(fr) == tools::md5sum(fn))
    r <- fread(fr, colClasses = c(chr = "character"))
    n <- fread(fn, colClasses = c(chr = "character"))
    tc <- setdiff(names(r), STRUCT)
    if (!identical(tc, setdiff(names(n), STRUCT)))
        stop(seed, " ", m, ": trait columns differ: ", paste(tc, collapse = ","), " vs ",
             paste(setdiff(names(n), STRUCT), collapse = ","))
    kr <- paste(r$chr, r$pos, sep = ":"); kn <- paste(n$chr, n$pos, sep = ":")
    same_keys <- identical(kr, kn)
    max_abs <- max_rel <- spear <- NA_real_; top_diff <- NA_integer_
    if (same_keys) {
        dp  <- unlist(lapply(tc, function(x) abs(n[[x]] - r[[x]])))
        pr  <- unlist(lapply(tc, function(x) r[[x]]))
        max_abs  <- max(dp, na.rm = TRUE)
        max_rel  <- max((dp / pr)[is.finite(pr) & pr > 0], na.rm = TRUE)
        spear    <- min(vapply(tc, function(x) cor(r[[x]], n[[x]], method = "spearman",
                                                   use = "complete.obs"), numeric(1)))
        top_diff <- sum(vapply(tc, function(x) {
            a <- kr[order(r[[x]])][1:100]; b <- kn[order(n[[x]])][1:100]
            length(setdiff(a, b)) + length(setdiff(b, a))
        }, integer(1)))
    }
    data.table(seed = seed, method = m, identical = identical_bytes, same_keys = same_keys,
               n_snps_ref = nrow(r), n_snps_new = nrow(n), max_abs_dp = max_abs,
               max_rel_dp = max_rel, n_top100_diff = top_diff, spearman_min = spear)
}

res <- rbindlist(lapply(SEEDS, function(s) rbindlist(lapply(METHODS, function(m) compare_one(s, m)))))
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
fwrite(res, OUT, sep = "\t")
print(res)
message("wrote ", OUT)
if (!all(res$identical)) {
    message("NOT REPRODUCED: ", sum(!res$identical), " of ", nrow(res), " (seed, method) tables differ")
    quit(status = 1)
}
message("REPRODUCED: all ", nrow(res), " tables byte-identical")
