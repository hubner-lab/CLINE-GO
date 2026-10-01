#!/usr/bin/env Rscript
# =============================================================================
# mvp_onefit_diag_600.R -- gather the Phase 2b side tables of all 600 SS-Clines replicates
# into long tables journal 18 can read.
#
# WHY. The one-fit re-run (mvp_run_sweep.sh HARVEST_DIAG=1) copied three things per replicate
# that the frozen arm never kept, and journal 17 left two questions open on exactly them:
#   - LD decay (Structure/tables/ld_decay_half_distances.tsv): the agreement window has to be
#     set against measured LD, which journal 17 had for ONE replicate;
#   - the pipeline's own PCA (pca.eigenvalues + pca.tracywidom): the leading-PC variance share,
#     the fourth K-ladder signal journal 17 could not measure after the 09-23 cleanup;
#   - RDA diagnostics / anova of the one fit (gif_lambda, axes, adj R2, VIF, K_floored).
# Nothing is recomputed from genotypes: this script only reads and reshapes harvested files.
#
# WRITES (OUT_DIR):
#   ld_decay.tsv    every row of every replicate's ld_decay_half_distances.tsv, + seed. Rows the
#                   pipeline could not fit keep their NA half-distance; `method` is kept verbatim
#                   (a non-converged nls can still read hill_weir -- known, unfixed; journal 18
#                   reports the method split beside every number, never the number alone)
#   pca.tsv         the first N_PC components: eigenvalue (full precision, pca.eigenvalues),
#                   share of total variance (eigenvalue / sum of all eigenvalues), Tracy-Widom
#                   statistic and p (pca.tracywidom, LEA's rounding)
#   rda_diag.tsv    RDA_diagnostics.tsv of every replicate, long (seed, key, value)
#   rda_anova.tsv   RDA_anova.tsv of every replicate, + seed
#
# ASSERTIONS (each stops the run): every file present for all 600; exactly one `All` /
# genome_wide LD row per replicate; the TW component index runs 1..n with n equal to the
# eigenvalue count or one short of it (LEA drops the null last component on some replicates);
# gif_lambda present in every diagnostics table.
#
#   PARAMS_DIR  default benchmarks/mvp_eval/params_onefit
#   OUT_DIR     default benchmarks/mvp_eval/remeasure600/onefit_diag
#   CELL default c1;  N_PC default 20;  NCORES default 16
#   MVP_ARM / MVP_ADDED / MVP_N_EXPECT   via benchmarks/mvp_arm.R (mandatory for SS-Clines)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(parallel) })

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
EVAL <- file.path(ROOT, "benchmarks/mvp_eval")
PDIR <- Sys.getenv("PARAMS_DIR", file.path(EVAL, "params_onefit"))
OUT  <- Sys.getenv("OUT_DIR", file.path(EVAL, "remeasure600", "onefit_diag"))
CELL <- Sys.getenv("CELL", "c1")
NPC  <- as.integer(Sys.getenv("N_PC", "20"))
NCOR <- as.integer(Sys.getenv("NCORES", "16"))
source(file.path(ROOT, "benchmarks/mvp_arm.R"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# pca.tracywidom is LEA's fixed-width-ish text: the header and the rows carry DIFFERENT numbers
# of empty tab fields, so fread misaligns it. Split on whitespace runs instead.
read_tw <- function(f) {
    x <- strsplit(trimws(readLines(f)), "[[:space:]]+")
    hdr <- x[[1]]
    if (!identical(hdr, c("N", "eigenvalues", "twstats", "pvalues", "effectn", "percentage")))
        stop("unexpected tracywidom header in ", f, ": ", paste(hdr, collapse = " "))
    m <- do.call(rbind, x[-1])
    if (ncol(m) != 6L) stop("tracywidom rows do not split into 6 fields: ", f)
    data.table(pc = as.integer(m[, 1]), twstat = suppressWarnings(as.numeric(m[, 3])),
               tw_p = suppressWarnings(as.numeric(m[, 4])))
}

one <- function(seed) {
    d <- file.path(PDIR, paste0("MVP", seed), CELL)
    need <- file.path(d, c("ld_decay_half_distances.tsv", "pca.eigenvalues", "pca.tracywidom",
                           "RDA_diagnostics.tsv", "RDA_anova.tsv"))
    miss <- need[!file.exists(need)]
    if (length(miss)) stop("missing: ", paste(miss, collapse = ", "))

    ld <- fread(need[1])
    if (nrow(ld[group == "All" & scope == "genome_wide"]) != 1L)
        stop("seed ", seed, ": expected exactly one All/genome_wide LD row")
    ld[, seed := seed]

    ev <- scan(need[2], quiet = TRUE)
    tw <- read_tw(need[3])
    # LEA omits the LAST component's Tracy-Widom row on some replicates (161 of 600): its
    # eigenvalue is numerically zero (~1e-10). Only that exact shortfall is accepted.
    if (!identical(tw$pc, seq_len(nrow(tw))) || !nrow(tw) %in% (length(ev) - 0:1) || nrow(tw) < NPC)
        stop("seed ", seed, ": ", length(ev), " eigenvalues vs ", nrow(tw), " Tracy-Widom rows")
    k <- seq_len(min(NPC, length(ev)))
    pca <- data.table(seed = seed, pc = k, eigenvalue = ev[k], share = ev[k] / sum(ev),
                      twstat = tw$twstat[k], tw_p = tw$tw_p[k], n_components = length(ev),
                      n_tw_rows = nrow(tw))

    dg <- fread(need[4], colClasses = "character", sep = "\t")
    if (!identical(names(dg), c("key", "value")) || !"gif_lambda" %in% dg$key)
        stop("seed ", seed, ": RDA_diagnostics.tsv lacks key/value or gif_lambda")
    dg[, seed := seed]
    an <- fread(need[5])[, seed := seed]
    list(ld = ld, pca = pca, dg = dg, an = an)
}

man  <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
PRIM <- mvp_prim(man)
stopifnot(nrow(PRIM) == mvp_n_expect(), uniqueN(PRIM$seed) == nrow(PRIM))
SEEDS <- sort(PRIM$seed)
message(sprintf("INFO: %d replicate(s) from %s, corpus [%s]", length(SEEDS), PDIR, mvp_arm_label()))

res <- mclapply(SEEDS, function(s) tryCatch(one(s), error = function(e)
    list(error = data.table(seed = s, error = conditionMessage(e)))), mc.cores = NCOR)
bad <- vapply(res, function(r) !is.list(r) || !is.null(r$error), logical(1))
if (any(bad)) {
    errs <- rbindlist(lapply(res[bad], function(r)
        if (is.list(r)) r$error else data.table(seed = NA, error = as.character(r))))
    print(errs); fwrite(errs, file.path(OUT, "errors.tsv"), sep = "\t")
    stop(nrow(errs), " replicate(s) failed (see errors.tsv)")
}
grab <- function(k) rbindlist(lapply(res, `[[`, k), use.names = TRUE, fill = TRUE)
LD <- grab("ld"); PCA <- grab("pca"); DG <- grab("dg"); AN <- grab("an")
stopifnot(uniqueN(LD$seed) == length(SEEDS), uniqueN(PCA$seed) == length(SEEDS),
          uniqueN(DG$seed) == length(SEEDS), uniqueN(AN$seed) == length(SEEDS))

w <- function(d, f) { fwrite(d, file.path(OUT, f), sep = "\t"); message("wrote ", f, " (", nrow(d), " rows)") }
setcolorder(LD, "seed"); setcolorder(DG, "seed"); setcolorder(AN, "seed")
w(LD, "ld_decay.tsv"); w(PCA, "pca.tsv"); w(DG, "rda_diag.tsv"); w(AN, "rda_anova.tsv")
writeLines(c(sprintf("params_dir\t%s", PDIR), sprintf("cell\t%s", CELL), sprintf("n_pc\t%d", NPC),
             sprintf("corpus\t%s", mvp_arm_label()), sprintf("n_replicates\t%d", length(SEEDS)),
             sprintf("run_at\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))),
           file.path(OUT, "provenance.tsv"))
message(sprintf("OK: %d replicates -> %s", length(SEEDS), OUT))
