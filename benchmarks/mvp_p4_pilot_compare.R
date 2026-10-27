#!/usr/bin/env Rscript
# mvp_p4_pilot_compare.R -- LD-window plan, Phase 4 pilot GATE: do fresh framework runs reproduce
# the 30 Sep p-value tables SNP for SNP?
#
#   MVP<seed>hd    LFMM, EMMAX   vs benchmarks/mvp_eval/params_onefit/MVP<seed>/c1/<M>_pvalues.tsv
#                  RDA (uncorr.) vs benchmarks/mvp_eval/params_rdaunc/MVP<seed>/c1/RDA_pvalues.tsv
#   MVP<seed>hdrc  RDA (cond.)   vs benchmarks/mvp_eval/params_onefit/MVP<seed>/c1/RDA_pvalues.tsv
#
# Per method and seed: SNPs old / new / compared, max |dp|, max |d log10 p|, and whether the top
# sets are identical (top 0.25 % per predictor, RDA's single column 0.5 %, inclusive ties: the
# rule of the declared sets). Also: every causal locus (data/mvp/MVP<seed>/truth_any.tsv,
# category causal) is among the tested SNPs of every new table, and the seed's LD distance
# (build log of gea.best) equals out600's decay_critical_bp.
#
# GATE: pass when every table has identical SNPs, max |dp| == 0 or < 1e-12, and identical top
# sets; written to gate.tsv (one row) next to the tables.
#
# Usage: Rscript benchmarks/mvp_p4_pilot_compare.R <seeds_csv> <out600_dir> <outdir>

suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 3)
seeds <- strsplit(a[1], ",")[[1]]; OUT600 <- a[2]; OUTDIR <- a[3]
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
EV <- "benchmarks/mvp_eval"

rd <- function(f) fread(f, colClasses = c(SNPID = "character", chr = "character"))
top_ids <- function(d, share) {
    cols <- setdiff(names(d), c("SNPID", "chr", "pos"))
    N <- ceiling(share * nrow(d) - 1e-9)
    unique(unlist(lapply(cols, function(cl) { cut <- sort(d[[cl]])[N]; d[get(cl) <= cut, SNPID] })))
}
kbest <- function(s) {
    y <- readLines(sprintf("config_MVP%shd.yaml", s))
    i <- grep("^sNMF:", y); as.integer(sub(".*k_best: *", "", grep("k_best:", y[i:length(y)], value = TRUE)[1]))
}

rows <- list(); causal <- list(); ld <- list()
for (s in seeds) {
    k <- kbest(s)
    pairs <- list(
        list(m = "LFMM",  proj = "hd",   old = file.path(EV, "params_onefit", paste0("MVP", s), "c1", "LFMM_pvalues.tsv")),
        list(m = "EMMAX", proj = "hd",   old = file.path(EV, "params_onefit", paste0("MVP", s), "c1", "EMMAX_pvalues.tsv")),
        list(m = "RDA",   proj = "hd",   old = file.path(EV, "params_rdaunc", paste0("MVP", s), "c1", "RDA_pvalues.tsv")),
        list(m = "RDA",   proj = "hdrc", old = file.path(EV, "params_onefit", paste0("MVP", s), "c1", "RDA_pvalues.tsv")))
    tr <- fread(sprintf("data/mvp/MVP%s/truth_any.tsv", s), colClasses = c(chr = "character"))
    cz <- tr[category == "causal", paste(chr, pos, sep = ":")]
    for (p in pairs) {
        newf <- sprintf("MVP%s%s_results/GEA/tables/methods/%s/%s_pvalues_K%d.tsv", s, p$proj, p$m, p$m, k)
        if (!file.exists(newf)) {
            rows[[length(rows) + 1]] <- data.table(seed = s, project = p$proj, method = p$m, status = "missing_new")
            next
        }
        o <- rd(p$old); n <- rd(newf)
        cols <- setdiff(names(o), c("SNPID", "chr", "pos"))
        stopifnot(identical(sort(cols), sort(setdiff(names(n), c("SNPID", "chr", "pos")))))
        m <- merge(o, n, by = "SNPID", suffixes = c(".o", ".n"))
        dp <- unlist(lapply(cols, function(cl) abs(m[[paste0(cl, ".o")]] - m[[paste0(cl, ".n")]])))
        dl <- unlist(lapply(cols, function(cl) abs(log10(m[[paste0(cl, ".o")]]) - log10(m[[paste0(cl, ".n")]]))))
        share <- if (p$m == "RDA") 0.005 else 0.0025
        to <- top_ids(o, share); tn <- top_ids(n, share)
        rows[[length(rows) + 1]] <- data.table(
            seed = s, project = p$proj, method = p$m, status = "ok", old_file = p$old,
            n_old = nrow(o), n_new = nrow(n), n_compared = nrow(m),
            same_snps = setequal(o$SNPID, n$SNPID),
            max_abs_dp = max(dp, na.rm = TRUE), max_abs_dlog10p = max(dl[is.finite(dl)], na.rm = TRUE),
            n_na_mismatch = sum(is.na(dp)),
            top_n_old = length(to), top_n_new = length(tn), top_identical = setequal(to, tn))
        causal[[length(causal) + 1]] <- data.table(seed = s, project = p$proj, method = p$m,
            n_causal = length(cz), n_causal_tested = sum(cz %in% n$SNPID),
            all_tested = all(cz %in% n$SNPID))
    }
    lg <- sprintf("MVP%shd_logs/snp_sets/build_gea.best.log", s)
    bp <- if (file.exists(lg)) {
        x <- grep("AUTO genome-wide distance \\(critical\\)", readLines(lg), value = TRUE)[1]
        as.integer(sub(".*\\(critical\\): ([0-9]+) bp.*", "\\1", x))
    } else NA_integer_
    t6 <- fread(file.path(OUT600, paste0("MVP", s), "table.tsv"))[group == "All" & scope == "genome_wide"]
    ld[[length(ld) + 1]] <- data.table(seed = s, resolved_bp = bp, out600_decay_critical_bp = t6$decay_critical_bp,
                                       equal = isTRUE(bp == round(t6$decay_critical_bp)))
}
rows <- rbindlist(rows, fill = TRUE); causal <- rbindlist(causal); ld <- rbindlist(ld)
fwrite(rows, file.path(OUTDIR, "pilot_pvalue_identity.tsv"), sep = "\t")
fwrite(causal, file.path(OUTDIR, "pilot_causal_tested.tsv"), sep = "\t")
fwrite(ld, file.path(OUTDIR, "pilot_ld_distances.tsv"), sep = "\t")
ok <- rows[, all(status == "ok")] &&
      rows[, all(same_snps & n_na_mismatch == 0 & (max_abs_dp == 0 | max_abs_dp < 1e-12) & top_identical)]
gate <- data.table(gate = if (isTRUE(ok)) "PASS" else "FAIL", n_tables = nrow(rows),
                   all_causal_tested = all(causal$all_tested), ld_equal = sum(ld$equal), n_seeds = nrow(ld))
fwrite(gate, file.path(OUTDIR, "gate.tsv"), sep = "\t")
print(rows[, .(seed, project, method, n_compared, same_snps, max_abs_dp, max_abs_dlog10p, top_identical)])
print(gate)
