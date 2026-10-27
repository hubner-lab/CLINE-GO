#!/usr/bin/env Rscript
# mvp_p4_detection_b.R -- LD-window plan, Phase 4: detection of the candidate sets under option B.
#
# Option B (user, 2026-10-26, plan D2/D2b): a selected SNP is a HIT when it is a causal locus, or a
# (neutral) SNP within the replicate's own LD distance d of a causal locus on the same chromosome;
# every other selected SNP is a MISS. d is the replicate's critical-r2 distance, read from its own
# LD-decay table (group All, genome_wide, decay_critical_bp rounded -- the value the framework's
# auto_genome_wide resolves to). Many-to-one: any number of SNPs may count against one causal locus.
#
# The old definition is reported next to it for old -> new: TP = causal + the deposit's
# "linked_neutral" category (any neutral SNP on a causal-bearing linkage group, at any distance),
# FP = the rest (mvp_j16_gallery.R:199).
#
# Inputs per seed:
#   truth   data/mvp/MVP<seed>/truth_any.tsv (chr, pos, category over every SNP)
#   tested  the universe of the new tables, MVP<seed><proj>_results/GEA/tables/methods/LFMM/*_pvalues_K*.tsv
#   d       MVP<seed><proj>_results/Structure/tables/ld_decay_half_distances.tsv
#   sets    new: MVP<seed><proj>_results/_intermediate/snp_sets/gea.<set>/selected_snps.tsv
#           old: MVP<seed>_results/_intermediate/snp_sets/<set>/selected_snps.tsv (offset13, 5 kb)
#
# Writes into OUTDIR:
#   detection_b_per_seed.tsv  seed, source (new|old), set, n, n_causal, n_linked_b, n_miss_b,
#                             prec_b, n_linked_cat, prec_old_def
#   chance_per_seed.tsv       seed, d_bp, n_tested, n_causal, n_linked_b, chance_b
#                             (= share of tested SNPs that are hits), chance_old_def
#   summary.tsv               per source x set: replicates, median n / TP / FP / precision (both
#                             definitions), empty and under-3 counts
#
# Usage: Rscript benchmarks/mvp_p4_detection_b.R <seeds_csv|seed file> <project_suffix> <outdir> [ncores]

suppressPackageStartupMessages({ library(data.table); library(parallel) })
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 3)
seeds <- if (file.exists(a[1])) trimws(readLines(a[1])) else strsplit(a[1], ",")[[1]]
seeds <- seeds[nzchar(seeds)]
SFX <- a[2]; OUTDIR <- a[3]; NC <- if (length(a) >= 4) as.integer(a[4]) else 1L
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
SETS <- c("solo_lfmm", "solo_emmax", "solo_rda", "union", "best", "intersect3")

# distance of each query SNP to the nearest causal locus on its chromosome (Inf when none)
near_causal <- function(q, cz) {
    out <- rep(Inf, nrow(q))
    for (ch in unique(q$chr)) {
        p <- sort(cz[chr == ch, pos]); i <- which(q$chr == ch)
        if (!length(p) || !length(i)) next
        x <- q$pos[i]; j <- findInterval(x, p)
        lo <- ifelse(j >= 1, x - p[pmax(j, 1)], Inf)
        hi <- ifelse(j < length(p), p[pmin(j + 1, length(p))] - x, Inf)
        out[i] <- pmin(abs(lo), abs(hi))
    }
    out
}
read_set <- function(f) {
    if (!file.exists(f)) return(data.table(chr = character(), pos = integer()))  # old builder: empty set = no file
    fread(f, colClasses = c(SNPID = "character", chr = "character"))[, .(chr, pos)]
}

one_seed <- function(s) {
    proj <- paste0("MVP", s, SFX, "_results")
    tr <- fread(sprintf("data/mvp/MVP%s/truth_any.tsv", s), colClasses = c(chr = "character"))
    lt <- fread(file.path(proj, "Structure/tables/ld_decay_half_distances.tsv"))[group == "All" & scope == "genome_wide"]
    stopifnot(nrow(lt) == 1)
    d <- round(lt$decay_critical_bp)
    lf <- Sys.glob(file.path(proj, "GEA/tables/methods/LFMM/LFMM_pvalues_K*.tsv"))
    lf <- lf[!grepl("_sig_|_wza", lf)]; stopifnot(length(lf) == 1)
    tested <- fread(lf, select = c("chr", "pos"), colClasses = c(chr = "character"))
    cz <- tr[category == "causal", .(chr, pos)]
    tt <- merge(tested, tr, by = c("chr", "pos"), all.x = TRUE)
    stopifnot(!anyNA(tt$category))
    tt[, dist := near_causal(tt, cz)]
    score <- function(src, set, f) {
        x <- merge(read_set(f), tt, by = c("chr", "pos"), all.x = TRUE)
        if (anyNA(x$category)) stop(sprintf("seed %s %s/%s: selected SNP outside the tested universe", s, src, set))
        n <- nrow(x); nc <- sum(x$category == "causal")
        nb <- sum(x$category != "causal" & x$dist <= d); nl <- sum(x$category == "linked_neutral")
        data.table(seed = s, source = src, set = set, d_bp = d, n = n, n_causal = nc, n_linked_b = nb,
                   n_miss_b = n - nc - nb, prec_b = if (n) (nc + nb) / n else NA_real_,
                   n_linked_cat = nl, prec_old_def = if (n) (nc + nl) / n else NA_real_)
    }
    det <- rbindlist(c(
        lapply(SETS, function(st) score("new", st, file.path(proj, "_intermediate/snp_sets", paste0("gea.", st), "selected_snps.tsv"))),
        lapply(SETS, function(st) score("old", st, sprintf("MVP%s_results/_intermediate/snp_sets/%s/selected_snps.tsv", s, st)))))
    ch <- data.table(seed = s, d_bp = d, n_tested = nrow(tt), n_causal = sum(tt$category == "causal"),
                     n_linked_b = sum(tt$category != "causal" & tt$dist <= d),
                     chance_b = mean(tt$category == "causal" | tt$dist <= d),
                     chance_old_def = mean(tt$category %in% c("causal", "linked_neutral")))
    list(det = det, ch = ch)
}

res <- mclapply(seeds, function(s) tryCatch(one_seed(s), error = function(e) conditionMessage(e)), mc.cores = NC)
bad <- vapply(res, is.character, logical(1))
if (any(bad)) for (i in which(bad)) message("FAILED seed ", seeds[i], ": ", res[[i]])
det <- rbindlist(lapply(res[!bad], `[[`, "det")); ch <- rbindlist(lapply(res[!bad], `[[`, "ch"))
fwrite(det, file.path(OUTDIR, "detection_b_per_seed.tsv"), sep = "\t")
fwrite(ch, file.path(OUTDIR, "chance_per_seed.tsv"), sep = "\t")
sm <- det[, .(replicates = .N, median_n = median(n), median_TP_b = median(n_causal + n_linked_b),
              median_FP_b = median(n_miss_b), median_prec_b = median(prec_b, na.rm = TRUE),
              median_TP_old = median(n_causal + n_linked_cat), median_FP_old = median(n - n_causal - n_linked_cat),
              median_prec_old = median(prec_old_def, na.rm = TRUE),
              n_with_prec = sum(!is.na(prec_b)), n_empty = sum(n == 0), n_under3 = sum(n < 3)),
          by = .(source, set)]
fwrite(sm, file.path(OUTDIR, "summary.tsv"), sep = "\t")
message(sprintf("seeds scored: %d, failed: %d; median chance_b %.4f, chance_old_def %.4f",
                nrow(ch), sum(bad), median(ch$chance_b), median(ch$chance_old_def)))
if (any(bad)) quit(status = 1)
