# =============================================================================
# lib_detection.R — shared scoring core for the causal-locus detection benchmark.
#
# Sourced by benchmarks/eval_detection.R (single-config scoring) and
# benchmarks/sweep_thresholds.R (grid + combine sweep) so both apply exactly the
# same TP/FP accounting. Do not duplicate this logic in a caller.
#
# Requires: data.table, and the caller to have sourced the pipeline's
# scripts/R/utils/pval_threshold.R when threshold-based calling is used.
# =============================================================================

# ------------------------------------------------------------------ arg parsing
parse_kv_args <- function(argv) {
    out <- list()
    for (a in argv) {
        if (!grepl("^--[A-Za-z0-9._-]+=", a)) stop("Malformed argument (expected --key=value): ", a)
        out[[gsub("-", "_", sub("^--([^=]+)=.*$", "\\1", a))]] <- sub("^--[^=]+=", "", a)
    }
    out
}

# ------------------------------------------------------------------ truth set
# Returns a list with the causal table, the linked/background key vectors, and
# the count of simulated causal loci. Keys are "chr:pos".
load_truth <- function(path) {
    truth <- fread(path, colClasses = c(chr = "character"))
    stopifnot(all(c("chr", "pos", "category") %in% names(truth)))
    truth[, key := paste(chr, pos, sep = ":")]
    list(
        causal      = truth[category == "causal"],
        linked_keys = truth[category == "linked_neutral", key],
        bg_keys     = truth[category == "background_neutral", key],
        n_simulated = truth[category == "causal", .N]
    )
}

# --------------------------------------------------------------- p-value table
# Structural columns are excluded from the trait list. n_snps/mean_maf are WZA
# metadata, dropped the same way find_significant_snps_per_trait() drops them.
STRUCTURAL_COLS <- c("SNPID", "chr", "pos", "key", "n_snps", "mean_maf")

load_pvalues <- function(path, traits = "all") {
    pv <- fread(path, colClasses = c(chr = "character"))
    stopifnot(all(c("chr", "pos") %in% names(pv)))
    setDT(pv)
    pv[, key := paste(chr, pos, sep = ":")]
    tc <- setdiff(names(pv), STRUCTURAL_COLS)
    if (!identical(traits, "all")) {
        want    <- trimws(strsplit(traits, ",", fixed = TRUE)[[1]])
        missing <- setdiff(want, tc)
        if (length(missing)) stop("Requested traits absent from ", path, ": ", paste(missing, collapse = ", "))
        tc <- want
    }
    if (!length(tc)) stop("No trait columns found in ", path)
    for (x in tc) pv[[x]] <- suppressWarnings(as.numeric(pv[[x]]))
    list(pv = pv, trait_cols = tc)
}

# Genomic-control recalibration, matching what LEA's lfmm2.test(genomic.control =
# TRUE) does internally: per trait, lambda = median(chisq)/qchisq(0.5, 1), then
# p_gc = pchisq(chisq / lambda, 1, lower = FALSE).
#
# PreGEA deliberately runs LFMM with genomic_control = FALSE so lambda stays
# informative across the K ladder, while production mode=gea uses TRUE. Applying
# this here lets both be scored from the same ladder output.
#
# Note it is a MONOTONE transform within a trait, so single-trait rankings are
# unchanged. It matters for (a) any threshold-based call and (b) multi-trait
# min-p ranking, where per-trait lambdas differ and so the cross-trait ordering
# does move.
apply_genomic_control <- function(pv, trait_cols, verbose = TRUE) {
    for (tc in trait_cols) {
        p  <- pv[[tc]]
        ok <- !is.na(p) & p > 0 & p <= 1
        if (!any(ok)) next
        chi    <- qchisq(p[ok], df = 1, lower.tail = FALSE)
        lambda <- median(chi) / qchisq(0.5, df = 1)
        if (!is.finite(lambda) || lambda <= 0) next
        newp     <- p
        newp[ok] <- pchisq(chi / lambda, df = 1, lower.tail = FALSE)
        pv[[tc]] <- newp
        if (verbose) message(sprintf("INFO: genomic control %s: lambda = %.4f", tc, lambda))
    }
    pv
}

# Causal loci actually present in this p-value table with a non-NA p-value for at
# least one scored trait — the honest denominator for "recall among testable".
testable_causal <- function(pv, trait_cols, truth) {
    any_tested <- Reduce(`|`, lapply(trait_cols, function(x) !is.na(pv[[x]])))
    intersect(truth$causal$key, pv$key[any_tested])
}

# ------------------------------------------------------------- window matching
# Exact match when window_bp == 0. Otherwise greedy nearest-first one-to-one
# assignment, so TP can never exceed either the called set or the causal set.
# Returns a list: $calls = called keys that matched, $causal = causal keys they
# matched to. The two have equal length (one-to-one), but they are DIFFERENT key
# sets whenever window_bp > 0, and the caller needs both: a call may match a
# causal locus that never entered the tested SNP set (e.g. one of the 74 Láruson
# QTNs dropped by MAF), which is a real detection via a linked marker but must
# not be counted against the "testable" denominator.
match_calls <- function(called_keys, window_bp, truth) {
    empty <- list(calls = character(0), causal = character(0))
    if (!length(called_keys)) return(empty)
    if (window_bp == 0) {
        hit <- intersect(called_keys, truth$causal$key)
        return(list(calls = hit, causal = hit))
    }

    # NB: the column is `snp`, not `key` — data.table(key = ...) is interpreted as
    # the constructor's key= argument and would silently create no column at all.
    calls <- data.table(snp = called_keys)
    calls[, c("chr", "pos") := tstrsplit(snp, ":", fixed = TRUE)]
    calls[, pos := as.numeric(pos)]

    pairs <- list()
    for (cc in unique(calls$chr)) {
        a <- calls[chr == cc]; b <- truth$causal[chr == cc]
        if (!nrow(b)) next
        for (i in seq_len(nrow(a))) {
            d <- abs(b$pos - a$pos[i]); j <- which(d <= window_bp)
            if (length(j)) pairs[[length(pairs) + 1L]] <-
                data.table(call_key = a$snp[i], causal_key = b$key[j], dist = d[j])
        }
    }
    if (!length(pairs)) return(empty)
    p <- rbindlist(pairs)[order(dist)]
    used_call <- character(0); used_causal <- character(0)
    for (r in seq_len(nrow(p))) {
        if (p$call_key[r] %in% used_call || p$causal_key[r] %in% used_causal) next
        used_call   <- c(used_call, p$call_key[r])
        used_causal <- c(used_causal, p$causal_key[r])
    }
    list(calls = used_call, causal = used_causal)
}

# ------------------------------------------------------------------- scoring
# `linked_neutral` hits are reported as expected_linked and are NEVER counted as
# false positives — they are physically linked to true QTNs. Only
# `background_neutral` hits are true FP.
# `testable_keys` is the set of causal loci actually present in the tested SNP
# table. Recall is reported against two denominators and they use DIFFERENT
# numerators, which only coincide at window 0:
#   recall_testable  = causal loci matched that were testable / |testable|
#   recall_simulated = ALL causal loci matched / 102
# Precision counts any matched causal locus as a true positive — detecting an
# untested QTN through a linked marker is a success, not a false positive.
score_calls <- function(called_keys, window_kb, truth, testable_keys) {
    called_keys <- unique(called_keys)
    n_testable  <- length(testable_keys)
    m           <- match_calls(called_keys, window_kb * 1000, truth)
    tp_any      <- length(m$causal)
    tp_testable <- sum(m$causal %in% testable_keys)
    non_tp      <- setdiff(called_keys, m$calls)
    expected_linked <- sum(non_tp %in% truth$linked_keys)
    fp_background   <- sum(non_tp %in% truth$bg_keys)

    precision_strict <- if (tp_any + fp_background > 0) tp_any / (tp_any + fp_background) else NA_real_
    precision_all    <- if (length(called_keys) > 0) tp_any / length(called_keys) else NA_real_
    recall_testable  <- if (n_testable > 0) tp_testable / n_testable else NA_real_
    recall_simulated <- tp_any / truth$n_simulated
    f1 <- if (!is.na(precision_strict) && !is.na(recall_testable) &&
              (precision_strict + recall_testable) > 0) {
        2 * precision_strict * recall_testable / (precision_strict + recall_testable)
    } else NA_real_

    # How many SNPs a user must inspect per real hit. 1 / precision_all, reported
    # separately because it is the quantity a user actually feels when reading a
    # Manhattan plot, and because it stays finite and interpretable when precision
    # is small (precision_all 0.02 reads as "1 in 50"; calls_per_tp 50 reads as
    # "fifty rows to find one"). Inf when a rule called SNPs but hit nothing.
    calls_per_tp <- if (length(called_keys) == 0) NA_real_
                    else if (tp_any > 0) length(called_keys) / tp_any else Inf

    list(n_called = length(called_keys), tp = tp_testable, tp_any_causal = tp_any,
         expected_linked = expected_linked, fp_background = fp_background,
         untracked = length(non_tp) - expected_linked - fp_background,
         fn_testable = n_testable - tp_testable,
         precision_strict = precision_strict, precision_all = precision_all,
         recall_testable = recall_testable, recall_simulated = recall_simulated, f1 = f1,
         calls_per_tp = calls_per_tp)
}

# --------------------------------------------------------- cross-method agreement
# For every SNP called by ANY method, how many DISTINCT methods placed a call within
# window_bp of it. Union / at_least_k / intersection are then thresholds on that count.
#
# Why a window at all: causal loci sit inside LD blocks. On MVP seed 1232548 the
# linked_neutral SNPs sit a median 2 351 bp from their nearest causal locus (p75 4 496,
# p90 6 819), so two methods that flag adjacent SNPs of one haplotype block are
# corroborating each other, not disagreeing. Exact-key set intersection scores that as a
# disagreement and systematically understates every combine rule.
#
# This mirrors what the PIPELINE already does — scripts/R/lib/combine_sigsnps.R's
# .overlap_cross_method() runs foverlaps over GEA.snp_clumping_distance (5 000 bp on this
# config) — so the harness stops being stricter than the code it benchmarks.
#
# At window_bp == 0 this reduces exactly to table(all_keys) over the pooled call vectors,
# which is the pre-existing behaviour and is asserted as a unit test by the caller.
#
# The RETAINED unit stays the SNP, not the window: at window > 0 an "intersection" keeps
# every SNP that has corroboration nearby, so each method contributes its own position.
# That is deliberate and matches the pipeline; it is not the exact key intersection.
combine_support <- function(parts, window_bp) {
    empty <- data.table(snp = character(0), n_methods = integer(0))
    parts <- parts[vapply(parts, length, integer(1)) > 0L]
    if (!length(parts)) return(empty)

    calls <- rbindlist(lapply(names(parts), function(m)
        data.table(method = m, snp = unique(parts[[m]]))))
    if (window_bp == 0) return(calls[, .(n_methods = uniqueN(method)), by = snp])

    calls[, c("chr", "pos") := tstrsplit(snp, ":", fixed = TRUE)]
    calls[, pos := as.numeric(pos)]
    # The query column is named qsnp, not snp: foverlaps only prefixes x-table columns
    # with "i." when the name COLLIDES with one in y, so a shared name would make the
    # output column name depend on what happens to be in y. Renaming makes it fixed.
    q <- unique(calls[, .(qsnp = snp, chr, pos)])[
        , `:=`(start = pos - window_bp, end = pos + window_bp)]
    s <- calls[, .(method, chr, start = pos, end = pos)]
    setkey(s, chr, start, end)
    ov <- foverlaps(q, s, by.x = c("chr", "start", "end"), nomatch = 0L)
    if (!nrow(ov)) return(empty)
    ov[, .(n_methods = uniqueN(method)), by = .(snp = qsnp)]
}

# --------------------------------------------------------------- calling rules
# Mirrors find_significant_snps_per_trait(): threshold computed PER TRAIT, mask
# is strict `<` with NA -> FALSE, positions unioned across traits.
call_by_threshold <- function(pv, trait_cols, adjust, value, quiet = TRUE) {
    called <- character(0); info <- list()
    for (tc in trait_cols) {
        p  <- pv[[tc]]
        th <- if (quiet) suppressMessages(compute_pval_threshold(p, adjust, as.numeric(value)))
              else compute_pval_threshold(p, adjust, as.numeric(value))
        info[[tc]] <- th
        if (identical(th$status, "ok") && !is.na(th$threshold)) {
            called <- c(called, pv$key[!is.na(p) & p < th$threshold])
        }
    }
    list(called = unique(called), info = info)
}

# Method-level ranking score = min p across scored traits. Anti-conservative as a
# calibrated p-value, but this is used only to ORDER SNPs when comparing ladder
# rungs independently of any threshold, which is exactly what it is valid for.
rank_keys <- function(pv, trait_cols) {
    score <- do.call(pmin, c(lapply(trait_cols, function(x) pv[[x]]), list(na.rm = TRUE)))
    ok    <- is.finite(score)
    list(keys = pv$key[ok][order(score[ok])], score = sort(score[ok]))
}

# ------------------------------------------------------------------- AUC-PR
# Threshold-free ranking quality, by the interpolation-free step rule: mean precision at
# each retrieved causal locus, with un-retrieved ones contributing zero. `background_neutral`
# hits are the only false positives, matching score_calls().
#
# This is a VERBATIM extraction of the block inlined in eval_detection.R's mode=rank
# (lines ~110-124). eval_detection.R deliberately still runs its own copy: rewriting it to
# call this would put journal 05's published AUC-PR numbers behind a refactor with nothing
# forcing the two to agree. New callers use this one.
auc_pr_from_rank <- function(ranked_keys, truth, n_testable) {
    if (!length(ranked_keys) || n_testable == 0) return(NA_real_)
    is_causal <- ranked_keys %in% truth$causal$key
    is_bg     <- ranked_keys %in% truth$bg_keys
    tp_cum <- cumsum(is_causal); fp_cum <- cumsum(is_bg)
    prec   <- tp_cum / (tp_cum + fp_cum); prec[is.nan(prec)] <- 1
    sum(prec[is_causal]) / n_testable
}

# ------------------------------------------------------------------- plotting
# scale_color_clinego() wraps CLINEGO_CATEGORICAL, the Okabe-Ito colourblind-safe palette,
# which has exactly EIGHT entries. Journal 06 plots up to 13 series (11 methods plus two
# combine pseudo-methods), and ggplot2 aborts rather than recycling:
#   "Insufficient values in manual scale. 11 needed but only 8 provided."
# So use the pipeline palette while it is sufficient -- keeping colour meaning consistent
# with every other figure in the project -- and fall back to evenly-spaced hues past eight.
# The fallback is NOT colourblind-safe at 11+ levels; nothing is, which is why these figures
# also carry direct labels or a table beside them.
#
# Caller passes the number of levels rather than the data, so it works the same for a factor,
# a character column, or a dcast'd wide frame.
scale_colour_n <- function(n) {
    if (n <= length(CLINEGO_CATEGORICAL)) scale_color_clinego()
    else ggplot2::scale_colour_hue(h = c(0, 330), c = 90, l = 60)
}
scale_fill_n <- function(n) {
    if (n <= length(CLINEGO_CATEGORICAL)) scale_fill_clinego()
    else ggplot2::scale_fill_hue(h = c(0, 330), c = 90, l = 60)
}

# ------------------------------------------------------------------- emitting
new_emitter <- function() {
    rows <- list()
    list(
        emit = function(config, category, metric, value) {
            rows[[length(rows) + 1L]] <<- data.table(
                axis = "detection", config = config, category = category,
                metric = metric, value = as.numeric(value))
        },
        emit_scored = function(config, s, window_kb) {
            cat0 <- paste0("window_", window_kb, "kb")
            for (nm in names(s)) rows[[length(rows) + 1L]] <<- data.table(
                axis = "detection", config = config, category = cat0,
                metric = nm, value = as.numeric(s[[nm]]))
        },
        collect = function() rbindlist(rows)
    )
}

append_eval <- function(dt, out_tsv) {
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
    exists_already <- file.exists(out_tsv)
    fwrite(dt, out_tsv, sep = "\t", append = exists_already, col.names = !exists_already)
    message("Appended ", nrow(dt), " rows to ", out_tsv)
}
