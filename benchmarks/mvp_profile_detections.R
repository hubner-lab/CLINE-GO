#!/usr/bin/env Rscript
# =============================================================================
# mvp_profile_detections.R -- WHICH causal loci get found, and which get missed.
#
# The threshold/combine sweep answers "how many". This answers "which ones", by
# joining every causal locus's detection status to the properties the SIMULATOR
# assigned it. That converts an aggregate recall number into a statement a user can
# act on: e.g. "half the causal loci on this replicate sit below MAF 0.05, and those
# are the ones every method misses" is actionable; "recall was 0.31" is not.
#
# It also replaces the MAF re-filter arm. Re-running the whole pipeline at
# Filter.maf 0.05 would answer "what is recall at MAF 0.05" with a new project tree,
# a fresh sNMF, and a truth set that has to be re-derived anyway. Stratifying the
# EXISTING detections by each locus's own MAF answers the same question -- which loci
# a MAF filter would have removed, and whether we were finding them at all -- at zero
# additional compute, and it answers the sharper version (is detection a FUNCTION of
# MAF, or is MAF just correlated with effect size?) that a second filtered run cannot.
#
# COVARIATES come from the deposit's own mutation table (67 columns), not from
# anything we computed: effect sizes, per-locus Va, allele-frequency clines, Fst.
# The join goes truth -> loci_freq.tsv (the converter's own chr/pos <-> mutname map)
# -> raw mutations. Going straight to the raw file would mean re-deriving the
# coordinate transform here, which is exactly the kind of duplicated mapping that
# silently drifts.
#
# Usage:
#   Rscript mvp_profile_detections.R --methods="EMMAX=/p/a.tsv,..." \
#       --truth=FILE --loci-freq=FILE --muts=FILE.gz --tag=STR [--outdir=DIR] \
#       [--points=EMMAX:bonf_0.05,...] [--combine-window-kb=5] [--near-kb=1]
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(qvalue)
    library(ggplot2)
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(PIPELINE_ROOT, "scripts/R/utils/pval_threshold.R"))
source(file.path(PIPELINE_ROOT, "scripts/R/utils/theme_clinego.R"))
source(file.path(PIPELINE_ROOT, "benchmarks/lib_detection.R"))

args <- parse_kv_args(commandArgs(trailingOnly = TRUE))
req  <- function(k) {
    if (is.null(args[[k]]) || !nzchar(args[[k]])) stop("Missing required --", gsub("-", "_", k))
    args[[k]]
}
opt  <- function(k, d) if (is.null(args[[k]]) || !nzchar(args[[k]])) d else args[[k]]

TAG        <- req("tag")
TRUTH_F    <- req("truth")
FREQ_F     <- req("loci_freq")
MUTS_F     <- req("muts")
OUTDIR     <- opt("outdir", file.path(PIPELINE_ROOT, "benchmarks/mvp_eval/profile"))
CW_KB      <- as.numeric(opt("combine_window_kb", "5"))
NEAR_KB    <- as.numeric(opt("near_kb", "1"))
LG_LEN     <- as.integer(opt("lg_len", "50000"))

spec <- strsplit(req("methods"), ",", fixed = TRUE)[[1]]
methods <- list()
for (s in spec) {
    kv <- strsplit(s, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) stop("Bad --methods entry (expected NAME=PATH): ", s)
    if (!file.exists(kv[2])) stop("Method p-value table not found: ", kv[2])
    methods[[kv[1]]] <- kv[2]
}

# Operating points, one per method. Default: the pipeline's own bonf 0.05 on every
# method -- the configuration a user gets without touching anything.
POINTS <- local({
    p <- setNames(rep("bonf_0.05", length(methods)), names(methods))
    if (!is.null(args$points) && nzchar(args$points)) {
        for (s in strsplit(args$points, ",", fixed = TRUE)[[1]]) {
            kv <- strsplit(s, ":", fixed = TRUE)[[1]]
            if (length(kv) != 2) stop("Bad --points entry (expected METHOD:adjust_value): ", s)
            if (!kv[1] %in% names(methods)) stop("--points names an unloaded method: ", kv[1])
            p[[kv[1]]] <- kv[2]
        }
    }
    p
})

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
truth <- load_truth(TRUTH_F)
causal <- copy(truth$causal)
message("INFO: ", nrow(causal), " causal loci in ", basename(TRUTH_F))

# ------------------------------------------------------------------ covariates
# Step 1: chr/pos -> mutname, via the converter's own table.
freq <- fread(FREQ_F, colClasses = c(chr = "character"))
stopifnot(all(c("chr", "pos", "mutname", "maf", "a_freq_subset") %in% names(freq)))
freq[, key := paste(chr, pos, sep = ":")]

# Step 2: mutname -> the deposit's 67-column mutation record.
muts <- fread(cmd = paste("zcat", shQuote(MUTS_F)))
setnames(muts, gsub('"', "", names(muts)))
if (!"mutname" %in% names(muts)) stop("mutations table has no `mutname` column: ", MUTS_F)
muts[, mutname := gsub('"', "", as.character(mutname))]

# a_freq_subset deliberately NOT taken from here: loci_freq.tsv already carries it, and
# pulling the same name from both tables would silently produce .x/.y columns that every
# downstream reference then has to guess between.
COVARS <- c("mutTempEffect", "mutSalEffect", "Va_temp", "Va_sal",
            "Va_temp_prop", "Va_sal_prop", "af_cor_temp", "af_slope_temp",
            "af_cor_sal", "af_slope_sal", "He_outflank", "Fst_outflank", "LG")
have <- intersect(COVARS, names(muts))
missing_cov <- setdiff(COVARS, have)
if (length(missing_cov)) message("WARN: covariates absent from mutations table: ",
                                 paste(missing_cov, collapse = ", "))

loci <- merge(causal[, .(key, chr, pos)], freq[, .(key, mutname, maf, a_freq_subset)],
              by = "key", all.x = TRUE)
# JOIN GATE. A causal locus with no covariate row would silently drop out of every
# stratified table and quietly shrink the denominator -- the same class of defect the
# truth-join gate guards against upstream. Fatal, not NA.
if (anyNA(loci$mutname))
    stop("FATAL: ", sum(is.na(loci$mutname)), " of ", nrow(loci),
         " causal loci have no row in ", basename(FREQ_F), " -- coordinate mapping drifted.")
n_before <- nrow(loci)
loci <- merge(loci, muts[, c("mutname", have), with = FALSE], by = "mutname", all.x = TRUE)
if (nrow(loci) != n_before)
    stop("FATAL: covariate join changed the causal-locus count: ", n_before, " -> ", nrow(loci),
         " -- mutname is not unique in ", basename(MUTS_F))
n_nocov <- sum(is.na(loci[[have[1]]]))
if (n_nocov) message("WARN: ", n_nocov, " causal loci have NA covariates after the mutations join")

# LD context, derived from the truth table itself rather than from any LD statistic:
# how isolated a causal locus is, and how much linked material surrounds it.
#
# Both loops pull the comparison values into LOCAL variables first. Writing
# `causal[chr == chr[i]]` would resolve BOTH sides to causal's own column and match every
# row -- the same data.table self-comparison trap mvp_write_sweep_configs.R documents for
# its methods_for() filter.
tr <- fread(TRUTH_F, colClasses = c(chr = "character"))
loci[, dist_nearest_causal := vapply(seq_len(.N), function(i) {
    cc <- chr[i]; pp <- pos[i]; kk <- key[i]
    o  <- causal[chr == cc & key != kk, pos]
    if (length(o)) min(abs(o - pp)) else NA_real_
}, numeric(1))]
loci[, n_linked_5kb := vapply(seq_len(.N), function(i) {
    cc <- chr[i]; pp <- pos[i]
    tr[category == "linked_neutral" & chr == cc & abs(pos - pp) <= 5000, .N]
}, integer(1))]

# ------------------------------------------------------ detections per method
loaded <- list(); called <- list()
for (m in names(methods)) {
    lp <- load_pvalues(methods[[m]], "all")
    loaded[[m]] <- lp
    ap <- strsplit(POINTS[[m]], "_", fixed = TRUE)[[1]]
    adj <- paste(head(ap, -1), collapse = "_"); val <- tail(ap, 1)
    called[[m]] <- call_by_threshold(lp$pv, lp$trait_cols, adj, val, quiet = TRUE)$called
    message(sprintf("INFO: %-8s %s -> %d calls", m, POINTS[[m]], length(called[[m]])))
}

# Combine rules enter as pseudo-methods, because the question "does combining rescue the
# loci no single method finds" cannot be asked of single methods alone.
support <- combine_support(called, CW_KB * 1000)
called[["COMBINE_union"]]      <- support$snp
called[["COMBINE_at_least_2"]] <- support[n_methods >= 2, snp]

# `detected` is exact-key. `detected_near` allows a call within NEAR_KB of the locus --
# reported separately and never as the headline, because on linkage groups 1-10 almost
# every non-causal SNP is linked_neutral, so a generous window turns nearly any call on
# those groups into a "detection".
# The column is `locus_key`, not `key`: data.table(key = ...) is swallowed by the
# constructor's key= argument and creates no column at all -- the same trap
# lib_detection.R:match_calls() and mvp_write_sweep_configs.R both carry a note about.
det <- rbindlist(lapply(names(called), function(m) {
    ck <- called[[m]]
    if (!length(ck)) return(data.table(method = m, locus_key = loci$key,
                                       detected = FALSE, detected_near = FALSE))
    cd <- data.table(snp = ck)
    cd[, c("chr", "pos") := tstrsplit(snp, ":", fixed = TRUE)][, pos := as.numeric(pos)]
    near <- vapply(seq_len(nrow(loci)), function(i)
        any(cd$chr == loci$chr[i] & abs(cd$pos - loci$pos[i]) <= NEAR_KB * 1000), logical(1))
    data.table(method = m, locus_key = loci$key,
               detected = loci$key %in% ck, detected_near = near)
}))

prof <- merge(det, loci, by.x = "locus_key", by.y = "key")
# Combine pseudo-methods have no single operating point; label them by the rule instead
# of leaving an NA that reads as "missing" rather than "not applicable".
prof[, point := ifelse(method %in% names(POINTS), POINTS[method],
                       paste0("combine@", CW_KB, "kb"))]
fwrite(prof, file.path(OUTDIR, paste0(TAG, "_locus_detection.tsv")), sep = "\t")
message("INFO: wrote ", TAG, "_locus_detection.tsv (", nrow(prof), " locus x method rows)")

# ------------------------------------------------------------- MAF stratification
# Bin edges are the MAF filters a user would plausibly set, not quantiles: the whole
# point is to read "what would Filter.maf = 0.05 have cost me" straight off the table.
BREAKS <- c(0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.5)
prof[, maf_bin := cut(maf, BREAKS, include.lowest = TRUE)]
maf_recall <- prof[, .(n_causal = .N, n_detected = sum(detected),
                       recall = mean(detected), recall_near = mean(detected_near)),
                   by = .(method, maf_bin)][order(method, maf_bin)]
fwrite(maf_recall, file.path(OUTDIR, paste0(TAG, "_maf_recall.tsv")), sep = "\t")

# The single number the MAF-filter question reduces to.
maf_cut <- prof[, .(n_causal = .N,
                    n_below_0.05 = sum(maf < 0.05),
                    detected_total = sum(detected),
                    detected_below_0.05 = sum(detected & maf < 0.05),
                    recall_below_0.05 = mean(detected[maf < 0.05]),
                    recall_above_0.05 = mean(detected[maf >= 0.05])), by = method]
fwrite(maf_cut, file.path(OUTDIR, paste0(TAG, "_maf_cut.tsv")), sep = "\t")
message("\n=== recall below vs above MAF 0.05 ===")
print(maf_cut)

# ------------------------------------------------------------ covariate models
# Univariate first (does each covariate matter on its own), then the joint fit (does MAF
# still matter once effect size is in the model). Reported as coefficients with SEs, not
# as p-values alone: with 70 loci the point estimate and its width are the honest output.
model_covars <- intersect(c("maf", "mutTempEffect", "mutSalEffect", "Va_temp_prop",
                            "Va_sal_prop", "af_cor_temp", "af_cor_sal", "Fst_outflank",
                            "dist_nearest_causal", "n_linked_5kb"), names(prof))
fit_rows <- list()
for (m in unique(prof$method)) {
    d <- prof[method == m]
    if (uniqueN(d$detected) < 2L) {
        message("INFO: ", m, " detected ", sum(d$detected), "/", nrow(d),
                " causal loci -- no variation, model skipped")
        next
    }
    for (cv in model_covars) {
        x <- suppressWarnings(as.numeric(d[[cv]]))
        if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) next
        fit <- try(glm(d$detected ~ x, family = binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) next
        co <- summary(fit)$coefficients
        if (nrow(co) < 2) next
        fit_rows[[length(fit_rows) + 1L]] <- data.table(
            method = m, model = "univariate", term = cv,
            estimate = co[2, 1], se = co[2, 2], z = co[2, 3], p = co[2, 4],
            n = sum(!is.na(x)), n_detected = sum(d$detected))
    }
    joint <- intersect(c("maf", "Va_temp_prop", "Va_sal_prop", "Fst_outflank"), model_covars)
    dj <- d[, c("detected", joint), with = FALSE]
    dj <- dj[complete.cases(dj)]
    if (nrow(dj) > 3 * length(joint) && uniqueN(dj$detected) == 2L) {
        fit <- try(glm(detected ~ ., data = dj, family = binomial()), silent = TRUE)
        if (!inherits(fit, "try-error")) {
            co <- summary(fit)$coefficients
            for (r in setdiff(rownames(co), "(Intercept)"))
                fit_rows[[length(fit_rows) + 1L]] <- data.table(
                    method = m, model = "joint", term = r,
                    estimate = co[r, 1], se = co[r, 2], z = co[r, 3], p = co[r, 4],
                    n = nrow(dj), n_detected = sum(dj$detected))
        }
    }
}
if (length(fit_rows)) {
    fits <- rbindlist(fit_rows)
    fwrite(fits, file.path(OUTDIR, paste0(TAG, "_detection_glm.tsv")), sep = "\t")
    message("\n=== univariate logistic slopes, |z| >= 2 ===")
    print(fits[model == "univariate" & abs(z) >= 2][order(-abs(z)),
        .(method, term, estimate = round(estimate, 3), z = round(z, 2), p = signif(p, 3))])
} else message("WARN: no logistic models could be fitted")

# ------------------------------------------------------------------- figures
p1 <- ggplot(maf_recall[!is.na(maf_bin)], aes(maf_bin, recall, group = method, colour = method)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    scale_colour_n(uniqueN(maf_recall$method)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste0("Causal-locus recall by MAF (", TAG, ")"),
         x = "MAF bin", y = "Recall (detected / causal in bin)", colour = NULL) +
    theme_clinego()
ggsave(file.path(OUTDIR, paste0(TAG, "_maf_recall.png")), p1, width = 9, height = 5.5, dpi = 150)
ggsave(file.path(OUTDIR, paste0(TAG, "_maf_recall.svg")), p1, width = 9, height = 5.5)

# |effect| across both axes. pmax(na.rm=TRUE) yields -Inf when BOTH are NA, so it is
# normalised back to NA rather than silently sorting into the bottom quartile.
eff <- copy(prof)
eff[, effect := pmax(abs(suppressWarnings(as.numeric(mutTempEffect))),
                     abs(suppressWarnings(as.numeric(mutSalEffect))), na.rm = TRUE)]
eff[!is.finite(effect), effect := NA_real_]
if (any(is.finite(eff$effect))) {
    eff[, eff_bin := cut(effect, quantile(effect, probs = seq(0, 1, 0.25), na.rm = TRUE),
                         include.lowest = TRUE)]
    er <- eff[!is.na(eff_bin), .(recall = mean(detected), n = .N), by = .(method, eff_bin)]
    fwrite(er, file.path(OUTDIR, paste0(TAG, "_effect_recall.tsv")), sep = "\t")
    p2 <- ggplot(er, aes(eff_bin, recall, group = method, colour = method)) +
        geom_line(linewidth = 0.7) + geom_point(size = 2) +
        scale_colour_n(uniqueN(er$method)) + scale_y_continuous(limits = c(0, 1)) +
        labs(title = paste0("Causal-locus recall by effect size (", TAG, ")"),
             x = "|effect| quartile", y = "Recall", colour = NULL) +
        theme_clinego()
    ggsave(file.path(OUTDIR, paste0(TAG, "_effect_recall.png")), p2, width = 9, height = 5.5, dpi = 150)
    ggsave(file.path(OUTDIR, paste0(TAG, "_effect_recall.svg")), p2, width = 9, height = 5.5)
}
message("INFO: profile outputs in ", OUTDIR)
