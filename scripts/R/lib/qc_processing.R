# qc_processing.R — the numbers behind Processing's QC tables, each taken from the step that
# made the decision rather than re-derived after the fact (release Phase 4c).
# Usage: source("/pipeline/scripts/R/lib/qc_processing.R")
# Requires: data.table
#
# Caller: plot_qc_processing.R
#   filtering_summary.tsv     filtering_attrition(): sample rows count the lists the sample
#                             filters wrote, SNP rows are plink's own report of filter_vcf's
#                             --geno and --maf steps (parse_plink_filter_log()), 7ccf66 + 7e11b6
#   het_vs_missingness.png    het_outlier_bounds(): the het_outlier_filter rule's formula, 7ccf66
#   sample_qc_by_site.tsv     site_qc_table(): missing calls and depth per sampling site, 787ffe

# IIDs of a plink keep/remove list: "FID IID" per line, space- or tab-separated (the sample
# filters write both). A line without a second field is an error, not a sample.
read_sample_list <- function(path) {
    lines <- trimws(readLines(path, warn = FALSE))
    lines <- lines[nzchar(lines)]
    fields <- strsplit(lines, "[ \t]+")
    bad <- lengths(fields) < 2L
    if (any(bad))
        stop(path, ": line(s) without an IID: ", paste(head(lines[bad], 3L), collapse = " | "))
    vapply(fields, `[`, character(1), 2L)
}

# plink 1.9's own report of a filtering run (its .log), as numbers.
#
# plink applies --geno before --maf (both after --keep), so the MAF count is over the variants
# that passed --geno. Its lines, as plink 1.90b7.2 writes them:
#   354 variants loaded from .bim file.
#   --keep: 47 samples remaining.
#   300 variants removed due to missing genotype data (--geno).
#   1 variant removed due to minor allele threshold(s)
#   53 variants and 47 samples pass filters and QC.
# At --geno 1 plink drops the missingness filter and prints no --geno line, so a missing
# --geno line reads as 0 removed only when `geno` >= 1. A missing --maf line is an error:
# filter_vcf always passes --maf, and plink refuses --maf 0.
#
# @param lines character, the log's lines
# @param geno  the --geno threshold that run used (Filter.snp_miss)
# @return list: n_loaded, n_removed_geno, n_removed_maf, n_pass, n_samples_pass
parse_plink_filter_log <- function(lines, geno) {
    lines <- trimws(lines)
    grab <- function(pattern, what, n_groups = 1L) {
        hit <- regmatches(lines, regexec(pattern, lines, perl = TRUE))
        hit <- hit[lengths(hit) > 0L]
        if (length(hit) == 0L) return(NULL)
        if (length(hit) > 1L) stop("plink log: more than one line reports ", what)
        as.numeric(hit[[1L]][seq_len(n_groups) + 1L])
    }
    loaded <- grab("^(\\d+) variants? loaded from \\.bim file\\.$", "variants loaded")
    geno_rm <- grab("^(\\d+) variants? removed due to missing genotype data \\(--geno\\)\\.$",
                    "the --geno removals")
    maf_rm <- grab("^(\\d+) variants? removed due to minor allele threshold\\(s\\)",
                   "the --maf removals")
    pass <- grab("^(\\d+) variants? and (\\d+) (?:samples?|people|person) pass filters and QC\\.$",
                 "variants passing", 2L)
    if (is.null(loaded) || is.null(pass))
        stop("plink log: no 'variants loaded' or 'pass filters and QC' line; is this plink 1.9's log of a filtering run?")
    if (is.null(geno_rm)) {
        if (as.numeric(geno) < 1)
            stop("plink log: no --geno line although --geno ", geno, " filters (plink omits it only at 1)")
        geno_rm <- 0
    }
    if (is.null(maf_rm)) stop("plink log: no --maf line; filter_vcf always passes --maf")
    if (loaded - geno_rm - maf_rm != pass[1L])
        stop(sprintf("plink log does not add up: %s loaded - %s (--geno) - %s (--maf) != %s passing",
                     loaded, geno_rm, maf_rm, pass[1L]))
    list(n_loaded = loaded, n_removed_geno = geno_rm, n_removed_maf = maf_rm,
         n_pass = pass[1L], n_samples_pass = pass[2L])
}

# filtering_summary.tsv. Every row is a step the pipeline ran, in the order it ran:
#   Raw VCF                        the keep-list's samples, bcftools' SNP count
#   sample_steps                   named counts, one per sample filter that ran, in order
#                                  (each the length of the list that filter wrote)
#   After SNP missingness filter   plink's --geno step on the final samples
#   After MAF + SNP missingness    plink's --maf step after it: filter_vcf's output
#   After LD pruning               ld_prune's output
# A sample filter removes samples only, so its rows keep the raw SNP count.
#
# @param plink  parse_plink_filter_log() of the filter_vcf_steps replay
# @return data.table stage, n_samples, n_snps, pct_samples, pct_snps
filtering_attrition <- function(n_raw_samples, n_raw_snps, sample_steps, plink, n_ld_snps) {
    n_final <- plink$n_samples_pass
    out <- data.table::data.table(
        stage = c("Raw VCF", names(sample_steps),
                  "After SNP missingness filter", "After MAF + SNP missingness filter",
                  "After LD pruning"),
        n_samples = c(n_raw_samples, unlist(sample_steps, use.names = FALSE),
                      n_final, n_final, n_final),
        n_snps = c(n_raw_snps, rep(n_raw_snps, length(sample_steps)),
                   plink$n_loaded - plink$n_removed_geno, plink$n_pass, n_ld_snps)
    )
    out[, pct_samples := round(100 * n_samples / n_raw_samples, 1)]
    out[, pct_snps    := round(100 * n_snps / n_raw_snps, 1)]
    out[, stage := factor(stage, levels = stage)]
    out[]
}

# The het_outlier_filter rule's bounds (processing.smk): observed heterozygosity
# (N_NM - O_HOM) / N_NM, mean +/- k POPULATION standard deviations (divided by n, as the
# rule's awk does), over the samples that rule read.
het_outlier_bounds <- function(obs_het, k) {
    m <- mean(obs_het)
    s <- sqrt(sum((obs_het - m)^2) / length(obs_het))
    c(mean = m, sd = s, lo = m - k * s, hi = m + k * s)
}

# Missing calls and depth per sampling site, over every sample of the keep-list (the samples
# the missingness filter judged), so a site whose samples failed it is still listed.
#
# @param meta    data.table with site, sample: the input metadata (Input.metadata)
# @param smiss   data.table with IID, F_MISS: calculate_sample_missing's per-sample rates
#                on the raw VCF
# @param passed_missingness IIDs the sample-missingness filter kept
# @param retained           IIDs of the final cohort (after every sample filter)
# @param depth   NULL, or data.table with INDV, MEAN_DEPTH (vcftools --depth on the raw VCF);
#                a sample's non-finite or missing depth is left out of its site's depth
# @return data.table site, n_samples, n_failed_missingness, n_retained, f_miss_mean,
#   f_miss_median, f_miss_min, f_miss_max, depth_mean, depth_median; sites in the
#   metadata's order
site_qc_table <- function(meta, smiss, passed_missingness, retained, depth = NULL) {
    site_of <- stats::setNames(as.character(meta$site), as.character(meta$sample))
    ids <- as.character(smiss$IID)
    unknown <- setdiff(ids, names(site_of))
    if (length(unknown) > 0L)
        stop(length(unknown), " sample(s) of the missingness table have no site in the metadata: ",
             paste(head(unknown, 5L), collapse = ", "))
    dt <- data.table::data.table(site = unname(site_of[ids]), sample = ids,
                                 f_miss = as.numeric(smiss$F_MISS))
    dt[, depth := NA_real_]
    if (!is.null(depth)) {
        d <- stats::setNames(as.numeric(depth$MEAN_DEPTH), as.character(depth$INDV))
        dt[, depth := unname(d[sample])]
    }
    dt[!is.finite(depth), depth := NA_real_]
    na_or <- function(x, f) if (all(is.na(x))) NA_real_ else f(x[!is.na(x)])
    out <- dt[, .(
        n_samples            = .N,
        n_failed_missingness = sum(!sample %in% passed_missingness),
        n_retained           = sum(sample %in% retained),
        f_miss_mean          = round(mean(f_miss), 4),
        f_miss_median        = round(stats::median(f_miss), 4),
        f_miss_min           = round(min(f_miss), 4),
        f_miss_max           = round(max(f_miss), 4),
        depth_mean           = round(na_or(depth, mean), 2),
        depth_median         = round(na_or(depth, stats::median), 2)
    ), by = site]
    out[order(match(site, unique(as.character(meta$site))))]
}
