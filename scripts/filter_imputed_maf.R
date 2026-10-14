# Filter the sNMF-imputed genotype set by minor-allele frequency (sNMF.imputed_maf).
#
# Filter.maf runs as plink --maf on OBSERVED calls (filter_vcf). LEA impute(method = "mode")
# then fills every missing call, mostly with the major allele, so a SNP whose minor allele is
# carried by one small site can pass the observed filter and be tested far below it. This
# script recomputes MAF on the imputed matrix -- the one every GEA/GWAS method tests -- and
# keeps SNPs with MAF >= THRESHOLD (same >= rule as plink --maf; THRESHOLD 0 keeps all).
#
# MAF is computed over the full post-filter cohort the imputation ran on. Climate-valid or
# per-trait sample subsets taken downstream can sit slightly below the threshold.
#
# Writes ONE consistent set, so column i of OUT_LFMM is row i of OUT_VCFSNP and body row i of
# OUT_VCF:
#   OUT_LFMM    kept columns of the imputed .lfmm (samples x SNPs, space-separated, 0/1/2)
#   OUT_VCFSNP  kept lines of LEA's .vcfsnp, byte-identical
#   OUT_VCF     kept body rows of the imputed VCF; header and REF/ALT copied verbatim
#   OUT_REPORT  one row: threshold, n_samples, n_snps_before, n_snps_after, n_snps_removed
#
# Usage: Rscript filter_imputed_maf.R LFMM_IMP VCFSNP VCF_IMP THRESHOLD \
#            OUT_LFMM OUT_VCFSNP OUT_VCF OUT_REPORT
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    stop("filter_imputed_maf.R: expected 8 arguments (LFMM_IMP VCFSNP VCF_IMP THRESHOLD ",
         "OUT_LFMM OUT_VCFSNP OUT_VCF OUT_REPORT), got ", length(args))
}
LFMM_IMP   <- args[1]
VCFSNP     <- args[2]
VCF_IMP    <- args[3]
THRESHOLD  <- as.numeric(args[4])
OUT_LFMM   <- args[5]
OUT_VCFSNP <- args[6]
OUT_VCF    <- args[7]
OUT_REPORT <- args[8]

if (is.na(THRESHOLD) || THRESHOLD < 0 || THRESHOLD >= 0.5) {
    stop("filter_imputed_maf.R: THRESHOLD must be a number in [0, 0.5), got '", args[4], "'")
}

# ---- pass 1: allele counts per SNP, one sample row at a time ---------------------------------
con <- file(LFMM_IMP, "r")
ac <- NULL
n_ind <- 0L
repeat {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    g <- as.integer(strsplit(trimws(line), "[ \t]+")[[1]])
    if (is.null(ac)) ac <- numeric(length(g))
    if (length(g) != length(ac)) {
        stop(sprintf("filter_imputed_maf.R: row %d of %s has %d SNPs, row 1 has %d",
                     n_ind + 1L, LFMM_IMP, length(g), length(ac)))
    }
    if (anyNA(g) || any(g == 9L)) {
        stop("filter_imputed_maf.R: ", LFMM_IMP, " still has missing calls (9) at row ", n_ind + 1L,
             " -- not an imputed matrix")
    }
    if (any(g < 0L | g > 2L)) {
        stop("filter_imputed_maf.R: ", LFMM_IMP, " has a dosage outside 0/1/2 at row ", n_ind + 1L,
             " (diploid only)")
    }
    ac <- ac + g
    n_ind <- n_ind + 1L
}
close(con)
if (n_ind == 0L) stop("filter_imputed_maf.R: ", LFMM_IMP, " is empty")
n_snp <- length(ac)

# Compare allele COUNTS, not frequencies: 1 - 0.9 is 0.0999..., which would drop a SNP that sits
# exactly on the threshold. The epsilon absorbs THRESHOLD * n2 rounding the same way.
n2 <- 2 * n_ind
mac <- pmin(ac, n2 - ac)
keep <- mac >= THRESHOLD * n2 - 1e-9
n_keep <- sum(keep)
message(sprintf("INFO: %d samples, %d imputed SNPs; MAF >= %s keeps %d, removes %d (%.2f%%)",
                n_ind, n_snp, format(THRESHOLD), n_keep, n_snp - n_keep,
                100 * (n_snp - n_keep) / n_snp))
message("INFO: MAF computed on the imputed matrix over the full post-filter cohort")
if (n_keep == 0L) {
    stop(sprintf(paste0("filter_imputed_maf.R: no SNP has imputed MAF >= %s (max %.4f). ",
                        "Lower sNMF.imputed_maf, or set it to 0 to disable the filter."),
                 format(THRESHOLD), max(mac) / n2))
}

# ---- vcfsnp: must label every lfmm column ------------------------------------------------------
snp_lines <- readLines(VCFSNP)
if (length(snp_lines) != n_snp) {
    stop(sprintf("filter_imputed_maf.R: %s has %d lines but %s has %d SNP columns",
                 VCFSNP, length(snp_lines), LFMM_IMP, n_snp))
}
snp_fields <- strsplit(trimws(snp_lines), "[ \t]+")
snp_chr <- vapply(snp_fields, `[`, "", 1)
snp_pos <- as.numeric(vapply(snp_fields, `[`, "", 2))

# ---- VCF: header verbatim, body rows checked against vcfsnp and filtered -----------------------
con_in <- file(VCF_IMP, "r")
con_out <- file(OUT_VCF, "w")
# A VCF error is found mid-stream, after OUT_VCF has been opened: remove the partial file so
# a failed run leaves no output behind (nothing else is written before this block finishes).
fail_vcf <- function(msg) {
    close(con_in)
    close(con_out)
    unlink(OUT_VCF)
    stop(msg, call. = FALSE)
}
row_i <- 0L
repeat {
    chunk <- readLines(con_in, n = 100000L)
    if (length(chunk) == 0) break
    is_head <- startsWith(chunk, "#")
    if (any(is_head)) {
        if (row_i > 0L || any(is_head[-seq_len(sum(cumprod(is_head)))])) {
            fail_vcf(paste0("filter_imputed_maf.R: header line after a body row in ", VCF_IMP))
        }
        writeLines(chunk[is_head], con_out)
    }
    body <- chunk[!is_head]
    if (length(body) == 0) next
    idx <- row_i + seq_along(body)
    if (max(idx) > n_snp) {
        fail_vcf(sprintf("filter_imputed_maf.R: %s has more body rows than the %d SNPs in %s",
                         VCF_IMP, n_snp, LFMM_IMP))
    }
    head2 <- strsplit(substr(body, 1, 200), "\t", fixed = TRUE)
    b_chr <- vapply(head2, `[`, "", 1)
    b_pos <- as.numeric(vapply(head2, `[`, "", 2))
    bad <- which(b_chr != snp_chr[idx] | b_pos != snp_pos[idx])
    if (length(bad)) {
        j <- idx[bad[1]]
        fail_vcf(sprintf("filter_imputed_maf.R: VCF row %d is %s:%s but vcfsnp row %d is %s:%s -- %s and %s are not row-aligned",
                         j, b_chr[bad[1]], format(b_pos[bad[1]], scientific = FALSE),
                         j, snp_chr[j], format(snp_pos[j], scientific = FALSE), VCF_IMP, VCFSNP))
    }
    writeLines(body[keep[idx]], con_out)
    row_i <- max(idx)
}
if (row_i != n_snp) {
    fail_vcf(sprintf("filter_imputed_maf.R: %s has %d body rows but %s has %d SNP columns",
                     VCF_IMP, row_i, LFMM_IMP, n_snp))
}
close(con_in)
close(con_out)

writeLines(snp_lines[keep], OUT_VCFSNP)

# ---- pass 2: write the kept lfmm columns -------------------------------------------------------
con_in <- file(LFMM_IMP, "r")
con_out <- file(OUT_LFMM, "w")
repeat {
    line <- readLines(con_in, n = 1)
    if (length(line) == 0) break
    g <- strsplit(trimws(line), "[ \t]+")[[1]]
    writeLines(paste(g[keep], collapse = " "), con_out)
}
close(con_in)
close(con_out)

write.table(data.frame(threshold = THRESHOLD, n_samples = n_ind, n_snps_before = n_snp,
                       n_snps_after = n_keep, n_snps_removed = n_snp - n_keep),
            OUT_REPORT, sep = "\t", quote = FALSE, row.names = FALSE)
message("INFO: wrote ", OUT_LFMM, ", ", OUT_VCFSNP, ", ", OUT_VCF, ", ", OUT_REPORT)
