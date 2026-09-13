# sig_snps.R — per-trait significance filtering and cross-trait overlap annotation
# Usage: source("/pipeline/scripts/R/lib/sig_snps.R")
# Requires: data.table, dplyr, qvalue, parallel, pval_threshold.R, io_selected_snps.R

# Find significant SNPs for each trait in a p-value table.
#
# @param pvalues_dt  data.table with columns SNPID, chr, pos, [trait_cols...]
# @param adjustment  character: "bonf" | "qval" | "top"
# @param value       numeric: threshold value for the adjustment
# @param cpu         integer: cores for mclapply
# @param exclude_cols character vector of non-trait columns to skip
# @param is_wza      logical: WZA mode (some windows have NA p per trait — already handled
#                    by compute_pval_threshold but logged differently)
# @return named list:
#   sig_snps         data.table (SNPID, chr, pos, pvalue, pval_threshold, trait)
#                    or empty data.table if nothing significant
#   diagnostics      data.table (trait, status, threshold, n_tested, n_na_dropped)
find_significant_snps_per_trait <- function(pvalues_dt, adjustment, value, cpu,
                                            exclude_cols = character(0),
                                            is_wza = FALSE) {
    pvalues_dt <- data.table::copy(pvalues_dt)
    pvalues_dt$chr <- as.character(pvalues_dt$chr)

    traits <- setdiff(
        colnames(pvalues_dt),
        c("SNPID", "chr", "pos", exclude_cols)
    )

    message(paste0("INFO: Traits: ", paste(traits, collapse = ", ")))
    if (is_wza) message("INFO: WZA mode — bonf denominator = non-NA windows per trait")

    diag_list <- vector("list", length(traits))
    names(diag_list) <- traits

    snps_per_trait <- parallel::mclapply(traits, function(trait) {
        result <- compute_pval_threshold(pvalues_dt[[trait]], adjustment, value)

        message(paste0("INFO: Trait=", trait,
                       " status=",       result$status,
                       " threshold=",    result$threshold,
                       " n_tested=",     result$n_tested,
                       " n_na_dropped=", result$n_na_dropped))

        diag_list[[trait]] <<- data.table::data.table(
            trait        = trait,
            status       = result$status,
            threshold    = result$threshold,
            n_tested     = result$n_tested,
            n_na_dropped = result$n_na_dropped
        )

        if (result$status != "ok") {
            message(paste0("WARNING: No threshold for trait ", trait,
                           " (status=", result$status, "). Writing empty result."))
            return(NULL)
        }

        # INCLUSIVE (<=). compute_pval_threshold() returns a cutoff that is
        # itself a member of the intended call set for two of the four rules:
        #   top  N -> max_pvalue_top() = the N-th smallest p, so `<` returns N-1
        #   qval f -> max_pvalue_fdr() = the LARGEST p whose q < f, so `<` drops
        #             the boundary SNP and every SNP tied with it
        # bonf/custom are unaffected in practice (exact equality on continuous p)
        # but `<=` is also the convention there, and it matches the Manhattan's
        # `log10p >= threshold_log10` (plot_manhattan.R, io_pvalues.R) so the
        # table and the plot agree on the boundary SNP.
        mask <- pvalues_dt[[trait]] <= result$threshold
        mask[is.na(mask)] <- FALSE

        pvalues_dt[mask, list(SNPID, chr, pos,
                              pvalue          = .SD[[trait]],
                              pval_threshold  = result$threshold,
                              trait           = trait),
                   .SDcols = trait]
    }, mc.cores = as.integer(cpu))

    diagnostics <- data.table::rbindlist(Filter(Negate(is.null), diag_list))

    sig_snps <- data.table::rbindlist(Filter(Negate(is.null), snps_per_trait))

    if (is.null(sig_snps) || nrow(sig_snps) == 0) {
        message("WARNING: No significant SNPs found")
        sig_snps <- data.table::data.table(
            SNPID          = character(),
            chr            = character(),
            pos            = numeric(),
            pvalue         = numeric(),
            pval_threshold = numeric(),
            trait          = character()
        )
    } else {
        message(paste0("INFO: Found ", nrow(sig_snps), " significant SNPs across ", length(traits), " traits"))
    }

    list(sig_snps = sig_snps, diagnostics = diagnostics)
}

# Annotate significant SNPs with cross-trait overlap information.
#
# For each sig SNP, finds other-trait SNPs within `distance` bp and records
# which traits they belong to.
#
# @param sig_snps  data.table (SNPID, chr, pos, pvalue, pval_threshold, trait)
# @param distance  integer: search radius in bp (= snp_clumping_distance)
# @return sig_snps with three extra columns:
#   overlap_traits   comma-sep trait names within distance (NA if none)
#   overlap_snps     comma-sep "chr:pos" strings (NA if none)
#   overlap_distance integer: the distance used
annotate_cross_trait_overlaps <- function(sig_snps, distance) {
    distance <- as.integer(distance)

    empty_cols <- function(n) {
        list(
            overlap_traits   = rep(NA_character_, n),
            overlap_snps     = rep(NA_character_, n),
            overlap_distance = rep(distance, n)
        )
    }

    if (is.null(sig_snps) || nrow(sig_snps) == 0) {
        sig_snps[, overlap_traits   := NA_character_]
        sig_snps[, overlap_snps     := NA_character_]
        sig_snps[, overlap_distance := distance]
        return(sig_snps)
    }

    sig_snps <- data.table::copy(sig_snps)
    sig_snps[, chr := as.character(chr)]
    sig_snps[, pos := as.integer(pos)]

    n <- nrow(sig_snps)
    ot <- character(n)
    os <- character(n)

    # Build a lookup table: extended windows for all SNPs
    # query: each sig SNP extended by distance on both sides
    # subject: all sig SNPs as point ranges
    query <- data.table::copy(sig_snps[, .(row_i = .I, SNPID, chr,
                                            s = pmax(1L, pos - distance),
                                            e = pos + distance,
                                            trait)])
    subj  <- data.table::copy(sig_snps[, .(SNPID, chr, s = pos, e = pos, trait)])

    data.table::setkey(subj, chr, s, e)
    data.table::setkey(query, chr, s, e)

    ov <- data.table::foverlaps(query, subj,
        by.x = c("chr", "s", "e"),
        by.y = c("chr", "s", "e"),
        type = "any", mult = "all", nomatch = NA)

    # After foverlaps(query, subj): plain cols = subj (the NEIGHBOUR), i.* = query
    # (the source row). nomatch = NA leaves the SUBJECT side NA, so the no-match
    # guard must test the plain column.
    # Filter: only other-trait matches (same SNPID but different trait is NOT other-trait)
    ov <- ov[!is.na(SNPID) & trait != i.trait]

    if (nrow(ov) > 0) {
        # Aggregate per source row
        ov_agg <- ov[, .(
            overlap_traits = paste(sort(unique(trait)),  collapse = ","),
            overlap_snps   = paste(sort(unique(paste0(chr, ":", s))), collapse = ",")
        ), by = "row_i"]

        for (idx in seq_len(nrow(ov_agg))) {
            ri <- ov_agg$row_i[idx]
            ot[ri] <- ov_agg$overlap_traits[idx]
            os[ri] <- ov_agg$overlap_snps[idx]
        }
    }

    ot[ot == ""] <- NA_character_
    os[os == ""] <- NA_character_

    sig_snps[, overlap_traits   := ot]
    sig_snps[, overlap_snps     := os]
    sig_snps[, overlap_distance := distance]

    sig_snps
}
