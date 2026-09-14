# combine_sigsnps.R — combine per-method significant SNPs into a single table
# Usage: source("/pipeline/scripts/R/lib/combine_sigsnps.R")
# Requires: data.table, dplyr, stringr
#
# data.table port of combine_selected_snps.R.
# Cross-method overlap detection uses foverlaps (instead of GenomicRanges::findOverlaps).

# Normalise legacy strategy names to canonical names
.normalise_strategy <- function(method) {
    switch(method,
        "All"           = "Union",
        "Sum"           = "Union",
        "Overlap"       = "Cross-method",
        "MethodOverlap" = "Cross-method per-trait",
        "PairOverlap"   = "Cross-method per-trait",
        method  # pass through canonical names and single-method names
    )
}

# Combine per-method significant-SNPs tables into a single selected_snps table.
#
# @param sig_snps_list  named list of data.tables (name = method name), each with
#                       columns: SNPID, chr, pos, trait, method, pvalue.
#                       Empty-row tables are allowed.
# @param strategy       character: "Union" | "Cross-method" | "Cross-method per-trait"
#                       or a single method name for passthrough.
#                       Legacy aliases: "All"/"Sum" → Union, "Overlap" → Cross-method,
#                       "MethodOverlap"/"PairOverlap" → Cross-method per-trait.
# @param clumping_distance  integer: maxgap for overlap detection (bp)
# @param predictors     character vector: traits to include (filter step)
# @return data.table with columns: SNPID, chr, pos, <method_cols>, min_pvalue
#   Each method column holds comma-sep trait names for SNPs significant in that method.
#   min_pvalue is the minimum over EVERY (method, trait) test of the SNP — the right
#   statistic for the trait-agnostic combined regions, and the wrong one for a per-trait
#   region (audit 2026-09-13 SC3): use combine_sigsnps_with_traits() for those.
#   Returns empty table if no SNPs found.
combine_sigsnps <- function(sig_snps_list, strategy, clumping_distance, predictors) {
    .combine_core(sig_snps_list, strategy, clumping_distance, predictors)$snps
}

# Same selection, plus the per-(SNP, trait) minimum p-value the wide table collapses away.
#
# @return list(snps = <the combine_sigsnps() table>,
#              trait_pvalues = data.table SNPID, chr, pos, trait, min_pvalue — one row per
#              selected SNP x trait it is significant for, min over the METHODS that tested
#              that trait). build_per_trait_regions() takes it so a bio_3 region reports
#              bio_3's own best p, not a co-located bio_2 hit's.
combine_sigsnps_with_traits <- function(sig_snps_list, strategy, clumping_distance, predictors) {
    .combine_core(sig_snps_list, strategy, clumping_distance, predictors)
}

.empty_trait_pvalues_dt <- function() {
    data.table::data.table(SNPID = character(), chr = character(), pos = integer(),
                           trait = character(), min_pvalue = numeric())
}

.combine_core <- function(sig_snps_list, strategy, clumping_distance, predictors) {
    strategy        <- .normalise_strategy(strategy)
    clumping_distance <- as.integer(clumping_distance)
    methods_vec     <- names(sig_snps_list)

    message(paste0("INFO: Combining SNPs. Strategy: ", strategy,
                   "  Clumping distance: ", clumping_distance))

    # Filter by predictors (climate traits)
    sig_filtered <- lapply(sig_snps_list, function(dt) {
        if (is.null(dt) || nrow(dt) == 0) {
            return(data.table::data.table(
                SNPID = character(), chr = character(), pos = integer(),
                trait = character(), method = character(), pvalue = numeric()
            ))
        }
        dt[, chr := as.character(chr)]
        dt[trait %in% predictors, .(SNPID, chr, pos, trait, method, pvalue)]
    }) %>% setNames(methods_vec)

    total_snps <- sum(vapply(sig_filtered, nrow, integer(1L)))
    if (total_snps == 0L) {
        message("WARNING: No significant SNPs in any method")
        return(list(snps = .empty_combined_snps_dt(methods_vec),
                    trait_pvalues = .empty_trait_pvalues_dt()))
    }

    # Select SNP set according to strategy
    if (strategy %in% methods_vec) {
        # Single-method passthrough
        snp_rows <- sig_filtered[[strategy]]

    } else if (strategy == "Union") {
        snp_rows <- data.table::rbindlist(sig_filtered)

    } else if (strategy == "Cross-method") {
        snp_rows <- .overlap_cross_method(sig_filtered, methods_vec, clumping_distance,
                                          per_trait = FALSE)

    } else if (strategy == "Cross-method per-trait") {
        snp_rows <- .overlap_cross_method(sig_filtered, methods_vec, clumping_distance,
                                          per_trait = TRUE)

    } else {
        stop(paste0("Unknown strategy: '", strategy,
                    "'. Use: Union, Cross-method, Cross-method per-trait, or a single method name."))
    }

    if (is.null(snp_rows) || nrow(snp_rows) == 0L) {
        message("WARNING: No SNPs remain after applying strategy")
        return(list(snps = .empty_combined_snps_dt(methods_vec),
                    trait_pvalues = .empty_trait_pvalues_dt()))
    }

    # Per-(SNP, trait) minimum: snp_rows is long (one row per method x trait test that
    # selected the SNP), so this is the trait-scoped statistic before the collapse below.
    trait_pvalues <- snp_rows[, .(min_pvalue = min(pvalue, na.rm = TRUE)),
                              by = c("SNPID", "chr", "pos", "trait")]
    data.table::setorder(trait_pvalues, chr, pos, trait)

    # Build output: one row per SNPID, method columns hold comma-sep traits
    selected_ids <- unique(snp_rows$SNPID)

    method_traits <- lapply(methods_vec, function(m) {
        dt <- sig_filtered[[m]]
        if (nrow(dt) == 0L) {
            return(data.table::data.table(SNPID = character(), traits = character()))
        }
        dt[SNPID %in% selected_ids,
           .(traits = paste(sort(unique(trait)), collapse = ",")),
           by = "SNPID"]
    }) %>% setNames(methods_vec)

    result <- snp_rows[, .(
        min_pvalue = min(pvalue, na.rm = TRUE)
    ), by = c("SNPID", "chr", "pos")]
    data.table::setorder(result, chr, pos)

    for (m in methods_vec) {
        mt <- method_traits[[m]]
        if (nrow(mt) > 0L) {
            result <- merge(result, mt[, .(SNPID, traits)],
                            by = "SNPID", all.x = TRUE)
            data.table::setnames(result, "traits", m)
        } else {
            result[[m]] <- NA_character_
        }
        result[[m]][is.na(result[[m]])] <- ""
    }

    col_order <- c("SNPID", "chr", "pos", methods_vec, "min_pvalue")
    result <- result[, ..col_order]

    message(paste0("INFO: ", nrow(result), " unique SNPs selected"))
    for (m in methods_vec) {
        n_m <- sum(result[[m]] != "", na.rm = TRUE)
        message(paste0("INFO: ", m, ": ", n_m, " SNPs with significant traits"))
    }

    list(snps = result, trait_pvalues = trait_pvalues)
}

# --- internal helpers ---

# Find SNPs significant in ≥2 methods with overlap within clumping_distance.
# per_trait=FALSE: any trait combo; per_trait=TRUE: must share the same trait.
.overlap_cross_method <- function(sig_filtered, methods_vec, clumping_distance,
                                  per_trait = FALSE) {
    results <- list()

    # UNORDERED pairs only. .snp_overlap() returns the matched rows from BOTH
    # tables and its match condition (|pos1 - pos2| <= clumping_distance) is
    # symmetric, so (n2, n1) reproduces exactly the rows (n1, n2) already
    # produced — the old ordered double loop ran every pairwise foverlaps twice
    # and threw the second copy away in unique(rbindlist(results)).
    for (i in seq_along(methods_vec)) {
        for (j in seq_len(i - 1L)) {
            n1 <- methods_vec[i]
            n2 <- methods_vec[j]
            dt1 <- sig_filtered[[n1]]
            dt2 <- sig_filtered[[n2]]
            if (nrow(dt1) == 0L || nrow(dt2) == 0L) next

            if (per_trait) {
                shared_traits <- intersect(dt1$trait, dt2$trait)
                for (bio in shared_traits) {
                    g1 <- dt1[trait == bio]
                    g2 <- dt2[trait == bio]
                    if (nrow(g1) == 0L || nrow(g2) == 0L) next
                    ov <- .snp_overlap(g1, g2, clumping_distance)
                    if (!is.null(ov)) results <- c(results, list(ov))
                }
            } else {
                ov <- .snp_overlap(dt1, dt2, clumping_distance)
                if (!is.null(ov)) results <- c(results, list(ov))
            }
        }
    }

    if (length(results) == 0L) return(NULL)
    unique(data.table::rbindlist(results))
}

# foverlaps-based overlap: find rows in dt2 within clumping_distance of any row in dt1.
# Returns rbind of matched rows from both tables.
.snp_overlap <- function(dt1, dt2, clumping_distance) {
    # Extend dt1 by clumping_distance on both sides
    q1 <- data.table::copy(dt1)
    q1[, s := pmax(1L, as.integer(pos) - clumping_distance)]
    q1[, e := as.integer(pos) + clumping_distance]

    s2 <- data.table::copy(dt2)
    s2[, s := as.integer(pos)]
    s2[, e := as.integer(pos)]

    data.table::setkey(s2, chr, s, e)
    data.table::setkey(q1, chr, s, e)

    ov <- data.table::foverlaps(q1, s2,
        by.x = c("chr", "s", "e"),
        by.y = c("chr", "s", "e"),
        type = "any", mult = "all", nomatch = NA)

    ov <- ov[!is.na(SNPID)]
    if (nrow(ov) == 0L) return(NULL)

    # Recover original dt1 rows that had matches
    q1_hits <- unique(ov[, .(SNPID = i.SNPID, chr = chr, pos = i.pos,
                              trait = i.trait, method = i.method, pvalue = i.pvalue)])
    dt2_hits <- unique(ov[, .(SNPID,  chr, pos = s,   # s = original pos
                               trait, method, pvalue)])

    rbind(q1_hits, dt2_hits)
}

.empty_combined_snps_dt <- function(methods_vec) {
    dt <- data.table::data.table(
        SNPID = character(),
        chr   = character(),
        pos   = integer()
    )
    for (m in methods_vec) dt[[m]] <- character()
    dt[, min_pvalue := numeric()]
    dt
}
