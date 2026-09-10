# On-the-fly GEA × GWAS region overlap computation for the Overlapping tab.
# All computation is done in Shiny (data.table + foverlaps); no GenomicRanges needed.

# ── Core overlap computation ───────────────────────────────────────────────────

#' Compute pairwise overlaps between GEA and GWAS region sets
#'
#' Returns one row per overlapping (gea_region, gwas_region) pair.
#' Pairs are produced "all-pairs separately" — no transitive merging.
#'
#' @param gea_regions  data.table with columns: region_id, chr, start, end
#' @param gwas_regions data.table with columns: region_id, chr, start, end
#' @return data.table with columns:
#'   gea_region_id, gwas_region_id, chr,
#'   gea_start, gea_end, gwas_start, gwas_end,
#'   overlap_start, overlap_end, union_start, union_end, overlap_pct
#'   (empty data.table if no overlaps)
#' @noRd
compute_region_overlaps <- function(gea_regions, gwas_regions) {
    if (is.null(gea_regions) || nrow(gea_regions) == 0 ||
        is.null(gwas_regions) || nrow(gwas_regions) == 0) {
        return(.empty_overlap_pairs())
    }

    # Build foverlaps-compatible tables
    gea_ov <- data.table::data.table(
        region_id = as.character(gea_regions$region_id),
        chr       = as.character(gea_regions$chr),
        start     = as.integer(gea_regions$start),
        end       = as.integer(gea_regions$end)
    )
    gwas_ov <- data.table::data.table(
        region_id = as.character(gwas_regions$region_id),
        chr       = as.character(gwas_regions$chr),
        start     = as.integer(gwas_regions$start),
        end       = as.integer(gwas_regions$end)
    )

    data.table::setkey(gea_ov,  chr, start, end)
    data.table::setkey(gwas_ov, chr, start, end)

    hits <- data.table::foverlaps(gwas_ov, gea_ov,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type     = "any",
        mult     = "all",
        nomatch  = NULL
    )

    if (nrow(hits) == 0) return(.empty_overlap_pairs())

    # hits columns: region_id (GEA), chr, start (GEA), end (GEA),
    #               i.region_id (GWAS), i.start (GWAS), i.end (GWAS)
    result <- data.table::data.table(
        gea_region_id  = hits$region_id,
        gwas_region_id = hits$i.region_id,
        chr            = hits$chr,
        gea_start      = hits$start,
        gea_end        = hits$end,
        gwas_start     = hits$i.start,
        gwas_end       = hits$i.end
    )

    result[, overlap_start := pmax(gea_start, gwas_start)]
    result[, overlap_end   := pmin(gea_end,   gwas_end)]
    result[, union_start   := pmin(gea_start, gwas_start)]
    result[, union_end     := pmax(gea_end,   gwas_end)]

    # Overlap percentage relative to smaller region
    result[, min_len := pmin(gea_end - gea_start, gwas_end - gwas_start)]
    result[, overlap_pct := round(100 * (overlap_end - overlap_start) / pmax(min_len, 1L), 1)]
    result[, min_len := NULL]

    data.table::setorder(result, chr, gea_start)
    result
}

#' Empty overlap pairs skeleton
#' @noRd
.empty_overlap_pairs <- function() {
    data.table::data.table(
        gea_region_id  = character(),
        gwas_region_id = character(),
        chr            = character(),
        gea_start      = integer(),
        gea_end        = integer(),
        gwas_start     = integer(),
        gwas_end       = integer(),
        overlap_start  = integer(),
        overlap_end    = integer(),
        union_start    = integer(),
        union_end      = integer(),
        overlap_pct    = numeric()
    )
}

# ── Overlap region building ────────────────────────────────────────────────────

#' Apply an overlap bounds strategy to a single pair row, returning region bounds
#'
#' @param pair_row Single-row data.table from compute_region_overlaps()
#' @param strategy One of: "union", "intersection", "gea_only", "gwas_only"
#' @return Named list: chr, start, end
#' @noRd
.apply_bounds <- function(pair_row, strategy) {
    switch(strategy,
        union        = list(chr   = pair_row$chr,
                            start = pair_row$union_start,
                            end   = pair_row$union_end),
        intersection = list(chr   = pair_row$chr,
                            start = pair_row$overlap_start,
                            end   = pair_row$overlap_end),
        gea_only     = list(chr   = pair_row$chr,
                            start = pair_row$gea_start,
                            end   = pair_row$gea_end),
        gwas_only    = list(chr   = pair_row$chr,
                            start = pair_row$gwas_start,
                            end   = pair_row$gwas_end),
        # default: union
        list(chr   = pair_row$chr,
             start = pair_row$union_start,
             end   = pair_row$union_end)
    )
}

#' Compute all overlap regions for the region explorer
#'
#' Iterates over every GEA×GWAS overlapping pair and builds one overlap region
#' per pair using the specified bounds strategy.
#' Returns a data.table compatible with compute_all_regions() output, plus
#' extra columns gea_region_id and gwas_region_id for shape building.
#'
#' @param overlap_pairs  data.table from compute_region_overlaps()
#' @param gea_regions    data.table of GEA regions (from compute_all_regions)
#' @param gwas_regions   data.table of GWAS regions (from compute_all_regions)
#' @param strategy       Character: "union" | "intersection" | "gea_only" | "gwas_only"
#' @param gea_sig_snps   data.table of GEA sig SNPs (SNPID, chr, pos, pvalue, method, trait)
#' @param gwas_sig_snps  data.table of GWAS sig SNPs
#' @return data.table with: region_id, chr, start, end, length, snp_count,
#'   snp_ids, traits, methods, min_pvalue, gea_region_id, gwas_region_id
#' @noRd
compute_all_overlap_regions <- function(overlap_pairs, gea_regions, gwas_regions,
                                         strategy = "union",
                                         gea_sig_snps = NULL, gwas_sig_snps = NULL) {
    if (is.null(overlap_pairs) || nrow(overlap_pairs) == 0) {
        return(.empty_overlap_regions())
    }

    rows <- vector("list", nrow(overlap_pairs))
    for (i in seq_len(nrow(overlap_pairs))) {
        pair <- overlap_pairs[i]
        bounds <- .apply_bounds(pair, strategy)

        chr_val   <- bounds$chr
        start_val <- as.integer(bounds$start)
        end_val   <- as.integer(bounds$end)

        # Collect SNPs within bounds from both sources
        snps_in <- data.table::rbindlist(list(
            if (!is.null(gea_sig_snps) && nrow(gea_sig_snps) > 0) {
                gea_sig_snps[as.character(chr) == chr_val &
                             as.integer(pos) >= start_val &
                             as.integer(pos) <= end_val]
            } else data.table::data.table(),
            if (!is.null(gwas_sig_snps) && nrow(gwas_sig_snps) > 0) {
                gwas_sig_snps[as.character(chr) == chr_val &
                              as.integer(pos) >= start_val &
                              as.integer(pos) <= end_val]
            } else data.table::data.table()
        ), use.names = TRUE, fill = TRUE)

        traits_str  <- if (nrow(snps_in) > 0) paste(sort(unique(snps_in$trait)),  collapse = ",") else ""
        methods_str <- if (nrow(snps_in) > 0) paste(sort(unique(snps_in$method)), collapse = ",") else ""
        snp_count   <- if (nrow(snps_in) > 0) length(unique(snps_in$SNPID))       else 0L
        snp_ids     <- if (nrow(snps_in) > 0) paste(unique(snps_in$SNPID), collapse = ",") else ""
        min_p       <- if (nrow(snps_in) > 0) min(snps_in$pvalue, na.rm = TRUE)   else NA_real_

        region_id <- paste0(chr_val, "_", start_val, "-", end_val)

        rows[[i]] <- data.table::data.table(
            region_id      = region_id,
            chr            = chr_val,
            start          = start_val,
            end            = end_val,
            length         = end_val - start_val,
            snp_count      = snp_count,
            snp_ids        = snp_ids,
            traits         = traits_str,
            methods        = methods_str,
            min_pvalue     = min_p,
            gea_region_id  = pair$gea_region_id,
            gwas_region_id = pair$gwas_region_id
        )
    }

    result <- data.table::rbindlist(rows, use.names = TRUE)

    # Deduplicate by region_id (same bounds can arise from multiple pairs)
    if (nrow(result) > 1) {
        result <- result[, .SD[1], by = region_id]
    }

    data.table::setorder(result, chr, start)
    result
}

#' Empty overlap regions skeleton
#' @noRd
.empty_overlap_regions <- function() {
    data.table::data.table(
        region_id      = character(),
        chr            = character(),
        start          = integer(),
        end            = integer(),
        length         = integer(),
        snp_count      = integer(),
        snp_ids        = character(),
        traits         = character(),
        methods        = character(),
        min_pvalue     = numeric(),
        gea_region_id  = character(),
        gwas_region_id = character()
    )
}

# ── Unified sig SNP region assignment ─────────────────────────────────────────

#' Assign region_ids to sig SNPs based on a computed regions table
#'
#' Replaces the pipeline-file-based assign_region_ids() with an on-the-fly stamp
#' from live-computed regions (works for both the GEA/GWAS region explorer and the
#' overlapping tab's overlap regions). SNPs that fall within a region's bounds get
#' that region's region_id. SNPs outside all regions get NA.
#'
#' @param snps_dt data.table with columns SNPID, chr (character), pos (integer)
#' @param regions data.table with columns chr, start, end, region_id (e.g. from
#'   compute_all_regions() or compute_all_overlap_regions())
#' @return snps_dt with region_id column set/updated
#' @noRd
assign_region_ids_from_regions <- function(snps_dt, regions) {
    snps_dt[, region_id := NA_character_]
    if (is.null(regions) || nrow(regions) == 0) return(snps_dt)
    if (!all(c("chr", "start", "end", "region_id") %in% names(regions))) return(snps_dt)

    regs_ov <- data.table::data.table(
        chr       = as.character(regions$chr),
        start     = as.integer(regions$start),
        end       = as.integer(regions$end),
        region_id = as.character(regions$region_id)
    )

    snp_pos <- unique(snps_dt[, .(SNPID, chr = as.character(chr), pos = as.integer(pos))])
    snp_ov <- data.table::data.table(
        SNPID = snp_pos$SNPID,
        chr   = snp_pos$chr,
        start = snp_pos$pos,
        end   = snp_pos$pos
    )

    data.table::setkey(regs_ov, chr, start, end)
    data.table::setkey(snp_ov,  chr, start, end)

    ov <- data.table::foverlaps(snp_ov, regs_ov,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "within", mult = "first", nomatch = NA)

    rid_map <- ov[, .(SNPID, region_id)]
    snps_dt[rid_map, region_id := i.region_id, on = "SNPID"]
    snps_dt[, region_id := as.character(region_id)]
    snps_dt
}
