# regions.R — single-linkage region clustering and combined region builder
# Usage: source("/pipeline/scripts/R/lib/regions.R")
# Requires: data.table, stringr
#
# data.table port of the GRanges pipeline behavior in create_regions.R.
# Clustering criterion: two consecutive SNPs (sorted by pos) merge when
# their gap <= 2*distance, matching the GRanges "extend by dist + reduce" approach.
# Region boundaries: [min(pos)-distance, max(pos)+distance].

# Resolve per-chromosome distance from a scalar or named integer vector.
.get_dist <- function(dist_spec, chr_name) {
    if (length(dist_spec) == 1 && is.null(names(dist_spec))) {
        return(as.integer(dist_spec))
    }
    chr_name <- as.character(chr_name)
    if (chr_name %in% names(dist_spec)) {
        return(as.integer(dist_spec[[chr_name]]))
    }
    fallback <- max(dist_spec, na.rm = TRUE)
    message(paste0("WARNING: chr '", chr_name, "' not in per-chr distance map; using max=", fallback, " bp"))
    return(as.integer(fallback))
}

# Cluster significant SNPs into genomic regions using single-linkage clustering.
#
# Matches create_regions.R create_regions_from_snps() GRanges behavior:
#   extend each SNP by dist → reduce → equivalent to merging when gap <= 2*dist.
#
# @param sig_snps    data.table with cols SNPID, chr, pos, min_pvalue, and
#                    method columns (each method col holds comma-sep trait values
#                    as produced by combine_selected_snps.R)
# @param dist_spec   scalar integer OR named integer vector (per-chromosome distances)
#                    as returned by resolve_clumping_distance()
# @param trait_label character or NULL. When non-NULL, appended to region_id as
#                    e.g. "1_100-500_bio_1". When NULL: combined regions (no suffix).
# @return data.table:
#   region_id, trait, chr, start, end, length, snp_count, snp_ids, methods, min_pvalue
cluster_snps_to_regions <- function(sig_snps, dist_spec, trait_label = NULL) {
    if (is.null(sig_snps) || nrow(sig_snps) == 0) {
        return(.empty_region_dt(with_trait = TRUE))
    }

    sig_snps <- data.table::copy(sig_snps)
    sig_snps[, chr := as.character(chr)]
    sig_snps[, pos := as.integer(pos)]
    if (!"min_pvalue" %in% colnames(sig_snps)) {
        # Almost always the wide-vs-long shape trap: the LONG sig-SNP table has no
        # min_pvalue column, so every region would report an empty p-value. Say it
        # out loud — this used to surface only as min()'s -Inf warning.
        message("WARNING: sig_snps has no min_pvalue column; every region will report NA")
        sig_snps[, min_pvalue := NA_real_]
    }

    method_cols <- setdiff(colnames(sig_snps), c("SNPID", "chr", "pos", "min_pvalue"))

    chrs <- unique(sig_snps$chr)
    all_regions <- lapply(chrs, function(this_chr) {
        dist <- .get_dist(dist_spec, this_chr)

        chr_snps <- sig_snps[chr == this_chr]
        data.table::setorder(chr_snps, pos)

        # Single-linkage: split cluster when consecutive gap > 2*dist
        # (matches GRanges extend-by-dist + reduce: two SNPs merge if gap <= 2*dist)
        chr_snps[, gap     := c(Inf, diff(pos))]
        chr_snps[, cluster := cumsum(gap > 2L * dist)]
        chr_snps[, gap     := NULL]

        lapply(unique(chr_snps$cluster), function(cl) {
            grp <- chr_snps[cluster == cl]
            region_start <- pmax(1L, min(grp$pos) - dist)
            region_end   <- max(grp$pos) + dist

            # Active methods: method columns with non-empty values
            active_methods <- character(0)
            for (m in method_cols) {
                vals <- grp[[m]]
                if (any(!is.na(vals) & vals != "" & vals != '""')) {
                    active_methods <- c(active_methods, m)
                }
            }

            rid <- if (!is.null(trait_label)) {
                paste0(this_chr, "_", region_start, "-", region_end, "_", trait_label)
            } else {
                paste0(this_chr, "_", region_start, "-", region_end)
            }

            data.table::data.table(
                region_id  = rid,
                trait      = if (!is.null(trait_label)) trait_label else "",
                chr        = this_chr,
                start      = as.integer(region_start),
                end        = as.integer(region_end),
                length     = as.integer(region_end - region_start),
                snp_count  = nrow(grp),
                snp_ids    = paste(grp$SNPID, collapse = ","),
                methods    = paste(sort(unique(active_methods)), collapse = ","),
                min_pvalue = if (all(is.na(grp$min_pvalue))) NA_real_
                             else min(grp$min_pvalue, na.rm = TRUE)
            )
        })
    })

    regions <- data.table::rbindlist(
        Filter(Negate(is.null), unlist(all_regions, recursive = FALSE))
    )
    if (nrow(regions) > 0) data.table::setorder(regions, chr, start)
    regions
}

# Build per-trait regions with cross-trait evidence columns.
#
# @param sig_snps   combined selected-SNPs table (from combine_selected_snps.R):
#                   SNPID, chr, pos, min_pvalue, <method_cols>
# @param dist_spec  scalar or named-vector distance (from resolve_clumping_distance)
# @return data.table: per-trait regions with columns
#   region_id, trait, chr, start, end, length, snp_count, snp_ids, methods,
#   min_pvalue, other_traits, other_snp_count
build_per_trait_regions <- function(sig_snps, dist_spec) {
    if (is.null(sig_snps) || nrow(sig_snps) == 0) return(.empty_region_dt(with_trait = TRUE))

    sig_snps  <- data.table::copy(sig_snps)
    sig_snps[, chr := as.character(chr)]
    method_cols <- setdiff(colnames(sig_snps), c("SNPID", "chr", "pos", "min_pvalue"))

    # Extract all unique traits from method columns
    all_traits <- character(0)
    for (m in method_cols) {
        vals <- sig_snps[[m]]
        t_m  <- unique(unlist(strsplit(
            gsub('"', '', vals[!is.na(vals) & vals != ""]),
            ","
        )))
        t_m  <- t_m[t_m != ""]
        all_traits <- c(all_traits, t_m)
    }
    all_traits <- unique(all_traits)

    message(paste0("INFO: Found ", length(all_traits), " traits: ",
                   paste(all_traits, collapse = ", ")))

    per_trait_list <- lapply(all_traits, function(trait) {
        message(paste0("INFO:   Processing trait: ", trait))
        pattern <- paste0("(^|,)", trait, "($|,)")

        trait_snps <- sig_snps[0, ]
        for (m in method_cols) {
            rows <- sig_snps[grepl(pattern, sig_snps[[m]]), ]
            if (nrow(rows) > 0) trait_snps <- unique(rbind(trait_snps, rows))
        }

        if (nrow(trait_snps) == 0) {
            message(paste0("INFO:     No SNPs for ", trait))
            return(NULL)
        }

        regions <- cluster_snps_to_regions(trait_snps, dist_spec, trait_label = trait)
        if (!is.null(regions) && nrow(regions) > 0)
            message(paste0("INFO:     Created ", nrow(regions), " regions for ", trait))
        regions
    })

    per_trait_regions <- data.table::rbindlist(Filter(Negate(is.null), per_trait_list))

    if (is.null(per_trait_regions) || nrow(per_trait_regions) == 0) {
        return(.empty_region_dt(with_trait = TRUE))
    }

    # Cross-trait evidence: for each region, count SNPs from OTHER traits
    message("INFO: Computing cross-trait evidence...")
    n <- nrow(per_trait_regions)
    other_traits_col    <- character(n)
    other_snp_count_col <- integer(n)

    for (i in seq_len(n)) {
        region       <- per_trait_regions[i, ]
        region_trait <- region$trait
        region_chr   <- region$chr

        snps_in <- sig_snps[chr == region_chr &
                            pos >= region$start &
                            pos <= region$end, ]

        if (nrow(snps_in) == 0) {
            other_traits_col[i]    <- ""
            other_snp_count_col[i] <- 0L
            next
        }

        cross_traits  <- character(0)
        cross_snp_ids <- character(0)

        for (m in method_cols) {
            for (j in seq_len(nrow(snps_in))) {
                val <- snps_in[[m]][j]
                if (!is.na(val) && val != "" && val != '""') {
                    traits_m <- strsplit(gsub('"', "", val), ",")[[1]]
                    traits_m <- traits_m[traits_m != ""]
                    other_t  <- traits_m[traits_m != region_trait]
                    if (length(other_t) > 0) {
                        cross_traits  <- c(cross_traits, other_t)
                        cross_snp_ids <- c(cross_snp_ids, snps_in$SNPID[j])
                    }
                }
            }
        }

        other_traits_col[i]    <- paste(unique(cross_traits), collapse = ",")
        other_snp_count_col[i] <- length(unique(cross_snp_ids))
    }

    per_trait_regions[, other_traits    := other_traits_col]
    per_trait_regions[, other_snp_count := other_snp_count_col]
    per_trait_regions
}

# Build combined regions (all traits together) with per-trait and per-method SNP counts.
#
# @param sig_snps   combined selected-SNPs table
# @param dist_spec  scalar or named-vector distance
# @return data.table: combined regions with traits column + per-trait/method _snps columns
build_combined_regions <- function(sig_snps, dist_spec) {
    if (is.null(sig_snps) || nrow(sig_snps) == 0) return(.empty_combined_region_dt())

    combined <- cluster_snps_to_regions(sig_snps, dist_spec, trait_label = NULL)

    if (is.null(combined) || nrow(combined) == 0) return(.empty_combined_region_dt())

    sig_snps    <- data.table::copy(sig_snps)
    sig_snps[, chr := as.character(chr)]
    method_cols <- setdiff(colnames(sig_snps), c("SNPID", "chr", "pos", "min_pvalue"))

    n <- nrow(combined)
    traits_col <- character(n)
    row_data   <- vector("list", n)

    for (i in seq_len(n)) {
        snp_ids <- strsplit(combined$snp_ids[i], ",")[[1]]
        region_traits  <- character(0)
        trait_snp_map  <- list()
        method_snp_map <- list()

        for (snp_id in snp_ids) {
            snp_row <- sig_snps[SNPID == snp_id, ]
            for (m in method_cols) {
                val <- snp_row[[m]]
                if (!is.na(val) && val != "" && val != '""') {
                    traits_m <- strsplit(gsub('"', "", val), ",")[[1]]
                    traits_m <- traits_m[traits_m != ""]
                    region_traits <- c(region_traits, traits_m)
                    for (t in traits_m)
                        trait_snp_map[[t]]  <- c(trait_snp_map[[t]],  snp_id)
                    method_snp_map[[m]] <- c(method_snp_map[[m]], snp_id)
                }
            }
        }

        traits_col[i] <- paste(unique(region_traits[region_traits != ""]), collapse = ",")
        row_data[[i]] <- list(
            trait_counts  = vapply(trait_snp_map,  function(x) length(unique(x)), integer(1L)),
            method_counts = vapply(method_snp_map, function(x) length(unique(x)), integer(1L))
        )
    }

    combined[, traits := traits_col]
    combined[, trait  := NULL]

    all_trait_names  <- sort(unique(unlist(lapply(row_data, function(x) names(x$trait_counts)))))
    all_method_names <- sort(unique(unlist(lapply(row_data, function(x) names(x$method_counts)))))

    for (t in all_trait_names) {
        combined[[paste0(t, "_snps")]] <- vapply(row_data, function(x) {
            if (t %in% names(x$trait_counts)) x$trait_counts[[t]] else 0L
        }, integer(1L))
    }
    for (m in all_method_names) {
        combined[[paste0(m, "_snps")]] <- vapply(row_data, function(x) {
            if (m %in% names(x$method_counts)) x$method_counts[[m]] else 0L
        }, integer(1L))
    }

    combined
}

# --- empty table constructors ---

.empty_region_dt <- function(with_trait = TRUE) {
    data.table::data.table(
        region_id       = character(),
        trait           = character(),
        chr             = character(),
        start           = integer(),
        end             = integer(),
        length          = integer(),
        snp_count       = integer(),
        snp_ids         = character(),
        methods         = character(),
        min_pvalue      = numeric(),
        other_traits    = character(),
        other_snp_count = integer()
    )
}

.empty_combined_region_dt <- function() {
    data.table::data.table(
        region_id  = character(),
        chr        = character(),
        start      = integer(),
        end        = integer(),
        length     = integer(),
        snp_count  = integer(),
        snp_ids    = character(),
        traits     = character(),
        methods    = character(),
        min_pvalue = numeric()
    )
}
