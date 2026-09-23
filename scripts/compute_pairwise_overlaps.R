#!/usr/bin/env Rscript
# Compute pairwise trait overlaps from method-collapsed significant SNPs.
#
# Reads Selected_SNPs.tsv from GEA (association) and/or GWAS (phenotype_association).
# Each file has method columns with comma-separated trait names per SNP.
# Extracts per-trait SNP sets, then computes all pairwise overlaps.
#
# Outputs:
#   1. pairwise_collapsed_snps.tsv  - long format: one row per SNP per trait
#   2. pairwise_overlap_table.tsv   - one row per trait pair with overlap counts

library(data.table)
library(stringr)

# trait_in_cell(): exact (non-regex) trait membership on a comma-separated cell.
# Shared with build_per_trait_regions() so both sides of the pipeline decide
# "is this SNP significant for trait X" by the same rule.
source("/pipeline/scripts/R/lib/regions.R")

args <- commandArgs(trailingOnly = TRUE)
options(scipen = 99999)

GEA_SELECTED  <- args[1]   # association/tables/selected_snps.tsv or "NULL"
GWAS_SELECTED <- args[2]   # phenotype_association/tables/selected_snps.tsv or "NULL"
WINDOW_SIZE   <- as.integer(args[3])
MIN_SNPS      <- as.integer(args[4])
OUT_COLLAPSED <- args[5]
OUT_PAIRWISE  <- args[6]

message("INFO: Pairwise overlap parameters:")
message("  GEA_SELECTED  = ", GEA_SELECTED)
message("  GWAS_SELECTED = ", GWAS_SELECTED)
message("  WINDOW_SIZE   = ", WINDOW_SIZE, " bp")
message("  MIN_SNPS      = ", MIN_SNPS)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: Extract per-trait SNP sets from a Selected_SNPs wide-format table
# ─────────────────────────────────────────────────────────────────────────────
extract_trait_snps <- function(sel_path, source_label) {
    if (is.null(sel_path) || sel_path == "NULL" || !file.exists(sel_path)) return(NULL)
    dt <- data.table::fread(sel_path, colClasses = c(chr = "character", SNPID = "character"))
    if (nrow(dt) == 0) return(NULL)

    meta_cols   <- c("SNPID", "chr", "pos", "min_pvalue")
    method_cols <- setdiff(names(dt), meta_cols)
    message("INFO: Reading ", source_label, " selected SNPs (", nrow(dt), " SNPs, methods: ",
            paste(method_cols, collapse = ", "), ")")

    # Collect all unique traits across all method columns
    all_traits <- character(0)
    for (m in method_cols) {
        vals <- dt[[m]]
        vals <- gsub('"', '', vals)
        traits_m <- unique(unlist(strsplit(vals[vals != "" & !is.na(vals)], ",")))
        all_traits <- unique(c(all_traits, traits_m))
    }
    all_traits <- all_traits[all_traits != "" & !is.na(all_traits)]
    if (length(all_traits) == 0) return(NULL)

    # For each trait, collect SNPs and count confirming methods
    rows <- lapply(all_traits, function(trait) {
        method_hit <- vapply(method_cols, function(m) {
            trait_in_cell(dt[[m]], trait)
        }, logical(nrow(dt)))
        if (is.vector(method_hit)) method_hit <- matrix(method_hit, ncol = length(method_cols))
        any_hit <- rowSums(method_hit) > 0
        if (!any(any_hit)) return(NULL)
        snps <- dt[any_hit, .(
            SNPID      = SNPID,
            chr        = chr,
            pos        = as.integer(pos),
            min_pvalue = min_pvalue,
            n_methods  = as.integer(rowSums(method_hit[any_hit, , drop = FALSE])),
            source     = source_label,
            trait      = trait
        )]
        snps
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0) return(NULL)
    data.table::rbindlist(rows)
}

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────
gea_dt  <- extract_trait_snps(GEA_SELECTED,  "climate")
gwas_dt <- extract_trait_snps(GWAS_SELECTED, "phenotype")

combined <- data.table::rbindlist(list(gea_dt, gwas_dt), use.names = TRUE, fill = TRUE)
if (nrow(combined) == 0) {
    message("WARNING: No significant SNPs found in either input. Writing empty outputs.")
    empty_collapsed <- data.table::data.table(
        trait = character(), source = character(), SNPID = character(),
        chr = character(), pos = integer(), min_pvalue = numeric(), n_methods = integer()
    )
    empty_pairwise <- data.table::data.table(
        trait_a = character(), source_a = character(), n_snps_a = integer(),
        trait_b = character(), source_b = character(), n_snps_b = integer(),
        exact_matches = integer(), window_overlaps = integer(),
        overlap_score = integer(), comparison_type = character()
    )
    data.table::fwrite(empty_collapsed, OUT_COLLAPSED, sep = "\t", quote = FALSE)
    data.table::fwrite(empty_pairwise,  OUT_PAIRWISE,  sep = "\t", quote = FALSE)
    quit(save = "no", status = 0)
}

# Deduplicate: if same SNPID appears for the same trait from both sources, keep min pvalue
combined <- combined[, .SD[which.min(min_pvalue)], by = .(trait, SNPID)]

# Resolve source: if trait appears in both, keep more specific (first occurrence wins, already deduped)
# Re-derive source per trait: if all SNPs are climate → climate; all phenotype → phenotype; mixed → both
trait_source <- combined[, .(
    source = {
        srcs <- unique(source)
        if (length(srcs) == 1) srcs else "both"
    }
), by = trait]
combined <- combined[, source := NULL]
combined <- trait_source[combined, on = "trait"]

data.table::setcolorder(combined, c("trait", "source", "SNPID", "chr", "pos", "min_pvalue", "n_methods"))
data.table::setorder(combined, trait, chr, pos)

message("INFO: Traits found: ", paste(sort(unique(combined$trait)), collapse = ", "))
message("INFO: Total trait-SNP rows: ", nrow(combined))

data.table::fwrite(combined, OUT_COLLAPSED, sep = "\t", quote = FALSE)
message("INFO: Saved pairwise_collapsed_snps.tsv (", nrow(combined), " rows)")

# ─────────────────────────────────────────────────────────────────────────────
# Single-linkage clustering for window overlap detection
# ─────────────────────────────────────────────────────────────────────────────
.cluster_snps <- function(pos, window_size) {
    n <- length(pos)
    if (n == 0L) return(integer(0))
    cluster <- integer(n)
    k <- 1L
    cluster[1L] <- k
    if (n > 1L) {
        for (i in 2L:n) {
            if (is.na(pos[i]) || is.na(pos[i - 1L]) || (pos[i] - pos[i - 1L]) > window_size) {
                k <- k + 1L
            }
            cluster[i] <- k
        }
    }
    cluster
}

count_window_overlaps <- function(pos_a, pos_b, chr_a, chr_b, window_size, min_snps) {
    dt_a <- data.table::data.table(chr = chr_a, pos = pos_a, grp = "A")
    dt_b <- data.table::data.table(chr = chr_b, pos = pos_b, grp = "B")
    all_pos <- data.table::rbindlist(list(dt_a, dt_b))
    data.table::setorder(all_pos, chr, pos)
    n_overlap <- 0L
    for (ch in unique(all_pos$chr)) {
        sub  <- all_pos[chr == ch]
        clus <- .cluster_snps(sub$pos, window_size)
        sub[, cluster := clus]
        cl_summary <- sub[, .(n_a = sum(grp == "A"), n_b = sum(grp == "B")), by = cluster]
        n_overlap <- n_overlap + nrow(cl_summary[n_a >= min_snps & n_b >= min_snps])
    }
    n_overlap
}

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise comparisons
# ─────────────────────────────────────────────────────────────────────────────
traits <- sort(unique(combined$trait))
n_traits <- length(traits)

if (n_traits < 2L) {
    message("INFO: Only ", n_traits, " trait(s) — no pairs to compare. Writing empty pairwise table.")
    empty_pairwise <- data.table::data.table(
        trait_a = character(), source_a = character(), n_snps_a = integer(),
        trait_b = character(), source_b = character(), n_snps_b = integer(),
        exact_matches = integer(), window_overlaps = integer(),
        overlap_score = integer(), comparison_type = character()
    )
    data.table::fwrite(empty_pairwise, OUT_PAIRWISE, sep = "\t", quote = FALSE)
    quit(save = "no", status = 0)
}

message("INFO: Computing pairwise overlaps for ", n_traits, " traits (",
        choose(n_traits, 2), " pairs)...")

# Pre-index SNPs per trait for speed
snps_by_trait <- split(combined, by = "trait")

pair_rows <- vector("list", choose(n_traits, 2))
idx <- 0L
for (i in 1L:(n_traits - 1L)) {
    for (j in (i + 1L):n_traits) {
        ta <- traits[i]
        tb <- traits[j]
        dt_a <- snps_by_trait[[ta]]
        dt_b <- snps_by_trait[[tb]]
        src_a <- dt_a$source[1L]
        src_b <- dt_b$source[1L]

        # Comparison type
        is_gea_a <- src_a %in% c("climate", "both")
        is_gea_b <- src_b %in% c("climate", "both")
        ctype <- if (is_gea_a && is_gea_b) "GEA-GEA" else
                 if (!is_gea_a && !is_gea_b) "GWAS-GWAS" else "GEA-GWAS"

        exact  <- length(intersect(dt_a$SNPID, dt_b$SNPID))
        wovlap <- count_window_overlaps(dt_a$pos, dt_b$pos,
                                        dt_a$chr, dt_b$chr,
                                        WINDOW_SIZE, MIN_SNPS)
        idx <- idx + 1L
        pair_rows[[idx]] <- data.table::data.table(
            trait_a         = ta,
            source_a        = src_a,
            n_snps_a        = nrow(dt_a),
            trait_b         = tb,
            source_b        = src_b,
            n_snps_b        = nrow(dt_b),
            exact_matches   = exact,
            window_overlaps = wovlap,
            overlap_score   = exact + wovlap,
            comparison_type = ctype
        )
    }
}

pairwise_dt <- data.table::rbindlist(pair_rows[seq_len(idx)])
data.table::setorder(pairwise_dt, -overlap_score, -exact_matches)

message("INFO: Pairs with any overlap: ", sum(pairwise_dt$overlap_score > 0), "/", nrow(pairwise_dt))

data.table::fwrite(pairwise_dt, OUT_PAIRWISE, sep = "\t", quote = FALSE)
message("INFO: Saved pairwise_overlap_table.tsv (", nrow(pairwise_dt), " pairs)")
