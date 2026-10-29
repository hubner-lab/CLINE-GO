#!/usr/bin/env Rscript
# On-demand pairwise trait overlap computation for the Shiny app.
#
# Unlike the pipeline compute_pairwise_overlaps.R (which reads pre-collapsed
# selected_snps.tsv at a fixed threshold), this script takes long-format
# sig-SNP TSVs written by Shiny from the user's currently-selected threshold
# and clumping-distance settings.
#
# Inputs (long-format, one row per SNP per method per trait):
#   gea_sig_tsv  — columns: SNPID, chr, pos, pvalue, method, trait
#   gwas_sig_tsv — same format; "NULL" or missing file = no GWAS sig SNPs
#
# Window sizes are per-source so the computation mirrors what the user sees
# in each filter bar:
#   gea_window   — clumping distance from the GEA filter bar (bp)
#   gwas_window  — clumping distance from the GWAS filter bar (bp)
#   For GEA-GWAS pairs: max(gea_window, gwas_window) is used.
#
# Output files:
#   out_collapsed — pairwise_collapsed_snps.tsv  (long: one row per trait-SNP)
#   out_table     — pairwise_overlap_table.tsv   (one row per trait pair)
#
# Schema of both output files matches the pipeline's compute_pairwise_overlaps.R
# exactly so the Shiny DT renderer in mod_pairwise_overlap.R is unchanged.

library(data.table)
options(scipen = 99999)

args <- commandArgs(trailingOnly = TRUE)

GEA_SIG_TSV  <- args[1]
GWAS_SIG_TSV <- args[2]
GEA_WINDOW   <- as.integer(args[3])
GWAS_WINDOW  <- as.integer(args[4])
MIN_SNPS     <- as.integer(args[5])
OUT_COLLAPSED <- args[6]
OUT_TABLE     <- args[7]

message("INFO: Pairwise on-demand parameters:")
message("  GEA_SIG_TSV  = ", GEA_SIG_TSV)
message("  GWAS_SIG_TSV = ", GWAS_SIG_TSV)
message("  GEA_WINDOW   = ", GEA_WINDOW, " bp")
message("  GWAS_WINDOW  = ", GWAS_WINDOW, " bp")
message("  MIN_SNPS     = ", MIN_SNPS)

# ─────────────────────────────────────────────────────────────────────────────
# Empty output skeletons
# ─────────────────────────────────────────────────────────────────────────────
empty_collapsed <- data.table::data.table(
    trait      = character(), source     = character(), SNPID      = character(),
    chr        = character(), pos        = integer(),   min_pvalue = numeric(),
    n_methods  = integer()
)
empty_table <- data.table::data.table(
    trait_a         = character(), source_a        = character(), n_snps_a     = integer(),
    trait_b         = character(), source_b        = character(), n_snps_b     = integer(),
    exact_matches   = integer(),   window_overlaps = integer(),   overlap_score = integer(),
    comparison_type = character()
)

write_empty <- function(msg) {
    message("INFO: ", msg)
    data.table::fwrite(empty_collapsed, OUT_COLLAPSED, sep = "\t", quote = FALSE)
    data.table::fwrite(empty_table,     OUT_TABLE,     sep = "\t", quote = FALSE)
    quit(save = "no", status = 0)
}

# ─────────────────────────────────────────────────────────────────────────────
# Read long-format sig-SNP TSVs
# ─────────────────────────────────────────────────────────────────────────────
read_sig_tsv <- function(path, source_label) {
    if (is.null(path) || path == "NULL" || !file.exists(path)) return(NULL)
    dt <- tryCatch(
        data.table::fread(path, sep = "\t",
                          colClasses = c(SNPID = "character", chr = "character",
                                         method = "character", trait = "character")),
        error = function(e) NULL
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    if (!all(c("SNPID","chr","pos","pvalue","method","trait") %in% names(dt))) {
        message("WARNING: unexpected columns in ", path, " — skipping ", source_label)
        return(NULL)
    }
    dt[, source := source_label]
    message("INFO: Read ", nrow(dt), " rows from ", source_label, " (",
            length(unique(dt$trait)), " traits, ", length(unique(dt$method)), " methods)")
    dt
}

gea_dt  <- read_sig_tsv(GEA_SIG_TSV,  "climate")
gwas_dt <- read_sig_tsv(GWAS_SIG_TSV, "phenotype")

combined_raw <- data.table::rbindlist(list(gea_dt, gwas_dt), use.names = TRUE, fill = TRUE)

if (nrow(combined_raw) == 0)
    write_empty("No significant SNPs in either input. Writing empty outputs.")

# ─────────────────────────────────────────────────────────────────────────────
# Collapse to per-(trait, SNPID): keep min pvalue, count methods
# ─────────────────────────────────────────────────────────────────────────────
# Ensure pvalue column is numeric
combined_raw[, pvalue := suppressWarnings(as.numeric(pvalue))]

combined <- combined_raw[
    !is.na(pvalue),
    .(
        chr       = chr[1],
        pos       = as.integer(pos[1]),
        min_pvalue = min(pvalue, na.rm = TRUE),
        n_methods  = data.table::uniqueN(method)
    ),
    by = .(trait, SNPID, source)
]

# If same SNPID appears for the same trait from both sources, keep min pvalue row
combined <- combined[, .SD[which.min(min_pvalue)], by = .(trait, SNPID)]

# Re-derive source per trait: if all SNPs are climate → climate;
# all phenotype → phenotype; mixed → both
trait_source <- combined[, .(
    source = {
        srcs <- unique(source)
        if (length(srcs) == 1L) srcs else "both"
    }
), by = trait]
combined[, source := NULL]
combined <- trait_source[combined, on = "trait"]

data.table::setcolorder(combined, c("trait","source","SNPID","chr","pos","min_pvalue","n_methods"))
data.table::setorder(combined, trait, chr, pos)

message("INFO: Traits found: ", paste(sort(unique(combined$trait)), collapse = ", "))
message("INFO: Total trait-SNP rows: ", nrow(combined))

data.table::fwrite(combined, OUT_COLLAPSED, sep = "\t", quote = FALSE)
message("INFO: Saved pairwise_collapsed_snps.tsv (", nrow(combined), " rows)")

# ─────────────────────────────────────────────────────────────────────────────
# Single-linkage clustering (copied verbatim from compute_pairwise_overlaps.R)
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

# Windows (single-linkage clusters of both traits' SNPs at window_size) holding >= min_snps
# SNPs of each trait. A window that holds an exact shared SNP is also counted in
# `with_exact`, so the overlap score can count that locus once (ad7fce).
count_window_overlaps <- function(pos_a, pos_b, chr_a, chr_b, window_size, min_snps,
                                  id_a, id_b) {
    shared <- intersect(id_a, id_b)
    dt_a <- data.table::data.table(chr = chr_a, pos = pos_a, grp = "A", exact = id_a %in% shared)
    dt_b <- data.table::data.table(chr = chr_b, pos = pos_b, grp = "B", exact = id_b %in% shared)
    all_pos <- data.table::rbindlist(list(dt_a, dt_b))
    data.table::setorder(all_pos, chr, pos)
    n_overlap <- 0L
    n_with_exact <- 0L
    for (ch in unique(all_pos$chr)) {
        sub  <- all_pos[chr == ch]
        clus <- .cluster_snps(sub$pos, window_size)
        sub[, cluster := clus]
        cl_summary <- sub[, .(n_a = sum(grp == "A"), n_b = sum(grp == "B"), has_exact = any(exact)),
                          by = cluster]
        hit <- cl_summary[n_a >= min_snps & n_b >= min_snps]
        n_overlap    <- n_overlap + nrow(hit)
        n_with_exact <- n_with_exact + sum(hit$has_exact)
    }
    list(windows = n_overlap, with_exact = n_with_exact)
}

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise comparisons
# ─────────────────────────────────────────────────────────────────────────────
traits   <- sort(unique(combined$trait))
n_traits <- length(traits)

if (n_traits < 2L) {
    message("INFO: Only ", n_traits, " trait(s) found — no pairs to compare. Writing empty pairwise table.")
    data.table::fwrite(empty_table, OUT_TABLE, sep = "\t", quote = FALSE)
    quit(save = "no", status = 0)
}

message("INFO: Computing pairwise overlaps for ", n_traits, " traits (",
        choose(n_traits, 2L), " pairs)...")

snps_by_trait <- split(combined, by = "trait")

pair_rows <- vector("list", choose(n_traits, 2L))
idx <- 0L
for (i in 1L:(n_traits - 1L)) {
    for (j in (i + 1L):n_traits) {
        ta <- traits[i]; tb <- traits[j]
        dt_a <- snps_by_trait[[ta]]; dt_b <- snps_by_trait[[tb]]
        src_a <- dt_a$source[1L]; src_b <- dt_b$source[1L]

        is_gea_a <- src_a %in% c("climate", "both")
        is_gea_b <- src_b %in% c("climate", "both")
        ctype <- if (is_gea_a && is_gea_b)  "GEA-GEA"  else
                 if (!is_gea_a && !is_gea_b) "GWAS-GWAS" else "GEA-GWAS"

        # Window size: per-source type to match filter-bar clumping
        win <- if (ctype == "GEA-GEA")   GEA_WINDOW  else
               if (ctype == "GWAS-GWAS") GWAS_WINDOW else
               max(GEA_WINDOW, GWAS_WINDOW)

        exact  <- length(intersect(dt_a$SNPID, dt_b$SNPID))
        wovlap <- count_window_overlaps(dt_a$pos, dt_b$pos,
                                        dt_a$chr, dt_b$chr,
                                        win, MIN_SNPS, dt_a$SNPID, dt_b$SNPID)
        idx <- idx + 1L
        pair_rows[[idx]] <- data.table::data.table(
            trait_a         = ta,
            source_a        = src_a,
            n_snps_a        = nrow(dt_a),
            trait_b         = tb,
            source_b        = src_b,
            n_snps_b        = nrow(dt_b),
            exact_matches   = exact,
            window_overlaps = wovlap$windows,
            # Each shared SNP once, plus each window that holds none: a window around
            # shared SNPs used to add a second count for them (ad7fce).
            overlap_score   = exact + wovlap$windows - wovlap$with_exact,
            comparison_type = ctype
        )
    }
}

pairwise_dt <- data.table::rbindlist(pair_rows[seq_len(idx)])
data.table::setorder(pairwise_dt, -overlap_score, -exact_matches)

message("INFO: Pairs with any overlap: ", sum(pairwise_dt$overlap_score > 0), "/", nrow(pairwise_dt))

data.table::fwrite(pairwise_dt, OUT_TABLE, sep = "\t", quote = FALSE)
message("INFO: Saved pairwise_overlap_table.tsv (", nrow(pairwise_dt), " pairs)")
