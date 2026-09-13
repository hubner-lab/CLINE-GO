# genes_in_regions.R — gene-region overlap and exon/promoter SNP counting
# Usage: source("/pipeline/scripts/R/lib/genes_in_regions.R")
# Requires: data.table, dplyr
# Assumes gff_parsing.R has been sourced (for extract_gene_id, read_gff, read_gff_exons)
#
# data.table port of find_genes_around_regions.R.
# Uses foverlaps instead of GenomicRanges::findOverlaps.

# Find genes overlapping a set of regions and count exon/promoter SNPs.
#
# @param regions_dt      data.table: region_id, chr, start, end (already extended)
# @param gff_dt          data.table from read_gff(): chr, start, end, gene_id, ...
# @param gff_path        character: path to GFF3 file (for loading exon features)
# @param promoter_length integer: bp upstream of gene start to count as promoter
# @param allsnps_dt      data.table: chr, pos (all VCF SNPs for counting)
# @return named list:
#   genes_per_region   data.table: one row per region-gene pair
#   genes_collapsed    data.table: one row per gene (across all regions)
find_genes_for_regions <- function(regions_dt, gff_dt, gff_path,
                                   promoter_length = 10000L,
                                   allsnps_dt = NULL) {
    promoter_length <- as.integer(promoter_length)

    if (nrow(regions_dt) == 0L) {
        message("WARNING: No regions to process")
        empty <- .empty_genes_dt()
        return(list(genes_per_region = empty, genes_collapsed = empty))
    }
    if (nrow(gff_dt) == 0L) {
        message("WARNING: No gene features in GFF")
        empty <- .empty_genes_dt()
        return(list(genes_per_region = empty, genes_collapsed = empty))
    }

    message(paste0("INFO: Processing ", nrow(regions_dt), " regions"))

    regions <- data.table::copy(regions_dt)
    genes   <- data.table::copy(gff_dt)
    regions[, chr := as.character(chr)]
    genes[,   chr := as.character(chr)]

    # --- gene-region overlap via foverlaps ---
    # x = regions (query), y = genes (keyed subject)
    # After foverlaps: plain cols = y (genes), i.* cols = x (regions)
    data.table::setkey(genes,   chr, start, end)

    reg_query <- regions[, .(region_id, chr, start, end)]
    data.table::setkey(reg_query, chr, start, end)

    ov <- data.table::foverlaps(reg_query, genes,
        by.x = c("chr", "start", "end"),
        by.y = c("chr", "start", "end"),
        type = "any", mult = "all", nomatch = NULL)

    if (is.null(ov) || nrow(ov) == 0L) {
        message("WARNING: No genes found overlapping any region")
        empty <- .empty_genes_dt()
        return(list(genes_per_region = empty, genes_collapsed = empty))
    }

    # Rename: plain start/end = gene coords (from y); i.* = region coords from x (drop)
    data.table::setnames(ov, c("start", "end"), c("gene_start", "gene_end"))
    drop_cols <- intersect(c("i.start", "i.end", "i.chr"), colnames(ov))
    if (length(drop_cols) > 0) ov[, (drop_cols) := NULL]

    message(paste0("INFO: Found ", nrow(ov), " gene-region pairs"))

    genes_per_region <- data.table::copy(ov)

    # --- exon/promoter SNP counting ---
    if (!is.null(allsnps_dt) && nrow(allsnps_dt) > 0L) {
        message("INFO: Counting exon/promoter SNPs")
        snps <- data.table::copy(allsnps_dt)
        snps[, chr := as.character(chr)]
        snps[, pos := as.integer(pos)]

        unique_gene_ids <- unique(genes_per_region$gene_id)

        # Exon SNPs
        exon_gff <- tryCatch(
            read_gff_exons(gff_path),
            error = function(e) {
                message("WARNING: Could not load exons: ", e$message)
                data.table::data.table(chr = character(), start = integer(),
                                       end = integer(), gene_id = character())
            }
        )
        exon_counts <- .count_snps_in_features(
            snps, exon_gff[gene_id %in% unique_gene_ids], col_name = "exon"
        )

        # Promoter SNPs: upstream of the TRANSCRIPTION START, which is gene_start
        # on the + strand and gene_end on the -. A gff_dt with no strand column,
        # or a '.' strand, is treated as + — the previous, strand-blind behaviour.
        gene_rows <- genes[gene_id %in% unique_gene_ids]
        g_strand <- if ("strand" %in% colnames(gene_rows)) {
            as.character(gene_rows$strand)
        } else {
            rep(NA_character_, nrow(gene_rows))
        }
        on_minus <- !is.na(g_strand) & g_strand == "-"
        promoter_dt <- data.table::data.table(
            chr     = gene_rows$chr,
            start   = data.table::fifelse(
                on_minus, as.integer(gene_rows$end),
                pmax(1L, as.integer(gene_rows$start) - promoter_length)),
            end     = data.table::fifelse(
                on_minus, as.integer(gene_rows$end) + promoter_length,
                as.integer(gene_rows$start)),
            gene_id = gene_rows$gene_id
        )
        promoter_counts <- .count_snps_in_features(snps, promoter_dt, col_name = "promoter")

        genes_per_region <- merge(genes_per_region, exon_counts,    by = "gene_id", all.x = TRUE)
        genes_per_region <- merge(genes_per_region, promoter_counts, by = "gene_id", all.x = TRUE)
    } else {
        genes_per_region[, exon_snps    := ""]
        genes_per_region[, promoter_snps := ""]
    }

    # Fill NAs and add count columns
    genes_per_region[is.na(exon_snps),    exon_snps    := ""]
    genes_per_region[is.na(promoter_snps), promoter_snps := ""]
    genes_per_region[, exon_snp_count := vapply(
        exon_snps, function(x) if (x == "") 0L else length(strsplit(x, ",")[[1]]), integer(1L)
    )]
    genes_per_region[, promoter_snp_count := vapply(
        promoter_snps, function(x) if (x == "") 0L else length(strsplit(x, ",")[[1]]), integer(1L)
    )]

    # Reorder columns: region_id first
    key_first <- c("region_id", "gene_id", "chr", "gene_start", "gene_end")
    rest_cols <- setdiff(colnames(genes_per_region), key_first)
    data.table::setcolorder(genes_per_region, c(key_first, rest_cols))

    # --- collapsed table: one row per gene ---
    collapse_cols <- setdiff(colnames(genes_per_region),
                             c("region_id", "gene_id", "chr", "gene_start", "gene_end",
                               "exon_snps", "promoter_snps",
                               "exon_snp_count", "promoter_snp_count"))

    genes_collapsed <- genes_per_region[, c(
        .(region_id          = paste(sort(unique(region_id)), collapse = ",")),
        .(exon_snps          = paste(sort(unique(exon_snps[exon_snps != ""])), collapse = ",")),
        .(promoter_snps      = paste(sort(unique(promoter_snps[promoter_snps != ""])), collapse = ",")),
        .(exon_snp_count     = max(exon_snp_count)),
        .(promoter_snp_count = max(promoter_snp_count)),
        lapply(.SD, function(x) paste(sort(unique(x[!is.na(x)])), collapse = ","))
    ), by = c("gene_id", "chr", "gene_start", "gene_end"), .SDcols = collapse_cols]

    data.table::setorder(genes_collapsed, chr, gene_start)

    n_genes <- data.table::uniqueN(genes_per_region$gene_id)
    message(paste0("INFO: ", n_genes, " unique genes across ", nrow(regions_dt), " regions"))

    list(genes_per_region = genes_per_region, genes_collapsed = genes_collapsed)
}

# --- internal helpers ---

# foverlaps-based count of SNPs overlapping feature ranges, grouped by gene.
# Returns data.table: gene_id, <col_name>_snps (comma-sep "chr:pos")
.count_snps_in_features <- function(snps_dt, feat_dt, col_name) {
    snps_col <- paste0(col_name, "_snps")
    empty_result <- data.table::data.table(gene_id = character())
    empty_result[[snps_col]] <- character()

    if (is.null(feat_dt) || nrow(feat_dt) == 0L) return(empty_result)

    snps  <- data.table::copy(snps_dt)
    feats <- data.table::copy(feat_dt)

    snps[,  s := as.integer(pos)][, e := as.integer(pos)]
    feats[, s := as.integer(start)][, e := as.integer(end)]

    data.table::setkey(feats, chr, s, e)
    data.table::setkey(snps,  chr, s, e)

    ov <- data.table::foverlaps(snps, feats,
        by.x = c("chr", "s", "e"),
        by.y = c("chr", "s", "e"),
        type = "within", mult = "all", nomatch = NULL)

    if (is.null(ov) || nrow(ov) == 0L) return(empty_result)

    # x = snps, y = feats, so after foverlaps the SNP position is i.s and the
    # plain s is the FEATURE start. The id must name the SNP.
    ov[, snp_id := paste0(chr, ":", i.s)]
    result <- ov[, .(snps = paste(sort(unique(snp_id)), collapse = ",")), by = "gene_id"]
    data.table::setnames(result, "snps", snps_col)
    result
}

.empty_genes_dt <- function() {
    data.table::data.table(
        region_id = character(), gene_id = character(), chr = character(),
        gene_start = integer(), gene_end = integer()
    )
}
