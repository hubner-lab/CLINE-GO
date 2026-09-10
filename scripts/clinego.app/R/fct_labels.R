#' Build region dropdown labels with enrichment summary
#' Format: "chr_start-end (N SNPs, XE YP, Z GO [HAP])"
#' @param regions data.table from load_regions_per_trait or load_regions_combined
#' @param genes   data.table from load_genes (optional, for exon/promoter counts)
#' @param enrich  data.table (optional, for GO term counts per region)
#' @param hap_regions character vector of region_ids that have haplotype data
#' @noRd
build_region_labels <- function(regions, genes = NULL, enrich = NULL,
                                  hap_regions = character(0)) {
    if (is.null(regions) || nrow(regions) == 0) return(character(0))

    labels <- vapply(seq_len(nrow(regions)), function(i) {
        rid <- regions$region_id[i]
        n_snps <- regions$snp_count[i] %||% 0L

        # Exon/promoter counts from genes table
        xe <- yp <- 0L
        if (!is.null(genes) && nrow(genes) > 0 && "region_id" %in% names(genes)) {
            rg <- genes[genes$region_id == rid, ]
            if (nrow(rg) > 0) {
                xe <- sum(rg$exon_snp_count,     na.rm = TRUE)
                yp <- sum(rg$promoter_snp_count, na.rm = TRUE)
            }
        }

        # GO term count
        n_go <- 0L
        if (!is.null(enrich) && nrow(enrich) > 0 && "region_id" %in% names(enrich)) {
            n_go <- sum(enrich$region_id == rid, na.rm = TRUE)
        }

        hap_tag <- if (rid %in% hap_regions) " [HAP]" else ""

        parts <- paste0(n_snps, " SNPs",
                        if (xe > 0 || yp > 0) paste0(", ", xe, "E ", yp, "P") else "",
                        if (n_go > 0) paste0(", ", n_go, " GO") else "",
                        hap_tag)

        paste0(rid, " (", parts, ")")
    }, character(1))

    setNames(as.list(regions$region_id), labels)
}

#' Format a region ID for display: chr_start-end → chr:start-end
#' @noRd
format_region_id <- function(region_id) {
    # Convert underscore-based IDs back to colon format for display
    gsub("^(\\S+?)_(\\d+)-(\\d+)$", "\\1:\\2-\\3", region_id)
}

#' Format a haplotype tag for the tag selector display
#' Format: "site_association" → "Site / Association"
#' @noRd
format_hap_tag <- function(tag) {
    parts <- strsplit(tag, "_")[[1]]
    paste(tools::toTitleCase(parts), collapse = " / ")
}

#' Format numbers with SI suffixes (k, M) for display in value boxes
#' @noRd
format_si <- function(x) {
    if (is.na(x)) return("—")
    if (x >= 1e6) sprintf("%.1fM", x / 1e6)
    else if (x >= 1e3) sprintf("%.0fk", x / 1e3)
    else as.character(x)
}
