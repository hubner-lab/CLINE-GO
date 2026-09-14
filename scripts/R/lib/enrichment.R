# enrichment.R — GO enrichment computation (clusterProfiler)
# Usage: source("/pipeline/scripts/R/lib/enrichment.R")
# Requires: data.table, dplyr, stringr, tidyr, clusterProfiler, qs
# Assumes gff_parsing.R sourced (for extract_gene_id)

# Build TERM2GENE and TERM2NAME from a GFF3 file.
#
# @param gff_path    character: path to GFF3
# @param gff_feature character: feature type to load (e.g. "mRNA")
# @param go_field    character: attribute field name with GO terms (e.g. "ontology")
# @return named list: term2gene (data.frame: term, gene), term2name (data.frame: term, name)
#   term2name is NULL when GO.db is unavailable.
build_term2gene_from_gff <- function(gff_path, gff_feature, go_field) {
    if (is.null(go_field) || go_field %in% c("NULL", "")) {
        message("WARNING: go_field not specified — skipping enrichment")
        return(NULL)
    }

    message("INFO: Loading GFF for background GO terms")
    gff_raw <- data.table::fread(
        cmd    = paste("grep -v '#'", shQuote(gff_path)),
        header = FALSE
    )
    gff_raw <- gff_raw[V3 == gff_feature, .(description = V9)]

    if (nrow(gff_raw) == 0L) {
        message(paste0("WARNING: No features '", gff_feature, "' found in GFF"))
        return(NULL)
    }

    gff_data <- gff_raw[, .(
        gene_id  = extract_gene_id(description),
        go_terms = stringr::str_extract(description,
                                         paste0("(?<=", go_field, "=)[^;]+"))
    )]
    gff_data <- gff_data[!is.na(go_terms) & go_terms != "" & go_terms != "NA"]

    if (nrow(gff_data) == 0L) {
        message(paste0("WARNING: No GO terms found in field '", go_field, "'"))
        return(NULL)
    }

    message(paste0("INFO: GO annotations for ", data.table::uniqueN(gff_data$gene_id),
                   " background genes"))

    term2gene <- tidyr::separate_rows(as.data.frame(gff_data), go_terms, sep = ",")
    term2gene <- term2gene[term2gene$go_terms != "" & !is.na(term2gene$go_terms), ]
    term2gene <- data.frame(
        term = term2gene$go_terms,
        gene = term2gene$gene_id,
        stringsAsFactors = FALSE
    )

    message(paste0("INFO: TERM2GENE: ", nrow(term2gene), " pairs, ",
                   length(unique(term2gene$term)), " unique GO terms"))

    term2name <- get_go_descriptions(unique(term2gene$term))

    list(term2gene = term2gene, term2name = term2name,
         all_genes_with_go = unique(term2gene$gene))
}

# Fetch GO term descriptions from GO.db (if available).
# @param go_ids character vector of GO IDs
# @return data.frame(term, name) or NULL when GO.db not installed
get_go_descriptions <- function(go_ids) {
    if (requireNamespace("GO.db", quietly = TRUE)) {
        message("INFO: Using GO.db for term descriptions")
        # GO.db::GO.db, NOT library(GO.db). The attach is global and permanent:
        # it also pulls AnnotationDbi, whose select() then MASKS dplyr::select()
        # for everything sourced afterwards — in a pipeline where every script is
        # sourced into one session and several use a bare select(). The test
        # suite had to undo it by hand (test-enrichment.R's
        # with_search_path_restored). A :: reference needs no attach at all.
        go_db <- GO.db::GO.db
        # select() HARD-ERRORS when none of the keys resolve, which killed the whole
        # enrichment rule for a GFF whose go_field held a typo or an obsolete term.
        # Validate first and let unresolvable ids be their own name — the same
        # degradation the "GO.db not available" branch below already performs.
        as_own_name <- function(ids) {
            data.frame(term = ids, name = ids, stringsAsFactors = FALSE)
        }
        valid <- intersect(go_ids, AnnotationDbi::keys(go_db, keytype = "GOID"))
        if (length(valid) == 0L) {
            message(paste0("WARNING: none of the ", length(go_ids),
                           " GO ids resolve in GO.db — using GO IDs as descriptions"))
            return(as_own_name(go_ids))
        }
        if (length(valid) < length(go_ids)) {
            message(paste0("WARNING: ", length(go_ids) - length(valid), " of ",
                           length(go_ids), " GO ids do not resolve in GO.db — ",
                           "those keep their GO ID as the description"))
        }
        go_terms <- AnnotationDbi::select(go_db,
            keys    = valid,
            columns = c("GOID", "TERM"),
            keytype = "GOID"
        )
        term2name <- data.frame(
            term = go_terms$GOID,
            name = go_terms$TERM,
            stringsAsFactors = FALSE
        )
        unresolved <- setdiff(go_ids, term2name$term)
        if (length(unresolved) > 0L) {
            term2name <- rbind(term2name, as_own_name(unresolved))
        }
        term2name$name[is.na(term2name$name)] <- term2name$term[is.na(term2name$name)]
        return(term2name)
    }
    message("WARNING: GO.db not available — descriptions will be GO IDs")
    NULL
}

# Run GO enrichment for a single region.
#
# @param region_id        character
# @param trait            character (used for logging only)
# @param genes_in_region  character vector of gene IDs in this region
# @param term2gene        data.frame(term, gene)
# @param term2name        data.frame(term, name) or NULL
# @param all_genes_with_go  character vector: background gene set
# @return named list(region_id, trait, enrich_obj) or NULL if no enrichment
run_enrichment_for_region <- function(region_id, trait, genes_in_region,
                                       term2gene, term2name, all_genes_with_go) {
    message(paste0("INFO: Running enrichment for region: ", region_id))

    genes_with_go <- intersect(genes_in_region, all_genes_with_go)

    if (length(genes_with_go) < 2L) {
        message(paste0("WARNING: Region ", region_id, " has <2 genes with GO (",
                       length(genes_with_go), ")"))
        return(NULL)
    }

    message(paste0("INFO: Region ", region_id, ": ", length(genes_in_region),
                   " genes, ", length(genes_with_go), " with GO"))

    tryCatch({
        enrich_result <- clusterProfiler::enricher(
            gene          = genes_with_go,
            TERM2GENE     = term2gene,
            TERM2NAME     = term2name,
            pvalueCutoff  = 0.05,
            qvalueCutoff  = 0.2,
            pAdjustMethod = "BH",
            minGSSize     = 2L,
            maxGSSize     = 500L
        )

        if (is.null(enrich_result) || nrow(enrich_result@result) == 0L) {
            message(paste0("WARNING: No significant enrichment for region ", region_id))
            return(NULL)
        }

        if (all(is.na(enrich_result@result$qvalue))) {
            message("INFO: qvalue all NA — substituting p.adjust")
            enrich_result@result$qvalue <- enrich_result@result$p.adjust
        }

        message(paste0("INFO: Region ", region_id, ": ",
                       nrow(enrich_result@result), " enriched terms"))

        list(region_id = region_id, trait = trait, enrich_obj = enrich_result)
    }, error = function(e) {
        message(paste0("ERROR in enrichment for region ", region_id, ": ", e$message))
        NULL
    })
}

# Save enrichment result (TSV + qs) to disk.
#
# @param result  list(region_id, trait, enrich_obj) from run_enrichment_for_region
# @param tables_dir      base tables directory (per-trait subdirs created inside)
# @param intermediate_dir base intermediate directory (per-trait subdirs created inside)
save_enrichment_result <- function(result, tables_dir, intermediate_dir) {
    region_id <- result$region_id
    trait     <- result$trait
    enrich_obj <- result$enrich_obj

    trait_tables <- file.path(tables_dir, trait)
    trait_inter  <- file.path(intermediate_dir, trait)
    dir.create(trait_tables, recursive = TRUE, showWarnings = FALSE)
    dir.create(trait_inter,  recursive = TRUE, showWarnings = FALSE)

    # Save qs
    qs_file <- file.path(trait_inter, paste0("Region_", region_id, "_enrichment.qs"))
    qs::qsave(enrich_obj, qs_file)
    message(paste0("INFO: Saved .qs to ", qs_file))

    # Save TSV
    result_dt <- as.data.table(enrich_obj@result)[, .(
        GO_id      = ID,
        description = Description,
        gene_ratio = GeneRatio,
        bg_ratio   = BgRatio,
        pvalue, p_adjust = p.adjust, qvalue,
        gene_ids   = geneID,
        gene_count = Count
    )]
    result_dt[, region_id := region_id]
    result_dt[, trait     := trait]
    data.table::setcolorder(result_dt,
        c("region_id", "trait", "GO_id", "description",
          "gene_ratio", "bg_ratio", "pvalue", "p_adjust", "qvalue",
          "gene_ids", "gene_count"))

    tsv_file <- file.path(trait_tables, paste0("Region_", region_id, "_enrichment.tsv"))
    data.table::fwrite(result_dt, tsv_file, sep = "\t")
    message(paste0("INFO: Saved TSV to ", tsv_file))

    invisible(result_dt)
}
