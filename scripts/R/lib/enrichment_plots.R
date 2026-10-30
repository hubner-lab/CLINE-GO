# enrichment_plots.R — enrichment visualization helpers (emapplot, cnetplot, dotplot)
# Usage: source("/pipeline/scripts/R/lib/enrichment_plots.R")
# Requires: data.table, dplyr, ggplot2, clusterProfiler, enrichplot, qs

# Build gene_id → label mapping from genes_per_region table.
#
# @param genes_file  path to Genes_per_region.tsv or "NULL"
# @param label_field character: column name to use as label (e.g. "Name", "description")
# @return named character vector (label_field values keyed by gene_id) or NULL
build_label_map <- function(genes_file, label_field) {
    if (is.null(genes_file) || genes_file == "NULL" ||
        !file.exists(genes_file) || label_field == "gene_id") {
        return(NULL)
    }

    genes <- data.table::fread(genes_file)
    if (!label_field %in% colnames(genes)) {
        message(paste0("WARNING: '", label_field, "' not in genes table. Available: ",
                       paste(colnames(genes), collapse = ", ")))
        message("INFO: Falling back to gene_id labels")
        return(NULL)
    }

    label_dt <- unique(genes[!is.na(get(label_field)) & get(label_field) != "",
                              .(gene_id, label = get(label_field))])
    mapping  <- setNames(label_dt$label, label_dt$gene_id)
    message(paste0("INFO: Label map: ", length(mapping), " gene_id → ", label_field))
    mapping
}

# Remap gene IDs in an enrichResult to a label field.
remap_enrich_genes <- function(enrich_obj, label_map) {
    if (is.null(label_map)) return(enrich_obj)

    enrich_obj@result$geneID <- vapply(
        enrich_obj@result$geneID,
        function(genes_str) {
            gene_ids <- strsplit(genes_str, "/")[[1]]
            remapped <- vapply(gene_ids, function(gid) {
                if (gid %in% names(label_map)) label_map[[gid]] else gid
            }, character(1L))
            paste(remapped, collapse = "/")
        },
        character(1L)
    )

    if (!is.null(enrich_obj@gene)) {
        enrich_obj@gene <- vapply(enrich_obj@gene, function(gid) {
            if (gid %in% names(label_map)) label_map[[gid]] else gid
        }, character(1L))
    }

    enrich_obj
}

# Get top N region IDs per trait by snp_count, sorted descending.
#
# @param regions_file  path to regions_per_trait.tsv or "NULL"
# @param top_n         integer: max regions per trait (0 = all)
# @return character vector of region_id strings, or NULL (keep all)
top_regions_for_plotting <- function(regions_file, top_n) {
    top_n <- as.integer(top_n)
    if (is.null(regions_file) || regions_file == "NULL" ||
        !file.exists(regions_file) || top_n <= 0L) {
        return(NULL)
    }

    regions <- data.table::fread(regions_file)
    if (nrow(regions) == 0L) return(NULL)

    top_ids <- regions[order(-snp_count, min_pvalue),
                       .SD[seq_len(min(.N, top_n))],
                       by = "trait"]$region_id

    message(paste0("INFO: Top ", top_n, " regions per trait: ", length(top_ids), " total"))
    top_ids
}

# Save a ggplot2 or enrichplot object as PNG + SVG + qs.
#
# A render that fails part-way (e.g. "Viewport has zero dimension(s)") leaves the bytes the device
# wrote before the error, and the app shows that file as a broken plot (file_ok() only checks it
# is non-empty): Trifolium86Main's emapplots left a 42 KB PNG and a 2 KB SVG (9700d6). So a failed
# save removes its file.
#
# @param p           ggplot/enrichResult-based plot
# @param output_base path without extension
# @param width       numeric: inches
# @param height      numeric: inches
save_enrichment_plot <- function(p, output_base, width, height) {
    png_path <- paste0(output_base, ".png")
    tryCatch({
        ggplot2::ggsave(png_path, plot = p,
                        width = width, height = height, dpi = 300)
        message(paste0("INFO: Saved PNG: ", png_path))
    }, error = function(e) {
        unlink(png_path)
        message(paste0("WARNING: PNG save failed, no PNG written: ", e$message))
    })

    svg_path <- paste0(output_base, ".svg")
    tryCatch({
        ggplot2::ggsave(svg_path, plot = p,
                        width = width, height = height)
        message(paste0("INFO: Saved SVG: ", svg_path))
    }, error = function(e) {
        unlink(svg_path)
        message(paste0("WARNING: SVG save failed, no SVG written: ", e$message))
    })

    tryCatch({
        qs::qsave(p, paste0(output_base, ".qs"))
        message(paste0("INFO: Saved .qs: ", output_base, ".qs"))
    }, error = function(e) message(paste0("WARNING: .qs save failed: ", e$message)))
}

# Create and save emapplot (GO term similarity network).
create_emapplot <- function(enrich_obj, n_terms, output_base, width, height) {
    tryCatch({
        enrich_obj <- enrichplot::pairwise_termsim(enrich_obj)
        p <- enrichplot::emapplot(enrich_obj, showCategory = n_terms) +
            ggplot2::ggtitle("GO Term Similarity Network") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
        save_enrichment_plot(p, output_base, width, height)
    }, error = function(e) {
        message(paste0("WARNING: emapplot failed: ", e$message))
    })
}

# Create and save cnetplot (gene-concept network).
create_cnetplot <- function(enrich_obj, n_terms, output_base, width, height,
                             label_map = NULL) {
    tryCatch({
        plot_obj <- remap_enrich_genes(enrich_obj, label_map)
        p <- enrichplot::cnetplot(plot_obj, showCategory = n_terms) +
            ggplot2::ggtitle("Gene-Concept Network") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
        save_enrichment_plot(p, output_base, width, height)
    }, error = function(e) {
        message(paste0("WARNING: cnetplot failed: ", e$message))
    })
}

# Create and save dotplot (GO enrichment dotplot).
create_dotplot <- function(enrich_obj, n_terms, output_base, width, height) {
    tryCatch({
        p <- enrichplot::dotplot(enrich_obj, showCategory = n_terms) +
            ggplot2::ggtitle("GO Enrichment") +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
        save_enrichment_plot(p, output_base, width, height)
    }, error = function(e) {
        message(paste0("WARNING: dotplot failed: ", e$message))
    })
}

# Process all per-region enrichment .qs files in a directory, generating plots.
#
# @param intermediate_dir  base .qs directory (per-trait subdirs inside)
# @param plots_dir         base plots directory
# @param top_terms         integer: max terms per plot
# @param width, height     plot dimensions in inches
# @param label_map         from build_label_map() or NULL
# @param top_region_ids    from top_regions_for_plotting() or NULL (keep all)
# @return integer: total plots created
process_all_enrichment <- function(intermediate_dir, plots_dir, top_terms,
                                   width, height, label_map = NULL,
                                   top_region_ids = NULL) {
    if (!dir.exists(intermediate_dir)) {
        message("WARNING: Enrichment .qs directory not found — skipping plotting")
        return(0L)
    }

    region_files <- list.files(intermediate_dir,
                                pattern   = "^Region_.*_enrichment\\.qs$",
                                full.names = TRUE,
                                recursive  = TRUE)

    if (length(region_files) == 0L) {
        message("WARNING: No per-region enrichment .qs files found")
        return(0L)
    }
    message(paste0("INFO: Found ", length(region_files), " enrichment .qs files"))

    plots_created <- 0L

    for (file in region_files) {
        rel_path  <- sub(paste0(intermediate_dir, "/?"), "", file)
        trait     <- dirname(rel_path)
        filename  <- basename(file)
        region_id <- sub("Region_", "", sub("_enrichment\\.qs$", "", filename))

        if (!is.null(top_region_ids) && !(region_id %in% top_region_ids)) {
            message(paste0("INFO: Skipping region ", region_id, " (not in top list)"))
            next
        }

        enrich_obj <- qs::qread(file)
        if (is.null(enrich_obj) || !inherits(enrich_obj, "enrichResult")) {
            message(paste0("WARNING: Invalid enrichment object in ", filename))
            next
        }

        n_terms  <- nrow(enrich_obj@result)
        n_plot   <- min(top_terms, n_terms)
        base_dir <- file.path(plots_dir, trait)
        dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
        base_name <- file.path(base_dir, paste0("Region_", region_id))

        create_dotplot(enrich_obj, n_plot, paste0(base_name, "_dotplot"), width, height)
        plots_created <- plots_created + 3L  # PNG + SVG + qs

        create_cnetplot(enrich_obj, n_plot, paste0(base_name, "_cnetplot"), width, height,
                        label_map = label_map)
        plots_created <- plots_created + 2L  # PNG + SVG (no qs for cnet)

        if (n_terms >= 2L) {
            create_emapplot(enrich_obj, n_plot, paste0(base_name, "_emapplot"), width, height)
            plots_created <- plots_created + 3L
        }
    }

    plots_created
}
