# gff_parsing.R — GFF3 loading and attribute parsing helpers
# Usage: source("/pipeline/scripts/R/lib/gff_parsing.R")
# Requires: data.table, stringr, tidyr

# Extract gene ID from GFF9 attributes string.
# Prefers Parent= (for exon/CDS/mRNA features); falls back to ID=.
# Strips isoform suffixes: .1, .1_exon_2, etc.
extract_gene_id <- function(attr) {
    id <- stringr::str_extract(attr, "(?<=Parent=)[^;]+")
    id <- ifelse(is.na(id), stringr::str_extract(attr, "(?<=ID=)[^;]+"), id)
    id <- stringr::str_remove(id, "\\.[0-9]+(_exon_[0-9]+)?$")
    id
}

# Remove "key=" prefix from a GFF attribute value string.
clean_attr_value <- function(x) {
    stringr::str_remove(x, "^[^=]+=")
}

# Load a GFF3 file and parse attribute fields into columns.
#
# @param gff_path   path to GFF3 file
# @param feature    GFF feature type to load (e.g., "mRNA", "gene")
# @return data.table with cols:
#   chr, start, end, gene_id, <all parsed attribute fields>
#   chr is character; start/end are integer.
read_gff <- function(gff_path, feature) {
    message(paste0("INFO: Loading GFF feature='", feature, "' from ", gff_path))

    raw <- data.table::fread(
        cmd    = paste("grep -v '#'", shQuote(gff_path)),
        header = FALSE
    )
    raw <- raw[V3 == feature]
    raw <- raw[, .(chr = as.character(V1),
                   start = as.integer(V4),
                   end   = as.integer(V5),
                   strand = as.character(V7),
                   description = V9)]

    if (nrow(raw) == 0L) {
        message(paste0("WARNING: No features of type '", feature, "' in GFF"))
        return(data.table::data.table(chr = character(), start = integer(),
                                       end = integer(), strand = character(),
                                       gene_id = character()))
    }

    raw[, gene_id := extract_gene_id(description)]

    # Infer attribute field names from first row
    sample_desc <- raw$description[1]
    fields <- unique(as.character(
        stringr::str_extract_all(sample_desc, "(?<=^|;)[^=;]+(?==)")[[1]]
    ))
    fields <- fields[nzchar(fields)]

    message(paste0("INFO: GFF attribute fields: ", paste(fields, collapse = ", ")))

    # separate() splits by POSITION, so a clashing name has to be renamed rather
    # than dropped. The core GFF3 field keeps the name 'strand'.
    into <- fields
    if ("strand" %in% into) {
        message("INFO: GFF attribute 'strand' shadows core GFF3 field 7; ",
                "attribute kept as 'strand_attr'")
        into[into == "strand"] <- "strand_attr"
    }

    raw2 <- tidyr::separate(
        raw, col = "description", into = into, sep = ";",
        fill = "right", extra = "drop"
    )

    for (f in into) {
        if (f %in% colnames(raw2)) {
            raw2[[f]] <- clean_attr_value(raw2[[f]])
            raw2[[f]][raw2[[f]] == "NA" | raw2[[f]] == ""] <- NA_character_
        }
    }

    data.table::as.data.table(raw2)
}

# Load exon (or CDS fallback) features from a GFF3 for exon-SNP counting.
#
# @param gff_path  path to GFF3
# @return data.table: chr, start, end, gene_id
read_gff_exons <- function(gff_path) {
    load_feature <- function(feat) {
        raw <- data.table::fread(
            cmd    = paste("grep -v '#'", shQuote(gff_path)),
            header = FALSE
        )
        raw <- raw[V3 == feat,
                   .(chr     = as.character(V1),
                     start   = as.integer(V4),
                     end     = as.integer(V5),
                     gene_id = extract_gene_id(V9))]
        raw
    }

    exons <- load_feature("exon")
    if (nrow(exons) > 0L) {
        message(paste0("INFO: Loaded ", nrow(exons), " exon features"))
        return(exons)
    }
    message("INFO: No exon features found; falling back to CDS")
    cds <- load_feature("CDS")
    message(paste0("INFO: Loaded ", nrow(cds), " CDS features"))
    cds
}
