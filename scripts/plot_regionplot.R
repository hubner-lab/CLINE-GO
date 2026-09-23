#!/usr/bin/env Rscript
# Regional Manhattan plots for top regions
# One plot per trait: overlays all methods (EMMAX, LFMM, etc.) for the same trait
# Optionally plots a custom user-specified region

library(topr)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(svglite)

args = commandArgs(trailingOnly=TRUE)
####################################
REGIONS_FILE = args[1]         # Regions.tsv
GFF_TOPR = args[2]             # topr gene annotation file
ASSOC_TABLES = args[3]         # comma-separated: method:adjust:filepath
TOP_REGIONS = args[4] %>% as.numeric  # Keep only top N regions (0 = all)
GENES_TO_HIGHLIGHT = args[5] %>% str_split('\\|') %>% unlist
PLOT_DIR = args[6]
CUSTOM_REGION = args[7]        # "NULL" or chr:start-end
CUSTOM_TRAITS = args[8]        # "NULL" or comma-separated trait names
CUSTOM_METHODS = args[9]       # "NULL" or comma-separated method names (subset of available)
REGIME = if (length(args) >= 10) tolower(args[10]) else "snp"  # "snp" or "wza"
####################################

is_wza <- REGIME == "wza"
if (is_wza) message('INFO: WZA regime — loading window-level p-values')
message('INFO: Regional Manhattan plot generation')

######################## Functions

# Extract numeric chromosome identifier from names like chr3H, Chr1, chr1H_1, 1
extract_chr_num <- function(chrom_names) {
  sapply(chrom_names, function(x) {
    digits <- paste0(unlist(regmatches(x, gregexpr("[0-9]+", x))), collapse = "")
    as.numeric(digits)
  }, USE.NAMES = FALSE)
}

# Calculate significance threshold from adjust string and association data
calc_threshold <- function(assoc, adjust_str) {
  method <- str_extract(adjust_str, '(.*)_', group = 1)
  signif <- str_extract(adjust_str, '_(.*)', group = 1) %>% as.numeric

  if (method == 'bonf') {
    return(signif / nrow(assoc))
  } else if (method == 'qval') {
    if (!requireNamespace("qvalue", quietly = TRUE)) {
      message('WARNING: qvalue package not available, using Bonferroni')
      return(signif / nrow(assoc))
    }
    qvals <- qvalue::qvalue(assoc$pvalue)$qvalues
    return(max(assoc$pvalue[qvals <= signif], na.rm = TRUE))
  } else if (method == 'custom') {
    return(signif)
  } else if (method == 'top') {
    return(sort(assoc$pvalue)[min(as.integer(signif), nrow(assoc))])
  } else {
    message(paste0('WARNING: Unknown threshold method: ', method, ', using Bonferroni'))
    return(signif / nrow(assoc))
  }
}

# Parse CUSTOM_TRAITS structured format: "bio_12:BLINK|EMMAX|LFMM,bio_7:BLINK|LFMM"
# Returns list of list(trait, methods) or NULL if plain single-trait format
parse_trait_method_pairs <- function(custom_traits_str) {
  if (!grepl(':', custom_traits_str)) return(NULL)
  parts <- strsplit(custom_traits_str, ',')[[1]]
  lapply(parts, function(p) {
    sp <- strsplit(p, ':')[[1]]
    list(trait = sp[1], methods = strsplit(sp[2], '\\|')[[1]])
  })
}

# Plot a region with ALL trait x method combos that have sig SNPs.
# One topr data frame per combo, legend label "trait (METHOD)".
plot_region_multi_trait <- function(region_str, trait_method_pairs, assoc_list, gff_table, filename_base) {
  topr_data       <- list()
  topr_labels     <- c()
  topr_thresholds <- c()

  for (pair in trait_method_pairs) {
    trait   <- pair$trait
    methods <- pair$methods
    for (md_name in methods) {
      md <- assoc_list[[md_name]]
      if (is.null(md) || !trait %in% md$trait_cols) next
      trait_data <- md$data %>%
        dplyr::select(SNPID, chr, pos, !!trait) %>%
        setNames(c('SNP', 'CHROM', 'POS', 'P')) %>%
        dplyr::mutate(CHROM = extract_chr_num(CHROM))
      topr_data       <- c(topr_data, list(trait_data))
      topr_labels     <- c(topr_labels, paste0(trait, ' (', md_name, ')'))
      topr_thresholds <- c(topr_thresholds, md$threshold)
    }
  }

  if (length(topr_data) == 0) {
    message('WARNING: No data for any trait/method pair, skipping')
    return(invisible(NULL))
  }

  sign_thresh <- min(topr_thresholds)
  rp_args <- list(
    build = gff_table,
    region = region_str,
    sign_thresh = sign_thresh,
    max.overlaps = 999,
    show_gene_names = TRUE,
    # show_genes left at topr's default (NULL = auto): exon structure below a
    # 1 Mb window, gene boxes above, where exons turn into unreadable noise.
    # TRUE here hid exon structure on every region, whatever the GFF carried.
    unit_overview = 1,
    unit_main = 2.5,
    unit_gene = 1,
    axis_text_size = 14,
    axis_title_size = 16,
    title_text_size = 16,
    legend_position = 'right',
    scale = 1,
    segment.size = 3,
    legend_labels = topr_labels,
    sign_thresh_label_size = 0,
    title = region_str
  )

  png(paste0(filename_base, '.png'), width = 14, height = 5, units = 'in', res = 300)
  do.call(regionplot, c(list(topr_data), rp_args, list(extract_plots = FALSE)))
  dev.off()
  message(paste0('INFO:   Saved multi-trait ', filename_base, '.png'))

  tryCatch({
    plot_list <- do.call(regionplot, c(list(topr_data), rp_args, list(extract_plots = TRUE)))
    if (!is.null(plot_list$main_plot))
      ggsave(paste0(filename_base, '_main.svg'), plot = plot_list$main_plot,
             device = svglite, width = 12, height = 6)
    if (!is.null(plot_list$overview_plot))
      ggsave(paste0(filename_base, '_overview.svg'), plot = plot_list$overview_plot,
             device = svglite, width = 12, height = 2)
    if (!is.null(plot_list$gene_plot))
      ggsave(paste0(filename_base, '_genes.svg'), plot = plot_list$gene_plot,
             device = svglite, width = 14, height = 2)
  }, error = function(e) {
    message(paste0('WARNING: Could not extract SVG plots: ', e$message))
  })
}

# Plot a single region for a single trait across multiple methods
# Returns nothing, saves PNG + SVG files
plot_region_trait <- function(region_str, trait, assoc_list, gff_table, filename_base) {
  topr_data <- list()
  topr_labels <- c()
  topr_thresholds <- c()

  for (md_name in names(assoc_list)) {
    md <- assoc_list[[md_name]]
    if (!trait %in% md$trait_cols) next

    trait_data <- md$data %>%
      dplyr::select(SNPID, chr, pos, !!trait) %>%
      setNames(c('SNP', 'CHROM', 'POS', 'P')) %>%
      dplyr::mutate(CHROM = extract_chr_num(CHROM))

    topr_data <- c(topr_data, list(trait_data))
    topr_labels <- c(topr_labels, md_name)
    topr_thresholds <- c(topr_thresholds, md$threshold)
  }

  if (length(topr_data) == 0) {
    message(paste0('WARNING: No method data for trait ', trait, ', skipping'))
    return(invisible(NULL))
  }

  sign_thresh <- min(topr_thresholds)

  # Common regionplot arguments
  rp_args <- list(
    build = gff_table,
    region = region_str,
    sign_thresh = sign_thresh,
    max.overlaps = 999,
    show_gene_names = TRUE,
    # show_genes left at topr's default (NULL = auto): exon structure below a
    # 1 Mb window, gene boxes above, where exons turn into unreadable noise.
    # TRUE here hid exon structure on every region, whatever the GFF carried.
    unit_overview = 1,
    unit_main = 2.5,
    unit_gene = 1,
    axis_text_size = 14,
    axis_title_size = 16,
    title_text_size = 16,
    legend_position = 'right',
    scale = 1,
    segment.size = 3,
    legend_labels = topr_labels,
    sign_thresh_label_size = 0,
    title = trait
  )

  # PNG
  png(paste0(filename_base, '.png'), width = 12, height = 5, units = 'in', res = 300)
  do.call(regionplot, c(list(topr_data), rp_args, list(extract_plots = FALSE)))
  dev.off()
  message(paste0('INFO:   Saved ', filename_base, '.png'))

  # SVG extracted components
  tryCatch({
    plot_list <- do.call(regionplot, c(list(topr_data), rp_args, list(extract_plots = TRUE)))

    if (!is.null(plot_list$main_plot)) {
      ggsave(paste0(filename_base, '_main.svg'), plot = plot_list$main_plot,
             device = svglite, width = 10, height = 6)
    }
    if (!is.null(plot_list$overview_plot)) {
      ggsave(paste0(filename_base, '_overview.svg'), plot = plot_list$overview_plot,
             device = svglite, width = 10, height = 2)
    }
    if (!is.null(plot_list$gene_plot)) {
      ggsave(paste0(filename_base, '_genes.svg'), plot = plot_list$gene_plot,
             device = svglite, width = 12, height = 2)
    }
  }, error = function(e) {
    message(paste0('WARNING: Could not extract SVG plots: ', e$message))
  })
}

######################## Main

# Parse association table specs: "METHOD:ADJUST:filepath,METHOD:ADJUST:filepath"
assoc_specs <- str_split(ASSOC_TABLES, ',') %>% unlist
method_data <- lapply(assoc_specs, function(spec) {
  parts <- str_split(spec, ':')[[1]]
  list(method = parts[1], adjust = parts[2], file = parts[3])
})

message(paste0('INFO: Methods: ', paste(sapply(method_data, `[[`, 'method'), collapse = ', ')))

# Load all association tables
wza_meta_cols <- c('n_snps', 'mean_maf')
assoc_list <- lapply(method_data, function(md) {
  message(paste0('INFO: Loading ', md$method, ' from ', md$file))
  assoc <- fread(md$file)
  assoc$chr <- as.character(assoc$chr)
  exclude <- if (is_wza) c('SNPID', 'chr', 'pos', wza_meta_cols) else c('SNPID', 'chr', 'pos')
  trait_cols <- setdiff(colnames(assoc), exclude)

  threshold <- calc_threshold(
    data.frame(pvalue = assoc[[trait_cols[1]]]),
    md$adjust
  )

  list(method = md$method, adjust = md$adjust, data = assoc,
       trait_cols = trait_cols, threshold = threshold)
})
names(assoc_list) <- sapply(assoc_list, `[[`, 'method')

# Load GFF annotation
# Exon columns are comma-joined lists that topr strsplit()s. With one exon per
# gene there is no comma and fread infers integer (all empty -> logical), which
# kills the exon-structure plot with "non-character argument".
gff_table <- fread(GFF_TOPR, colClasses = c(exon_chromstart = "character",
                                            exon_chromend   = "character")) %>%
  dplyr::mutate(chrom = extract_chr_num(chrom))

# Filter genes based on GENES_TO_HIGHLIGHT
if (!any(GENES_TO_HIGHLIGHT %in% c('nameless', 'all'))) {
  gff_table <- gff_table %>%
    dplyr::filter(gene_symbol %in% GENES_TO_HIGHLIGHT)
}
if (any(GENES_TO_HIGHLIGHT == 'nameless')) {
  gff_table <- gff_table %>%
    dplyr::mutate(gene_symbol = '')
}

if (nrow(gff_table) == 0) {
  message('WARNING: No genes found in GFF annotation')
}

# ===================== Part 1: Auto regions from Regions.tsv =====================

if (REGIONS_FILE != "NULL" && REGIONS_FILE != "" && file.exists(REGIONS_FILE)) {
regions <- fread(REGIONS_FILE)

if (nrow(regions) == 0) {
  message('WARNING: No regions in Regions.tsv')
} else {
  # Sort by snp_count descending, then by min_pvalue ascending
  regions <- regions %>%
    dplyr::arrange(desc(snp_count), min_pvalue)

  # Filter to top N regions
  if (TOP_REGIONS > 0 && nrow(regions) > TOP_REGIONS) {
    message(paste0('INFO: Filtering to top ', TOP_REGIONS, ' regions (from ', nrow(regions), ')'))
    regions <- regions %>% dplyr::slice_head(n = TOP_REGIONS)
  }

  message(paste0('INFO: Plotting ', nrow(regions), ' regions'))

  for (i in seq_len(nrow(regions))) {
    region <- regions[i, ]
    region_id <- region$region_id
    region_chr <- as.character(region$chr)
    region_start <- region$start
    region_end <- region$end
    region_traits <- if ('traits' %in% colnames(regions)) {
      str_split(region$traits, ',')[[1]]
    } else {
      region$trait
    }

    message(paste0('INFO: Region ', i, '/', nrow(regions), ': ', region_id))

    # Padding
    region_span <- region_end - region_start
    padding <- max(region_span * 0.2, 500000)
    plot_start <- max(1, region_start - padding)
    plot_end <- region_end + padding

    chr_num <- extract_chr_num(region_chr)
    region_str <- paste0(chr_num, ':', plot_start, '-', plot_end)
    region_safe <- gsub("[:\\-]", "_", region_id)

    # One plot per trait
    for (trait in region_traits) {
      message(paste0('INFO:   Trait: ', trait))
      filename_base <- paste0(PLOT_DIR, 'regionplot_', region_safe)
      plot_region_trait(region_str, trait, assoc_list, gff_table, filename_base)
    }
  }
} # end Part 1 regions loop

} # end if REGIONS_FILE exists

# ===================== Part 2: Custom region (if specified) =====================

if (!is.null(CUSTOM_REGION) && CUSTOM_REGION != 'NULL' && CUSTOM_REGION != '') {
  message(paste0('INFO: Plotting custom region: ', CUSTOM_REGION))

  # Parse custom region
  custom_chr <- str_extract(CUSTOM_REGION, '(.*):', group = 1)
  custom_startend <- str_extract(CUSTOM_REGION, ':(.*)', group = 1)
  custom_start <- str_extract(custom_startend, '(.*)\\-', group = 1) %>% as.numeric
  custom_end <- str_extract(custom_startend, '\\-(.*)', group = 1) %>% as.numeric

  if (is.na(custom_start) || is.na(custom_end) || custom_start > custom_end) {
    stop(paste0('ERROR: Invalid custom region format: ', CUSTOM_REGION, '. Expected chr:start-end'))
  }

  chr_num <- extract_chr_num(custom_chr)
  region_str <- paste0(chr_num, ':', custom_start, '-', custom_end)
  region_safe <- gsub("[:\\-]", "_", CUSTOM_REGION)

  # Determine traits
  if (!is.null(CUSTOM_TRAITS) && CUSTOM_TRAITS != 'NULL' && CUSTOM_TRAITS != '') {
    custom_trait_list <- str_split(CUSTOM_TRAITS, ',')[[1]]
  } else {
    # Use all available traits from the first method
    custom_trait_list <- assoc_list[[1]]$trait_cols
  }

  # Filter methods if custom methods specified
  custom_assoc_list <- assoc_list
  if (!is.null(CUSTOM_METHODS) && CUSTOM_METHODS != 'NULL' && CUSTOM_METHODS != '') {
    custom_method_names <- str_split(CUSTOM_METHODS, ',')[[1]]
    custom_assoc_list <- assoc_list[intersect(names(assoc_list), custom_method_names)]
    if (length(custom_assoc_list) == 0) {
      stop(paste0('ERROR: None of custom methods found: ', CUSTOM_METHODS,
                   '. Available: ', paste(names(assoc_list), collapse = ', ')))
    }
  }

  # Detect multi-trait structured format ("bio_12:BLINK|EMMAX,bio_7:LFMM")
  trait_method_pairs <- parse_trait_method_pairs(CUSTOM_TRAITS)

  if (!is.null(trait_method_pairs)) {
    message(paste0('INFO:   Multi-trait mode: ',
      paste(sapply(trait_method_pairs, function(p) paste0(p$trait, '(', paste(p$methods, collapse=','), ')')),
            collapse = ', ')))
    filename_base <- paste0(PLOT_DIR, 'regionplot_custom_', region_safe)
    plot_region_multi_trait(region_str, trait_method_pairs, custom_assoc_list, gff_table, filename_base)
  } else {
    # Legacy single-trait path
    message(paste0('INFO:   Traits: ', paste(custom_trait_list, collapse = ', ')))
    message(paste0('INFO:   Methods: ', paste(names(custom_assoc_list), collapse = ', ')))
    for (trait in custom_trait_list) {
      message(paste0('INFO:   Custom trait: ', trait))
      filename_base <- paste0(PLOT_DIR, 'regionplot_custom_', region_safe, '_', trait)
      plot_region_trait(region_str, trait, custom_assoc_list, gff_table, filename_base)
    }
  }
}

message('INFO: Region plots complete')
