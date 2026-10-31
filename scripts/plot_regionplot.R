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
library(patchwork)

source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/factor_order.R")
source("/pipeline/scripts/R/lib/region_plot.R")

args = commandArgs(trailingOnly=TRUE)
####################################
REGIONS_FILE = args[1]         # Regions.tsv
GFF_TOPR = args[2]             # topr gene annotation file
ASSOC_TABLES = args[3]         # comma-separated: method:adjust:filepath
TOP_REGIONS = args[4] %>% as.numeric  # Keep only top N regions (0 = all)
GENES_TO_HIGHLIGHT = args[5] %>% str_split('\\|') %>% unlist
                               # gene symbols, 'all', 'nameless', or "NULL" = all (5080b0)
if (identical(GENES_TO_HIGHLIGHT, 'NULL')) GENES_TO_HIGHLIGHT <- 'all'
PLOT_DIR = args[6]
CUSTOM_REGION = args[7]        # "NULL" or chr:start-end
CUSTOM_TRAITS = args[8]        # "NULL" or comma-separated trait names
CUSTOM_METHODS = args[9]       # "NULL" or comma-separated method names (subset of available)
REGIME = if (length(args) >= 10) tolower(args[10]) else "snp"  # "snp" or "wza"
FLANK = if (length(args) >= 11 && args[11] != "NULL") suppressWarnings(as.numeric(args[11])) else 0
                               # bp drawn on either side of the custom region (ad9f15)
####################################
if (length(args) < 9 || length(args) > 11)
  stop('ERROR: plot_regionplot.R takes 9 to 11 arguments, got ', length(args))
if (is.na(FLANK) || FLANK < 0)
  stop('ERROR: FLANK (argument 11) must be a number of bp >= 0, got ', args[11])

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

# One series per trait x method: topr's data frame, its legend label, its method and its
# threshold, computed per (method, trait) under the method's own rule (c34204).
region_series <- function(md, trait, label, method) {
  list(data = md$data %>%
         dplyr::select(SNPID, chr, pos, !!trait) %>%
         setNames(c('SNP', 'CHROM', 'POS', 'P')) %>%
         dplyr::mutate(CHROM = extract_chr_num(CHROM)),
       label = label,
       method = method,
       threshold = series_threshold(md$data[[trait]], md$adjust))
}

# topr's ntop, its default, passed explicitly: series past the 10th are drawn below the axis, and
# every threshold line is then mirrored there as topr mirrored its own.
TOPR_NTOP <- 10

# Draw one region: topr's three panels (overview, main, genes) with every series' threshold line
# (c34204) in one red, a dash pattern and a legend entry per line (F7), the region's edges when
# the window is wider (ad9f15) and gene labels that do not overprint (e9466a), stacked in topr's
# unit ratios. The PNG is that stack; the three panels are also written as SVGs, as before.
draw_region <- function(series, region_str, title, gff_table, filename_base, vlines,
                        png_width, svg_width) {
  series_col <- rep_len(topr::get_topr_colors(), length(series))
  series_lab <- vapply(series, `[[`, character(1), 'label')
  lines <- threshold_lines(vapply(series, `[[`, numeric(1), 'threshold'), series_lab,
                           vapply(series, `[[`, character(1), 'method'))
  if (nrow(lines) == 0)
    message('WARNING: no series has a threshold under its rule; no significance line drawn')
  for (k in seq_len(nrow(lines)))
    message(sprintf('INFO:   Threshold line %d (%s): p = %.4g, %s', k, lines$linetype[k],
                    lines$thresh[k], gsub('\n', ' ', lines$label[k])))
  rp_args <- list(
    build = gff_table,
    region = region_str,
    color = series_col,
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
    legend_labels = series_lab,
    # topr draws no threshold line: add_threshold_lines() draws them below. numeric(0), not NULL,
    # which topr 2.0.4's overview stops on ("object 'i' not found" unless a global i exists), and
    # not topr's default of 5e-08, a line no rule set.
    sign_thresh = numeric(0),
    ntop = TOPR_NTOP,
    title = title
  )
  if (length(vlines) > 0)
    rp_args <- c(rp_args, list(vline = vlines, vline_color = CLINEGO_COL$muted,
                               vline_linetype = 'dotted'))

  plot_list <- do.call(regionplot, c(list(lapply(series, `[[`, 'data')), rp_args,
                                     list(extract_plots = TRUE)))
  # The lines on the main panel, named in the legend, and on the overview, thinner and unnamed,
  # at the widths topr drew its own (0.5, 0.2).
  mirror <- length(series) > TOPR_NTOP
  plot_list$main_plot <- add_threshold_lines(plot_list$main_plot, lines, linewidth = 0.5,
                                             mirror = mirror)
  if (!is.null(plot_list$overview_plot))
    plot_list$overview_plot <- add_threshold_lines(plot_list$overview_plot, lines,
                                                   linewidth = 0.2, legend = FALSE,
                                                   mirror = mirror)
  if (!is.null(plot_list$gene_plot)) {
    plot_list$gene_plot <- tryCatch(repel_gene_labels(plot_list$gene_plot, CLINEGO_COL$fg),
      error = function(e) {
        message('WARNING: gene labels left as topr draws them: ', conditionMessage(e))
        plot_list$gene_plot
      })
  }
  panels <- Filter(Negate(is.null), plot_list[c('overview_plot', 'main_plot', 'gene_plot')])
  units  <- c(overview_plot = rp_args$unit_overview, main_plot = rp_args$unit_main,
              gene_plot = rp_args$unit_gene)
  # One legend column beside the stack, as topr's own figure has it: per-panel legends cut a long
  # series legend at its panel's height.
  # Transparent legend boxes: the column is centred on the whole stack, so a legend can sit
  # beside the x axis, and a white box hid the end of the last tick label (A4 F7, d150 region 12:
  # "6000000" read "600000(").
  fig <- patchwork::wrap_plots(panels, ncol = 1, heights = unname(units[names(panels)]),
                               guides = "collect") &
    ggplot2::theme(legend.background = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank())
  ggsave(paste0(filename_base, '.png'), fig, width = png_width, height = 5, units = 'in',
         dpi = 300)
  message(paste0('INFO:   Saved ', filename_base, '.png'))

  tryCatch({
    if (!is.null(plot_list$main_plot))
      ggsave(paste0(filename_base, '_main.svg'), plot = plot_list$main_plot,
             device = svglite, width = svg_width, height = 6)
    if (!is.null(plot_list$overview_plot))
      ggsave(paste0(filename_base, '_overview.svg'), plot = plot_list$overview_plot,
             device = svglite, width = svg_width, height = 2)
    if (!is.null(plot_list$gene_plot))
      ggsave(paste0(filename_base, '_genes.svg'), plot = plot_list$gene_plot,
             device = svglite, width = svg_width + 2, height = 2)
  }, error = function(e) {
    message(paste0('WARNING: Could not extract SVG plots: ', e$message))
  })
}

# Plot a region with ALL trait x method combos that have sig SNPs.
# One topr data frame per combo, legend label "trait (METHOD)", in numeric bioclim order
# (9e5dd0: the calls list the traits text-sorted, bio_1, bio_10, bio_2, ...).
plot_region_multi_trait <- function(region_str, trait_method_pairs, assoc_list, gff_table,
                                    filename_base, title, vlines = NULL) {
  pair_traits <- vapply(trait_method_pairs, `[[`, character(1), 'trait')
  trait_method_pairs <- trait_method_pairs[numbered_order(pair_traits)]
  series <- list()
  for (pair in trait_method_pairs) {
    for (md_name in pair$methods) {
      md <- assoc_list[[md_name]]
      if (is.null(md) || !pair$trait %in% md$trait_cols) next
      series <- c(series, list(region_series(md, pair$trait,
                                             paste0(pair$trait, ' (', md_name, ')'), md_name)))
    }
  }

  if (length(series) == 0) {
    message('WARNING: No data for any trait/method pair, skipping')
    return(invisible(NULL))
  }
  draw_region(series, region_str, title, gff_table, filename_base, vlines,
              png_width = 14, svg_width = 12)
}

# Plot a single region for a single trait across multiple methods
# Returns nothing, saves PNG + SVG files
plot_region_trait <- function(region_str, trait, assoc_list, gff_table, filename_base,
                              vlines = NULL) {
  series <- list()
  for (md_name in names(assoc_list)) {
    md <- assoc_list[[md_name]]
    if (!trait %in% md$trait_cols) next
    series <- c(series, list(region_series(md, trait, md_name, md_name)))
  }

  if (length(series) == 0) {
    message(paste0('WARNING: No method data for trait ', trait, ', skipping'))
    return(invisible(NULL))
  }
  draw_region(series, region_str, trait, gff_table, filename_base, vlines,
              png_width = 12, svg_width = 10)
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

  list(method = md$method, adjust = md$adjust, data = assoc, trait_cols = trait_cols)
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

    # Padding, the region's edges marked inside it (ad9f15)
    region_span <- region_end - region_start
    padding <- max(region_span * 0.2, 500000)
    window <- region_window(region_start, region_end, padding)

    chr_num <- extract_chr_num(region_chr)
    region_str <- paste0(chr_num, ':', window[1], '-', window[2])
    region_safe <- gsub("[:\\-]", "_", region_id)
    edges <- region_edge_lines(region_start, region_end, window[1], window[2])

    # One plot per trait
    for (trait in region_traits) {
      message(paste0('INFO:   Trait: ', trait))
      filename_base <- file.path(PLOT_DIR, paste0('regionplot_', region_safe))
      plot_region_trait(region_str, trait, assoc_list, gff_table, filename_base, edges)
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

  # The window is the region plus FLANK bp each side; the region's edges are marked when it is
  # wider (ad9f15). Output names stay keyed by the region, which is what the app collects.
  chr_num <- extract_chr_num(custom_chr)
  region_label <- paste0(chr_num, ':', custom_start, '-', custom_end)
  window <- region_window(custom_start, custom_end, FLANK)
  region_str <- paste0(chr_num, ':', window[1], '-', window[2])
  edges <- region_edge_lines(custom_start, custom_end, window[1], window[2])
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
    filename_base <- file.path(PLOT_DIR, paste0('regionplot_custom_', region_safe))
    plot_region_multi_trait(region_str, trait_method_pairs, custom_assoc_list, gff_table,
                            filename_base, region_label, edges)
  } else {
    # Legacy single-trait path
    message(paste0('INFO:   Traits: ', paste(custom_trait_list, collapse = ', ')))
    message(paste0('INFO:   Methods: ', paste(names(custom_assoc_list), collapse = ', ')))
    for (trait in custom_trait_list) {
      message(paste0('INFO:   Custom trait: ', trait))
      filename_base <- file.path(PLOT_DIR, paste0('regionplot_custom_', region_safe, '_', trait))
      plot_region_trait(region_str, trait, custom_assoc_list, gff_table, filename_base, edges)
    }
  }
}

message('INFO: Region plots complete')
