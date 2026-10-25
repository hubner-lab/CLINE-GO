#!/usr/bin/env Rscript
# DERIVED COPY (second generation) for the Trifolium paper figures (plan ~/.claude/plans/fancy-scribbling-garden.md,
# Phase 3): regions/figures/plot_regionplot_thresholds.R, itself the pipeline's scripts/plot_regionplot.R plus every
# method threshold + region-edge vlines (its own header below). User 2026-10-16: "we used topr in the pipeline to
# draw the regions so I'd like to work with that. One figure per region". Changes here, all tagged `# DERIVED (Dn)`
# and listed in manuscript_fig3/pipeline_deviations.md:
#   D7 topr's own size arguments set for a 180 mm figure at 7/6 pt (axis/title/legend text, gene labels, points,
#      threshold and vline widths, label segments); the three topr panels (overview, main, genes) are taken with
#      extract_plots = TRUE and stacked with patchwork in topr's unit ratios, then written by ms_save() as one editable
#      SVG + 600 dpi PNG at SIZE_MM (arg 12); fonts checked there. No separate _main/_overview/_genes SVGs.
#   D1 minou colours: EMMAX series = the predictor colour (fig3_style.R FIG3_PRED_COL, all 16 predictors since
#      2026-10-17), RDA = sage; thresholds EMMAX red / RDA sage (FIG3_THR) instead of darkred /
#      steelblue4; region-edge vlines minou grey; one gene colour FIG3_GENE (topr gene_color).
#   D17 (2026-10-18) two-entry legend — WITHDRAWN the same day (user: "return to simple legend one color per track");
#      replaced by D21.
#   D21 (2026-10-18) topr's own legend, one key per series, with display labels: "BIO1 (EMMAX)" ... and
#      "multivariate (RDA)" instead of "bio_1 (EMMAX)" / "climate_multivariate (RDA)"; the "trait (METHOD)" labels still
#      pick each series' colour.
# Check with `diff regions/figures/plot_regionplot_thresholds.R <this file>`.
# DERIVED COPY of scripts/plot_regionplot.R (pipeline, 2026-10-13 state) for the Trifolium86Main deep-dive figures,
# kept in ANALYSIS/ because the pipeline is USE-ONLY here. ONE change, in plot_region_multi_trait():
#   the pipeline draws only min(topr_thresholds), so RDA hits (q 0.05, p <= 6.9e-5) sat under the Bonferroni line;
#   here every distinct method threshold is drawn (topr::regionplot accepts a vector), the lowest in darkred
#   (EMMAX / LFMM Bonferroni) and the others in steelblue4 (RDA q-value). Everything else is byte-identical:
#   check with `diff scripts/plot_regionplot.R <this file>`. Filed in docs/pipeline_improvement_requests.md.
# SECOND CHANGE (2026-10-15, ANALYSIS/Trifolium86MainK5/regionsER region-QC report): optional 11th argument VLINES,
#   comma-separated positions (the region edges) drawn as dotted grey vertical lines via topr::regionplot(vline=).
# Regional Manhattan plots for top regions
# One plot per trait: overlays all methods (EMMAX, LFMM, etc.) for the same trait
# Optionally plots a custom user-specified region

library(topr)
library(dplyr)
library(data.table)
library(stringr)
library(ggplot2)
library(svglite)
library(patchwork)  # DERIVED (D7)
source('/pipeline/ANALYSIS/Trifolium86MainK5ERpc5/fig3/fig3_style.R')  # DERIVED (D1, D7)

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
VLINES = if (length(args) >= 11 && args[11] != "NULL") as.numeric(str_split(args[11], ",")[[1]]) else NULL  # DERIVED
SIZE_MM = if (length(args) >= 12) as.numeric(str_split(args[12], "x")[[1]]) else c(180, 79)  # DERIVED (D7)
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

  # D8 withdrawn 2026-10-17 (user: "draw both, simpler than explaining"): bio_11 and bio_19 are drawn as their own
  # series in their own colours, although they equal bio_8 / bio_16 at every plant (exact overplot; caption says so)
  # DERIVED (D7, 2026-10-17): series (and legend) in numeric bioclim order, the multivariate RDA series last; the
  # calls list them as text-sorted (bio_1, bio_10, bio_2, ...)
  ord_key <- vapply(trait_method_pairs, function(p) {
    n <- suppressWarnings(as.numeric(sub('^bio_', '', p$trait))); if (is.na(n)) Inf else n }, 0)
  trait_method_pairs <- trait_method_pairs[order(ord_key)]
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

  sign_thresh <- sort(unique(topr_thresholds))   # DERIVED: all method thresholds, not only the lowest
  sign_thresh_color <- ifelse(sign_thresh == min(sign_thresh), FIG3_THR[['EMMAX']], FIG3_THR[['RDA']])  # DERIVED (D1)
  series_trait <- sub(' .*', '', topr_labels)   # first name of the series label                                                     # DERIVED (D1)
  series_col <- ifelse(grepl('(RDA)', topr_labels, fixed = TRUE), FIG3_THR[['RDA']],                  # DERIVED (D1)
                       unname(FIG3_PRED_COL[series_trait]))                                              # DERIVED (D1)
  stopifnot(!anyNA(series_col))                                                                          # DERIVED (D1)
  series_lab <- ifelse(grepl('(RDA)', topr_labels, fixed = TRUE), 'multivariate (RDA)',               # DERIVED (D21)
                       paste0(sub('^bio_', 'BIO', series_trait), sub('^[^ ]+', '', topr_labels)))         # DERIVED (D21)
  stopifnot(!any(grepl('bio_|climate_', series_lab)))                                                    # DERIVED (D21)
  message('INFO:   legend: ', paste(series_lab, collapse = ' | '))                                       # DERIVED (D21)
  rp_args <- list(
    build = gff_table,
    region = region_str,
    sign_thresh = sign_thresh,
    sign_thresh_color = sign_thresh_color,
    max.overlaps = 999,
    show_gene_names = TRUE,
    # show_genes left at topr's default (NULL = auto): exon structure below a
    # 1 Mb window, gene boxes above, where exons turn into unreadable noise.
    # TRUE here hid exon structure on every region, whatever the GFF carried.
    unit_overview = 1,
    unit_main = 2.5,
    unit_gene = 1,
    axis_text_size = MS_TICK,                 # DERIVED (D7) was 14
    axis_title_size = MS_BASE,                # DERIVED (D7) was 16
    title_text_size = MS_BASE,                # DERIVED (D7) was 16
    legend_text_size = MS_TICK,               # DERIVED (D7) topr default 11
    legend_title_size = MS_BASE,              # DERIVED (D7) topr default 12
    gene_label_size = MS_LAB,                 # DERIVED (D7) 6 pt gene names
    size = 0.7,                               # DERIVED (D7) topr default 2
    sign_thresh_size = 0.3,                   # DERIVED (D7) topr default 0.5
    color = series_col,                       # DERIVED (D1)
    gene_color = FIG3_GENE,                   # DERIVED (D1)
    legend_position = 'right',
    scale = 1,
    segment.size = 0.2,                       # DERIVED (D7) was 3
    legend_labels = series_lab,               # DERIVED (D21) display labels
    sign_thresh_label_size = 0,
    title = region_str
  )
  if (!is.null(VLINES)) rp_args <- c(rp_args, list(vline = VLINES, vline_color = MINOU[['grey']], vline_linetype = 'dotted',  # DERIVED (D1)
                                                   vline_size = 0.3))  # DERIVED (D7)

  # DERIVED (D7): one figure from topr's three panels, at the final mm size
  plot_list <- do.call(regionplot, c(list(topr_data), rp_args, list(extract_plots = TRUE)))
  # topr's y title is a plotmath expression, whose minus and brackets need the Symbol font (not in the journal font
  # set, ms_save() refuses it); same words as plain text, as in the paper's other Manhattan plots
  plot_list$main_plot <- plot_list$main_plot + labs(y = '-log10 P')
  # DERIVED (D7, 2026-10-17): topr writes gene names with plain geom_text at gene_start, so two short tandem copies
  # (RCC1 x 2 in 5_5557888-5883272) printed over each other; the same labels redrawn with ggrepel, kept on their row
  if (!is.null(plot_list$gene_plot)) {
    gp <- plot_list$gene_plot
    is_txt <- vapply(gp$layers, function(l) inherits(l$geom, 'GeomText'), TRUE)
    for (i in which(is_txt)) {
      l <- gp$layers[[i]]
      lab_data <- l$data[nzchar(l$data$gene_symbol), , drop = FALSE]   # unlabelled genes must not repel the labels
      # label starts just right of its gene box (topr starts it on the box), thin white halo over neighbouring boxes
      xr <- diff(range(c(l$data$gene_start, l$data$gene_end)))
      lab_data$gene_start <- lab_data$gene_end + 0.004 * xr
      gp$layers[[i]] <- ggrepel::geom_text_repel(data = lab_data, mapping = l$mapping, inherit.aes = FALSE,
          size = l$aes_params$size %||% MS_LAB, hjust = 0, vjust = 0, colour = MS_INK, direction = 'y',
          bg.colour = 'white', bg.r = 0.12,
          box.padding = 0.05, point.padding = 0, min.segment.length = Inf, max.overlaps = Inf, seed = 1)
    }
    plot_list$gene_plot <- gp
  }
  panels <- Filter(Negate(is.null), plot_list[c('overview_plot', 'main_plot', 'gene_plot')])
  units  <- c(overview_plot = rp_args$unit_overview, main_plot = rp_args$unit_main, gene_plot = rp_args$unit_gene)
  fig <- wrap_plots(panels, ncol = 1, heights = unname(units[names(panels)]))
  ms_save(filename_base, fig, SIZE_MM[1], SIZE_MM[2])
  message(paste0('INFO:   Saved multi-trait ', filename_base, '.svg/.png (', paste(names(panels), collapse = ', '), ')'))
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
