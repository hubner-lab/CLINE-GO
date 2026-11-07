#!/usr/bin/env Rscript
# region_check.R: the numbers the region stop turns on, beside score_regions.R's fixed table.
#
#   Rscript manuscript/arabidopsis/region_check.R <project> [n_draws] [seed]
#
# score_regions.R holds the scoring rule fixed before any region existed, so it stays as committed
# and its table stays the primary result. Everything here was added after that table was read, and
# is reported as such. It asks four things that table cannot:
#
#   Power.     Each scan's regions are padded by the clumping distance, so a random region can touch
#              a listed gene most of the time. Arms: the rule (bonf 0.05, qval 0.05) by the distance
#              (the run's own, 10 kb). Agreement sets are rebuilt at each distance by the pipeline's
#              own selection, because re-clustering a stored set keeps its old membership.
#   Recall.    score_regions.R reports the share of regions touching the list. Listed items cluster
#              (two to three GrENE-net blocks can sit under one region), so the share of list items
#              touched is reported too, over the same kind of random placement.
#   Density.   Regions per gene: a region set in gene-rich arms beats a null that also lands in
#              pericentromeres whatever the list holds.
#   Classics.  in_region says little when every SNP is padded, so the nearest called SNP is given.
#
# Also checks the union set against regions_combined.tsv row for row, scores WZA's significant
# windows, and, when rda_axes_sweep.R full has run, each RDA axis count.
#
# Read the qval rows with care. Random regions may overlap each other and real ones cannot, so the
# null covers less genome than the scan does: where regions cover most of the genome the null recall
# falls far below the scan's, and recall_p is an artefact there. Where they cover a small share the gap
# is small, but it runs the same way.
#
# Anchors, each a hard stop: bonf at the run's distance must reproduce the stored sets, the GWAS
# per-trait regions and every row of score_regions.R's scores.tsv.
#
# Reads <project>_results/, data/arabidopsis1135/lists/, ANALYSIS/Arabidopsis1135/regions/<project>/
# and, if present, ANALYSIS/<project>/pregea/rda_axes_full_hits.tsv.
# Writes ANALYSIS/Arabidopsis1135/regions/<project>/check/:
#   sets.tsv             arm scan n_snps n_regions span_mb genome_share median_region_kb n_genes
#                        share_of_union_snps share_of_union_regions
#   scores.tsv           arm scan list n_items n_regions n_with_item share null_mean excess empirical_p
#                        items_touched recall null_recall recall_excess recall_p
#   gene_density.tsv     arm scan n_regions genes_per_region null_genes_per_region ratio empirical_p
#   classic_nearest.tsv  arm scan gene agi in_region nearest_snp distance_bp min_pvalue
#                        (for WZA, nearest_snp holds the nearest window's id)
#   fdr_cells.tsv        arm pair type value threshold status n_sig
#
# Run from the repo root, in the image (the pipeline libraries are baked at /pipeline).

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(stringr)
  library(yaml)
})

source("/pipeline/scripts/R/lib/region_distance.R")
source("/pipeline/scripts/R/lib/regions.R")
source("/pipeline/scripts/R/utils/pval_threshold.R")
source("/pipeline/scripts/R/lib/snp_set_select.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("usage: region_check.R <project> [n_draws] [seed]")
project <- args[1]
n_draws <- if (length(args) >= 2L) as.integer(args[2]) else 10000L
seed <- if (length(args) >= 3L) as.integer(args[3]) else 42L

# score_scan() and its helpers from score_regions.R. Only its top-level function definitions and list
# names are evaluated: the file's body runs the scoring.
list_names <- c("gea_lists", "gwas_lists")
for (statement in parse(here("manuscript", "arabidopsis", "score_regions.R"))) {
  if (!is.call(statement) || !identical(statement[[1]], as.name("<-"))) next
  assigned <- statement[[3]]
  is_function <- is.call(assigned) && identical(assigned[[1]], as.name("function"))
  if (is.name(statement[[2]]) && (is_function || as.character(statement[[2]]) %in% list_names)) eval(statement)
}

results <- here(paste0(project, "_results"))
list_dir <- here("data", "arabidopsis1135", "lists")
region_dir <- here("ANALYSIS", "Arabidopsis1135", "regions", project)
out_dir <- file.path(region_dir, "check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

config <- yaml::read_yaml(here(paste0("config_", project, ".yaml")))
clumping <- resolve_clumping_distance_strict(
  config$GEA$snp_clumping_distance,
  file.path(results, "Structure", "tables", "ld_decay_half_distances.tsv"),
  config$GEA$clumping_ld_rule, "All"
)
k_best <- config$sNMF$k_best

sequence_regions <- fread(cmd = sprintf("awk '/^##sequence-region/ {print $2, $4}' %s",
                                        shQuote(here("data", "Arabidopsis_1001G.gff3"))),
                          header = FALSE, col.names = c("chr", "length"),
                          colClasses = c(V1 = "character"))
chrom_lengths <- setNames(as.numeric(sequence_regions$length), sequence_regions$chr)
genome_bp <- sum(chrom_lengths[as.character(1:5)])

# Lists -----------------------------------------------------------------------------------------
lists <- lapply(setNames(union(gea_lists, gwas_lists), union(gea_lists, gwas_lists)), read_list)
aragwas <- lists$aragwas_ft[above_permutation == TRUE]
lists$aragwas_ft <- NULL
classic <- read_list("classic_genes")

# Study 261 and 262 are the 1001 Genomes FT10 and FT16 values this project uses, on 1,163 and 1,123
# accessions. Studies 30 and 31 are the same traits on a separate panel of about 194.
aragwas_lists <- function(trait) {
  list(aragwas_ft = aragwas[label == trait],
       aragwas_1001g = aragwas[label == trait & study %in% c(261, 262)],
       aragwas_panel194 = aragwas[label == trait & study %in% c(30, 31)])
}

gene_table <- fread(cmd = sprintf("awk -F'\\t' '$3==\"%s\" {print $1, $4, $5}' %s", config$GFF$feature,
                                  shQuote(file.path(results, "_intermediate", "annotation", "normalized.gff3"))),
                    header = FALSE, col.names = c("chr", "start", "end"),
                    colClasses = c(V1 = "character"))
gene_table <- gene_table[chr %in% as.character(1:5)]

# Random placement, shared by recall and density ---------------------------------------------------
place_regions <- function(regions) {
  set.seed(seed)
  lapply(split(regions, by = "chr"), function(on_chr) {
    chr_length <- chrom_lengths[[on_chr$chr[1]]]
    widths <- pmin(on_chr$end, chr_length) - on_chr$start
    room <- pmax(chr_length - widths, 1)
    list(chr = on_chr$chr[1], widths = widths, from = on_chr$start,
         random_start = matrix(1 + floor(runif(nrow(on_chr) * n_draws) * room), nrow = nrow(on_chr)))
  })
}

# Items touched by the regions and by each random draw. Overlaps are counted from the sorted ends, so
# overlapping random regions need no merging.
recall_scan <- function(placed, items) {
  observed <- 0
  null_touched <- numeric(n_draws)
  for (part in placed) {
    on_chr <- items[chr == part$chr]
    if (nrow(on_chr) == 0L) next
    observed <- observed + sum(findInterval(on_chr$end, sort(part$from)) >
                                 findInterval(on_chr$start - 1, sort(part$from + part$widths)))
    null_touched <- null_touched + vapply(seq_len(n_draws), function(draw) {
      random_from <- part$random_start[, draw]
      sum(findInterval(on_chr$end, sort(random_from)) >
            findInterval(on_chr$start - 1, sort(random_from + part$widths)))
    }, numeric(1))
  }
  data.table(items_touched = observed, recall = observed / nrow(items),
             null_recall = mean(null_touched) / nrow(items),
             recall_excess = (observed - mean(null_touched)) / nrow(items),
             recall_p = (1 + sum(null_touched >= observed)) / (1 + n_draws))
}

density_scan <- function(placed, n_regions) {
  observed <- 0
  null_genes <- numeric(n_draws)
  for (part in placed) {
    on_chr <- gene_table[chr == part$chr]
    gene_starts <- sort(on_chr$start)
    gene_ends <- sort(on_chr$end)
    observed <- observed + sum(findInterval(part$from + part$widths, gene_starts) -
                                 findInterval(part$from - 1, gene_ends))
    random_genes <- findInterval(part$random_start + part$widths, gene_starts) -
      findInterval(part$random_start - 1, gene_ends)
    null_genes <- null_genes + colSums(matrix(random_genes, nrow = nrow(part$random_start)))
  }
  data.table(n_regions = n_regions, genes_per_region = observed / n_regions,
             null_genes_per_region = mean(null_genes) / n_regions,
             ratio = observed / mean(null_genes),
             empirical_p = (1 + sum(null_genes >= observed)) / (1 + n_draws))
}

# The nearest called SNP to each classic gene, 0 when it lies inside the gene.
nearest_snp <- function(snps, regions) {
  classic[, {
    on_chr <- snps[snps$chr == .BY$chr]
    regions_on_chr <- regions[regions$chr == .BY$chr]
    inside <- any(regions_on_chr$start <= end & regions_on_chr$end >= start)
    if (nrow(on_chr) == 0L) {
      .(in_region = inside, nearest_snp = NA_character_, distance_bp = NA_real_, min_pvalue = NA_real_)
    } else {
      gap <- pmax(start - on_chr$pos, on_chr$pos - end, 0)
      nearest <- which.min(gap)
      .(in_region = inside, nearest_snp = on_chr$SNPID[nearest], distance_bp = gap[nearest],
        min_pvalue = on_chr$min_pvalue[nearest])
    }
  }, by = .(chr, start, end, item, label)]
}

# The same columns for a scan of windows, which have no SNPs: the nearest window, 0 when it holds the gene.
nearest_window <- function(regions) {
  classic[, {
    windows <- regions[regions$chr == .BY$chr]
    if (nrow(windows) == 0L) {
      .(in_region = FALSE, nearest_snp = NA_character_, distance_bp = NA_real_, min_pvalue = NA_real_)
    } else {
      gap <- pmax(windows$start - end, start - windows$end, 0)
      nearest <- which.min(gap)
      .(in_region = gap[nearest] == 0, nearest_snp = windows$region_id[nearest], distance_bp = gap[nearest],
        min_pvalue = NA_real_)
    }
  }, by = .(chr, start, end, item, label)]
}

set_rows <- list()
score_rows <- list()
density_rows <- list()
nearest_rows <- list()
cell_rows <- list()

# One scan: its regions scored on its lists, with recall, density and classic genes.
add_scan <- function(arm, scan, regions, snps, scan_lists) {
  if (nrow(regions) == 0L) {
    set_rows[[length(set_rows) + 1L]] <<- data.table(arm = arm, scan = scan, n_snps = 0L, n_regions = 0L)
    return(invisible())
  }
  merged <- merge_intervals(regions)
  placed <- place_regions(regions)
  n_genes <- sum(vapply(split(gene_table, by = "chr"), function(on_chr) {
    on_regions <- merged[chr == on_chr$chr[1]]
    sum(touches(on_chr$start, on_chr$end, on_regions$start, on_regions$end))
  }, numeric(1)))
  set_rows[[length(set_rows) + 1L]] <<- data.table(
    arm = arm, scan = scan, n_snps = if (is.null(snps)) NA_integer_ else nrow(snps),
    n_regions = nrow(regions), span_mb = sum(merged$end - merged$start) / 1e6,
    genome_share = sum(merged$end - merged$start) / genome_bp,
    median_region_kb = median(regions$end - regions$start) / 1e3, n_genes = n_genes
  )
  density_rows[[length(density_rows) + 1L]] <<- cbind(data.table(arm = arm, scan = scan),
                                                      density_scan(placed, nrow(regions)))
  for (list_name in names(scan_lists)) {
    items <- scan_lists[[list_name]]
    if (nrow(items) == 0L) next
    score_rows[[length(score_rows) + 1L]] <<- cbind(
      data.table(arm = arm, scan = scan, list = list_name, n_items = nrow(items)),
      score_scan(regions, items, chrom_lengths, n_draws, seed),
      recall_scan(placed, items)
    )
  }
  nearest_rows[[length(nearest_rows) + 1L]] <<- cbind(
    data.table(arm = arm, scan = scan),
    if (is.null(snps)) nearest_window(regions) else nearest_snp(snps, regions)
  )
}

# Union set against the pipeline, row for row --------------------------------------------------
sorted_ids <- function(ids) vapply(strsplit(ids, ",", fixed = TRUE), function(x) paste(sort(x), collapse = ","), "")
ours <- fread(file.path(region_dir, "regions_gea.union.tsv"), colClasses = c(chr = "character"))
theirs <- fread(file.path(results, "GEA", "tables", "regions_combined.tsv"), colClasses = c(chr = "character"))
setorder(ours, chr, start)
setorder(theirs, chr, start)
same_rows <- nrow(ours) == nrow(theirs) &&
  identical(ours$region_id, theirs$region_id) &&
  identical(ours$chr, theirs$chr) && all(ours$start == theirs$start) && all(ours$end == theirs$end) &&
  all(ours$snp_count == theirs$snp_count) && identical(sorted_ids(ours$snp_ids), sorted_ids(theirs$snp_ids))
message(sprintf("Union set against regions_combined.tsv, row for row: %d rows, %s", nrow(ours),
                if (same_rows) "identical" else "DIFFERENT"))
if (!same_rows) stop("The union set's regions differ from regions_combined.tsv")

# Arms ------------------------------------------------------------------------------------------
read_method_pvalues <- function(module, method) {
  fread(file.path(results, module, "tables", "methods", method, sprintf("%s_pvalues_K%s.tsv", method, k_best)),
        colClasses = c(chr = "character", SNPID = "character"))
}
gea_pvalues <- lapply(setNames(nm = c("EMMAX", "LFMM", "RDA")), function(m) read_method_pvalues("GEA", m))
gwas_pvalues <- lapply(setNames(nm = c("EMMAX", "BLINK")), function(m) read_method_pvalues("GWAS", m))

set_dirs <- Sys.glob(file.path(results, "_intermediate", "snp_sets", "gea.*"))
declarations <- lapply(setNames(set_dirs, basename(set_dirs)),
                       function(d) read_snp_set_decl(file.path(d, "declaration.json")))

arms <- CJ(rule = c("bonf", "qval"), distance = c(as.integer(clumping), 10000L), sorted = FALSE)
for (row in seq_len(nrow(arms))) {
  rule <- arms$rule[row]
  distance <- arms$distance[row]
  arm <- sprintf("%s_%.0fkb", rule, distance / 1e3)
  primary <- rule == "bonf" && distance == clumping

  for (set_name in names(declarations)) {
    declaration <- declarations[[set_name]]
    cells <- lapply(declaration$cells, function(cell) list(type = rule, value = 0.05))
    selected <- select_snp_set(gea_pvalues, cells, declaration$pairs, strategy = declaration$strategy,
                               gap = distance)
    if (primary) {
      stored <- fread(file.path(results, "_intermediate", "snp_sets", set_name, "selected_snps.tsv"),
                      colClasses = c(chr = "character", SNPID = "character"))
      if (!setequal(stored$SNPID, selected$snps$SNPID)) stop("Selection does not reproduce ", set_name)
    }
    if (set_name == "gea.union") {
      cell_rows[[length(cell_rows) + 1L]] <- cbind(data.table(arm = arm), as.data.table(selected$cells))
    }
    add_scan(arm, set_name, cluster_snps_to_regions(selected$snps, distance), selected$snps,
             lists[gea_lists])
    message(sprintf("%s %s: %d SNPs", arm, set_name, nrow(selected$snps)))
  }

  for (trait in strsplit(config$GWAS$traits, ",", fixed = TRUE)[[1]]) {
    pairs <- paste0(trait, "::", names(gwas_pvalues))
    cells <- setNames(lapply(pairs, function(pair) list(type = rule, value = 0.05)), pairs)
    selected <- select_snp_set(gwas_pvalues, cells, pairs, strategy = "Union", gap = distance)
    regions <- cluster_snps_to_regions(selected$snps, distance)
    if (primary) {
      stored <- fread(file.path(region_dir, sprintf("regions_gwas.%s.tsv", trait)), colClasses = c(chr = "character"))
      if (!setequal(paste(stored$chr, stored$start, stored$end), paste(regions$chr, regions$start, regions$end))) {
        stop("Selection does not reproduce the GWAS regions of ", trait)
      }
    }
    cell_rows[[length(cell_rows) + 1L]] <- cbind(data.table(arm = arm), as.data.table(selected$cells))
    add_scan(arm, paste0("gwas.", trait), regions, selected$snps,
             c(lists["flor_id"], aragwas_lists(trait)))
    message(sprintf("%s gwas.%s: %d SNPs", arm, trait, nrow(selected$snps)))
  }
}

# WZA: significant windows, adjacent ones merged, at the threshold the run used --------------------
window_regions <- function(windows) {
  bounds <- str_match(windows$SNPID, ":([0-9]+)-([0-9]+)$")
  spans <- unique(data.table(chr = windows$chr, start = as.numeric(bounds[, 2]),
                             end = as.numeric(bounds[, 3]) + 1))
  merged <- merge_intervals(spans)
  merged[, `:=`(end = end - 1, region_id = sprintf("%s_%.0f-%.0f", chr, start, end - 1))][]
}

read_windows <- function(module, method) {
  fread(file.path(results, module, "tables", "methods", method,
                  sprintf("%s_wza_K%s_sig_windows_bonf_0.05.tsv", method, k_best)),
        colClasses = c(chr = "character", SNPID = "character"))
}

gea_windows <- lapply(setNames(nm = names(gea_pvalues)), function(m) read_windows("GEA", m))
for (method in names(gea_windows)) {
  if (nrow(gea_windows[[method]]) == 0L) next
  add_scan("wza_bonf", paste0("gea.", method), window_regions(gea_windows[[method]]), NULL, lists[gea_lists])
}
all_windows <- rbindlist(gea_windows)
add_scan("wza_bonf", "gea.union", window_regions(all_windows), NULL, lists[gea_lists])

gwas_windows <- rbindlist(lapply(names(gwas_pvalues), function(m) read_windows("GWAS", m)))
for (this_trait in unique(gwas_windows$trait)) {
  add_scan("wza_bonf", paste0("gwas.", this_trait), window_regions(gwas_windows[trait == this_trait]),
           NULL, c(lists["flor_id"], aragwas_lists(this_trait)))
}

# RDA by axis count, once rda_axes_sweep.R full has written its hits --------------------------------
axis_hits_file <- here("ANALYSIS", project, "pregea", "rda_axes_full_hits.tsv")
if (file.exists(axis_hits_file)) {
  axis_hits <- fread(axis_hits_file, colClasses = c(chr = "character", SNPID = "character"))
  setnames(axis_hits, "pvalue", "min_pvalue")
  for (k in sort(unique(axis_hits$K))) {
    hits <- axis_hits[K == k, .(SNPID, chr, pos, min_pvalue)]
    if (nrow(hits) > 0L) {
      add_scan("rda_axes_79kb", sprintf("rda.K%d", k), cluster_snps_to_regions(hits, clumping), hits,
               lists[gea_lists])
    }
  }
}

# Write -----------------------------------------------------------------------------------------
sets <- rbindlist(set_rows, fill = TRUE)
union_rows <- sets[scan == "gea.union", .(arm, union_snps = n_snps, union_regions = n_regions)]
sets <- merge(sets, union_rows, by = "arm", all.x = TRUE, sort = FALSE)
sets[, `:=`(share_of_union_snps = n_snps / union_snps, share_of_union_regions = n_regions / union_regions)]
sets[, c("union_snps", "union_regions") := NULL]
scores <- rbindlist(score_rows)
classic_rows <- rbindlist(nearest_rows)
classic_rows <- classic_rows[, .(arm, scan, gene = label, agi = item, in_region, nearest_snp, distance_bp,
                                 min_pvalue)]

# The fixed table is the anchor: the same regions and seed give the same hits, and the same null up to
# the order of a scan's rows (a GWAS scan's rows arrive in another order, which reassigns its draws).
fixed <- fread(file.path(region_dir, "scores.tsv"))
anchor <- merge(fixed, scores[arm == sprintf("bonf_%.0fkb", clumping / 1e3)], by = c("scan", "list"),
                suffixes = c("", ".check"))
expected_rows <- nrow(fixed[!startsWith(scan, "gxg.")])
largest_gap <- anchor[, max(abs(null_mean - null_mean.check))]
if (nrow(anchor) != expected_rows || anchor[, any(n_with_item != n_with_item.check)] || largest_gap > 0.02) {
  stop("The bonf arm does not reproduce score_regions.R's scores.tsv")
}
message(sprintf("Anchor: %d rows of scores.tsv reproduced, largest null difference %.4f", nrow(anchor),
                largest_gap))

fwrite(sets, file.path(out_dir, "sets.tsv"), sep = "\t")
fwrite(scores, file.path(out_dir, "scores.tsv"), sep = "\t")
fwrite(rbindlist(density_rows), file.path(out_dir, "gene_density.tsv"), sep = "\t")
fwrite(classic_rows, file.path(out_dir, "classic_nearest.tsv"), sep = "\t")
fwrite(rbindlist(cell_rows), file.path(out_dir, "fdr_cells.tsv"), sep = "\t")
message("Wrote ", out_dir)
