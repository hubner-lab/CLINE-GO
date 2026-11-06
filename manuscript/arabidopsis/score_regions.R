#!/usr/bin/env Rscript
# score_regions.R: how many association regions hold a known gene, and how many chance would give.
#
#   Rscript manuscript/arabidopsis/score_regions.R <project> [n_draws] [seed]
#   Rscript manuscript/arabidopsis/score_regions.R --self-test
#
# Real data has no complete truth: a region without a listed gene is not a false positive, so nothing
# here counts false positives. Each scan reports the share of its regions touching a listed item and
# that share's excess over random regions of the same lengths on the same chromosomes.
#
# Scans:
#   gea.<set>        each declared SNP set; regions from the pipeline's own clustering at the run's
#                    LD distance (a set's table has no method columns, so cluster_snps_to_regions()
#                    gives what build_combined_regions() would, without its per-SNP loop)
#   gwas.<trait>     GWAS/tables/regions_per_trait.tsv
#   gxg.<set>.<trait>  stretches shared by a GEA set's regions and a trait's GWAS regions. The
#                    pipeline writes only trait-level overlap counts, so this one is built here.
#
# Reads <project>_results/ and data/arabidopsis1135/lists/. Writes
# ANALYSIS/Arabidopsis1135/regions/<project>/: regions_<scan>.tsv, scores.tsv, region_hits.tsv,
# classic_genes.tsv, shared_regions.tsv.
#
# Run from the repo root, in the image (the pipeline libraries are baked at /pipeline).

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

source("/pipeline/scripts/R/lib/region_distance.R")
source("/pipeline/scripts/R/lib/regions.R")

# Which lists each kind of scan is scored against (design step 6).
gea_lists <- c("flor_id", "grenenet_blocks", "abiotic")
gwas_lists <- c("flor_id", "aragwas_ft")
shared_lists <- c("flor_id", "grenenet_blocks", "abiotic")

# Sorted, non-overlapping intervals per chromosome, so findInterval can count them.
merge_intervals <- function(items) {
  items <- items[, .(chr, start, end)]
  setorder(items, chr, start)
  items[, block := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = chr]
  items[, .(start = min(start), end = max(end)), by = .(chr, block)][, block := NULL][]
}

# TRUE for each [from, to] that touches a merged interval of one chromosome.
touches <- function(from, to, merged_start, merged_end) {
  findInterval(to, merged_start) > findInterval(from - 1, merged_end)
}

score_scan <- function(regions, items, chrom_lengths, n_draws, seed) {
  merged <- merge_intervals(items)
  n_regions <- nrow(regions)
  if (n_regions == 0L) {
    return(data.table(n_regions = 0L, n_with_item = 0L, share = NA_real_, null_mean = NA_real_,
                      null_q025 = NA_real_, null_q975 = NA_real_, excess = NA_real_,
                      empirical_p = NA_real_))
  }

  observed <- 0L
  null_hits <- numeric(n_draws)
  set.seed(seed)
  for (this_chr in sort(unique(regions$chr))) {
    on_chr <- regions[chr == this_chr]
    chr_length <- chrom_lengths[[this_chr]]
    widths <- pmin(on_chr$end, chr_length) - on_chr$start
    room <- pmax(chr_length - widths, 1)
    random_start <- 1 + floor(runif(nrow(on_chr) * n_draws) * room)
    item_rows <- merged[chr == this_chr]
    if (nrow(item_rows) == 0L) next
    observed <- observed + sum(touches(on_chr$start, on_chr$end, item_rows$start, item_rows$end))
    hit <- touches(random_start, random_start + widths, item_rows$start, item_rows$end)
    null_hits <- null_hits + colSums(matrix(hit, nrow = nrow(on_chr)))
  }

  null_share <- null_hits / n_regions
  data.table(n_regions = n_regions, n_with_item = observed, share = observed / n_regions,
             null_mean = mean(null_share), null_q025 = quantile(null_share, 0.025, names = FALSE),
             null_q975 = quantile(null_share, 0.975, names = FALSE),
             excess = observed / n_regions - mean(null_share),
             empirical_p = (1 + sum(null_hits >= observed)) / (1 + n_draws))
}

# Each region with the listed items it touches.
region_items <- function(regions, items) {
  if (nrow(regions) == 0L || nrow(items) == 0L) return(NULL)
  keyed <- copy(items)[, .(chr, start, end, item, label)]
  setkey(keyed, chr, start, end)
  found <- foverlaps(regions[, .(region_id, chr, start, end)], keyed, nomatch = NULL)
  if (nrow(found) == 0L) return(NULL)
  found[, .(items = paste(item, collapse = ","), labels = paste(unique(label), collapse = "; ")),
        by = region_id]
}

# Stretches where two region tables overlap.
shared_stretches <- function(a, b) {
  if (nrow(a) == 0L || nrow(b) == 0L) return(data.table())
  keyed <- copy(b)[, .(chr, b_start = start, b_end = end, b_region = region_id,
                       start = start, end = end)]
  setkey(keyed, chr, start, end)
  found <- foverlaps(a[, .(chr, start, end, region_id)], keyed, nomatch = NULL)
  found[, .(region_id = sprintf("%s_%d-%d", chr, pmax(i.start, b_start), pmin(i.end, b_end)),
            chr, start = pmax(i.start, b_start), end = pmin(i.end, b_end),
            gea_region = region_id, gwas_region = b_region)]
}

self_test <- function() {
  # One chromosome of 1,000 bp, two overlapping items that merge into 101-200.
  items <- data.table(chr = "1", start = c(101, 151), end = c(160, 200), item = c("a", "b"),
                      label = c("a", "b"))
  regions <- data.table(region_id = c("hit", "miss"), chr = "1", start = c(90, 500),
                        end = c(110, 520))
  merged <- merge_intervals(items)
  stopifnot(nrow(merged) == 1L, merged$start == 101, merged$end == 200)
  stopifnot(identical(touches(c(90, 500, 201, 200, 100), c(110, 520, 300, 300, 100), 101, 200),
                      c(TRUE, FALSE, FALSE, TRUE, FALSE)))

  score <- score_scan(regions, items, c("1" = 1000), n_draws = 20000L, seed = 42L)
  # A 20 bp wide region touches 101-200 when its start is in 81-200: 120 of the 980 possible starts.
  expected <- 120 / 980
  stopifnot(score$n_regions == 2L, score$n_with_item == 1L, score$share == 0.5,
            abs(score$null_mean - expected) < 0.003, score$empirical_p < 0.3, score$empirical_p > 0.1)

  stopifnot(identical(region_items(regions, items)$region_id, "hit"))
  stopifnot(identical(region_items(regions, items)$items, "a"))
  other <- data.table(region_id = "x", chr = "1", start = 100, end = 95 + 10)
  shared <- shared_stretches(regions, other)
  stopifnot(nrow(shared) == 1L, shared$start == 100, shared$end == 105)
  message("score_regions.R self-test passed (null mean ", round(score$null_mean, 4),
          ", expected ", round(expected, 4), ")")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop("usage: score_regions.R <project> [n_draws] [seed] | --self-test")
if (args[1] == "--self-test") {
  self_test()
  quit(save = "no", status = 0)
}

project <- args[1]
n_draws <- if (length(args) >= 2L) as.integer(args[2]) else 10000L
seed <- if (length(args) >= 3L) as.integer(args[3]) else 42L
results <- here(paste0(project, "_results"))
list_dir <- here("data", "arabidopsis1135", "lists")
out_dir <- here("ANALYSIS", "Arabidopsis1135", "regions", project)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_list <- function(name) {
  fread(file.path(list_dir, paste0(name, ".tsv")), colClasses = c(chr = "character", label = "character"))
}
lists <- lapply(setNames(union(gea_lists, gwas_lists), union(gea_lists, gwas_lists)), read_list)
lists$aragwas_ft <- lists$aragwas_ft[above_permutation == TRUE]
classic <- read_list("classic_genes")

sequence_regions <- fread(cmd = sprintf("awk '/^##sequence-region/ {print $2, $4}' %s",
                                        shQuote(here("data", "Arabidopsis_1001G.gff3"))),
                          header = FALSE, col.names = c("chr", "length"),
                          colClasses = c(V1 = "character"))
chrom_lengths <- setNames(as.numeric(sequence_regions$length), sequence_regions$chr)

# The run's own distance, from the config that made it: a one-chromosome pilot must declare bp.
gea_config <- yaml::read_yaml(here(paste0("config_", project, ".yaml")))$GEA
clumping <- resolve_clumping_distance_strict(
  gea_config$snp_clumping_distance,
  file.path(results, "Structure", "tables", "ld_decay_half_distances.tsv"),
  gea_config$clumping_ld_rule, "All"
)
message("Clumping distance: ", clumping, " bp")

# Regions of each scan ------------------------------------------------------------------------
set_dirs <- Sys.glob(file.path(results, "_intermediate", "snp_sets", "gea.*"))
gea_regions <- lapply(setNames(set_dirs, basename(set_dirs)), function(set_dir) {
  snps <- fread(file.path(set_dir, "selected_snps.tsv"), colClasses = c(chr = "character"))
  cluster_snps_to_regions(snps, clumping)
})

union_regions <- gea_regions[["gea.union"]]
pipeline_union <- fread(file.path(results, "GEA", "tables", "regions_combined.tsv"),
                        colClasses = c(chr = "character"))
union_matches <- setequal(union_regions$region_id, pipeline_union$region_id)
message(sprintf("Union set against regions_combined.tsv: %d and %d regions, %s", nrow(union_regions),
                nrow(pipeline_union), if (union_matches) "identical" else "DIFFERENT"))

per_trait <- fread(file.path(results, "GWAS", "tables", "regions_per_trait.tsv"),
                   colClasses = c(chr = "character"))
gwas_regions <- split(per_trait, by = "trait", keep.by = TRUE)
names(gwas_regions) <- paste0("gwas.", names(gwas_regions))

shared_regions <- list()
for (set_name in names(gea_regions)) {
  for (trait_name in names(gwas_regions)) {
    shared_regions[[sprintf("gxg.%s.%s", sub("^gea\\.", "", set_name), sub("^gwas\\.", "", trait_name))]] <-
      shared_stretches(gea_regions[[set_name]], gwas_regions[[trait_name]])
  }
}

# Score ---------------------------------------------------------------------------------------
score_rows <- list()
hit_rows <- list()
classic_rows <- list()

score_family <- function(scans, list_names) {
  for (scan_name in names(scans)) {
    regions <- scans[[scan_name]]
    fwrite(regions, file.path(out_dir, sprintf("regions_%s.tsv", scan_name)), sep = "\t")
    for (list_name in list_names) {
      items <- lists[[list_name]]
      if (list_name == "aragwas_ft") {
        items <- items[label == sub("^gwas\\.", "", scan_name)]
      }
      score_rows[[length(score_rows) + 1L]] <<- cbind(
        data.table(scan = scan_name, list = list_name, n_items = nrow(items)),
        score_scan(regions, items, chrom_lengths, n_draws, seed)
      )
      found <- region_items(regions, items)
      if (!is.null(found)) {
        hit_rows[[length(hit_rows) + 1L]] <<- cbind(data.table(scan = scan_name, list = list_name), found)
      }
    }
    if (nrow(regions) > 0L) {
      classic_rows[[length(classic_rows) + 1L]] <<- classic[, {
        on_chr <- regions[chr == .BY$chr]
        gap <- if (nrow(on_chr)) pmax(on_chr$start - end, start - on_chr$end, 0) else Inf
        .(gene = label, agi = item, in_region = any(gap == 0), distance_bp = min(gap))
      }, by = .(chr, start, end, item, label)][, .(scan = scan_name, gene, agi, in_region, distance_bp)]
    }
  }
}

score_family(gea_regions, gea_lists)
score_family(gwas_regions, gwas_lists)
score_family(shared_regions[vapply(shared_regions, nrow, integer(1)) > 0L], shared_lists)

scores <- rbindlist(score_rows)
fwrite(scores, file.path(out_dir, "scores.tsv"), sep = "\t")
fwrite(rbindlist(hit_rows), file.path(out_dir, "region_hits.tsv"), sep = "\t")
fwrite(rbindlist(classic_rows), file.path(out_dir, "classic_genes.tsv"), sep = "\t")
fwrite(rbindlist(shared_regions, idcol = "scan", fill = TRUE), file.path(out_dir, "shared_regions.tsv"),
       sep = "\t")

message(sprintf("Wrote %d scan-by-list rows to %s", nrow(scores), out_dir))
if (!union_matches) stop("The union set's regions differ from the pipeline's regions_combined.tsv")
