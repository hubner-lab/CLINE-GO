#!/usr/bin/env Rscript
# offset_similarity.R: how alike are the offsets of different SNP sets and methods? Post hoc.
#
#   Rscript manuscript/arabidopsis/offset_similarity.R <project>
#
# Written after offset_eval.R's first result on Arabidopsis1135 (Phase 14), to weigh whether random
# SNP sets (level 3) are worth running. Per garden, Spearman across the founders that have an offset
# between (a) the sets of one method, (b) the methods on one set, (c) each run and the BIO1 climate
# distance (BIO1 home - BIO1 garden)^2; then the median over gardens. Founders and home BIO1 as in
# offset_eval.R: founder_sites.tsv joined to climate_present.tsv.
#
# Reads data/arabidopsis1135/ and <project>_results/Maladaptation/tables/. Writes
# ANALYSIS/Arabidopsis1135/offset/<project>/similarity/: pairs.tsv (every run pair, median rho) and
# climate_distance.tsv (every run, median rho). Run from the repo root, in the image.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: offset_similarity.R <project>")
project <- args[1]
data_dir <- here("data", "arabidopsis1135")
tables <- here(paste0(project, "_results"), "Maladaptation", "tables")
out_dir <- here("ANALYSIS", "Arabidopsis1135", "offset", project, "similarity")
methods <- c("gradient_forest", "geometric_offset", "rda_offset")

founders <- fread(file.path(data_dir, "founder_sites.tsv"), select = c("sample", "in_metadata", "site"))
climate <- fread(file.path(data_dir, "climate_present.tsv"), select = c("site", "bio_1"))
home <- founders[in_metadata == TRUE][climate, on = "site", nomatch = NULL][, .(sample, bio_home = bio_1)]
gardens <- fread(file.path(data_dir, "garden_sites.tsv"), select = c("garden", "bio_1"))
setnames(gardens, "bio_1", "bio_garden")

paths <- Sys.glob(file.path(tables, "*", "*", "garden_*", "genetic_offset_site.tsv"))
if (length(paths) == 0L) stop("No garden offsets under ", tables)
offsets <- rbindlist(lapply(paths, function(path) {
  run_dir <- dirname(dirname(path))
  x <- fread(path, select = c("sample", "genetic_offset"), colClasses = c(sample = "character"))
  x[, `:=`(method = basename(dirname(run_dir)),
           set = sub("_nospatial$", "", basename(run_dir)),
           garden = as.integer(sub("^garden_", "", basename(dirname(path)))))]
}))
offsets <- offsets[home, on = "sample", nomatch = NULL][gardens, on = "garden", nomatch = NULL]
offsets[, climate_distance := (bio_home - bio_garden)^2]
offsets[, run := paste(method, set, sep = "/")]
stopifnot(offsets[, .N, by = .(garden, sample, run)][, all(N == 1L)])

wide <- dcast(offsets, garden + sample + climate_distance ~ run, value.var = "genetic_offset")
runs <- setdiff(names(wide), c("garden", "sample", "climate_distance"))
garden_ids <- sort(unique(wide$garden))
rho <- array(NA_real_, c(length(runs), length(runs), length(garden_ids)))
rho_climate <- matrix(NA_real_, length(runs), length(garden_ids))
n_founders <- integer(length(garden_ids))
for (k in seq_along(garden_ids)) {
  one <- wide[garden == garden_ids[k]]
  m <- as.matrix(one[, ..runs])
  if (anyNA(m)) stop("A run lacks founders that another run has, in garden ", garden_ids[k])
  n_founders[k] <- nrow(m)
  rho[, , k] <- cor(m, method = "spearman")
  rho_climate[, k] <- cor(m, one$climate_distance, method = "spearman")
}
median_rho <- apply(rho, c(1, 2), median)
upper <- which(upper.tri(median_rho), arr.ind = TRUE)
split_run <- function(x) tstrsplit(x, "/", fixed = TRUE)
pairs <- data.table(run_a = runs[upper[, 1]], run_b = runs[upper[, 2]], median_rho = median_rho[upper])
pairs[, c("method_a", "set_a") := split_run(run_a)][, c("method_b", "set_b") := split_run(run_b)]
pairs[, kind := fifelse(method_a == method_b, "same method", fifelse(set_a == set_b, "same set", "other"))]
by_climate <- data.table(run = runs, median_rho = apply(rho_climate, 1, median))
by_climate[, c("method", "set") := split_run(run)]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(pairs[, .(kind, method_a, set_a, method_b, set_b, median_rho)], file.path(out_dir, "pairs.tsv"), sep = "\t")
fwrite(by_climate[, .(method, set, median_rho)], file.path(out_dir, "climate_distance.tsv"), sep = "\t")

message("Founders per garden: ", paste(unique(n_founders), collapse = ", "), "; gardens: ", length(garden_ids))
message("(a) sets within a method, range of the median rho over set pairs:")
print(pairs[kind == "same method", .(min = round(min(median_rho), 3), median = round(median(median_rho), 3),
                                     max = round(max(median_rho), 3)), by = .(method = method_a)])
message("(b) methods on one set, median rho:")
print(dcast(pairs[kind == "same set", .(set = set_a, pair = paste(method_a, method_b, sep = " ~ "),
                                        rho = round(median_rho, 3))], set ~ pair, value.var = "rho"))
message("(c) each run against the BIO1 climate distance, median rho:")
print(dcast(by_climate[, .(set, method, rho = round(median_rho, 3))], set ~ method, value.var = "rho"))
message("Wrote ", out_dir)
