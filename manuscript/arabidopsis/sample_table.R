#!/usr/bin/env Rscript
# sample_table.R: what became of each of the 1,135 accessions in the VCF.
#
# Reads the 1,135 accessions of data/Arabidopsis_1001G_metadata.tsv.bak, the project's inputs under
# data/arabidopsis1135/, and what mode=processing removed under Arabidopsis1135_results/. Writes
# ANALYSIS/Arabidopsis1135/sample_table.tsv, one row per accession:
#   sample accession country latitude longitude site FT10 FT16 F_MISS in_grenenet p0 status reason ibs
#
# latitude and longitude are the accession's own, empty where the source has 0, 0. site is the
# WorldClim 30" cell, for accessions that entered the project. FT10 and FT16 are the project's where
# it has them, else the source table's. in_grenenet covers the 225 founders in the VCF. status is
# kept or excluded, and reason is north_america, no_coordinates or no_worldclim_value before the
# pipeline, sample_miss or relatedness inside it. ibs is the pair's IBS for a relatedness removal.
#
# The counts inside the pipeline are not hard-coded: they must match filtering_summary.tsv.
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -e HOME=/tmp -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/sample_table.R

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

project <- "Arabidopsis1135"
input_dir <- here("data", "arabidopsis1135")
results <- here(paste0(project, "_results"))
samples_dir <- file.path(results, "_intermediate", "samples")
tables_dir <- file.path(results, "Processing", "tables")
out_file <- here("ANALYSIS", project, "sample_table.tsv")

# What build_inputs.R measured; anything else stops the script.
n_accessions <- 1135L
n_before_pipeline <- c(north_america = 125L, no_coordinates = 4L, no_worldclim_value = 37L)

#-------------------------------------------------------------------------- inputs
accessions <- fread(here("data", "Arabidopsis_1001G_metadata.tsv.bak"),
                    colClasses = c(site = "character", sample = "character"))
setnames(accessions, c("site", "sample"), c("country", "accession"))
accessions[, sample := paste0("S", accession)]
stopifnot(nrow(accessions) == n_accessions, !anyDuplicated(accessions$sample))

excluded <- fread(file.path(input_dir, "excluded_before_pipeline.tsv"),
                  select = c("sample", "reason"), colClasses = c(sample = "character"))
project_samples <- fread(file.path(input_dir, "Arabidopsis1135_metadata.tsv"),
                         colClasses = c(site = "character", sample = "character"))
founders <- fread(file.path(input_dir, "founder_sites.tsv"),
                  select = c("sample", "p0"), colClasses = c(sample = "character"))

missing_stats <- fread(file.path(samples_dir, "samples_missing_stats.tsv"),
                       select = c("IID", "F_MISS"), colClasses = c(IID = "character"))
removed_sample_miss <- fread(file.path(samples_dir, "samples_removed.list"), header = FALSE,
                             colClasses = "character")$V2
removed_relatedness <- fread(file.path(tables_dir, "relatedness_removed.tsv"),
                             colClasses = c(IID = "character"))
kept_by_pipeline <- fread(file.path(tables_dir, "metadata.tsv"),
                          select = "sample", colClasses = c(sample = "character"))$sample
filtering <- fread(file.path(tables_dir, "filtering_summary.tsv"))

stopifnot(!anyDuplicated(c(excluded$sample, removed_sample_miss, removed_relatedness$IID)))

#-------------------------------------------------------------------------- one row per accession
accessions[latitude == 0 & longitude == 0, c("latitude", "longitude") := NA_real_]
accessions[, `:=`(site = NA_character_, F_MISS = NA_real_, in_grenenet = FALSE, p0 = NA_real_,
                  reason = NA_character_, ibs = NA_real_)]

accessions[project_samples, on = "sample", `:=`(site = i.site, FT10 = i.FT10, FT16 = i.FT16)]
accessions[missing_stats, on = c(sample = "IID"), F_MISS := i.F_MISS]
accessions[founders, on = "sample", `:=`(in_grenenet = TRUE, p0 = i.p0)]

accessions[excluded, on = "sample", reason := i.reason]
accessions[sample %chin% removed_sample_miss, reason := "sample_miss"]
accessions[removed_relatedness, on = c(sample = "IID"), `:=`(reason = "relatedness", ibs = i.IBS)]
accessions[, status := fifelse(is.na(reason), "kept", "excluded")]

#-------------------------------------------------------------------------- the counts must close
samples_at <- function(step) filtering[stage == step, n_samples]
count_reason <- function(why) accessions[reason == why, .N]

stopifnot(
  vapply(names(n_before_pipeline), count_reason, integer(1)) == n_before_pipeline,
  n_accessions - sum(n_before_pipeline) == samples_at("Raw VCF"),
  samples_at("Raw VCF") - samples_at("After sample missingness filter") == count_reason("sample_miss"),
  samples_at("After sample missingness filter") - samples_at("After relatedness filter") ==
    count_reason("relatedness"),
  accessions[status == "kept", .N] == samples_at("After relatedness filter"),
  setequal(accessions[status == "kept", sample], kept_by_pipeline),
  # F_MISS exists exactly for the accessions that entered the pipeline.
  setequal(accessions[!is.na(F_MISS), sample], project_samples$sample)
)

#-------------------------------------------------------------------------- write
setcolorder(accessions, c("sample", "accession", "country", "latitude", "longitude", "site", "FT10",
                          "FT16", "F_MISS", "in_grenenet", "p0", "status", "reason", "ibs"))
fwrite(accessions, out_file, sep = "\t", na = "")

print(accessions[, .N, by = .(status, reason)][order(status, reason)])
print(accessions[in_grenenet == TRUE, .N, by = .(status, reason)][order(status, reason)])
message("Wrote ", out_file)
