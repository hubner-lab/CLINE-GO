#!/usr/bin/env Rscript
# build_inputs.R: the input tables of the Arabidopsis1135 project.
#
# Reads the 1,135 accessions of data/Arabidopsis_1001G_metadata.tsv.bak, WorldClim 2.1 at 30"
# (data/wc2.1_30s/), AraPheno FT10 and FT16, and the GrENE-net tables copied to
# data/arabidopsis1135/sources/. Writes only under data/arabidopsis1135/:
#   Arabidopsis1135_metadata.tsv   site sample latitude longitude FT10 FT16
#   excluded_before_pipeline.tsv   sample, country, coordinates, reason
#   climate_present.tsv            site + bio_1..bio_19, one row per site
#   gardens/garden_<id>.tsv        site + bio_1..bio_19, the garden's values on every row
#   garden_sites.tsv               the 30 gardens, their cell, climate and range check
#   founder_sites.tsv              GrENE-net founders: sample, p0, in VCF, in metadata, site
#
# site = the WorldClim 30" cell an accession falls in, coordinates = that cell's centre: the custom
# climate route needs every site in its own cell, and inside a cell WorldClim has one value anyway.
#
# Run from the repo root, in the image:
#   docker run --rm --user $(id -u):$(id -g) -e USER=pipeline -v $PWD:/data -w /data \
#     cline-go:6cdedcc Rscript manuscript/arabidopsis/build_inputs.R

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(readxl)
  library(terra)
})

# What the design session measured on these files; anything else stops the script.
n_accessions <- 1135L
n_north_america <- 125L
n_no_coordinates <- 4L
n_no_worldclim <- 37L
n_samples_kept <- 969L
n_sites_kept <- 703L
n_founders <- 231L
n_gardens <- 30L
n_replicates <- 325L

excluded_garden <- 33L             # one replicate, left out by GrENE-net as well
staging_resolution <- 0.008        # Climate.custom.grid_resolution
staging_buffer <- 0.5              # stage_custom_climate.R: max(GRID_RES, 0.5)
bio_names <- paste0("bio_", 1:19)
offset_predictors <- paste0("bio_", c(1, 2, 4, 8, 12, 15, 18))

out_dir <- here("data", "arabidopsis1135")
src_dir <- file.path(out_dir, "sources")
grene_dir <- file.path(src_dir, "dl_grene")

site_id <- function(cell) sprintf("c%.0f", cell)   # %.0f: no 1e+08 for a round cell number

write_tsv <- function(x, path) fwrite(x, path, sep = "\t", na = "")

bio_at_cells <- function(cells) {
  values <- as.data.table(terra::extract(bio_rast, cells))
  setnames(values, bio_names)
  values
}

bio_rast <- rast(here("data", "wc2.1_30s", sprintf("wc2.1_30s_bio_%d.tif", 1:19)))
names(bio_rast) <- bio_names

#-------------------------------------------------------------------------- accessions
accessions <- fread(here("data", "Arabidopsis_1001G_metadata.tsv.bak"),
                    colClasses = c(site = "character", sample = "character"))
setnames(accessions, "site", "country")
stopifnot(nrow(accessions) == n_accessions, !anyDuplicated(accessions$sample))

accessions[, cell := cellFromXY(bio_rast, as.matrix(accessions[, .(longitude, latitude)]))]
accessions[, has_worldclim := complete.cases(bio_at_cells(cell))]
accessions[, reason := fcase(
  country %in% c("US", "Canada"), "north_america",
  latitude == 0 & longitude == 0, "no_coordinates",
  !has_worldclim, "no_worldclim_value",
  default = NA_character_
)]

reason_counts <- accessions[, .N, by = reason]
print(reason_counts)
stopifnot(
  reason_counts[reason == "north_america", N] == n_north_america,
  reason_counts[reason == "no_coordinates", N] == n_no_coordinates,
  reason_counts[reason == "no_worldclim_value", N] == n_no_worldclim
)

#-------------------------------------------------------------------------- phenotypes
# AraPheno 261 and 262 are the phenotypes the old metadata came from: the two must agree. One
# accession (8424, Kas-2) has two observations, so an accession takes the mean of its own, and the
# old value must match one of them.
read_trait <- function(phenotype_id, trait) {
  values <- fread(file.path(src_dir, "arapheno_ft", sprintf("phenotype_%d.csv", phenotype_id)),
                  select = c("accession_id", "phenotype_value"),
                  colClasses = c(accession_id = "character"))
  values <- values[, .(mean_value = mean(phenotype_value), low = min(phenotype_value),
                       high = max(phenotype_value)), by = accession_id]
  setnames(values, c("sample", paste0(trait, c("_arapheno", "_low", "_high"))))
  values
}

accessions <- read_trait(261L, "FT10")[accessions, on = "sample"]
accessions <- read_trait(262L, "FT16")[accessions, on = "sample"]

for (trait in c("FT10", "FT16")) {
  old <- accessions[[trait]]
  new <- accessions[[paste0(trait, "_arapheno")]]
  n_differ <- sum(abs(old - accessions[[paste0(trait, "_low")]]) > 1e-6 &
                    abs(old - accessions[[paste0(trait, "_high")]]) > 1e-6, na.rm = TRUE)
  n_lost <- sum(!is.na(old) & is.na(new))
  n_added <- sum(is.na(old) & !is.na(new))
  message(sprintf("%s: %d values in both, %d differ, %d only in the old file, %d only in AraPheno",
                  trait, sum(!is.na(old) & !is.na(new)), n_differ, n_lost, n_added))
  stopifnot(n_differ == 0L, n_lost == 0L, n_added == 0L)
}

#-------------------------------------------------------------------------- metadata, climate
kept <- accessions[is.na(reason)]
stopifnot(nrow(kept) == n_samples_kept, uniqueN(kept$cell) == n_sites_kept)

centres <- xyFromCell(bio_rast, kept$cell)
metadata <- kept[, .(site = site_id(cell), sample = paste0("S", sample),
                     latitude = centres[, "y"], longitude = centres[, "x"],
                     FT10 = FT10_arapheno, FT16 = FT16_arapheno)]

excluded <- accessions[!is.na(reason),
                       .(sample = paste0("S", sample), country, latitude, longitude, reason)]

site_cells <- sort(unique(kept$cell))
climate_present <- cbind(data.table(site = site_id(site_cells)), bio_at_cells(site_cells))
stopifnot(!anyNA(climate_present))

# stage_custom_climate.R stops when two sites share a cell of its grid. Rebuild that grid here.
site_xy <- unique(metadata[, .(site, longitude, latitude)])
stopifnot(!anyDuplicated(site_xy$site))
staging_grid <- rast(
  ext(min(site_xy$longitude) - staging_buffer, max(site_xy$longitude) + staging_buffer,
      min(site_xy$latitude) - staging_buffer, max(site_xy$latitude) + staging_buffer),
  resolution = staging_resolution, crs = "EPSG:4326"
)
staging_cells <- cellFromXY(staging_grid, as.matrix(site_xy[, .(longitude, latitude)]))
stopifnot(!anyNA(staging_cells), uniqueN(staging_cells) == n_sites_kept)
message(sprintf("Staging grid at %g degrees: %d x %d = %.3g cells, %d sites in distinct cells",
                staging_resolution, nrow(staging_grid), ncol(staging_grid), ncell(staging_grid),
                uniqueN(staging_cells)))

#-------------------------------------------------------------------------- gardens
s9_columns <- names(fread(file.path(grene_dir, "ds", "Dataset_S9_evovled_accession_frequency.txt"),
                          nrows = 0L))
year1_columns <- grep("^[0-9]+_1_[0-9]+$", s9_columns, value = TRUE)
year1_gardens <- as.integer(sub("_.*", "", year1_columns))
year1_columns <- year1_columns[year1_gardens != excluded_garden]
garden_ids <- setdiff(unique(year1_gardens), excluded_garden)
stopifnot(length(garden_ids) == n_gardens, length(year1_columns) == n_replicates)

# Only the columns needed: Table S1 also holds contributors' names and email addresses.
gardens <- as.data.table(read_excel(
  file.path(grene_dir, "Supplemental_tables", "Table_S1_site_information.xlsx"), sheet = 1
))[, .(garden = as.integer(SITE_CODE), name = SITE_NAME, country = COUNTRY,
       longitude = as.numeric(LONGITUDE), latitude = as.numeric(LATITUDE))]

# Two typos in Table S1: the sign of garden 5 (Madrid), the decimal point of garden 12 (New York).
stopifnot(gardens[garden == 5L, longitude] == 3.689585,
          gardens[garden == 12L, longitude] == -73873887)
gardens[garden == 5L, longitude := -3.689585]
gardens[garden == 12L, longitude := -73.873887]

gardens <- gardens[garden %in% garden_ids][order(garden)]
stopifnot(nrow(gardens) == n_gardens)

gardens[, cell := cellFromXY(bio_rast, as.matrix(gardens[, .(longitude, latitude)]))]
garden_bio <- bio_at_cells(gardens$cell)
if (anyNA(garden_bio)) {
  stop("No WorldClim value at garden(s): ",
       paste(gardens$garden[!complete.cases(garden_bio)], collapse = ", "))
}
gardens <- cbind(gardens, garden_bio)

# Which predictors put a garden outside the climate the models train on. The offset models see only
# offset_predictors; the design's four out-of-range gardens were counted over all 19.
outside_range <- function(predictors) {
  training_range <- climate_present[, lapply(.SD, range), .SDcols = predictors]
  vapply(seq_len(nrow(gardens)), function(i) {
    outside <- predictors[vapply(predictors, function(p) {
      value <- gardens[[p]][i]
      value < training_range[[p]][1] || value > training_range[[p]][2]
    }, logical(1))]
    paste(outside, collapse = ",")
  }, character(1))
}
gardens[, outside_offset_predictors := outside_range(offset_predictors)]
gardens[, outside_all_bio := outside_range(bio_names)]

#-------------------------------------------------------------------------- founders
vcf_connection <- file(here("data", "Arabidopsis_1001G.vcf"), open = "r")
repeat {
  vcf_line <- readLines(vcf_connection, n = 1L)
  stopifnot(length(vcf_line) == 1L)
  if (startsWith(vcf_line, "#CHROM")) break
}
close(vcf_connection)
vcf_samples <- strsplit(vcf_line, "\t", fixed = TRUE)[[1]][-(1:9)]
stopifnot(length(vcf_samples) == n_accessions,
          all(paste0("S", accessions$sample) %in% vcf_samples))

founders <- fread(file.path(grene_dir, "ds", "Dataset_S5_founder_accession_frequency.txt"),
                  header = FALSE, col.names = c("founder", "p0"), colClasses = c(V1 = "character"))
stopifnot(nrow(founders) == n_founders, !anyDuplicated(founders$founder))
founders <- accessions[, .(founder = sample, country, cell, reason)][founders, on = "founder"]
founder_sites <- founders[, .(
  founder, sample = paste0("S", founder), p0, country,
  in_vcf = paste0("S", founder) %in% vcf_samples,
  in_metadata = paste0("S", founder) %in% metadata$sample,
  site = fifelse(paste0("S", founder) %in% metadata$sample, site_id(cell), NA_character_),
  excluded_because = reason
)]
print(founder_sites[, .N, by = .(in_vcf, in_metadata, excluded_because)])

#-------------------------------------------------------------------------- write
dir.create(file.path(out_dir, "gardens"), recursive = TRUE, showWarnings = FALSE)
write_tsv(metadata, file.path(out_dir, "Arabidopsis1135_metadata.tsv"))
write_tsv(excluded, file.path(out_dir, "excluded_before_pipeline.tsv"))
write_tsv(climate_present, file.path(out_dir, "climate_present.tsv"))
write_tsv(gardens, file.path(out_dir, "garden_sites.tsv"))
write_tsv(founder_sites, file.path(out_dir, "founder_sites.tsv"))

site_column <- climate_present[, .(site)]
for (i in seq_len(nrow(gardens))) {
  garden_table <- cbind(site_column, garden_bio[i])   # data.table recycles the garden's row
  stopifnot(nrow(garden_table) == n_sites_kept)
  write_tsv(garden_table, file.path(out_dir, "gardens", sprintf("garden_%d.tsv", gardens$garden[i])))
}

message(sprintf("Wrote %d samples, %d sites, %d gardens, %d founders in metadata",
                nrow(metadata), nrow(climate_present), nrow(gardens),
                founder_sites[in_metadata == TRUE, .N]))
print(gardens[outside_all_bio != "", .(garden, country, outside_offset_predictors, outside_all_bio)])
