#!/usr/bin/env Rscript
# offset_eval.R: do genomic offsets predict first-year fitness in the GrENE-net gardens?
#
#   Rscript manuscript/arabidopsis/offset_eval.R <project> [--drop-gardens 10,13,26,60]
#   Rscript manuscript/arabidopsis/offset_eval.R --reproduce
#   Rscript manuscript/arabidopsis/offset_eval.R --self-test
#
# Metric, from GrENE-net (GENOMIC_OFFSET_climate_distance_loo_prediction.R): for each replicate of
# generation 1, garden 33 left out, 325 in all, the Spearman correlation across founders between
# p1/p0 (Dataset S9 over Dataset S5) and a predictor of fitness. A higher offset should mean lower
# fitness, so the predictor is the reversed offset. The baseline is climate distance,
# -(BIO1 home - BIO1 garden)^2.
#
# Rules, fixed here before any offset exists:
#   unit     All plots of a garden share one offset column, so the garden is the independent unit.
#            Its value is the mean rho of its replicates, and tests run over gardens, never over
#            replicates.
#   pairing  Each run's baseline uses exactly the founders that have an offset in that garden, with
#            home BIO1 from climate_present.tsv, the table the offsets use.
#   spread   Mean, median and quartiles over the replicates, as Table S12 allows.
#   level 1  Mean rho over the replicates above 0, and garden mean rho above 0 in most gardens.
#   level 2  Mean rho above the run's own baseline, and a one-sided paired Wilcoxon signed-rank test
#            over the gardens below 0.05. level_2_holm repeats it after Holm over the sets scored.
#   level 3  The same test against the mean of the set's random draws, which are runs named
#            random.<set>.<draw>, for instance random.gea.emmax.7. With no draws it stays NA.
# Table S12's own means (0.181 climate distance, 0.340 offset) are references, not thresholds.
#
# Dataset S9 has no founder column: its rows are Dataset S5's by position. Nothing states that
# except the result, since misaligned rows would give r near 0 against Table S12. So --reproduce
# stops unless it returns the design run's mean 0.175 and r 0.90 on the same 302 replicates.
#
# Reads data/arabidopsis1135/ and <project>_results/Maladaptation/tables/. Writes
# ANALYSIS/Arabidopsis1135/offset/: baseline_reproduction.tsv, and per project summary.tsv,
# garden_means.tsv, replicate_rho.tsv, coverage.tsv. --drop-gardens writes to a subfolder.
#
# Run from the repo root, in the image.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})
setDTthreads(4L)

excluded_garden <- 33L             # one replicate, left out by GrENE-net as well
n_replicates_all <- 325L
alpha <- 0.05
offset_methods <- c("gradient_forest", "geometric_offset", "rda_offset")
draw_pattern <- "^random\\.(.+)\\.[0-9]+$"

data_dir <- here("data", "arabidopsis1135")
thin_dir <- file.path(data_dir, "sources", "thin")
grene_dir <- file.path(data_dir, "sources", "dl_grene")
out_root <- here("ANALYSIS", "Arabidopsis1135", "offset")

# Inputs --------------------------------------------------------------------------------------

# p1/p0: founders in rows, in Dataset S5's order, replicates in columns, named <garden>_1_<plot>.
read_fitness <- function() {
  founders <- fread(file.path(grene_dir, "ds", "Dataset_S5_founder_accession_frequency.txt"),
                    header = FALSE, col.names = c("founder", "p0"), colClasses = c(V1 = "character"))
  frequencies <- fread(file.path(grene_dir, "ds", "Dataset_S9_evovled_accession_frequency.txt"))
  columns <- grep("^[0-9]+_1_[0-9]+$", names(frequencies), value = TRUE)
  columns <- columns[!startsWith(columns, paste0(excluded_garden, "_"))]
  stopifnot(nrow(frequencies) == nrow(founders), length(columns) == n_replicates_all)
  fitness <- as.matrix(frequencies[, ..columns]) / founders$p0
  rownames(fitness) <- paste0("S", founders$founder)
  stopifnot(!anyNA(fitness))
  fitness
}

read_table_s12 <- function() {
  s12 <- fread(file.path(grene_dir, "Supplemental_tables", "Table_S12_LOO_accuracy.csv"))
  s12 <- s12[source %in% c("climate_distance", "genomic offset"),
             .(column = sprintf("%d_1_%d", site, plot), source, sp_r)]
  s12 <- dcast(s12, column ~ source, value.var = "sp_r")
  setnames(s12, c("climate_distance", "genomic offset"), c("s12_climate", "s12_offset"))
  s12
}

# Home BIO1 of the founders in the project, from the table the offsets use.
read_home_climate <- function() {
  founders <- fread(file.path(data_dir, "founder_sites.tsv"), select = c("sample", "in_metadata", "site"))
  climate <- fread(file.path(data_dir, "climate_present.tsv"), select = c("site", "bio_1"))
  founders[in_metadata == TRUE][climate, on = "site", nomatch = NULL][, .(sample, bio_home = bio_1)]
}

read_inputs <- function() {
  fitness <- read_fitness()
  replicates <- data.table(column = colnames(fitness),
                           garden = as.integer(sub("_.*", "", colnames(fitness))),
                           plot = as.integer(sub(".*_", "", colnames(fitness))))
  s12 <- read_table_s12()
  gardens <- fread(file.path(data_dir, "garden_sites.tsv"), select = c("garden", "bio_1"))
  setnames(gardens, "bio_1", "bio_garden")
  stopifnot(setequal(s12$column, replicates$column), setequal(gardens$garden, replicates$garden))
  list(fitness = fitness, replicates = replicates, s12 = s12, gardens = gardens,
       home = read_home_climate(),
       samples = fread(file.path(data_dir, "Arabidopsis1135_metadata.tsv"), select = "sample")$sample)
}

# Scoring -------------------------------------------------------------------------------------

# rho of each given replicate, over the samples in the given order. predictor matches samples.
garden_rho <- function(fitness, columns, samples, predictor) {
  as.vector(cor(fitness[samples, columns, drop = FALSE], predictor, method = "spearman"))
}

# One-sided: does the offset beat what it is paired with? Exact for 30 untied gardens; a garden
# with a tied or zero difference moves it to the normal approximation, with a warning we expect.
paired_p <- function(differences) {
  if (all(differences == 0)) return(1)
  suppressWarnings(wilcox.test(differences, alternative = "greater"))$p.value
}

spread <- function(rho, prefix) {
  quartiles <- quantile(rho, c(0.25, 0.75), names = FALSE)
  setNames(list(mean(rho), median(rho), quartiles[1], quartiles[2]),
           paste0(prefix, c("_mean", "_median", "_q1", "_q3")))
}

discover_runs <- function(results) {
  paths <- Sys.glob(file.path(results, "Maladaptation", "tables", "*", "*", "garden_*",
                              "genetic_offset_site.tsv"))
  if (length(paths) == 0L) stop("No garden offsets under ", results)
  run_dirs <- dirname(dirname(paths))
  data.table(method = basename(dirname(run_dirs)), set = sub("_nospatial$", "", basename(run_dirs)),
             garden = as.integer(sub("^garden_", "", basename(dirname(paths)))), path = paths)
}

# One method and set: rho of the offset and of its paired baseline, per replicate. Offsets join the
# founders by sample id, so an id in another spelling (the 0_ of the LD-pruned VCFs) would lose every
# founder without an error: each row must name a sample of the metadata.
score_run <- function(files, inputs, replicates) {
  fitness <- inputs$fitness
  offsets <- rbindlist(lapply(seq_len(nrow(files)), function(i) {
    offset_table <- fread(files$path[i], select = c("sample", "genetic_offset"),
                          colClasses = c(sample = "character"))
    unknown <- setdiff(offset_table$sample, inputs$samples)
    if (length(unknown) > 0L) {
      stop(length(unknown), " sample ids not in the metadata (", unknown[1], ") in ", files$path[i])
    }
    stopifnot(!anyDuplicated(offset_table$sample))
    offset_table[is.finite(genetic_offset)][, garden := files$garden[i]][]
  }))
  founders <- offsets[inputs$home, on = "sample", nomatch = NULL]
  founders <- founders[inputs$gardens, on = "garden", nomatch = NULL]
  if (nrow(founders) == 0L) stop("No founder to score in ", dirname(dirname(files$path[1])))
  founders[, climate_distance := -(bio_home - bio_garden)^2]
  columns_of <- split(replicates$column, replicates$garden)
  founders[, {
    columns <- columns_of[[as.character(.BY$garden)]]
    .(column = columns, n_founders = .N,
      rho_offset = garden_rho(fitness, columns, sample, -genetic_offset),
      rho_climate = garden_rho(fitness, columns, sample, climate_distance))
  }, by = garden]
}

# Every stored set against every method. A set the pipeline used but with no offsets is reported;
# a run missing some gardens would make its gardens unpaired, so that stops.
check_coverage <- function(results, files, n_gardens) {
  stored <- fread(file.path(results, "Maladaptation", "tables", "snp_sets_summary.tsv"),
                  select = c("set", "n_snps", "used", "reason"))
  stopifnot(all(files$set %in% stored$set))
  found <- files[, .(gardens_found = uniqueN(garden)), by = .(method, set)]
  coverage <- found[CJ(method = offset_methods, set = stored$set), on = .(method, set)]
  coverage <- stored[coverage, on = "set"]
  coverage[is.na(gardens_found), gardens_found := 0L]

  partial <- coverage[gardens_found > 0L & gardens_found != n_gardens]
  if (nrow(partial) > 0L) {
    stop("Runs without every garden: ", partial[, paste(method, set, gardens_found, collapse = "; ")])
  }
  missing <- coverage[used == TRUE & gardens_found == 0L]
  if (nrow(missing) > 0L) {
    message("MISSING, used by the pipeline but no offsets: ",
            missing[, paste(method, set, collapse = "; ")])
  }
  unused <- stored[used == FALSE]
  if (nrow(unused) > 0L) {
    message("Sets the pipeline left out: ",
            unused[, paste0(set, " (", reason, ", ", n_snps, " SNPs)", collapse = "; ")])
  }
  coverage
}

# Paired test of each real set against the mean of its random draws, per garden. Empty without draws.
compare_to_draws <- function(per_replicate, garden_means) {
  draw_means <- per_replicate[grepl(draw_pattern, set), .(rho_mean = mean(rho_offset)), by = .(method, set)]
  draw_means[, target := sub(draw_pattern, "\\1", set)]
  draw_gardens <- garden_means[grepl(draw_pattern, set), .(garden, method, rho_offset,
                                                           target = sub(draw_pattern, "\\1", set))]
  draw_gardens <- draw_gardens[, .(rho_random = mean(rho_offset)), by = .(method, target, garden)]
  real_gardens <- garden_means[!grepl(draw_pattern, set), .(method, target = set, garden, rho_offset)]
  paired <- draw_gardens[real_gardens, on = .(method, target, garden), nomatch = NULL]
  tests <- paired[, .(gardens_above_random = sum(rho_offset > rho_random),
                      random_p = paired_p(rho_offset - rho_random)), by = .(method, target)]
  draws <- draw_means[, .(random_draws = .N, random_mean = mean(rho_mean),
                          random_q95 = quantile(rho_mean, 0.95, names = FALSE)),
                      by = .(method, target)]
  merge(draws, tests, by = c("method", "target"))
}

evaluate_tree <- function(results, out_dir, inputs, drop_gardens = integer()) {
  replicates <- inputs$replicates[!garden %in% drop_gardens]
  files <- discover_runs(results)[!garden %in% drop_gardens]
  coverage <- check_coverage(results, files, uniqueN(replicates$garden))

  per_replicate <- files[, score_run(.SD, inputs, replicates), by = .(method, set)]
  if (anyNA(per_replicate$rho_offset) || anyNA(per_replicate$rho_climate)) {
    stop("No rho (a constant offset?) in: ",
         per_replicate[is.na(rho_offset) | is.na(rho_climate), paste(unique(paste(method, set)), collapse = "; ")])
  }
  per_replicate <- replicates[per_replicate, on = .(column, garden)]
  per_replicate <- inputs$s12[per_replicate, on = "column"][, column := NULL]
  setcolorder(per_replicate, c("method", "set", "garden", "plot", "n_founders", "rho_offset", "rho_climate"))

  garden_means <- per_replicate[, .(n_replicates = .N, n_founders = n_founders[1],
                                    rho_offset = mean(rho_offset), rho_climate = mean(rho_climate),
                                    s12_climate = mean(s12_climate), s12_offset = mean(s12_offset)),
                                by = .(method, set, garden)]
  garden_means[, difference := rho_offset - rho_climate]

  real <- per_replicate[!grepl(draw_pattern, set)]
  run_summary <- real[, c(.(n_founders_min = min(n_founders), n_founders_max = max(n_founders),
                            n_replicates = .N),
                          spread(rho_offset, "rho"), spread(rho_climate, "climate"),
                          .(s12_climate = mean(s12_climate), s12_offset = mean(s12_offset))),
                      by = .(method, set)]
  garden_tests <- garden_means[!grepl(draw_pattern, set),
                               .(n_gardens = .N, gardens_positive = sum(rho_offset > 0),
                                 gardens_above_climate = sum(difference > 0),
                                 climate_p = paired_p(difference)),
                               by = .(method, set)]
  run_summary <- garden_tests[run_summary, on = .(method, set)]
  run_summary[, climate_p_holm := p.adjust(climate_p, "holm")]

  draws <- compare_to_draws(per_replicate, garden_means)
  run_summary <- merge(run_summary, draws, by.x = c("method", "set"), by.y = c("method", "target"),
                       all.x = TRUE, sort = FALSE)
  run_summary[is.na(random_draws), random_draws := 0L]

  run_summary[, `:=`(
    above_s12_climate = rho_mean > s12_climate,
    level_1 = rho_mean > 0 & gardens_positive > n_gardens / 2,
    level_2 = rho_mean > climate_mean & climate_p < alpha,
    level_2_holm = rho_mean > climate_mean & climate_p_holm < alpha,
    level_3 = rho_mean > random_mean & random_p < alpha
  )]
  setorder(run_summary, -rho_mean)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(run_summary, file.path(out_dir, "summary.tsv"), sep = "\t")
  fwrite(garden_means, file.path(out_dir, "garden_means.tsv"), sep = "\t")
  fwrite(per_replicate, file.path(out_dir, "replicate_rho.tsv"), sep = "\t")
  fwrite(coverage, file.path(out_dir, "coverage.tsv"), sep = "\t")
  list(summary = run_summary, garden_means = garden_means, per_replicate = per_replicate,
       coverage = coverage)
}

print_summary <- function(run_summary) {
  shown <- run_summary[, .(method, set, n = n_founders_max, rho = round(rho_mean, 3),
                           q1 = round(rho_q1, 3), q3 = round(rho_q3, 3),
                           climate = round(climate_mean, 3), gardens_positive, gardens_above_climate,
                           p = signif(climate_p, 2), p_holm = signif(climate_p_holm, 2),
                           level_1, level_2, level_2_holm, level_3)]
  print(shown, nrows = nrow(shown))
}

# Baseline reproduction -----------------------------------------------------------------------

# Climate distance over founders with a home BIO1, for gardens with a BIO1.
baseline_rho <- function(inputs, home, gardens) {
  rbindlist(lapply(gardens[!is.na(bio_garden), garden], function(this_garden) {
    columns <- inputs$replicates[garden == this_garden, column]
    distance <- -(home$bio_home - gardens[garden == this_garden, bio_garden])^2
    data.table(column = columns, rho = garden_rho(inputs$fitness, columns, home$sample, distance))
  }))
}

reproduce <- function(inputs) {
  suppressPackageStartupMessages(library(terra))
  bio_1 <- rast(here("data", "wc2.1_30s", "wc2.1_30s_bio_1.tif"))
  bio_at <- function(longitude, latitude) {
    terra::extract(bio_1, vect(cbind(longitude, latitude), crs = "EPSG:4326"), ID = FALSE)[[1]]
  }

  score <- function(step, home, gardens) {
    rho <- inputs$s12[baseline_rho(inputs, home, gardens), on = "column"]
    data.table(step = step, founders = nrow(home), replicates = nrow(rho), mean_rho = mean(rho$rho),
               s12_mean_same_replicates = mean(rho$s12_climate), r_against_s12 = cor(rho$rho, rho$s12_climate))
  }

  # The design run: founders at their rounded coordinates, gardens as printed in Table S1.
  rounded <- fread(file.path(thin_dir, "grene_founders.tsv"), colClasses = c(ecotype = "character"))
  rounded[, sample := paste0("S", ecotype)]
  rounded[, bio_home := bio_at(lon, lat)]
  stopifnot(setequal(rounded$sample, rownames(inputs$fitness)))
  rounded <- rounded[!is.na(bio_home), .(sample, bio_home)]

  printed <- fread(file.path(thin_dir, "grene_sites.tsv"), select = c("site", "lon", "lat"))
  printed <- printed[site %in% inputs$gardens$garden, .(garden = site, lon, lat)]
  stopifnot(nrow(printed) == nrow(inputs$gardens))
  printed[, bio_garden := bio_at(lon, lat)]
  stopifnot(identical(printed[is.na(bio_garden), sort(garden)], c(5L, 12L)))

  # Table S1's two typos fixed: the only change among the 30 gardens.
  fixed <- inputs$gardens[printed[, .(garden, printed_bio = bio_garden)], on = "garden"]
  stopifnot(fixed[!garden %in% c(5L, 12L), all.equal(bio_garden, printed_bio)] == TRUE)

  steps <- rbind(
    score("design: rounded coordinates, 2 gardens missing", rounded, printed),
    score("rounded coordinates, gardens 5 and 12 fixed", rounded, inputs$gardens),
    score("project: founders at their 30s cells", inputs$home, inputs$gardens)
  )
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  fwrite(steps, file.path(out_root, "baseline_reproduction.tsv"), sep = "\t")
  print(steps[, .(step, founders, replicates, mean_rho = round(mean_rho, 3),
                  s12 = round(s12_mean_same_replicates, 3), r = round(r_against_s12, 3))])

  design <- steps[1]
  if (design$replicates != 302L || round(design$mean_rho, 3) != 0.175 ||
      round(design$r_against_s12, 2) != 0.90) {
    stop("The baseline does not reproduce the design run: an input changed, or S9 and S5 are misaligned")
  }
  message("Baseline reproduced: mean 0.175, r 0.90, 302 replicates")
}

# Self-test -----------------------------------------------------------------------------------

# A fake results tree on the real fitness and climate, with answers known in advance.
self_test <- function(inputs) {
  # Six gardens all above: one sign pattern of 2^6, so exactly 1/64.
  stopifnot(isTRUE(all.equal(paired_p(1:6), 1 / 64)), paired_p(rep(0, 6)) == 1, paired_p(-(1:6)) > 0.9)
  toy <- matrix(1:4, ncol = 1, dimnames = list(c("a", "b", "c", "d"), "r"))
  stopifnot(isTRUE(all.equal(garden_rho(toy, "r", c("d", "b", "c", "a"), c(4, 3, 2, 1)), 0.8)))

  results <- file.path(tempfile("offset_eval_"), "Fake_results")
  tables <- file.path(results, "Maladaptation", "tables")
  metadata <- fread(file.path(data_dir, "Arabidopsis1135_metadata.tsv"), select = c("site", "sample"))
  climate <- fread(file.path(data_dir, "climate_present.tsv"), select = c("site", "bio_1"))
  samples <- climate[metadata, on = "site"]
  stopifnot(!anyNA(samples$bio_1), nrow(samples) == 969L)

  columns_of <- split(inputs$replicates$column, inputs$replicates$garden)
  set.seed(42L)
  write_run <- function(method, set, offset_of, dropped = character(), blanked = character(),
                        id_prefix = "") {
    for (this_garden in inputs$gardens$garden) {
      garden_bio <- inputs$gardens[garden == this_garden, bio_garden]
      rows <- samples[!sample %in% dropped]
      rows[, genetic_offset := offset_of(this_garden, garden_bio, sample, bio_1)]
      rows[sample %in% blanked, genetic_offset := NA_real_]
      rows[, sample := paste0(id_prefix, sample)]
      # All 17 digits: fwrite's 15 turn a near-tie in garden 49 into a tie, and rho moves by 8e-5.
      rows[, genetic_offset := sprintf("%.17g", genetic_offset)]
      rows <- rows[sample.int(nrow(rows)), .(site, sample, genetic_offset)]   # order must not matter
      garden_dir <- file.path(tables, method, paste0(set, "_nospatial"), paste0("garden_", this_garden))
      dir.create(garden_dir, recursive = TRUE)
      fwrite(rows, file.path(garden_dir, "genetic_offset_site.tsv"), sep = "\t")
    }
  }

  squared_distance <- function(garden, garden_bio, ids, home_bio) (home_bio - garden_bio)^2
  flipped <- function(garden, garden_bio, ids, home_bio) -(home_bio - garden_bio)^2
  # High offset where the founder did badly in that garden: rho far above any baseline.
  oracle <- function(garden, garden_bio, ids, home_bio) {
    founders <- ids %in% rownames(inputs$fitness)
    offset <- rep(0, length(ids))
    offset[founders] <- -rowMeans(inputs$fitness[ids[founders], columns_of[[as.character(garden)]],
                                                 drop = FALSE])
    offset
  }

  dropped <- inputs$home$sample[1:3]
  blanked <- inputs$home$sample[4:5]
  write_run("gradient_forest", "gea.same", squared_distance, dropped, blanked)
  write_run("rda_offset", "gea.flipped", flipped, dropped, blanked)
  write_run("geometric_offset", "gea.oracle", oracle)
  for (draw in 1:3) write_run("geometric_offset", paste0("random.gea.oracle.", draw), squared_distance)

  write_store <- function(sets) {
    fwrite(rbind(data.table(set = sets, n_snps = 100L, used = TRUE, reason = "used"),
                 data.table(set = "gea.tiny", n_snps = 2L, used = FALSE, reason = "below_min_offset_snps")),
           file.path(tables, "snp_sets_summary.tsv"), sep = "\t")
  }
  sets <- c("gea.same", "gea.flipped", "gea.oracle", paste0("random.gea.oracle.", 1:3))
  write_store(sets)

  out_dir <- file.path(dirname(results), "out")
  full <- evaluate_tree(results, out_dir, inputs)
  print_summary(full$summary)
  n_home <- nrow(inputs$home)
  rows <- full$per_replicate
  stopifnot(nrow(rows) == 6L * n_replicates_all)

  # An offset equal to the squared distance must return the baseline on every replicate, despite
  # shuffled rows, and on the founders it has: five gone, three as rows and two as NA.
  same <- rows[set == "gea.same"]
  stopifnot(all(same$n_founders == n_home - 5L), identical(same$rho_offset, same$rho_climate))
  # The offset is reversed once: flipping it flips every rho.
  flip <- rows[set == "gea.flipped"]
  stopifnot(all(flip$n_founders == n_home - 5L), isTRUE(all.equal(flip$rho_offset, -flip$rho_climate)))
  stopifnot(all(rows[set == "gea.oracle", n_founders] == n_home))

  # Table S12 stays attached to the right replicate: the column key, not row order.
  s12_check <- inputs$s12[column == "1_1_1", s12_climate]
  stopifnot(identical(same[garden == 1L & plot == 1L, s12_climate], s12_check))

  verdict <- full$summary
  same_row <- verdict[set == "gea.same"]
  flip_row <- verdict[set == "gea.flipped"]
  oracle_row <- verdict[set == "gea.oracle"]
  stopifnot(
    nrow(verdict) == 3L, verdict$n_gardens == 30L, verdict$n_replicates == n_replicates_all,
    same_row$climate_p == 1, same_row$gardens_above_climate == 0L, !same_row$level_2,
    is.na(same_row$level_3), same_row$random_draws == 0L,
    flip_row$rho_mean < 0, !flip_row$level_1, !flip_row$level_2,
    oracle_row$level_1, oracle_row$level_2, oracle_row$level_2_holm, oracle_row$level_3,
    oracle_row$random_draws == 3L
  )
  # Draws are scored but never ranked as sets, and the pairing is by garden.
  stopifnot(nrow(full$garden_means[grepl(draw_pattern, set)]) == 3L * 30L)

  covered <- full$coverage
  stopifnot(nrow(covered) == 3L * 7L, covered[method == "geometric_offset" & set == "gea.oracle", gardens_found] == 30L,
            covered[set == "gea.tiny", all(gardens_found == 0L)])

  dropped_gardens <- c(10L, 13L, 26L, 60L)
  subset_run <- evaluate_tree(results, file.path(out_dir, "subset"), inputs, dropped_gardens)
  n_kept <- inputs$replicates[!garden %in% dropped_gardens, .N]
  stopifnot(subset_run$summary$n_gardens == 26L, subset_run$summary$n_replicates == n_kept,
            !any(subset_run$per_replicate$garden %in% dropped_gardens))

  stops_with <- function(pattern) {
    tryCatch({ evaluate_tree(results, out_dir, inputs); FALSE }, error = function(e) {
      grepl(pattern, conditionMessage(e))
    })
  }

  # Sample ids in another spelling would lose every founder without an error.
  write_store(c(sets, "gea.prefixed"))
  write_run("rda_offset", "gea.prefixed", squared_distance, id_prefix = "0_")
  stopifnot(stops_with("not in the metadata"))
  unlink(file.path(tables, "rda_offset", "gea.prefixed_nospatial"), recursive = TRUE)

  # So would a run with no founder left: the join returns nothing and the run would just vanish.
  write_store(c(sets, "gea.empty"))
  write_run("rda_offset", "gea.empty", squared_distance, blanked = inputs$home$sample)
  stopifnot(stops_with("No founder to score"))
  unlink(file.path(tables, "rda_offset", "gea.empty_nospatial"), recursive = TRUE)

  # An incomplete run breaks the pairing.
  write_store(sets)
  unlink(file.path(tables, "rda_offset", "gea.flipped_nospatial", "garden_1"), recursive = TRUE)
  stopifnot(stops_with("without every garden"))

  unlink(dirname(results), recursive = TRUE)
  message("offset_eval.R self-test passed")
}

# Main ----------------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  stop("usage: offset_eval.R <project> [--drop-gardens a,b,c] | --reproduce | --self-test")
}
inputs <- read_inputs()

if (args[1] == "--self-test") {
  self_test(inputs)
} else if (args[1] == "--reproduce") {
  reproduce(inputs)
} else {
  project <- args[1]
  drop_gardens <- integer()
  out_dir <- file.path(out_root, project)
  if (length(args) == 3L && args[2] == "--drop-gardens") {
    drop_gardens <- as.integer(strsplit(args[3], ",", fixed = TRUE)[[1]])
    stopifnot(!anyNA(drop_gardens), all(drop_gardens %in% inputs$gardens$garden))
    out_dir <- file.path(out_dir, paste0("without_", paste(drop_gardens, collapse = "_")))
  } else if (length(args) != 1L) {
    stop("usage: offset_eval.R <project> [--drop-gardens a,b,c] | --reproduce | --self-test")
  }
  outcome <- evaluate_tree(here(paste0(project, "_results")), out_dir, inputs, drop_gardens)
  print_summary(outcome$summary)
  message("Wrote ", out_dir)
}
