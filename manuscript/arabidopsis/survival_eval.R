#!/usr/bin/env Rscript
# survival_eval.R: does the founder population's offset at a GrENE-net garden predict how many of the
# garden's plots survive? A garden-level companion to offset_eval.R, after GrENE-net's Fig. S66.
#
#   Rscript manuscript/arabidopsis/survival_eval.R <project>
#
# Rule, written into ANALYSIS/Arabidopsis1135/paper_notes.md ("Garden survival test: the rule",
# 2026-11-15 17:35 IST) before any survival number existed:
#   unit       The garden: the 30 fitness gardens. Every garden was sown with one founder mix, so a
#              garden's value depends only on its climate.
#   predictor  The founder population's offset at the garden: the founders' offsets averaged with their
#              starting frequency p0 (Dataset S5) as weights, over the founders with a finite offset there.
#              The same for each random draw, and for the climate distance: the p0-weighted mean of
#              (BIO1 home - BIO1 garden)^2, over the founders of the Gradient Forest union run, with BIO1
#              read as offset_eval.R reads it. Higher means worse for all of them.
#   response   Primary: the share of a garden's fitness plots alive in year 3. Fitness plots: the 325 plots
#              of the year-1 evaluation (Dataset S9's <garden>_1_<plot> columns, garden 33 out), the set
#              of Fig. 5E's year-1 count. Alive: a generation-3 row of Table S4 with
#              totalplantnumber_complete > 0. A garden whose Table S1 SURVIVALYEARS is below 3 counts every
#              plot as dead; one with SURVIVALYEARS 3 and no generation-3 row for any of its plots is
#              excluded and listed. Where the two tables disagree (plants in a year Table S1 says the
#              garden did not reach), Table S4's observed plants count and the garden is listed.
#              Correction (attempt 28, after the first run): the first rule divided by the garden's plots
#              with a generation-1 row, which gave shares of 1.43 (garden 48) and 1.08 (garden 52), where
#              Table S4 records later plots with no generation-1 row, and counted extra trays (gardens 2
#              and 28). The first run's tables are kept in offset/<project>/survival_first_run/.
#              Secondary: SURVIVALYEARS itself, and the year-1 and year-2 shares built the same way.
#              SURVIVALYEARS never exceeds 3 here while Table S4 runs to generation 6, so years 4 and 5
#              are described (gardens with rows, plots with plants), not tested.
#   test       One-sided Spearman over gardens, a negative rho expected, ties by the asymptotic test: the
#              Gradient Forest, geometric and RDA offsets of the union set (role primary) and the climate
#              distance (role baseline). Holm over the three offset methods on the primary response. The
#              union's rank among itself and its draws per method (rank 1 = the most negative rho). The
#              two-sided Spearman between each union offset and the climate distance over gardens (role
#              axis), to show whether they are one axis. Sensitivity: every test again without gardens 26
#              and 60, the two outside the fitted range on the offset predictors.
#              Descriptive, outside the Holm family: the other five sets (role descriptive) and every draw
#              (role draw).
#   counts     The paper's counts are written beside ours, not forced: Fig. 5E evaluates 325, 243 and 139
#              plots in years 1, 3 and 5.
#   stop       Fewer than 20 gardens with a year-3 status stops the run, to ask the user; the next source
#              would be Dataset S1 on Zenodo, fetched only with their yes.
#
# Reads <project>_results/Maladaptation/tables/ (offset_eval.R's run discovery), Dataset S5, and GrENE-net
# Tables S1 and S4. Of Table S1 only columns B:C (SITE_CODE, SURVIVALYEARS) are read: the sheet also holds
# the participants' names and emails. Writes ANALYSIS/Arabidopsis1135/offset/<project>/survival/:
#   garden_survival.tsv  one row per garden: climate, Table S1 and S4 status per year, shares, the union
#                        offsets and the climate distance
#   garden_offsets.tsv   the founder-population offset of every run and draw at every garden (long)
#   survival_tests.tsv   one row per predictor x response x garden set: n, rho, one-sided p, Holm, rank
#   counts.tsv           per year: gardens with rows, plots with plants, gardens and plots with a status,
#                        and the paper's count
#
# Run from the repo root, in the image.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(readxl)
})
source(here("manuscript", "arabidopsis", "offset_eval.R"))   # readers only: its main block does not run

usage <- "usage: Rscript manuscript/arabidopsis/survival_eval.R <project>"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || startsWith(args[1], "-")) stop(usage)
project <- args[1]
results <- here(paste0(project, "_results"))
out_dir <- file.path(out_root, project, "survival")
tables_dir <- file.path(grene_dir, "Supplemental_tables")

primary_year <- 3L
tested_years <- 1:3
described_years <- 4:5
min_gardens <- 20L
out_of_range <- c(26L, 60L)
paper_counts <- c(`1` = 325L, `3` = 243L, `5` = 139L)
union_set <- "gea.union"

# Inputs ---------------------------------------------------------------------------------------------

gardens <- fread(file.path(data_dir, "garden_sites.tsv"),
                 select = c("garden", "name", "country", "bio_1", "outside_offset_predictors"))
setnames(gardens, "bio_1", "bio_garden")
fitness_gardens <- sort(gardens$garden)
stopifnot(length(fitness_gardens) == 30L)

home <- read_home_climate()
p0 <- fread(file.path(grene_dir, "ds", "Dataset_S5_founder_accession_frequency.txt"),
            header = FALSE, col.names = c("founder", "p0"), colClasses = c(V1 = "character"))
p0[, sample := paste0("S", founder)]
founders <- home[p0, on = "sample", nomatch = NULL][, .(sample, bio_home, p0)]
stopifnot(nrow(founders) > 0L, all(founders$p0 > 0))
samples <- fread(file.path(data_dir, "Arabidopsis1135_metadata.tsv"), select = "sample")$sample

s4 <- fread(file.path(tables_dir, "Table_S4_population_size_estimations.csv"),
            select = c("site", "plot", "generation", "totalplantnumber_complete", "predicted"))
s4 <- s4[site %in% fitness_gardens]
setnames(s4, "site", "garden")

s1 <- as.data.table(read_excel(file.path(tables_dir, "Table_S1_site_information.xlsx"),
                               range = cell_cols("B:C")))
if (!identical(names(s1), c("SITE_CODE", "SURVIVALYEARS"))) {
  stop("Table S1 columns B:C are not SITE_CODE, SURVIVALYEARS: ", paste(names(s1), collapse = ", "))
}
s1 <- s1[!is.na(SITE_CODE), .(garden = as.integer(SITE_CODE), survival_years = as.integer(SURVIVALYEARS))]
s1 <- s1[garden %in% fitness_gardens]
if (!setequal(s1$garden, fitness_gardens) || anyDuplicated(s1$garden)) {
  stop("Table S1 does not give one SURVIVALYEARS per fitness garden")
}

# Survival per garden and year -----------------------------------------------------------------------

s9_columns <- names(fread(file.path(grene_dir, "ds", "Dataset_S9_evovled_accession_frequency.txt"), nrows = 0L))
s9_columns <- grep("^[0-9]+_1_[0-9]+$", s9_columns, value = TRUE)
fitness_plots <- data.table(garden = as.integer(sub("_.*", "", s9_columns)),
                            plot = as.integer(sub(".*_", "", s9_columns)))[garden %in% fitness_gardens]
if (nrow(fitness_plots) != n_replicates_all || anyDuplicated(fitness_plots)) {
  stop("expected the ", n_replicates_all, " fitness plots of the year-1 evaluation")
}
started <- fitness_plots[, .(plots_started = .N), by = garden]
is_alive <- function(n) !is.na(n) & n > 0
fitness_rows <- s4[fitness_plots, on = .(garden, plot), nomatch = NULL]
if (fitness_rows[generation == 1L & is_alive(totalplantnumber_complete), uniqueN(paste(garden, plot))] !=
    n_replicates_all) {
  message("NOTE: not every fitness plot has a generation-1 row with plants in Table S4")
}

survival <- CJ(garden = fitness_gardens, year = c(tested_years, described_years))
reporting <- s4[, .(rows = .N, plots_with_rows = uniqueN(plot),
                    plots_with_plants = uniqueN(plot[is_alive(totalplantnumber_complete)]),
                    na_rows = sum(is.na(totalplantnumber_complete)),
                    predicted_rows = sum(predicted %in% c(TRUE, "True", "TRUE"))),
                by = .(garden, year = generation)]
alive_fitness <- fitness_rows[is_alive(totalplantnumber_complete), .(alive = uniqueN(plot)),
                              by = .(garden, year = generation)]
survival <- reporting[survival, on = .(garden, year)]
survival <- alive_fitness[survival, on = .(garden, year)]
for (column in c("rows", "plots_with_rows", "plots_with_plants", "na_rows", "predicted_rows", "alive")) {
  set(survival, which(is.na(survival[[column]])), column, 0L)
}
survival <- started[survival, on = "garden"]
survival <- s1[survival, on = "garden"]
if (anyNA(survival$plots_started)) stop("a fitness garden has no fitness plot")

survival[, status := fifelse(year > max(tested_years), "described",
                     fifelse(survival_years < year & alive > 0L, "tables_disagree",
                     fifelse(survival_years < year, "dead_by_table_s1",
                     fifelse(rows == 0L, "excluded_no_rows", "observed"))))]
survival[, share := fifelse(status %in% c("observed", "tables_disagree"), alive / plots_started,
                    fifelse(status == "dead_by_table_s1", 0, NA_real_))]

known_primary <- survival[year == primary_year & !is.na(share), uniqueN(garden)]
if (known_primary < min_gardens) {
  stop(sprintf("only %d gardens have a year-%d status (rule: at least %d); stop and ask the user",
               known_primary, primary_year, min_gardens))
}

# Founder-population offsets ---------------------------------------------------------------------------

files <- discover_runs(results)
files <- files[garden %in% fitness_gardens]
garden_offsets <- files[, {
  offset_table <- fread(path, select = c("sample", "genetic_offset"), colClasses = c(sample = "character"))
  unknown <- setdiff(offset_table$sample, samples)
  if (length(unknown) > 0L) stop(length(unknown), " sample ids not in the metadata in ", path)
  scored <- offset_table[is.finite(genetic_offset)][founders, on = "sample", nomatch = NULL]
  .(n_founders = nrow(scored), offset = sum(scored$p0 * scored$genetic_offset) / sum(scored$p0),
    founder_ids = list(sort(scored$sample)))
}, by = .(method, set, garden)]
runs <- unique(garden_offsets[, .(method, set)])
if (garden_offsets[, uniqueN(garden), by = .(method, set)][, any(V1 != 30L)]) {
  stop("a run lacks some of the 30 gardens")
}

reference <- garden_offsets[method == "gradient_forest" & set == union_set]
if (nrow(reference) != 30L) stop("no Gradient Forest union run with 30 gardens")
reference_ids <- setNames(reference$founder_ids, reference$garden)
same_founders <- garden_offsets[, .(same = identical(founder_ids[[1]], reference_ids[[as.character(garden)]])),
                                by = .(method, set, garden)][, all(same)]
if (!same_founders) message("NOTE: not every run scores the same founders as the Gradient Forest union run")

climate <- reference[, {
  ids <- founder_ids[[1]]
  f <- founders[sample %in% ids]
  g <- gardens[garden == .BY$garden, bio_garden]
  .(climate_distance = sum(f$p0 * (f$bio_home - g)^2) / sum(f$p0), n_founders_climate = nrow(f))
}, by = garden]

# Tests ----------------------------------------------------------------------------------------------

spearman_less <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", alternative = "less", exact = FALSE))
  list(n = sum(ok), rho = unname(test$estimate), p = test$p.value)
}

responses <- rbind(
  survival[year %in% tested_years, .(garden, response = sprintf("share_alive_year_%d", year), value = share)],
  unique(survival[, .(garden, response = "survival_years", value = as.numeric(survival_years))]))
primary_response <- sprintf("share_alive_year_%d", primary_year)

predictors <- rbind(
  garden_offsets[, .(predictor = sprintf("ours.%s.%s", method, set), method, set,
                     role = fifelse(grepl(draw_pattern, set), "draw",
                            fifelse(set == union_set, "primary", "descriptive")),
                     garden, x = offset)],
  climate[, .(predictor = "ours.climate_distance", method = NA_character_, set = NA_character_,
              role = "baseline", garden, x = climate_distance)])

garden_sets <- list(all = fitness_gardens, without_26_60 = setdiff(fitness_gardens, out_of_range))
tests <- rbindlist(lapply(names(garden_sets), function(gs) {
  keep <- garden_sets[[gs]]
  joined <- predictors[garden %in% keep][responses[garden %in% keep], on = "garden", allow.cartesian = TRUE]
  joined[, c(list(gardens = gs), spearman_less(x, value)), by = .(predictor, method, set, role, response)]
}))
setnames(tests, c("rho", "p"), c("rho", "p_one_sided"))

tests[, p_holm := NA_real_]
tests[role == "primary" & response == primary_response, p_holm := p.adjust(p_one_sided, "holm"), by = gardens]

tests[, union_rank := NA_integer_]
tests[, n_compared := NA_integer_]
for (m in unique(na.omit(tests$method))) {
  pool <- tests[method == m & (set == union_set | grepl(paste0("^random\\.", union_set, "\\.[0-9]+$"), set))]
  ranks <- pool[, .(union_rank = if (any(set == union_set) && !all(is.na(rho)))
                      as.integer(rank(rho, ties.method = "min", na.last = "keep")[set == union_set])
                    else NA_integer_,
                    n_compared = sum(!is.na(rho))), by = .(response, gardens)]
  tests[ranks, on = .(response, gardens), `:=`(union_rank = fifelse(method == m & set == union_set,
                                                                    i.union_rank, union_rank),
                                               n_compared = fifelse(method == m & set == union_set,
                                                                    i.n_compared, n_compared))]
}

axis <- rbindlist(lapply(names(garden_sets), function(gs) {
  keep <- garden_sets[[gs]]
  union_runs <- garden_offsets[set == union_set & garden %in% keep]
  union_runs[climate, on = "garden", nomatch = NULL][, {
    test <- suppressWarnings(cor.test(offset, climate_distance, method = "spearman", exact = FALSE))
    .(predictor = sprintf("ours.%s.%s", .BY$method, union_set), set = union_set, role = "axis",
      response = "garden_climate_distance", gardens = gs, n = .N, rho = unname(test$estimate),
      p_two_sided = test$p.value)
  }, by = method]
}))
tests <- rbind(tests, axis, fill = TRUE)
tests[, project := project]
setcolorder(tests, c("project", "predictor", "method", "set", "role", "response", "gardens", "n", "rho",
                     "p_one_sided", "p_two_sided", "p_holm", "union_rank", "n_compared"))
setorder(tests, gardens, response, role, predictor)

# Counts beside the paper's -------------------------------------------------------------------------

counts <- survival[, .(gardens_with_rows = sum(rows > 0L), plots_with_plants = sum(plots_with_plants),
                       fitness_plots_alive = sum(alive),
                       gardens_with_status = sum(!is.na(share)),
                       plots_evaluated = sum(plots_started[!is.na(share)]),
                       gardens_dead_by_table_s1 = sum(status == "dead_by_table_s1"),
                       gardens_excluded = sum(status == "excluded_no_rows"),
                       gardens_tables_disagree = sum(status == "tables_disagree")), by = year]
counts[year > max(tested_years), `:=`(gardens_with_status = NA_integer_, plots_evaluated = NA_integer_)]
counts[, paper_plots_evaluated := paper_counts[as.character(year)]]
counts[, project := project]
setcolorder(counts, "project")

# Write ----------------------------------------------------------------------------------------------

wide <- dcast(survival[year %in% c(tested_years, described_years)], garden ~ year,
              value.var = c("rows", "plots_with_plants", "alive", "status", "share"))
union_wide <- dcast(garden_offsets[set == union_set], garden ~ method, value.var = "offset")
setnames(union_wide, setdiff(names(union_wide), "garden"),
         paste0("offset_union_", setdiff(names(union_wide), "garden")))
garden_survival <- Reduce(function(a, b) merge(a, b, by = "garden", all.x = TRUE),
                          list(gardens, s1, started, wide, climate, union_wide,
                               reference[, .(garden, n_founders)]))
garden_survival[, project := project]
setcolorder(garden_survival, "project")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(garden_survival, file.path(out_dir, "garden_survival.tsv"), sep = "\t")
fwrite(garden_offsets[, .(project = project, method, set, garden, n_founders, offset)],
       file.path(out_dir, "garden_offsets.tsv"), sep = "\t")
fwrite(tests, file.path(out_dir, "survival_tests.tsv"), sep = "\t")
fwrite(counts, file.path(out_dir, "counts.tsv"), sep = "\t")

# Report ---------------------------------------------------------------------------------------------

listed <- survival[status %in% c("excluded_no_rows", "tables_disagree")]
if (nrow(listed) > 0L) {
  message("Listed gardens: ", listed[, paste0(garden, " (year ", year, ", ", status, ")", collapse = "; ")])
}
message(sprintf("%s: %d runs at 30 gardens, %d gardens with a year-%d status, %d founders per garden in the reference run",
                project, nrow(runs), known_primary, primary_year, reference$n_founders[1]))
print(counts[, .(year, gardens_with_rows, plots_with_plants, fitness_plots_alive, gardens_with_status, plots_evaluated,
                 paper_plots_evaluated)])
print(tests[role %in% c("primary", "baseline", "axis") & response %in% c(primary_response, "survival_years",
                                                                          "garden_climate_distance"),
            .(gardens, response, predictor, n, rho = round(rho, 3), p_one_sided = signif(p_one_sided, 3),
              p_two_sided = signif(p_two_sided, 3), p_holm = signif(p_holm, 3), union_rank, n_compared)])
message("Wrote ", out_dir)
