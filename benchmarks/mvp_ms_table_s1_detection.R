#!/usr/bin/env Rscript
# Table S1 (detection), built from Data S2 only.
#   (a) per candidate set, median (min-max) over the 600 SS-Clines replicates: SNPs identified,
#       true positives (causal or linked SNPs), false positives (neutral SNPs), precision,
#       causal loci identified (count, and % of the replicate's causal loci), linked SNPs.
#   (b) the replicates: tested SNPs, causal loci by genetic architecture, causal loci as % of
#       tested SNPs.
# Values are those of Figure 2a (0.25% cutoff, 5 kb window for the 2/3 and 3/3 sets). An empty
# set (the 3/3 set in 2 replicates) counts as 0 in every column except precision, which is
# undefined there and left out (as in detection_per_seed.tsv, where it is blank); column n says
# how many replicates each value is taken over.
#
# Usage: Rscript benchmarks/mvp_ms_table_s1_detection.R <Data_S2.zip> <out_dir>
# Writes <out_dir>/Table_S1_detection.tsv (long) and <out_dir>/Table_S1_detection.md.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: mvp_ms_table_s1_detection.R <Data_S2.zip> <out_dir>")
zip <- args[1]
out_dir <- args[2]
if (!file.exists(zip)) stop("Data S2 not found: ", zip)
if (!dir.exists(out_dir)) stop("output directory not found: ", out_dir)

read_member <- function(member) {
  fread(cmd = paste("unzip -p", shQuote(zip), shQuote(member)), sep = "\t")
}
det <- read_member("figures_ssclines_gallery13/detection_per_seed.tsv")
tested <- read_member("remeasure600/rdaunc/rank_metrics.tsv")[method == "LFMM", .(seed, n_snps)]

sets <- c("LFMM" = "LFMM", "EMMAX" = "EMMAX", "RDA" = "RDA",
          "1/3 methods" = "1/3 set", "2/3 methods" = "2/3 set", "3/3 methods" = "3/3 set")
det <- det[label %in% names(sets)]
det[, set := factor(sets[label], levels = unname(sets))]
counts <- det[, .N, by = set]
if (nrow(counts) != 6 || any(counts$N != 600)) stop("expected 600 replicates for each of 6 candidate sets")

# ---- (a) per candidate set ----
metrics_a <- c(n = "SNPs identified", TP = "True positives (causal or linked SNPs)",
               FP = "False positives (neutral SNPs)", prec_pct = "Precision, %",
               n_causal = "Causal loci identified", rec_pct = "Causal loci identified, % of the replicate's",
               n_linked = "Linked SNPs identified")
det[, `:=`(prec_pct = 100 * prec_cl, rec_pct = 100 * rec_c)]
det[, (names(metrics_a)) := lapply(.SD, as.numeric), .SDcols = names(metrics_a)]
long_a <- melt(det, id.vars = c("seed", "set"), measure.vars = names(metrics_a),
               variable.name = "metric", variable.factor = FALSE)
long_a <- long_a[!is.na(value)]
sum_a <- long_a[, .(median = median(value), min = min(value), max = max(value), n = .N),
                by = .(set, metric)]
setorder(sum_a, set)
sum_a[, `:=`(part = "a", stratum = as.character(set), label = unname(metrics_a[metric]))]

# ---- (b) the replicates ----
reps <- unique(det[, .(seed, arch, causal_total)])
if (nrow(reps) != 600 || anyDuplicated(reps$seed)) stop("expected one causal-locus count per replicate")
reps <- merge(reps, tested, by = "seed")
if (nrow(reps) != 600) stop("tested-SNP counts missing for some replicates")
reps[, causal_pct := 100 * causal_total / n_snps]
arch_levels <- c("oligogenic", "moderately polygenic", "highly polygenic")
if (!setequal(unique(reps$arch), arch_levels)) stop("unexpected genetic architectures: ",
                                                    paste(unique(reps$arch), collapse = ", "))
summ <- function(x) list(median = median(x), min = min(x), max = max(x), n = length(x))
sum_b <- rbind(
  reps[, c(stratum = "all", metric = "n_snps", label = "Tested SNPs", summ(n_snps))],
  reps[, c(stratum = "all", metric = "causal_total", label = "Causal loci", summ(causal_total))],
  reps[, c(metric = "causal_total", label = "Causal loci", summ(causal_total)), by = .(stratum = arch)],
  reps[, c(stratum = "all", metric = "causal_pct", label = "Causal loci, % of tested SNPs", summ(causal_pct))],
  use.names = TRUE)
sum_b[, part := "b"]

# ---- long TSV ----
cols <- c("part", "stratum", "metric", "label", "median", "min", "max", "n")
out <- rbind(sum_a[, ..cols], sum_b[, ..cols])
fwrite(out, file.path(out_dir, "Table_S1_detection.tsv"), sep = "\t")

# ---- markdown ----
fmt_count <- function(x) ifelse(x == round(x), format(round(x), big.mark = ",", trim = TRUE),
                                 formatC(x, format = "f", digits = 1, big.mark = ","))
fmt <- function(med, lo, hi, kind) {
  f <- switch(kind,
              count = fmt_count,
              pct1 = function(x) formatC(x, format = "f", digits = 1),
              pct2 = function(x) formatC(x, format = "f", digits = 2),
              int = function(x) format(round(x), big.mark = ",", trim = TRUE))
  sprintf("%s (%s–%s)", f(med), f(lo), f(hi))
}
kind_a <- c(n = "count", TP = "count", FP = "count", prec_pct = "pct1",
            n_causal = "count", rec_pct = "pct1", n_linked = "count")
sum_a[, cell := fmt(median, min, max, kind_a[metric]), by = seq_len(nrow(sum_a))]
wide_a <- dcast(sum_a, set ~ factor(metric, levels = names(metrics_a)), value.var = "cell")

md_row <- function(x) paste0("| ", paste(x, collapse = " | "), " |")
lines <- c(
  "# Table S1 — Detection by the six candidate sets over 600 simulated replicates",
  "",
  "SS-Clines, 600 replicates. Median (minimum–maximum) over replicates. Values: `Table_S1_detection.tsv` (long format); source in Data S2 (`detection_per_seed.tsv`, `rank_metrics.tsv`); built by `benchmarks/mvp_ms_table_s1_detection.R`.",
  "",
  "## (a) Detection by candidate set",
  "",
  md_row(c("Candidate set", unname(metrics_a))),
  md_row(rep("---", length(metrics_a) + 1)),
  vapply(seq_len(nrow(wide_a)), function(i) md_row(c(as.character(wide_a$set[i]),
                                                      unlist(wide_a[i, -1]))), character(1)),
  "",
  sprintf("Every value is taken over %s replicates, except the precision of the 3/3 set (%s replicates): the set is empty in %s replicates, where precision is undefined and every other column counts 0.",
          max(sum_a$n), sum_a[set == "3/3 set" & metric == "prec_pct", n],
          max(sum_a$n) - sum_a[set == "3/3 set" & metric == "prec_pct", n]),
  "",
  "## (b) The replicates",
  "",
  md_row(c("Quantity", "Replicates", "Median (minimum–maximum)")),
  md_row(rep("---", 3))
)
row_b <- function(label, d, kind) md_row(c(label, d$n, fmt(d$median, d$min, d$max, kind)))
lines <- c(lines,
  row_b("Tested SNPs", sum_b[metric == "n_snps"], "int"),
  row_b("Causal loci, all replicates", sum_b[metric == "causal_total" & stratum == "all"], "int"),
  vapply(arch_levels, function(a) row_b(paste0("Causal loci, ", a), sum_b[metric == "causal_total" & stratum == a], "int"),
         character(1)),
  row_b("Causal loci, % of tested SNPs", sum_b[metric == "causal_pct"], "pct2")
)
writeLines(lines, file.path(out_dir, "Table_S1_detection.md"))
message("wrote Table_S1_detection.{tsv,md} to ", out_dir)
