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

# =============================================================================
# [added 2026-10-25, supplement build A2] Three slug tables, written beside Table_S1_detection.*,
# which stay byte-identical (this block only appends; nothing above changed). Same long format as
# Table_S1_detection.tsv without its `part` and `metric` columns: a delivered table carries display
# labels only (`label`; the data keys stay in Table_S1_detection.tsv). Markdown titled by slug, no
# table number, no caption text.
#   t_detection       = part (a)
#   t_detection_arch  = the quantities of (a), separately per genetic architecture (200 replicates
#                       each; column n says where an empty set leaves precision undefined)
#   t_replicates      = part (b)
# Asserts: both slug TSVs read back equal the written Table_S1_detection.tsv parts, value for value;
# the pooled, per-architecture and replicate values the supplement text quotes.
# =============================================================================
sum_arch <- long_a[, .(median = median(value), min = min(value), max = max(value), n = .N),
                   by = .(set, metric, arch = det$arch[match(paste(seed, set), paste(det$seed, det$set))])]
if (anyNA(sum_arch$arch) || !setequal(unique(sum_arch$arch), arch_levels)) stop("t_detection_arch: architecture lookup failed")
sum_arch[, arch := factor(arch, levels = arch_levels)]
setorder(sum_arch, arch, set)
sum_arch[, `:=`(stratum = as.character(set), label = unname(metrics_a[metric]))]
if (sum_arch[metric != "prec_pct", any(n != 200)]) stop("t_detection_arch: expected 200 replicates per architecture")

cols_a <- c("stratum", "label", "median", "min", "max", "n")
fwrite(setnames(sum_a[, ..cols_a], "stratum", "set"), file.path(out_dir, "t_detection.tsv"), sep = "\t")
fwrite(setnames(sum_arch[, c("arch", cols_a), with = FALSE], c("arch", "stratum"), c("architecture", "set")),
       file.path(out_dir, "t_detection_arch.tsv"), sep = "\t")
fwrite(sum_b[, setdiff(cols, c("part", "metric")), with = FALSE], file.path(out_dir, "t_replicates.tsv"), sep = "\t")

# ---- asserts on the files as written ----
s1 <- fread(file.path(out_dir, "Table_S1_detection.tsv"), sep = "\t")
td <- fread(file.path(out_dir, "t_detection.tsv"), sep = "\t")
tr <- fread(file.path(out_dir, "t_replicates.tsv"), sep = "\t")
ta <- fread(file.path(out_dir, "t_detection_arch.tsv"), sep = "\t")
same <- function(x, y) nrow(x) == nrow(y) && isTRUE(all.equal(as.data.frame(x), as.data.frame(y), tolerance = 0, check.attributes = FALSE))
if (!same(td, setnames(s1[part == "a", -c("part", "metric")], "stratum", "set"))) stop("t_detection differs from Table_S1_detection.tsv part (a)")
if (!same(tr, s1[part == "b", -c("part", "metric")])) stop("t_replicates differs from Table_S1_detection.tsv part (b)")
# the asserts below address rows by data key: re-attach it from the in-memory tables, in written order
stopifnot(identical(td$label, sum_a$label), identical(tr$label, sum_b$label))
td[, metric := sum_a$metric]; tr[, metric := sum_b$metric]
near <- function(x, want, tol) length(x) == 1L && abs(x - want) < tol
v <- function(D, s, m, col = "median") D[set == s & metric == m][[col]]
chk <- c(
  lfmm_n       = near(v(td, "LFMM", "n"), 127, 1e-9) && near(v(td, "LFMM", "n", "min"), 26, 1e-9) && near(v(td, "LFMM", "n", "max"), 162, 1e-9),
  lfmm_tp_fp   = near(v(td, "LFMM", "TP"), 109, 1e-9) && near(v(td, "LFMM", "FP"), 16.5, 1e-9),
  lfmm_prec    = near(v(td, "LFMM", "prec_pct"), 85.2, 0.05),
  lfmm_rec     = near(v(td, "LFMM", "rec_pct"), 22.3, 0.05),
  emmax        = near(v(td, "EMMAX", "TP"), 71, 1e-9) && near(v(td, "EMMAX", "FP"), 62, 1e-9) && near(v(td, "EMMAX", "prec_pct"), 54.4, 0.05),
  set13        = near(v(td, "1/3 set", "TP"), 213, 1e-9) && near(v(td, "1/3 set", "FP"), 89, 1e-9) && near(v(td, "1/3 set", "rec_pct"), 25.9, 0.05),
  set23        = near(v(td, "2/3 set", "TP"), 178, 1e-9) && near(v(td, "2/3 set", "FP"), 36, 1e-9),
  set33_prec_n = near(v(td, "3/3 set", "prec_pct", "n"), 598, 0.5),
  # the median is 25,321.5 (600 replicates, even n); the text quotes the cell as this script prints it
  tested       = tr[metric == "n_snps", fmt(median, min, max, "int")] == "25,322 (5,141–32,130)",
  causal_arch  = all(mapply(function(a, md, lo, hi) { x <- tr[metric == "causal_total" & stratum == a]
                   near(x$median, md, 1e-9) && near(x$min, lo, 1e-9) && near(x$max, hi, 1e-9) },
                   arch_levels, c(7, 43, 460), c(3, 17, 231), c(14, 129, 771))),
  causal_pct   = near(tr[metric == "causal_pct", median], 0.25, 0.005))
want_arch <- data.table(architecture = rep(arch_levels, each = 3), set = rep(c("1/3 set", "2/3 set", "LFMM"), 3),
                        want = c(85.71, 83.33, 75.00, 25.90, 24.53, 22.29, 3.56, 3.15, 2.28))
got_arch <- merge(want_arch, ta[label == metrics_a[["rec_pct"]], .(architecture, set, median)], by = c("architecture", "set"))
chk["arch_rec"] <- nrow(got_arch) == 9L && all(abs(got_arch$median - got_arch$want) < 0.005)
chk["arch_n"] <- ta[label != metrics_a[["prec_pct"]], all(n == 200)] && ta[label == metrics_a[["prec_pct"]], sum(200 - n)] == 2L
print(got_arch[, .(architecture, set, median, want)])
if (!all(chk)) stop("A2 assertion failed: ", paste(names(chk)[!chk], collapse = ", "))
message("A2 asserts passed: ", paste(names(chk), collapse = ", "))

# ---- markdown, one file per slug ----
md_table_a <- function(S, first_col, first_vals) {
  S <- copy(S)
  S[, cell := fmt(median, min, max, kind_a[metric]), by = seq_len(nrow(S))]
  W <- dcast(S, set ~ factor(metric, levels = names(metrics_a)), value.var = "cell")
  c(md_row(c(first_col, unname(metrics_a))), md_row(rep("---", length(metrics_a) + 1)),
    vapply(seq_len(nrow(W)), function(i) md_row(c(as.character(W$set[i]), unlist(W[i, -1]))), character(1)))
}
prec_note <- function(S, what) {
  e <- S[metric == "prec_pct" & n < max(S$n)]
  if (!nrow(e)) return(sprintf("Every value is taken over %s replicates%s.", max(S$n), what))
  sprintf("Every value is taken over %s replicates%s, except the precision of the %s (%s replicates): the set is empty in %s replicates, where precision is undefined and every other column counts 0.",
          max(S$n), what, paste(e$set, collapse = ", "), paste(e$n, collapse = ", "), paste(max(S$n) - e$n, collapse = ", "))
}
src_line <- "Cells: median (minimum–maximum) over the replicates. Values: `%s.tsv`; source: `Data_S2.zip` (`detection_per_seed.tsv`%s); built by `benchmarks/mvp_ms_table_s1_detection.R`."
writeLines(c("# t_detection", "",
             sprintf(src_line, "t_detection", ""), "",
             md_table_a(sum_a, "Candidate set"), "", prec_note(sum_a, "")),
           file.path(out_dir, "t_detection.md"))
arch_md <- unlist(lapply(arch_levels, function(a) {
  S <- sum_arch[arch == a]
  c(sprintf("## %s (%d replicates)", paste0(toupper(substring(a, 1, 1)), substring(a, 2)), max(S$n)), "",
    md_table_a(S, "Candidate set"), "", prec_note(S, " of this architecture"), "")
}))
writeLines(c("# t_detection_arch", "",
             sprintf(src_line, "t_detection_arch", ""), "", arch_md),
           file.path(out_dir, "t_detection_arch.md"))
writeLines(c("# t_replicates", "",
             sprintf(src_line, "t_replicates", ", `rank_metrics.tsv`"), "",
             md_row(c("Quantity", "Replicates", "Median (minimum–maximum)")), md_row(rep("---", 3)),
             row_b("Tested SNPs", sum_b[metric == "n_snps"], "int"),
             row_b("Causal loci, all replicates", sum_b[metric == "causal_total" & stratum == "all"], "int"),
             vapply(arch_levels, function(a) row_b(paste0("Causal loci, ", a), sum_b[metric == "causal_total" & stratum == a], "int"),
                    character(1)),
             row_b("Causal loci, % of tested SNPs", sum_b[metric == "causal_pct"], "pct2")),
           file.path(out_dir, "t_replicates.md"))
message("wrote t_detection, t_detection_arch, t_replicates {tsv,md} to ", out_dir)
