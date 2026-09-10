#!/usr/bin/env Rscript
# mvp_gate_compare.R -- Phase B validation gate.
#
# Compares genetic_offset_site.tsv produced by the SCENARIO-aware pipeline against
# the harvested reference produced by the old per-garden driver, for one seed.
#
# Two verdicts are reported SEPARATELY and they are not the same thing:
#   values_equal : every site's genetic_offset agrees within --tol. This is the
#                  hard gate -- a difference here means the refactor changed the
#                  math and the run must stop.
#   bytes_equal  : the files are byte-identical. A value-equal/byte-different
#                  result is a PASS WITH A NOTE: the scenario loop can legitimately
#                  reorder tied rows or shift fwrite's digit rendering. Do not
#                  chase a byte diff that carries no value diff.
#
# Rows are matched on `sample`, never on position: the site tables are written
# sorted by descending offset, so ties can permute rows without any value changing.
#
# Reference layout (old driver):  {ref}/{garden}/{method}__{panel}.tsv
# New layout (scenario pipeline): {results}/Maladaptation/tables/{method}/{panel}_nospatial/{garden}/genetic_offset_site.tsv
#
# The reference calls the Condition()-corrected RDA `rda_corrected`; the pipeline
# now calls it `rda_offset_corrected` (its own registry entry). Mapped below.
#
# Usage:
#   Rscript benchmarks/mvp_gate_compare.R --results=MVP1232218p1_results \
#       --ref=benchmarks/mvp_eval/offset09/_gate_reference_1232218 \
#       [--out=benchmarks/mvp_eval/gate/gate_compare.tsv] [--tol=1e-9]

suppressPackageStartupMessages({
    library(data.table)
    library(tools)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/mnt/data/eugene/CLINE-GO")
args <- commandArgs(trailingOnly = TRUE)
kv <- function(k, d = NULL) {
    hit <- grep(paste0("^--", k, "="), args, value = TRUE)
    if (!length(hit)) return(d)
    sub(paste0("^--", k, "="), "", hit[1])
}

RESULTS <- kv("results")
REF     <- kv("ref")
OUT     <- kv("out", file.path(ROOT, "benchmarks/mvp_eval/gate/gate_compare.tsv"))
TOL     <- as.numeric(kv("tol", "1e-9"))
if (is.null(RESULTS) || is.null(REF)) stop("--results= and --ref= are required")
if (!grepl("^/", RESULTS)) RESULTS <- file.path(ROOT, RESULTS)
if (!grepl("^/", REF))     REF     <- file.path(ROOT, REF)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# reference method token -> pipeline method directory
METHOD_MAP <- c(gradient_forest   = "gradient_forest",
                geometric_offset  = "geometric_offset",
                rda_offset        = "rda_offset",
                rda_corrected     = "rda_offset_corrected")

read_offset <- function(p) {
    fread(p, colClasses = c("site" = "character", "sample" = "character"))
}

gardens <- sort(list.dirs(REF, full.names = FALSE, recursive = FALSE))
gardens <- gardens[nzchar(gardens)]
message(sprintf("== gate compare: %d reference gardens ==", length(gardens)))

rows <- list()
for (g in gardens) {
    ref_files <- list.files(file.path(REF, g), pattern = "\\.tsv$", full.names = FALSE)
    for (rf in ref_files) {
        token <- sub("__.*$", "", rf)
        panel <- sub("\\.tsv$", "", sub("^.*__", "", rf))
        method <- METHOD_MAP[[token]]
        if (is.null(method) || is.na(method)) {
            rows[[length(rows) + 1]] <- data.table(
                garden = g, method = token, panel = panel,
                status = "unmapped_method", n = NA_integer_,
                max_abs_diff = NA_real_, values_equal = NA, bytes_equal = NA)
            next
        }
        ref_p <- file.path(REF, g, rf)
        new_p <- file.path(RESULTS, "Maladaptation", "tables", method,
                           paste0(panel, "_nospatial"), g, "genetic_offset_site.tsv")
        if (!file.exists(new_p)) {
            rows[[length(rows) + 1]] <- data.table(
                garden = g, method = method, panel = panel,
                status = "missing_new", n = NA_integer_,
                max_abs_diff = NA_real_, values_equal = NA, bytes_equal = NA)
            next
        }
        a <- read_offset(ref_p); b <- read_offset(new_p)
        m <- merge(a[, .(sample, ref = genetic_offset)],
                   b[, .(sample, new = genetic_offset)],
                   by = "sample", all = TRUE)
        if (nrow(m) != nrow(a) || nrow(m) != nrow(b) || anyNA(m$ref) || anyNA(m$new)) {
            rows[[length(rows) + 1]] <- data.table(
                garden = g, method = method, panel = panel,
                status = "sample_set_mismatch", n = nrow(m),
                max_abs_diff = NA_real_, values_equal = FALSE, bytes_equal = FALSE)
            next
        }
        d <- max(abs(m$ref - m$new))
        rows[[length(rows) + 1]] <- data.table(
            garden = g, method = method, panel = panel,
            status = "compared", n = nrow(m),
            max_abs_diff = d,
            values_equal = d <= TOL,
            bytes_equal  = identical(unname(md5sum(ref_p)), unname(md5sum(new_p))))
    }
}

res <- rbindlist(rows)
fwrite(res, OUT, sep = "\t")

cmp <- res[status == "compared"]
message("\n=== gate result ===")
message(sprintf("compared            : %d", nrow(cmp)))
message(sprintf("missing in new tree : %d", sum(res$status == "missing_new")))
message(sprintf("sample-set mismatch : %d", sum(res$status == "sample_set_mismatch")))
message(sprintf("unmapped method     : %d", sum(res$status == "unmapped_method")))
if (nrow(cmp)) {
    message(sprintf("values equal (tol %g): %d / %d", TOL, sum(cmp$values_equal), nrow(cmp)))
    message(sprintf("bytes  equal         : %d / %d", sum(cmp$bytes_equal), nrow(cmp)))
    message(sprintf("max abs diff overall : %.3e", max(cmp$max_abs_diff)))
    message("\nper method:")
    print(cmp[, .(n = .N,
                  values_equal = sum(values_equal),
                  bytes_equal  = sum(bytes_equal),
                  max_abs_diff = max(max_abs_diff)), by = method])
}
message(sprintf("\nwrote %s", OUT))

failed <- res[status == "sample_set_mismatch" | (status == "compared" & !values_equal)]
if (nrow(failed)) {
    message("\n!! VALUE DIFFERENCES -- GATE FAILS. First 10:")
    print(head(failed, 10))
    quit(status = 1)
}
if (any(res$status == "missing_new")) {
    message("\n!! Some reference files have no counterpart in the new tree (see status=missing_new).")
    message("   Expected for panels the sweep no longer runs; investigate before treating as a pass.")
}
message("\nGate: no value differences.")
