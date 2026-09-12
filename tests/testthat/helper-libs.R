# Sources the pipeline libraries under test.
#
# None of scripts/R/lib/*.R or scripts/R/utils/*.R calls library() itself — every
# one of them assumes the sourcing wrapper has already attached its packages (see
# the header of scripts/find_sig_snps.R for the production pattern). Several use a
# BARE %>% or a bare qvalue()/covRob(), so those packages must be ATTACHED, not
# merely installed.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)      # bare %>% in combine_sigsnps.R, io_pvalues.R, manhattan_utils.R
    library(stringr)
    library(tidyr)
    library(qvalue)     # bare qvalue() in pval_threshold.R, rdadapt.R
    library(robust)     # bare covRob() in rdadapt.R
    library(parallel)
    library(ggplot2)    # theme_clinego.R / manhattan_utils.R return theme + layer objects
})

.clinego_root <- getOption(
    "clinego.repo_root",
    normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
)
.clinego_R <- file.path(.clinego_root, "scripts", "R")

# Order matters: no lib sources another, the wrapper always does the ordering.
#   pval_threshold.R  <- sig_snps.R, io_pvalues.R
#   gff_parsing.R     <- genes_in_regions.R
#   gff_parsing.R     <- enrichment.R   (extract_gene_id(), declared enrichment.R:4)
for (.f in c("utils/pval_threshold.R",
             "utils/io_pvalues.R",
             "utils/io_selected_snps.R",
             "utils/manhattan_utils.R",   # parse_assoc_files_str(), used by load_assoc_data()
             "lib/gff_parsing.R",
             "lib/regions.R",
             "lib/region_distance.R",
             "lib/sig_snps.R",
             "lib/combine_sigsnps.R",
             "lib/genes_in_regions.R",
             "lib/rdadapt.R",
             "lib/invariants.R",
             "utils/theme_clinego.R",
             "lib/enrichment.R",
             "lib/enrichment_plots.R")) {
    source(file.path(.clinego_R, .f))
}
rm(.f)

# Deliberately NOT sourced here — and the reason is NOT "unsourceable". Measured
# 2026-09-12: all of scripts/R/lib and scripts/R/utils sys.source() cleanly in the
# image with zero symbol collisions against the vector above. Two exceptions
# remain, for different reasons:
#
#   utils/emmax_core.R  — sourced LOCALLY inside tests/testthat/test-emmax_core.R
#                         instead. EMMAX_BIN/EMMAX_RUN (:5-6) are plain rebindable
#                         string globals, and that test rebinds EMMAX_BIN to
#                         /bin/false to exercise run_emmax()'s stop-on-nonzero
#                         contract. test_dir() shares ONE parent environment
#                         across every test file, so a rebinding here would leak
#                         into every alphabetically-later file.
#   utils/logging.R     — dead code. log_info/log_warn/log_error are one-line
#                         message() wrappers with zero callers repo-wide and
#                         nothing sources the file. Filed for deletion; not
#                         tested, because a test would be the only caller.

# Shared test helper. Every lib under test logs progress with message(), so almost
# every call needs wrapping. Defined ONCE here: test_dir() sources all test files
# into a shared parent environment, so a per-file copy in each of six files would
# mean the last one sourced silently wins for all of them.
quiet <- function(expr) suppressMessages(expr)

# Build a minimal enrichResult S4 object for the enrichment tests, WITHOUT running
# an enrichment and without GO.db. The class is defined in DOSE (installed), but
# methods::new() resolves S4 classes from the calling topenv, so the class must be
# fetched from DOSE's namespace explicitly — attaching clusterProfiler just for a
# class definition would also drag AnnotationDbi's select() onto the search path
# and mask dplyr::select() for every later file in the run.
new_enrich_result <- function(result, gene) {
    methods::new(methods::getClass("enrichResult", where = asNamespace("DOSE")),
                 result = result, gene = gene)
}
