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
})

.clinego_root <- getOption(
    "clinego.repo_root",
    normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
)
.clinego_R <- file.path(.clinego_root, "scripts", "R")

# Order matters: no lib sources another, the wrapper always does the ordering.
#   pval_threshold.R  <- sig_snps.R, io_pvalues.R
#   gff_parsing.R     <- genes_in_regions.R
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
             "lib/rdadapt.R")) {
    source(file.path(.clinego_R, .f))
}
rm(.f)

# Deliberately NOT sourced:
#   utils/emmax_core.R    — sets EMMAX_BIN/EMMAX_RUN to /pipeline paths at source
#                           time and shells out to a binary that refuses to run on
#                           an emulated host
#   utils/theme_clinego.R — seven global constants at source time; plotting only
#   lib/enrichment*.R     — need GO.db / clusterProfiler / enrichplot; annotation
#                           and plotting, not the region/SNP science Tier 1 covers

# Shared test helper. Every lib under test logs progress with message(), so almost
# every call needs wrapping. Defined ONCE here: test_dir() sources all test files
# into a shared parent environment, so a per-file copy in each of six files would
# mean the last one sourced silently wins for all of them.
quiet <- function(expr) suppressMessages(expr)
