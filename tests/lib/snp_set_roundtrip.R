# snp_set_roundtrip.R — the fixture of the Shiny <-> config SNP-set round trip.
#
# A SNP set saved from the GEA tab writes a record (snp_set.yaml) and shows a
# config block (GEA.snp_sets). Three things must agree, and each is pinned by a
# different suite against ONE static golden pair in
# tests/fixtures/snp_set_roundtrip/:
#
#   app save   -> record + snippet   == golden record.json + snippet.yaml
#                 (tests/testthat/test-equivalence-app-pipeline.R, Block 7b)
#   record     -> build_snp_set.R    == the app's own selected_snps.tsv, byte for byte
#                 (same block)
#   snippet    -> common.smk resolve == golden record.json's pairs / cells
#                 (tests/python/test_snp_set_declarations.py, TestRoundTrip)
#
# Python cannot run R and the R suite cannot load common.smk, so the goldens are
# the hand-off; neither test writes into the repository.
#
# The universe is the python harness's BASE_CONFIG (tests/python/_smk_harness.py):
# predictors bio_1..bio_3, methods EMMAX / LFMM / RDA, RDA's single column
# climate_multivariate, RDA pinned to its registry rule (bonf 0.01) whenever a
# threshold is declared. The selection is deliberately awkward:
#   - non-rectangular pairs (bio_3 and bio_2::LFMM left out), so the snippet
#     takes the exact `pairs:` form;
#   - bio_2::LFMM carries the smallest p in the fixture, so dropping the pair
#     filter changes the TSV;
#   - overrides: top 3 on bio_2::EMMAX, a 9-significant-digit custom cutoff on
#     bio_1::LFMM (a SNP sits exactly on it — inclusive <=), and one on the
#     UNselected bio_3::LFMM, which the snippet must leave out (common.smk
#     refuses an override key outside the selected pairs);
#   - Cross-method at 50 kb, so EMMAX's 1:80000 is dropped for having no partner
#     in another method within 50 kb.
#
# Regenerate the goldens (only after a deliberate change to the record or
# snippet format; review the diff):
#   docker run ... -v $PWD:/pipeline cline-go:latest \
#     bash /pipeline/tests/with_app_from_mount.sh \
#     Rscript /pipeline/tests/lib/snp_set_roundtrip.R --write-golden

RT_K          <- 3L
RT_NAME       <- "rt.set"
RT_MASTER     <- list(type = "custom", value = 1e-4)
RT_STRATEGY   <- "Cross-method"
RT_GAP        <- 50000L
RT_SELECTED   <- c("bio_1::EMMAX", "bio_2::EMMAX", "bio_1::LFMM", "climate_multivariate::RDA")
RT_OVERRIDES  <- list(
    "bio_2::EMMAX" = list(type = "top",    value = 3),
    "bio_1::LFMM"  = list(type = "custom", value = 1.23456789e-3),
    "bio_3::LFMM"  = list(type = "custom", value = 0.5)
)
RT_DECL_FIELDS <- c("source_module", "regime", "strategy", "clumping_distance",
                    "pairs", "cells")

rt_golden_dir <- function(root) file.path(root, "tests", "fixtures", "snp_set_roundtrip")

# Wide p-value tables, method -> data.table(SNPID, chr, pos, <trait columns>).
# Filler p-values are >= 0.01 and spread 100 kb apart, so only the planted SNPs
# can pass any rule here and no filler SNP is a Cross-method partner.
rt_pvalue_tables <- function() {
    planted <- data.table::data.table(
        chr = c("1", "1", "1", "1", "1", "2", "2", "2", "2", "3", "3", "3"),
        pos = c(1000L, 1500L, 5000L, 5100L, 80000L, 3000L, 20000L, 40000L, 70000L,
                9000L, 9100L, 50000L))
    filler <- data.table::data.table(
        chr = rep(c("1", "2", "3"), each = 10L),
        pos = rep(200000L + 100000L * (0:9), 3L))
    snps <- rbind(planted, filler)
    snps[, SNPID := paste0(chr, ":", pos)]
    data.table::setcolorder(snps, c("SNPID", "chr", "pos"))
    fill <- function(j) ((seq_len(nrow(snps)) * 37L + j * 11L) %% 97L + 1L) / 100
    plant <- function(p, hits) {
        p[match(names(hits), snps$SNPID)] <- unname(hits)
        p
    }
    emmax <- data.table::copy(snps)
    emmax[, bio_1 := plant(fill(1L), c("1:1000" = 1e-6, "1:5000" = 2e-5, "2:3000" = 5e-5))]
    emmax[, bio_2 := plant(fill(2L), c("3:9100" = 1e-7, "2:40000" = 2e-7, "1:80000" = 3e-7))]
    emmax[, bio_3 := plant(fill(3L), c("3:50000" = 1e-8))]
    lfmm <- data.table::copy(snps)
    lfmm[, bio_1 := plant(fill(4L), c("1:1500" = 1e-3, "3:9000" = 1.2e-3,
                                      "2:20000" = 1.23456789e-3, "2:70000" = 1.1e-3))]
    lfmm[, bio_2 := plant(fill(5L), c("1:1000" = 1e-9))]
    lfmm[, bio_3 := fill(6L)]
    rda <- data.table::copy(snps)
    rda[, climate_multivariate := plant(fill(7L), c("1:5100" = 1e-5, "2:70000" = 1e-4))]
    list(EMMAX = emmax, LFMM = lfmm, RDA = rda)
}

# Write the tables where the app's loader finds them; returns the file paths.
rt_write_tables <- function(project) {
    tabs <- rt_pvalue_tables()
    vapply(names(tabs), function(m) {
        f <- clinego.app:::mod_path(project, "GEA", "tables", "methods", m,
                                    paste0(m, "_pvalues_K", RT_K, ".tsv"))
        dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
        data.table::fwrite(tabs[[m]], f, sep = "\t")
        f
    }, character(1))
}

# Run the GEA tab's save path end to end on the fixture, through the app's own
# functions in the order mod_gea.R calls them. The caller sets
# options(clinego.pipeline_path) to a scratch directory first.
rt_app_save <- function(project, registry_root, selection = RT_SELECTED,
                        name = RT_NAME, overrides = RT_OVERRIDES) {
    files <- rt_write_tables(project)
    reg <- clinego.app:::gea_method_significance_defaults(
        clinego.app:::gea_method_registry(pipeline_path = registry_root))
    pv <- clinego.app:::load_all_method_pvalues(project, "GEA", RT_K)
    sig <- clinego.app:::compute_method_sigsnps_cached(
        pv, RT_MASTER$type, RT_MASTER$value, k = RT_K, regime = "snp",
        project = project, module = "GEA",
        overrides = overrides, registry_defaults = reg)
    # combo_counts as mod_gea.R builds it: n > 0 cells only
    counts <- list()
    for (m in names(sig)) {
        tc <- sig[[m]][, .(n = .N), by = "trait"]
        for (i in seq_len(nrow(tc))) counts[[paste0(tc$trait[i], "::", m)]] <- tc$n[i]
    }
    traits <- sort(unique(unlist(lapply(pv, function(dt) setdiff(names(dt), c("SNPID", "chr", "pos"))))))
    json <- if (is.null(selection)) NULL else as.character(jsonlite::toJSON(selection))
    dt <- clinego.app:::compute_interactive_sigsnps(
        sig, json, counts, traits, RT_STRATEGY, RT_GAP, project, "GEA")
    pairs <- clinego.app:::snp_set_pairs_for_save(json, counts, traits, pv)
    record <- clinego.app:::snp_set_record_for_save(
        name = name, source_module = "GEA", pairs = pairs, overrides = overrides,
        master_type = RT_MASTER$type, master_value = RT_MASTER$value,
        registry_defaults = reg, strategy = RT_STRATEGY, clumping_distance = RT_GAP,
        regime = "snp", k_best = RT_K)
    n <- clinego.app:::save_snp_set(project, name, dt, list(source_module = "GEA"),
                                    record = record)
    list(files = files, registry = reg, pvalues = pv, counts = counts, record = record,
         n = n, snippet = clinego.app:::snp_set_config_snippet(record, reg))
}

# The declaration fields of a record, as written to record.json.
rt_decl_json <- function(record) {
    out <- record[RT_DECL_FIELDS]
    out$pairs <- I(as.character(unlist(out$pairs)))
    jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = NA)
}

if (sys.nframe() == 0L && "--write-golden" %in% commandArgs(trailingOnly = TRUE)) {
    suppressPackageStartupMessages(library(data.table))
    root <- "/pipeline"
    options(clinego.pipeline_path = tempfile("rt_golden_"))
    dir.create(getOption("clinego.pipeline_path"))
    res <- rt_app_save("RT_GOLDEN", registry_root = root)
    out <- rt_golden_dir(root)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines(res$snippet$yaml, file.path(out, "snippet.yaml"), sep = "")
    writeLines(rt_decl_json(res$record), file.path(out, "record.json"))
    message("Wrote ", out, " (", res$n, " SNPs)")
}
