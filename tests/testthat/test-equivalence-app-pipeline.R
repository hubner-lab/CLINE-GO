# Shiny app <-> Snakemake pipeline equivalence (test plan Tier 5).
#
# ── What this file is for ─────────────────────────────────────────────────────
# CLINE-GO computes the same science twice: the pipeline (scripts/*.R wrappers
# over scripts/R/lib/) and the Shiny app (scripts/clinego.app/R/). zzz.R:39-42
# shares only TWO libs — regions.R and pval_threshold.R. Gene finding and
# combining are REIMPLEMENTED in the app, and region_distance.R is never loaded
# by it. Tiers 1 and 2 pinned each side's own behaviour; nothing until now
# asserted the two AGREE.
#
# There are FIVE divergence surfaces (the dossier's count of four is stale —
# region_distance-by-omission is the fifth, see the last block).
#
# ── This file MEASURES divergence; it does not fix it ─────────────────────────
# Where the two sides genuinely differ, the divergence is filed in
# docs/pipeline_improvement_requests.md and the correct-behaviour assertion is
# skip()ped naming the filing — the convention tests/testthat/test-known-bugs.R
# established. Never weaken an assertion here to make the output green.
#
# ── The two environments, and why they do not collide ────────────────────────
# zzz.R:57-60 sources CLINEGO_SHARED_LIBS into asNamespace("clinego.app");
# helper-libs.R:28-41 sources the pipeline libs into the test environment. They
# are distinct environments, so `library(clinego.app)` cannot clobber the
# functions the helper loaded (verified: cluster_snps_to_regions is identical()
# across the two, resolve_clumping_distance exists in only one).
#
# NAME COLLISION: combine_sigsnps exists on BOTH sides with DIFFERENT arity —
# pipeline (sig_snps_list, strategy, clumping_distance, predictors), app
# (sigsnps_list, strategy, gap). ALWAYS write clinego.app:::combine_sigsnps for
# the app's. cluster_snps_to_regions / compute_pval_threshold need no
# qualification: they are the same object, which is the whole point.
#
# helper-libs.R must NOT gain library(clinego.app) — test_dir() sources every
# helper for the directory, so the entire pipeline suite would fail whenever the
# app is not installed. The skip lives here instead.

app_available <- function() {
    requireNamespace("clinego.app", quietly = TRUE)
}

skip_without_app <- function() {
    skip_if_not(app_available(),
                "clinego.app not installed (needs the image's own package + -v $PWD:/pipeline)")
}

# --- fixtures ---------------------------------------------------------------

# BUILDERS, not shared objects. The pipeline's combine_sigsnps() does
# `dt[, chr := as.character(chr)]` (combine_sigsnps.R:50) — data.table `:=`
# mutates by reference, so a shared fixture would be modified out from under a
# later test. Same reason test-regions.R:10-18 uses a builder.

# The LONG per-method significant-SNP table. This is the single source of truth
# for every block below: it is the shape BOTH sides accept
# (pipeline combine_sigsnps.R:51, app fct_combine.R:63 / fct_regions.R:23).
sig_long <- function(SNPID, chr, pos, trait, method, pvalue) {
    data.table::data.table(SNPID = SNPID, chr = as.character(chr),
                           pos = as.integer(pos), trait = trait,
                           method = method, pvalue = as.numeric(pvalue))
}

# Two methods, two traits, two chromosomes, with a deliberate cross-method
# overlap (s2 is significant in both) so Union and Cross-method differ.
two_method_fixture <- function() {
    list(
        EMMAX = sig_long(c("1:1000", "1:2000", "2:5000"),
                         c("1", "1", "2"),
                         c(1000L, 2000L, 5000L),
                         c("bio_1", "bio_1", "bio_2"),
                         "EMMAX",
                         c(1e-8, 1e-7, 1e-6)),
        LFMM  = sig_long(c("1:2000", "2:9000"),
                         c("1", "2"),
                         c(2000L, 9000L),
                         c("bio_2", "bio_1"),
                         "LFMM",
                         c(1e-9, 1e-5))
    )
}

PREDICTORS <- c("bio_1", "bio_2")

# The WIDE table the pipeline's region builder consumes is NOT hand-projected —
# it is derived with the pipeline's own combine_sigsnps(), which is exactly what
# the production DAG does (find_sig_snps -> combine_selected_snps ->
# create_regions). Hand-shaping a second table would make a projection bug read
# as a scientific divergence.
wide_from_long <- function(sig_list, distance, strategy = "Union") {
    quiet(combine_sigsnps(sig_list, strategy, distance, PREDICTORS))
}

# GFF fixture written to disk so read_gff() (the PIPELINE's own loader) produces
# the table handed to BOTH sides. gene_id derivation differs between the sides
# (the app additionally strips isoform suffixes, fct_data_loading.R:543), so
# gene-set comparisons key on (chr, start, end) and gene_id is asserted apart.
write_gff <- function(lines, envir = parent.frame()) {
    path <- withr::local_tempfile(fileext = ".gff3", .local_envir = envir)
    writeLines(c("##gff-version 3", lines), path)
    path
}

gff_line <- function(chr, start, end, attrs, feature = "gene") {
    paste(chr, "test", feature, start, end, ".", "+", ".", attrs, sep = "\t")
}

# --- 0. the gate: are we really holding TWO implementations? ----------------
#
# A suite that has not proven it holds two distinct functions is not evidence of
# agreement — it can compare an implementation against itself and go green. Each
# assertion below records a MEASURED fact about the sharing boundary, so a future
# change to zzz.R's CLINEGO_SHARED_LIBS turns one of them red.

test_that("region clustering is genuinely ONE shared implementation", {
    skip_without_app()
    app_fn <- get("cluster_snps_to_regions", envir = asNamespace("clinego.app"))
    # Same source file (zzz.R:40 -> lib/regions.R), so the bodies must match.
    # If this goes red, the app stopped sharing regions.R and every "strict
    # equality" assertion in block 1 has to be re-justified as cross-impl.
    expect_identical(body(app_fn), body(cluster_snps_to_regions))
})

test_that("the p-value threshold core is genuinely ONE shared implementation", {
    skip_without_app()
    app_fn <- get("compute_pval_threshold", envir = asNamespace("clinego.app"))
    expect_identical(body(app_fn), body(compute_pval_threshold))
})

test_that("gene finding and combining are TWO implementations, not one", {
    skip_without_app()
    ns <- asNamespace("clinego.app")

    # Different names AND different shapes: many-regions vs one-region.
    expect_true(exists("find_genes_in_region", envir = ns, inherits = FALSE))
    expect_named(formals(get("find_genes_in_region", envir = ns)),
                 c("gff_genes", "region_row", "promoter_length"))
    expect_named(formals(find_genes_for_regions),
                 c("regions_dt", "gff_dt", "gff_path", "promoter_length", "allsnps_dt"))

    # Same NAME, different arity — this is the collision that makes ::: mandatory.
    expect_named(formals(clinego.app:::combine_sigsnps),
                 c("sigsnps_list", "strategy", "gap"))
    expect_named(formals(combine_sigsnps),
                 c("sig_snps_list", "strategy", "clumping_distance", "predictors"))
    expect_false(identical(body(clinego.app:::combine_sigsnps), body(combine_sigsnps)))
})

# --- 1. region clustering ---------------------------------------------------
#
# Shared core (cluster_snps_to_regions), divergent wrappers. Geometry goes
# through one code path, so it must match EXACTLY. traits/methods are two
# independent mechanisms — the pipeline parses the wide method columns
# (regions.R:253-261), the app DELETES the shared output (fct_regions.R:32) and
# re-derives it from the long rows via foverlaps (:35-51) — so those are set
# comparisons.

DIST <- 500L

regions_both_sides <- function(distance = DIST) {
    sig_list <- two_method_fixture()
    wide <- wide_from_long(sig_list, distance)
    pipe <- quiet(build_combined_regions(wide, distance))

    long <- data.table::rbindlist(two_method_fixture())
    app  <- quiet(clinego.app:::compute_all_regions(long, distance))

    list(pipe = pipe, app = app)
}

test_that("region GEOMETRY is identical (the shared code path)", {
    skip_without_app()
    r <- regions_both_sides()

    geom <- c("chr", "start", "end", "length", "snp_count", "snp_ids")
    # Sorted copies: the app's join does not preserve row order (see the next
    # test), so a positional comparison would fail for the wrong reason.
    p <- data.table::setorder(r$pipe[, ..geom], chr, start)
    a <- data.table::setorder(r$app[, ..geom], chr, start)

    expect_identical(nrow(p), nrow(a))
    expect_equal(as.data.frame(p), as.data.frame(a))
})

test_that("region_id is identical (both use the combined {chr}_{start}-{end} form)", {
    skip_without_app()
    r <- regions_both_sides()
    expect_setequal(r$pipe$region_id, r$app$region_id)
})

test_that("traits and methods agree as SETS per region (two independent mechanisms)", {
    skip_without_app()
    r <- regions_both_sides()

    split_sorted <- function(x) lapply(strsplit(x, ","), function(v) sort(v[nzchar(v)]))
    key <- function(dt, col) setNames(split_sorted(dt[[col]]), dt$region_id)

    pt <- key(r$pipe, "traits");  at <- key(r$app, "traits")
    pm <- key(r$pipe, "methods"); am <- key(r$app, "methods")

    expect_setequal(names(pt), names(at))
    for (rid in names(pt)) {
        expect_identical(at[[rid]], pt[[rid]],
                         info = paste("traits disagree for region", rid))
        expect_identical(am[[rid]], pm[[rid]],
                         info = paste("methods disagree for region", rid))
    }
})

test_that("min_pvalue per region agrees", {
    skip_without_app()
    r <- regions_both_sides()
    p <- r$pipe[, .(region_id, min_pvalue)]
    a <- r$app[,  .(region_id, min_pvalue)]
    m <- merge(p, a, by = "region_id", suffixes = c(".pipe", ".app"))
    expect_identical(nrow(m), nrow(p))
    expect_equal(m$min_pvalue.app, m$min_pvalue.pipe)
})

test_that("the pipeline adds per-trait/per-method SNP-count columns the app has not", {
    skip_without_app()
    # Not a defect — a deliberate shape difference (regions.R:266-275). Pinned so
    # that comparing on the intersection of column names stays justified.
    r <- regions_both_sides()
    expect_true(any(grepl("_snps$", names(r$pipe))))
    expect_false(any(grepl("_snps$", names(r$app))))
})

test_that("KNOWN DIVERGENCE: the app's join breaks the shared chr/start ordering", {
    # CORRECT behaviour, deliberately not made to pass.
    # cluster_snps_to_regions() ends with setorder(regions, chr, start)
    # (regions.R:102). fct_regions.R:51 then does
    # `regions <- trait_methods[regions, on = "region_id"]`, a data.table right
    # join keyed by region_id, which does not preserve that order. The pipeline
    # path keeps it. Anything reading the app's table positionally (a "first
    # region" pick, a head()) gets a different region than the pipeline would.
    # Filed 2026-09-10 in docs/pipeline_improvement_requests.md.
    skip("known divergence: app compute_all_regions loses chr/start order — filed 2026-09-10")

    skip_without_app()
    r <- regions_both_sides()
    expect_identical(r$app$chr,   r$pipe$chr)
    expect_identical(r$app$start, r$pipe$start)
})

test_that("both sides return their empty skeleton for zero input", {
    skip_without_app()
    empty <- sig_long(character(), character(), integer(),
                      character(), character(), numeric())
    a <- quiet(clinego.app:::compute_all_regions(empty, DIST))
    p <- quiet(build_combined_regions(NULL, DIST))
    expect_identical(nrow(a), 0L)
    expect_identical(nrow(p), 0L)
    # Same ten column names, same order — verified, and worth pinning because
    # downstream Shiny code rbinds onto this skeleton.
    common <- c("region_id", "chr", "start", "end", "length",
                "snp_count", "snp_ids", "traits", "methods", "min_pvalue")
    expect_identical(names(a), common)
    expect_identical(intersect(names(p), common), common)
})

# --- 2. gene finding -------------------------------------------------------
#
# TWO implementations. The pipeline uses region bounds AS-IS and spends
# promoter_length only on the promoter-SNP counting window
# (genes_in_regions.R:95); the app spends it on the QUERY WINDOW
# (fct_regions.R:73), changing which genes come back. Feeding both the same
# non-zero promoter_length is therefore guaranteed to differ for a reason that
# has nothing to do with overlap logic.

GFF_LINES <- c(
    gff_line("1",  1000,  2000, "ID=g1;Name=alpha;biotype=protein_coding"),
    gff_line("1",  4000,  5000, "ID=g2;Name=beta;biotype=protein_coding"),
    gff_line("1",  9000, 10000, "ID=g3;Name=gamma;biotype=protein_coding"),
    gff_line("2",  1000,  2000, "ID=g4;Name=delta;biotype=protein_coding")
)

# App loops one region at a time; union the per-region results to compare
# against the pipeline's all-regions call.
app_genes <- function(gff_dt, regions_dt, promoter_length) {
    out <- lapply(seq_len(nrow(regions_dt)), function(i) {
        g <- clinego.app:::find_genes_in_region(gff_dt, regions_dt[i], promoter_length)
        if (is.null(g) || nrow(g) == 0) return(NULL)
        g[, .(chr = as.character(chr), start, end)]
    })
    unique(data.table::rbindlist(out))
}

pipe_genes <- function(gff_dt, regions_dt, gff_path, promoter_length) {
    r <- quiet(find_genes_for_regions(regions_dt, gff_dt, gff_path,
                                      promoter_length = promoter_length))
    unique(r$genes_per_region[, .(chr = as.character(chr),
                                  start = gene_start, end = gene_end)])
}

gene_regions <- function() {
    data.table::data.table(region_id = c("R1", "R2"),
                           chr = c("1", "2"),
                           start = c(4500L, 500L),
                           end   = c(9500L, 1500L))
}

test_that("gene SETS agree when promoter_length is normalised to 0", {
    skip_without_app()
    path <- write_gff(GFF_LINES)
    gff  <- quiet(read_gff(path, "gene"))
    regs <- gene_regions()

    a <- app_genes(gff, regs, 0L)
    p <- pipe_genes(gff, regs, path, 0L)

    data.table::setorder(a, chr, start); data.table::setorder(p, chr, start)
    expect_equal(as.data.frame(a), as.data.frame(p))
    # Guard against a vacuous pass: the fixture must actually find genes.
    expect_gt(nrow(a), 0L)
})

test_that("the promoter_length divergence is EXACTLY a query-window shift", {
    skip_without_app()
    # The paired half of the normalisation above. Without this, the suite would
    # sidestep the divergence and never exercise the app's production behaviour
    # (a non-zero promoter_length). Widening the PIPELINE's region start by PL
    # by hand must reproduce the app's own extension at that same PL.
    # If this fails, the divergence is bigger than a window shift and the filing
    # understates it.
    PL   <- 3000L
    path <- write_gff(GFF_LINES)
    gff  <- quiet(read_gff(path, "gene"))
    regs <- gene_regions()

    widened <- data.table::copy(regs)[, start := pmax(1L, start - PL)]

    a <- app_genes(gff, regs, PL)            # app extends internally
    p <- pipe_genes(gff, widened, path, PL)  # pipeline needs it pre-extended

    data.table::setorder(a, chr, start); data.table::setorder(p, chr, start)
    expect_equal(as.data.frame(a), as.data.frame(p))
    # The shift must actually change the answer, or this proves nothing.
    expect_gt(nrow(a), nrow(app_genes(gff, regs, 0L)))
})

test_that("the app extends the gene window upstream only, like the pipeline's promoter window", {
    skip_without_app()
    path <- write_gff(GFF_LINES)
    gff  <- quiet(read_gff(path, "gene"))
    # A region ending just before g3 (9000-10000). Extending DOWNSTREAM would
    # pull g3 in; neither side may do that.
    reg <- data.table::data.table(region_id = "R1", chr = "1",
                                  start = 4500L, end = 8000L)
    a <- app_genes(gff, reg, 3000L)
    expect_false(9000L %in% a$start)
})

test_that("KNOWN DIVERGENCE: the app's gene_id fallback is dead code", {
    # CORRECT behaviour, deliberately not made to pass.
    # fct_data_loading.R:533-544 builds gene_id with
    # regmatches(regexpr("(?<=Parent=)...")), which DROPS non-matching elements
    # instead of yielding "". So `na_idx <- !nzchar(id)` (:536) is never TRUE and
    # the ID= fallback at :538-540 can never execute. On a GFF whose selected
    # feature rows mix Parent= and ID=-only, the length mismatch throws inside
    # the tryCatch at :552 and load_gff_genes() returns an EMPTY table — every
    # region silently finds zero genes. The pipeline's extract_gene_id
    # (gff_parsing.R:8-13) uses length-preserving str_extract + ifelse and is
    # correct on all three cases (all-Parent, mixed, no-Parent).
    # Filed 2026-09-10 in docs/pipeline_improvement_requests.md.
    skip("known divergence: app load_gff_genes gene_id fallback unreachable — filed 2026-09-10")

    skip_without_app()
    mixed <- c(
        gff_line("1", 1000, 2000, "ID=g1;Name=alpha"),              # ID= only
        gff_line("1", 4000, 5000, "Parent=t2;ID=g2;Name=beta")      # has Parent=
    )
    path <- write_gff(mixed)
    genes <- quiet(clinego.app:::load_gff_genes(path, "gene"))
    expect_identical(nrow(genes), 2L)
    expect_setequal(genes$gene_id, c("g1", "g2"))
})

test_that("the pipeline's gene_id extraction handles all three attribute shapes", {
    # The correct-behaviour reference the app is measured against above. Not a
    # skip: this must hold today.
    path <- write_gff(c(
        gff_line("1", 1000, 2000, "ID=g1;Name=alpha"),
        gff_line("1", 4000, 5000, "Parent=t2;ID=g2;Name=beta")
    ))
    gff <- quiet(read_gff(path, "gene"))
    expect_identical(nrow(gff), 2L)
    expect_false(any(is.na(gff$gene_id)))
})

# --- 3. combining ----------------------------------------------------------
#
# TWO implementations sharing the name combine_sigsnps. Output shapes differ by
# design: pipeline WIDE one row per SNPID (combine_sigsnps.R:116-117), app LONG
# one row per (SNPID, method, trait) (fct_combine.R:131-135). So membership is
# compared first and min_pvalue separately.

app_combine <- function(sig_list, strategy, gap) {
    quiet(clinego.app:::combine_sigsnps(sig_list, strategy, gap))
}

test_that("the strategy alias tables are identical on both sides", {
    skip_without_app()
    # Asserted rather than assumed: a silent drift here would send the two sides
    # down different branches for the same user-visible label.
    for (alias in c("All", "Sum", "Overlap", "MethodOverlap", "PairOverlap",
                    "Union", "Cross-method", "Cross-method per-trait")) {
        expect_identical(clinego.app:::.normalize_strategy(alias),
                         .normalise_strategy(alias),
                         info = paste("alias", alias))
    }
})

test_that("both sides reject an unknown strategy rather than defaulting", {
    skip_without_app()
    sig <- two_method_fixture()
    expect_error(app_combine(sig, "NoSuchStrategy", 200L))
    expect_error(quiet(combine_sigsnps(two_method_fixture(), "NoSuchStrategy",
                                       200L, PREDICTORS)))
})

test_that("Union selects the same SNP set on both sides", {
    skip_without_app()
    gap <- 200L
    a <- app_combine(two_method_fixture(), "All", gap)
    p <- wide_from_long(two_method_fixture(), gap, "Union")
    expect_setequal(unique(a$SNPID), unique(p$SNPID))
})

test_that("Union agrees on min_pvalue per SNP", {
    skip_without_app()
    gap <- 200L
    a <- app_combine(two_method_fixture(), "All", gap)
    p <- wide_from_long(two_method_fixture(), gap, "Union")
    am <- unique(a[, .(SNPID, min_pvalue)])
    m  <- merge(am, p[, .(SNPID, min_pvalue)], by = "SNPID",
                suffixes = c(".app", ".pipe"))
    expect_identical(nrow(m), nrow(p))
    expect_equal(m$min_pvalue.app, m$min_pvalue.pipe)
})

test_that("single-method passthrough selects the same SNP set on both sides", {
    skip_without_app()
    a <- app_combine(two_method_fixture(), "EMMAX", 200L)
    p <- wide_from_long(two_method_fixture(), 200L, "EMMAX")
    expect_setequal(unique(a$SNPID), unique(p$SNPID))
})

test_that("Cross-method selects the same SNP set on both sides", {
    skip_without_app()
    gap <- 200L
    a <- app_combine(two_method_fixture(), "Overlap", gap)
    p <- wide_from_long(two_method_fixture(), gap, "Cross-method")
    expect_setequal(unique(a$SNPID), unique(p$SNPID))
    # Guard against a vacuous pass: the consensus must actually select something,
    # and strictly less than the union.
    expect_gt(length(unique(a$SNPID)), 0L)
    expect_lt(length(unique(a$SNPID)),
              length(unique(app_combine(two_method_fixture(), "All", gap)$SNPID)))
})

test_that("KNOWN DIVERGENCE: app Cross-method keeps rows that never matched", {
    # CORRECT behaviour, deliberately not made to pass.
    # After selecting consensus SNP ids, fct_combine.R:96-97 returns
    # all_snps[SNPID %in% selected_ids] — EVERY row for a selected SNP,
    # including (method, trait) rows that were not part of the match. The
    # pipeline returns only the matched rows (combine_sigsnps.R:193-198). So the
    # app's min_pvalue for a consensus SNP can be smaller than the pipeline's
    # whenever an unmatched row for that SNP carries a lower p-value.
    # Filed 2026-09-10 in docs/pipeline_improvement_requests.md.
    skip("known divergence: app Cross-method row set is a superset — filed 2026-09-10")

    skip_without_app()
    gap <- 200L
    a <- app_combine(two_method_fixture(), "Overlap", gap)
    p <- wide_from_long(two_method_fixture(), gap, "Cross-method")
    am <- unique(a[, .(min_pvalue = min(pvalue)), by = "SNPID"])
    m  <- merge(am, p[, .(SNPID, min_pvalue)], by = "SNPID",
                suffixes = c(".app", ".pipe"))
    expect_equal(m$min_pvalue.app, m$min_pvalue.pipe)
})

test_that("both sides return an empty result for no significant SNPs", {
    skip_without_app()
    empty <- list(EMMAX = sig_long(character(), character(), integer(),
                                   character(), character(), numeric()))
    a <- app_combine(empty, "All", 200L)
    # suppressWarnings: on an empty table combine_sigsnps.R:110 assigns a method
    # column with `result[[m]] <- NA_character_`, which data.table flags as a
    # shallow copy. That is the pipeline's own code path, not something this test
    # controls — kept quiet so the suite stays warning-free.
    p <- suppressWarnings(quiet(combine_sigsnps(
        list(EMMAX = sig_long(character(), character(), integer(),
                              character(), character(), numeric())),
        "Union", 200L, PREDICTORS)))
    expect_identical(nrow(a), 0L)
    expect_identical(nrow(p), 0L)
})

# --- 4. thresholds ---------------------------------------------------------
#
# Shared core: the app's compute_method_thresholds() dispatches through
# compute_pval_threshold() (fct_data_loading.R:1055), the same function
# find_sig_snps.R uses. On top of it the app adds a rule layer
# (effective_rule_for, fct_threshold_rules.R:57-73) with NO pipeline analogue.
# Both sides select with <= (sig_snps.R:67-68, fct_data_loading.R:990).

pv_fixture <- function(n = 200L, seed = 42L) {
    set.seed(seed)
    data.table::data.table(
        SNPID = paste0("s", seq_len(n)),
        chr   = "1",
        pos   = seq_len(n),
        bio_1 = c(runif(n - 5L, 0, 1), c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5))
    )
}

# Master rule only: overrides/registry_defaults empty, so the app's extra rule
# layer is bypassed and the comparison is core-vs-core.
app_threshold <- function(pv, type, value) {
    suppressMessages(clinego.app:::compute_method_thresholds(
        list(EMMAX = pv), type, value,
        overrides = list(), registry_defaults = list()
    ))[["bio_1::EMMAX"]]
}

pipe_threshold <- function(pv, type, value) {
    quiet(compute_pval_threshold(pv$bio_1, type, value))$threshold
}

test_that("bonf gives the same cutoff on both sides", {
    skip_without_app()
    pv <- pv_fixture()
    expect_equal(app_threshold(pv, "bonf", 0.05), pipe_threshold(pv, "bonf", 0.05))
    expect_equal(app_threshold(pv, "bonf", 0.05), 0.05 / 200)
})

test_that("top gives the same cutoff on both sides", {
    skip_without_app()
    pv <- pv_fixture()
    for (n in c(1, 5, 20)) {
        expect_equal(app_threshold(pv, "top", n), pipe_threshold(pv, "top", n),
                     info = paste("top", n))
    }
})

test_that("custom gives the same cutoff on both sides", {
    skip_without_app()
    pv <- pv_fixture()
    expect_equal(app_threshold(pv, "custom", 1e-6), pipe_threshold(pv, "custom", 1e-6))
})

test_that("the same cutoff selects the same SNPs on both sides (inclusive <=)", {
    skip_without_app()
    pv <- pv_fixture()
    cut <- app_threshold(pv, "top", 5)
    expect_identical(sum(pv$bio_1 <= cut), 5L)
    expect_equal(cut, pipe_threshold(pv, "top", 5))
})

test_that("too-few-tests refusal reaches the app as NA rather than a bogus cutoff", {
    skip_without_app()
    pv <- data.table::data.table(SNPID = paste0("s", 1:5), chr = "1", pos = 1:5,
                                 bio_1 = c(1e-4, 0.01, 0.2, 0.4, 0.8))
    # compute_pval_threshold refuses qval below 10 tests (pval_threshold.R:68-74)
    # and the app maps a non-"ok" status to NA_real_ (fct_data_loading.R:1058-1059).
    expect_identical(quiet(compute_pval_threshold(pv$bio_1, "qval", 0.1))$status,
                     "too_few_tests")
    expect_true(is.na(app_threshold(pv, "qval", 0.1)))
})

test_that("qvalue is missing from the app's declared Imports (the qval defect)", {
    skip_without_app()
    # STRUCTURAL and attachment-independent, on purpose. helper-libs.R:14
    # attaches qvalue, and namespace lookup falls through to the search path, so
    # the app's qval branch WORKS under this suite while returning NA in
    # production. A live app-vs-pipeline qval assertion would therefore pass and
    # prove nothing. Assert the cause instead: qvalue is absent from
    # DESCRIPTION Imports, so the bare qvalue() call in max_pvalue_fdr()
    # (pval_threshold.R:10) cannot resolve from the package namespace and
    # compute_method_thresholds() catches it into NA (fct_data_loading.R:1054-1057).
    # When someone adds the Import, THIS goes red — that is the prompt to delete
    # the skip() in the next test.
    # Filed 2026-09-10 in docs/pipeline_improvement_requests.md.
    imports <- packageDescription("clinego.app", fields = "Imports")
    expect_false(grepl("\\bqvalue\\b", imports),
                 info = "qvalue now declared — delete the skip() in the next test")
})

test_that("KNOWN DIVERGENCE: qval agrees only because the harness attaches qvalue", {
    # CORRECT behaviour, deliberately not made to pass. In production both app
    # run paths return NA for every qval cell; under this suite they return the
    # right number. Fixing the Import means deleting this skip() line AND
    # inverting the structural assertion above.
    #
    # UN-SKIPPING ALONE IS NOT ENOUGH, and this is the same masking trap one
    # level down: once un-skipped, this assertion still runs in a session where
    # helper-libs.R:14 has attached qvalue, so it passes whether or not the
    # Import fix actually resolves in production. Verify a fix in a session with
    # qvalue NOT attached (or via the package path with a clean search path) —
    # the structural assertion above is what proves the declared dependency, and
    # this one only proves the numbers agree once it resolves.
    # Filed 2026-09-10 in docs/pipeline_improvement_requests.md.
    skip("known divergence: app qval works only under the test harness — filed 2026-09-10")

    skip_without_app()
    pv <- pv_fixture()
    expect_equal(app_threshold(pv, "qval", 0.1), pipe_threshold(pv, "qval", 0.1))
})

test_that("the app's per-cell rule layer has no pipeline analogue", {
    skip_without_app()
    # Divergence by design, pinned so the boundary is explicit: the app can pin
    # one (trait, method) cell to its own rule (fct_threshold_rules.R:57-73)
    # while find_sig_snps.R:35-36 takes a single ADJUST string per invocation
    # and has no per-trait granularity. Any equivalence claim holds ONLY for an
    # empty overrides/registry_defaults.
    pv <- pv_fixture()
    over <- list("bio_1::EMMAX" = list(type = "custom", value = 1e-3))
    got <- suppressMessages(clinego.app:::compute_method_thresholds(
        list(EMMAX = pv), "bonf", 0.05,
        overrides = over, registry_defaults = list()
    ))[["bio_1::EMMAX"]]
    expect_equal(got, 1e-3)
    expect_false(isTRUE(all.equal(got, pipe_threshold(pv, "bonf", 0.05))))
})

# --- 5. region_distance: divergence by omission -----------------------------

test_that("the app cannot express a per-chromosome clumping distance", {
    skip_without_app()
    ns <- asNamespace("clinego.app")
    # region_distance.R is NOT in CLINEGO_SHARED_LIBS (zzz.R:39-42) and the app
    # never calls into it, so `auto_per_chromosome` / `auto_genome_wide` are
    # inexpressible there: every app path passes a scalar
    # (utils_helpers.R:132, mod_gea.R:177-184). This is why every fixture above
    # pins the pipeline to a scalar distance — not a convenience, a
    # precondition for comparability. region_distance.R itself is covered
    # standalone by test-region_distance.R.
    for (f in c("resolve_clumping_distance", "resolve_row_distance",
                "hill_weir", "invert_hill_weir", "loess_fallback")) {
        expect_false(exists(f, envir = ns, inherits = FALSE),
                     info = paste(f, "unexpectedly present in the app namespace"))
    }
    # And the pipeline does have them, so the assertion above is about sharing,
    # not about the functions having been renamed away.
    expect_true(exists("resolve_clumping_distance", mode = "function"))
})

test_that("a scalar distance reaches the shared clusterer identically", {
    skip_without_app()
    # The narrow thing that IS shared: .get_dist's scalar branch (regions.R:12-14).
    # The named-vector branch (:16-21) is unreachable from the app.
    expect_identical(.get_dist(500L, "1"), 500L)
    r <- regions_both_sides(500L)
    expect_identical(nrow(r$pipe), nrow(r$app))
})
