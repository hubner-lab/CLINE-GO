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
    #
    # fct_data_loading.R:533-544 builds gene_id with
    # regmatches(regexpr("(?<=Parent=)...")), which returns ONLY THE ELEMENTS THAT
    # MATCHED rather than one entry per row. So `id` is shorter than nrow(dt),
    # `na_idx <- !nzchar(id)` (:536) never aligns with the rows it is meant to
    # index, and the ID= fallback at :538-540 is unreachable. The pipeline's
    # extract_gene_id (gff_parsing.R:8-13) uses length-preserving str_extract +
    # ifelse and is correct on every shape.
    #
    # SYMPTOM CORRECTED 2026-09-12, measured in-container. The 2026-09-10 filing
    # said the mixed case "throws inside the tryCatch at :552 and returns an EMPTY
    # table". That is only one of THREE outcomes, and not the worst:
    #
    #   all rows have Parent=        -> gene_id is the PARENT value, never the ID=
    #   mixed, exactly ONE Parent=   -> data.table RECYCLES that single id across
    #                                   every row. Right row count, wrong ids, no
    #                                   error, nothing stale-looking. The bad one.
    #   mixed, k>1 Parent=, k!=nrow  -> assignment refused -> empty table (as filed)
    #   no row has Parent=           -> gene_id is NA for every row
    #
    # The recycling case is the reason this block asserts all four shapes rather
    # than just the empty-table one: a test that only checked nrow would PASS on
    # the silent-corruption case.
    # Filed 2026-09-10, corrected 2026-09-12, docs/pipeline_improvement_requests.md.
    skip("known divergence: app load_gff_genes gene_id fallback unreachable — filed 2026-09-10")

    skip_without_app()
    # CALL SHAPE — corrected 2026-09-12, a TEST-code fix. This block previously
    # called load_gff_genes(path, "gene"), but the real signature is
    # load_gff_genes(project, config) (fct_data_loading.R:502): it resolves the GFF
    # as file.path(get_pipeline_path(), Input$dir, Input$gff) and reads GFF$feature
    # from the config. With the old call, deleting the skip() above would have
    # produced an arity/argument error rather than the empty table the defect
    # actually causes — i.e. the skip was hiding a broken test, and "un-skip and
    # watch it fail" would have proved nothing about the defect.
    #
    # The unique project name is mandatory, not hygiene: load_cached's key here is
    # "gff_genes_<project>" with NO fingerprint (:503-504), so it is sticky for the
    # whole session.
    root <- withr::local_tempdir()
    withr::local_options(clinego.pipeline_path = root)
    dir.create(file.path(root, "data"), recursive = TRUE, showWarnings = FALSE)

    # A UNIQUE project per call is mandatory, not hygiene: load_gff_genes' cache key
    # is "gff_genes_<project>" with NO fingerprint (fct_data_loading.R:503, filed),
    # so it is sticky for the whole session.
    app_gene_ids <- function(lines, tag) {
        writeLines(c("##gff-version 3", lines),
                   file.path(root, "data", paste0(tag, ".gff3")))
        cfg <- list(Input = list(dir = "data", gff = paste0(tag, ".gff3")),
                    GFF   = list(feature = "gene"))
        quiet(clinego.app:::load_gff_genes(
            paste0("equiv_gff_", tag, "_",
                   as.integer(stats::runif(1, 1, 1e6))), cfg))
    }

    # 1. Mixed, one Parent=: currently recycles "t2" into both rows.
    mixed1 <- app_gene_ids(c(
        gff_line("1", 1000, 2000, "ID=g1;Name=alpha"),
        gff_line("1", 4000, 5000, "Parent=t2;ID=g2;Name=beta")), "mixed1")
    expect_identical(nrow(mixed1), 2L)
    expect_setequal(mixed1$gene_id, c("g1", "g2"))

    # 2. Mixed, two Parent= among three rows: currently an empty table.
    mixed2 <- app_gene_ids(c(
        gff_line("1", 1000, 2000, "ID=g1"),
        gff_line("1", 3000, 4000, "Parent=t2;ID=g2"),
        gff_line("1", 5000, 6000, "Parent=t3;ID=g3")), "mixed2")
    expect_identical(nrow(mixed2), 3L)
    expect_setequal(mixed2$gene_id, c("g1", "g2", "g3"))

    # 3. No Parent= anywhere: currently NA for every row.
    noparent <- app_gene_ids(c(
        gff_line("1", 1000, 2000, "ID=g1"),
        gff_line("1", 3000, 4000, "ID=g2")), "noparent")
    expect_identical(nrow(noparent), 2L)
    expect_false(any(is.na(noparent$gene_id)))
    expect_setequal(noparent$gene_id, c("g1", "g2"))
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

# suppressWarnings, not just quiet(), on BOTH sides. When a strategy selects
# NOTHING, each implementation's rbindlist over an empty list produces a table
# that data.table then flags with its "shallow copy ... := can add or remove
# columns by reference" notice. The wart is symmetric — pipeline
# combine_sigsnps.R:166 and app fct_combine.R:121 — belongs to the empty path
# rather than to anything these blocks measure, and this suite's warning
# baseline is zero, so letting it through would leave standing noise for the
# next real warning to hide in.
app_combine <- function(sig_list, strategy, gap) {
    suppressWarnings(quiet(clinego.app:::combine_sigsnps(sig_list, strategy, gap)))
}

pipe_combine <- function(sig_list, strategy, distance, predictors = PREDICTORS) {
    suppressWarnings(quiet(combine_sigsnps(sig_list, strategy, distance, predictors)))
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

# --- 6. Cross-method per-trait ----------------------------------------------
#
# Block 3 above compares Union, single-method and plain Cross-method. The
# per-trait branch — pipeline .overlap_cross_method(per_trait = TRUE)
# (combine_sigsnps.R:132-167) against the app's loop (fct_combine.R:99-121) —
# has never been compared, only its ALIAS was (block 3's first test).
#
# SCOPE NOTE, because the filing is easy to over-read: the superset defect
# skip()ped at the Cross-method block is the plain Cross-method branch ONLY
# (fct_combine.R:97 returns all_snps[SNPID %in% selected_ids], every row for a
# selected SNP). The per-trait branch returns unique(rbindlist(selected_rows)) —
# matched rows only — and is structurally correct. Do not read that filing as
# covering this branch.

# The shared builders use PREDICTORS = c("bio_1", "bio_2"), which covers their
# whole trait set — so every existing block is structurally blind to the
# pipeline's `predictors` whitelist (combine_sigsnps.R:51). This builder adds a
# trait OUTSIDE that vector. It is deliberately a NEW builder rather than a
# widening of two_method_fixture(): the shared one is used by blocks 1-4, and
# changing it would silently change what they assert.
fixture_with_extra_trait <- function() {
    list(
        EMMAX = sig_long(c("1:1000", "1:2000", "3:7000"),
                         c("1", "1", "3"),
                         c(1000L, 2000L, 7000L),
                         c("bio_1", "bio_1", "bio_99"),   # bio_99 is NOT in PREDICTORS
                         "EMMAX",
                         c(1e-8, 1e-7, 1e-9)),
        LFMM  = sig_long(c("1:2000", "3:7000"),
                         c("1", "3"),
                         c(2000L, 7000L),
                         c("bio_1", "bio_99"),
                         "LFMM",
                         c(1e-9, 1e-10))
    )
}

# two_method_fixture() yields ZERO rows under the per-trait rule, and that is
# correct: its only cross-method co-location (1:2000) is bio_1 in EMMAX and
# bio_2 in LFMM, so no trait is shared by both methods anywhere. The vacuity
# guard caught that on the first run — "both sides agree" over an empty result
# would have asserted nothing. Hence a builder with a genuine SAME-TRAIT,
# both-methods, within-gap overlap.
fixture_per_trait_overlap <- function() {
    list(
        EMMAX = sig_long(c("1:1000", "1:2000", "2:5000"),
                         c("1", "1", "2"),
                         c(1000L, 2000L, 5000L),
                         c("bio_1", "bio_1", "bio_2"),
                         "EMMAX",
                         c(1e-8, 1e-7, 1e-6)),
        LFMM  = sig_long(c("1:1500", "2:5200"),
                         c("1", "2"),
                         c(1500L, 5200L),
                         c("bio_1", "bio_2"),      # same traits as EMMAX's, nearby
                         "LFMM",
                         c(1e-9, 1e-5))
    )
}

test_that("Cross-method per-trait selects the same SNPs on both sides", {
    skip_without_app()
    gap <- 1000L
    pipe <- pipe_combine(fixture_per_trait_overlap(), "Cross-method per-trait", gap)
    app  <- app_combine(fixture_per_trait_overlap(), "Cross-method per-trait", gap)

    expect_gt(nrow(app), 0L)   # never green on nothing
    expect_gt(nrow(pipe), 0L)
    expect_setequal(unique(app$SNPID), unique(pipe$SNPID))
})

test_that("Cross-method per-trait returns nothing when no trait is shared across methods", {
    skip_without_app()
    # The zero case, asserted deliberately rather than left as an accident of a
    # fixture: both sides must agree that an empty result is the right answer.
    gap  <- 1000L
    pipe <- pipe_combine(two_method_fixture(), "Cross-method per-trait", gap)
    app  <- app_combine(two_method_fixture(), "Cross-method per-trait", gap)
    expect_identical(nrow(pipe), 0L)
    expect_identical(nrow(app), 0L)
})

test_that("Cross-method per-trait requires the SAME trait on both methods", {
    skip_without_app()
    # This is what separates per-trait from plain Cross-method. In the shared
    # fixture 1:2000 is EMMAX/bio_1 and LFMM/bio_2 — a cross-method overlap at
    # the same position but under DIFFERENT traits, so the per-trait rule must
    # reject it while plain Cross-method accepts it.
    gap <- 1000L
    per_trait <- pipe_combine(two_method_fixture(), "Cross-method per-trait", gap)
    plain     <- pipe_combine(two_method_fixture(), "Cross-method", gap)
    expect_false("1:2000" %in% per_trait$SNPID)
    expect_true("1:2000" %in% plain$SNPID)

    app_pt <- app_combine(two_method_fixture(), "Cross-method per-trait", gap)
    expect_false("1:2000" %in% app_pt$SNPID)
})

test_that("Cross-method per-trait agrees on min_pvalue for every shared SNP", {
    skip_without_app()
    gap  <- 1000L
    pipe <- pipe_combine(fixture_with_extra_trait(), "Cross-method per-trait", gap,
                         c(PREDICTORS, "bio_99"))
    app  <- app_combine(fixture_with_extra_trait(), "Cross-method per-trait", gap)

    app_min <- unique(app[, .(SNPID, min_pvalue)])
    m <- merge(pipe[, .(SNPID, min_pvalue)], app_min, by = "SNPID",
               suffixes = c("_pipe", "_app"))
    expect_gt(nrow(m), 0L)
    expect_equal(m$min_pvalue_pipe, m$min_pvalue_app)
})

test_that("the pipeline's `predictors` whitelist has NO app analogue", {
    skip_without_app()
    # combine_sigsnps.R:51 filters to `predictors`; the app filters nothing. With
    # bio_99 excluded pipeline-side, the pipeline drops 3:7000 and the app keeps
    # it. Every other block is blind to this because their PREDICTORS covers the
    # whole trait set.
    gap <- 1000L
    pipe <- pipe_combine(fixture_with_extra_trait(), "Cross-method per-trait", gap)  # bio_99 NOT allowed
    app  <- app_combine(fixture_with_extra_trait(), "Cross-method per-trait", gap)

    expect_false("3:7000" %in% pipe$SNPID)
    expect_true("3:7000" %in% app$SNPID)
})

test_that("both sides iterate unordered method pairs to the same result set", {
    skip_without_app()
    # The pipeline fixed its double loop to unordered pairs
    # (combine_sigsnps.R:136-142); the app still runs the ordered loop AND calls
    # both directions inside it (fct_combine.R:104-112). Same answer, 4x the
    # foverlaps calls for two methods. Asserted as agreement so that porting the
    # optimisation cannot silently change the result.
    gap <- 1000L
    swapped <- rev(two_method_fixture())
    a <- app_combine(two_method_fixture(), "Cross-method per-trait", gap)
    b <- app_combine(swapped, "Cross-method per-trait", gap)
    expect_setequal(unique(a$SNPID), unique(b$SNPID))
})

# --- 7. SNP-set store: a CONTRACT, not two implementations -------------------
#
# The app WRITES the store and Snakemake READS it. Nothing computes the same
# thing twice here, so this block pins an interface: the four TSV columns, and
# the glob that has to find them. Both sides currently agree — these assertions
# exist so a one-sided change cannot pass unnoticed.
#
# Structural assertions over source text, in the manner of the qvalue/DESCRIPTION
# assertion in block 4: the pipeline halves are a shell script and a .smk file,
# neither loadable into R.

SNP_SET_COLS <- c("SNPID", "chr", "pos", "min_pvalue")

test_that("the app writes exactly the four columns promote_snp_set.R requires", {
    skip_without_app()
    withr::local_options(clinego.pipeline_path = withr::local_tempdir())
    project <- basename(tempfile("EQ_SS_"))

    dt <- data.table::data.table(
        SNPID = c("1:100", "1:100", "2:300"), chr = c("1", "1", "2"),
        pos = c(100L, 100L, 300L), pvalue = c(1e-8, 3e-9, 5e-7),
        method = c("EMMAX", "LFMM", "EMMAX"), trait = "bio_1",
        region_id = NA_character_)
    dt[, min_pvalue := min(pvalue, na.rm = TRUE), by = "SNPID"]

    n <- clinego.app:::save_snp_set(project, "setA", dt[], list(source_module = "GEA"))
    expect_identical(n, 2L)

    written <- data.table::fread(clinego.app:::snp_set_path(project, "setA"), sep = "\t")
    expect_identical(names(written), SNP_SET_COLS)

    # The reader's own declared requirement, read from its source.
    src <- readLines(file.path(.clinego_root, "scripts", "promote_snp_set.R"))
    req <- grep("^required_cols <- ", src, value = TRUE)
    expect_length(req, 1L)
    for (col in SNP_SET_COLS) expect_true(grepl(paste0('"', col, '"'), req), info = col)
})

test_that("the app's store layout matches the glob common.smk resolves sets with", {
    skip_without_app()
    withr::local_options(clinego.pipeline_path = withr::local_tempdir())
    project <- basename(tempfile("EQ_SS2_"))
    dt <- data.table::data.table(SNPID = "1:100", chr = "1", pos = 100L,
                                 min_pvalue = 1e-8)
    clinego.app:::save_snp_set(project, "setA", dt, list())

    written <- clinego.app:::snp_set_path(project, "setA")
    # Snakemake globs "{INTER}snp_sets/*/selected_snps.tsv" and recovers the set
    # name with basename(dirname(p)) — common.smk:1423-1428.
    expect_identical(basename(written), "selected_snps.tsv")
    expect_identical(basename(dirname(written)), "setA")
    expect_identical(basename(dirname(dirname(written))), "snp_sets")

    smk <- readLines(file.path(.clinego_root, "workflow", "rules", "common.smk"))
    expect_true(any(grepl('snp_sets/', smk, fixed = TRUE)))
    expect_true(any(grepl('*/selected_snps.tsv', smk, fixed = TRUE)))
})

test_that("the manifest cannot be mistaken for a set by that glob", {
    skip_without_app()
    withr::local_options(clinego.pipeline_path = withr::local_tempdir())
    project <- basename(tempfile("EQ_SS3_"))
    clinego.app:::save_snp_set(project, "setA",
        data.table::data.table(SNPID = "1:100", chr = "1", pos = 100L,
                               min_pvalue = 1e-8), list())
    # manifest.json is a FILE at the store root, so "*/selected_snps.tsv" cannot
    # reach it. Pinned because moving it into a subdirectory would silently
    # register a bogus set named after that directory.
    man <- clinego.app:::snp_sets_manifest_path(project)
    expect_identical(basename(dirname(man)), "snp_sets")
    expect_identical(basename(man), "manifest.json")
})

test_that("delete_snp_set's method list matches the maladaptation registry", {
    skip_without_app()
    # fct_snp_sets.R:112 hardcodes the methods whose result dirs it cleans. A new
    # pipeline method added to workflow/methods/maladaptation.py and not here
    # silently leaks result directories on every delete.
    src <- readLines(file.path(.clinego_root, "scripts", "clinego.app", "R",
                               "fct_snp_sets.R"))
    line <- grep("all_methods <- ", src, value = TRUE)
    expect_length(line, 1L)

    py <- readLines(file.path(.clinego_root, "workflow", "methods", "maladaptation.py"))
    registry <- sub("^\\s*[\"']([A-Za-z_]+)[\"']\\s*:\\s*\\{.*$", "\\1",
                    grep("^\\s*[\"'][A-Za-z_]+[\"']\\s*:\\s*\\{", py, value = TRUE))
    registry <- unique(registry[nzchar(registry)])
    expect_gt(length(registry), 0L)
    # CONFIRMED DRIFT, 2026-09-12: maladaptation.py registers FOUR methods
    # (gradient_forest, geometric_offset, rda_offset, rda_offset_corrected) and
    # fct_snp_sets.R:112 lists three. Deleting a SNP set therefore orphans every
    # rda_offset_corrected result directory for that set. Filed; the assertion
    # below is the correct behaviour and is skip()ped until the list is fixed.
    skip("known bug: rda_offset_corrected missing from delete_snp_set's all_methods — filed 2026-09-12")
    for (m in registry) {
        expect_true(grepl(m, line, fixed = TRUE),
                    info = paste("maladaptation method", m,
                                 "is not in delete_snp_set's all_methods"))
    }
})

# --- 8. coords JSON: the Manhattan hand-off ---------------------------------
#
# The pipeline renders a static background PNG and writes the axis geometry to a
# *_coords.json beside it; the app reads that JSON and positions a plotly overlay
# on top. Zero shared function names and two different plotting stacks, so
# expect_identical(body(...)) does not apply — the assertion has to be that the
# NUMBERS survive the round trip.
#
# The JSON is written and read here exactly as production does it
# (plot_manhattan.R:246-262 / fct_manhattan.R:5-15), through a real file, so the
# setNames -> toJSON -> unlist name-preservation chain is actually covered.

manhattan_fixture <- function() {
    data.table::data.table(
        SNPID = c("1:100", "1:900", "2:50", "2:800"),
        chr   = c("1", "1", "2", "2"),
        pos   = c(100L, 900L, 50L, 800L),
        pvalue = c(1e-8, 1e-3, 1e-6, 1e-2))
}

write_coords_like_pipeline <- function(chr_info, path) {
    # Verbatim the shape plot_manhattan.R:247-249 builds. as.list is load
    # bearing: a bare named numeric vector under auto_unbox = TRUE serialises as
    # an unnamed array and the chromosome keys are lost.
    coords <- list(
        chr_offsets  = stats::setNames(as.list(chr_info$tot), as.character(chr_info$chr_f)),
        chr_lengths  = stats::setNames(as.list(chr_info$chr_len), as.character(chr_info$chr_f)),
        gap_fraction = 0.02,
        x_range      = c(0, max(chr_info$tot + chr_info$chr_len)),
        y_range      = c(0, 10))
    jsonlite::write_json(coords, path, auto_unbox = TRUE, digits = 6)
    path
}

test_that("chr_offsets survive the JSON round trip keyed by chromosome", {
    skip_without_app()
    prep <- prepare_manhattan_data(manhattan_fixture(), pval_col = "pvalue")
    path <- withr::local_tempfile(fileext = ".json")
    write_coords_like_pipeline(prep$chr_info, path)

    coords <- clinego.app:::safe_read_json(path)
    offsets <- unlist(coords$chr_offsets)
    expect_identical(names(offsets), as.character(prep$chr_info$chr_f))
    expect_equal(unname(offsets), prep$chr_info$tot, tolerance = 1e-6)
})

test_that("the app's chr_midpoints reproduces the pipeline's `center`", {
    skip_without_app()
    # center is computed in chr_info (manhattan_utils.R:41) but never serialised,
    # so the app RE-DERIVES it from the two transported vectors
    # (fct_manhattan.R:14). This is the one quantity the app recomputes rather
    # than reads, which is exactly why it needs an equivalence assertion.
    prep <- prepare_manhattan_data(manhattan_fixture(), pval_col = "pvalue")
    path <- withr::local_tempfile(fileext = ".json")
    write_coords_like_pipeline(prep$chr_info, path)

    coords <- clinego.app:::safe_read_json(path)
    expect_equal(unname(clinego.app:::chr_midpoints(coords)),
                 prep$chr_info$center, tolerance = 1e-6)
})

test_that("the app's add_cum_pos reproduces the pipeline's pos_cum", {
    skip_without_app()
    prep <- prepare_manhattan_data(manhattan_fixture(), pval_col = "pvalue")
    path <- withr::local_tempfile(fileext = ".json")
    write_coords_like_pipeline(prep$chr_info, path)
    coords <- clinego.app:::safe_read_json(path)

    sig <- data.table::data.table(
        SNPID = prep$data$SNPID, chr = prep$data$chr,
        pos = prep$data$pos, pvalue = prep$data$pvalue)
    # add_cum_pos COPIES and returns (fct_manhattan.R:490) — it does not mutate
    # by reference, despite using `:=` internally. The result must be assigned.
    sig <- clinego.app:::add_cum_pos(sig, coords)

    m <- merge(sig[, .(SNPID, cum_pos, log10p)],
               data.table::as.data.table(prep$data)[, .(SNPID, pos_cum, log10p_pipe = log10p)],
               by = "SNPID")
    expect_identical(nrow(m), 4L)
    expect_equal(m$cum_pos, m$pos_cum, tolerance = 1e-6)
    expect_equal(m$log10p, m$log10p_pipe)
})

test_that("a chromosome missing from the JSON keys yields NA, silently", {
    skip_without_app()
    # chr_offsets[chr] is a NAME lookup, so an unknown chromosome gives NA and
    # the point simply vanishes from the plotly overlay. Nothing in the app
    # asserts key-set containment; pinned here as the live behaviour, with the
    # missing guard filed.
    prep <- prepare_manhattan_data(manhattan_fixture(), pval_col = "pvalue")
    path <- withr::local_tempfile(fileext = ".json")
    write_coords_like_pipeline(prep$chr_info, path)
    coords <- clinego.app:::safe_read_json(path)

    sig <- data.table::data.table(SNPID = "9:1", chr = "9", pos = 1L, pvalue = 1e-5)
    sig <- clinego.app:::add_cum_pos(sig, coords)
    expect_true(is.na(sig$cum_pos))
})

test_that("safe_read_json degrades to NULL rather than erroring", {
    skip_without_app()
    expect_null(clinego.app:::safe_read_json(file.path(tempdir(), "no-such-coords.json")))
    bad <- withr::local_tempfile(fileext = ".json")
    writeLines("{ not json", bad)
    expect_null(clinego.app:::safe_read_json(bad))
})

# --- 9. assign_region_ids: the app recomputes what the pipeline wrote --------
#
# fct_data_loading.R:166-209 does NOT read the region_id column out of
# regions_combined.tsv — it re-derives membership with its own
# foverlaps(type = "within", mult = "first") over the region bounds. So the two
# sides agree only as long as that recomputation matches the clustering that
# produced the file.

test_that("regions_combined regions provably cannot overlap", {
    # The precondition that makes mult = "first" safe, asserted rather than
    # assumed. Clusters split when pB - pA > 2*dist and bounds are
    # [min-dist, max+dist], so consecutive regions are separated by
    # (pB - pA) - 2*dist > 0; the pmax(1L, ...) clamp only raises a start.
    dist <- 1000L
    # cluster_snps_to_regions reads min_pvalue (regions.R:94), which lives on the
    # WIDE table, so the input is derived with the pipeline's own combine step —
    # the same reason wide_from_long() exists for block 1. Handing it the long
    # shape yields Inf and a warning instead.
    snps <- sig_long(c("1:1000", "1:1500", "1:9000", "1:9500"),
                     "1", c(1000L, 1500L, 9000L, 9500L), "bio_1", "EMMAX",
                     c(1e-8, 1e-7, 1e-6, 1e-5))
    regs <- quiet(cluster_snps_to_regions(wide_from_long(list(EMMAX = snps), dist), dist))
    data.table::setorder(regs, chr, start)
    expect_gt(nrow(regs), 1L)
    expect_true(all(regs$start[-1] > regs$end[-nrow(regs)]))
})

test_that("the app's recomputation reproduces the pipeline's region_id", {
    skip_without_app()
    withr::local_options(clinego.pipeline_path = withr::local_tempdir())
    project <- basename(tempfile("EQ_RID_"))
    dist <- 1000L

    snps <- sig_long(c("1:1000", "1:1500", "1:9000", "2:400"),
                     c("1", "1", "1", "2"),
                     c(1000L, 1500L, 9000L, 400L),
                     "bio_1", "EMMAX", c(1e-8, 1e-7, 1e-6, 1e-9))
    regs <- quiet(cluster_snps_to_regions(wide_from_long(list(EMMAX = snps), dist), dist))

    out <- clinego.app:::regions_combined_path(project, "GEA")
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(regs, out, sep = "\t")

    # The pipeline's own answer: which region each SNP was clustered into.
    truth <- regs[, .(region_id, ids = strsplit(snp_ids, ",", fixed = TRUE))][
        , .(SNPID = unlist(ids)), by = region_id]

    got <- data.table::copy(snps)
    clinego.app:::assign_region_ids(got, project, "GEA")

    m <- merge(truth, got[, .(SNPID, app_region = region_id)], by = "SNPID")
    expect_identical(nrow(m), nrow(truth))
    expect_identical(m$app_region, m$region_id)
})

test_that("a SNP outside every region gets NA, not a neighbouring region", {
    skip_without_app()
    withr::local_options(clinego.pipeline_path = withr::local_tempdir())
    project <- basename(tempfile("EQ_RID2_"))

    snps <- sig_long(c("1:1000", "1:1500"), "1", c(1000L, 1500L),
                     "bio_1", "EMMAX", c(1e-8, 1e-7))
    regs <- quiet(cluster_snps_to_regions(wide_from_long(list(EMMAX = snps), 1000L), 1000L))
    out <- clinego.app:::regions_combined_path(project, "GEA")
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(regs, out, sep = "\t")

    far <- sig_long("1:99999", "1", 99999L, "bio_1", "EMMAX", 1e-8)
    clinego.app:::assign_region_ids(far, project, "GEA")
    expect_true(is.na(far$region_id))
})

test_that("assign_region_ids and its fct_overlap twin are the same code", {
    skip_without_app()
    # fct_data_loading.R:199-207 and fct_overlap.R:263-271 are byte-equivalent
    # setkey / foverlaps(type="within", mult="first") / i.region_id blocks with
    # nothing keeping them in sync. Assert the shared shape so a fix to one that
    # is not applied to the other shows up here.
    ns <- asNamespace("clinego.app")
    expect_true(exists("assign_region_ids", envir = ns, inherits = FALSE))
    expect_true(exists("assign_region_ids_from_regions", envir = ns, inherits = FALSE))
    for (f in c("assign_region_ids", "assign_region_ids_from_regions")) {
        src <- paste(deparse(body(get(f, envir = ns))), collapse = " ")
        expect_true(grepl('type = "within"', src, fixed = TRUE), info = f)
        expect_true(grepl('mult = "first"', src, fixed = TRUE), info = f)
        expect_true(grepl("i.region_id", src, fixed = TRUE), info = f)
    }
})
