# scripts/R/lib/snp_set_select.R — the ONE selection implementation behind a
# curated SNP set, shared by the Shiny app (save from the GEA tab) and the
# pipeline (GEA.snp_sets / GWAS.snp_sets built by scripts/build_snp_set.R).
#
# That the APP really runs this code is asserted elsewhere: the app's zzz.R gate
# (test-zzz.R) checks the namespace holds these functions, and the app suites'
# combine tests (test-fct_combine.R, Tier 5 Blocks 3/7 in
# test-equivalence-app-pipeline.R) still pass against the moved code. This file
# pins the lib's own contract: inclusive thresholds, per-cell rules, the pair
# filter, the strategies, the collapse, and the declaration/record round trip.
#
# Fixture (custom 1e-4 in every cell selects the starred values):
#              EMMAX bio_1  EMMAX bio_2  LFMM bio_1  LFMM bio_2
#   1:100        1e-6 *       0.3          1e-5 *      0.6
#   1:200        0.5          1e-7 *       0.4         1e-6 *
#   2:300        1e-4 * (=)   0.2          0.8         1e-5 *
#   2:5000       0.9          1e-3         0.7         0.9

ss_pvals <- function() {
    snp <- data.table::data.table(SNPID = c("1:100", "1:200", "2:300", "2:5000"),
                                  chr = c("1", "1", "2", "2"),
                                  pos = c(100L, 200L, 300L, 5000L))
    list(
        EMMAX = cbind(snp, bio_1 = c(1e-6, 0.5, 1e-4, 0.9),
                           bio_2 = c(0.3, 1e-7, 0.2, 1e-3)),
        LFMM  = cbind(snp, bio_1 = c(1e-5, 0.4, 0.8, 0.7),
                           bio_2 = c(0.6, 1e-6, 1e-5, 0.9))
    )
}

ALL_PAIRS <- c("bio_1::EMMAX", "bio_2::EMMAX", "bio_1::LFMM", "bio_2::LFMM")

ss_cells <- function(pairs = ALL_PAIRS, type = "custom", value = 1e-4) {
    stats::setNames(rep(list(list(type = type, value = value)), length(pairs)), pairs)
}

ss_select <- function(pairs = ALL_PAIRS, cells = ss_cells(pairs),
                      strategy = "Union", gap = 0L, pvals = ss_pvals()) {
    quiet(select_snp_set(pvals, cells, pairs, strategy = strategy, gap = gap))
}

# ---------------------------------------------------------------- thresholds

test_that("the cutoff is inclusive: a p equal to the custom value is selected", {
    # One cell only: 2:300 is also selected via LFMM bio_2, which would mask a
    # strict '<' here.
    r <- ss_select(pairs = "bio_1::EMMAX", cells = ss_cells("bio_1::EMMAX"))
    expect_identical(r$snps$SNPID, c("1:100", "2:300"))   # 2:300 p = 1e-4 exactly
    expect_identical(r$cells$n_sig, 2L)
})

test_that("each cell uses its own rule", {
    cells <- ss_cells()
    cells[["bio_2::EMMAX"]] <- list(type = "bonf", value = 0.05)   # 0.05 / 4 = 0.0125
    r <- ss_select(cells = cells)
    expect_true("2:5000" %in% r$snps$SNPID)      # EMMAX bio_2 p = 1e-3
    row <- r$cells[r$cells$pair == "bio_2::EMMAX", ]
    expect_equal(row$threshold, 0.0125)
    expect_identical(row$n_sig, 2L)
})

test_that("top N keeps the N smallest p of the cell", {
    r <- ss_select(pairs = "bio_1::EMMAX",
                   cells = ss_cells("bio_1::EMMAX", type = "top", value = 2))
    expect_identical(r$snps$SNPID, c("1:100", "2:300"))
})

test_that("a rule that cannot run yields no SNPs and says why, without stopping", {
    cells <- ss_cells()
    cells[["bio_1::EMMAX"]] <- list(type = "qval", value = 0.05)   # 4 tests < 10
    r <- ss_select(cells = cells)
    expect_identical(r$cells$status[r$cells$pair == "bio_1::EMMAX"], "too_few_tests")
    expect_identical(r$cells$n_sig[r$cells$pair == "bio_1::EMMAX"], 0L)
    # The other cells still select: 1:100 via LFMM, 2:300 via LFMM bio_2.
    expect_identical(r$snps$SNPID, c("1:100", "1:200", "2:300"))
})

# ---------------------------------------------------------------- pairs

test_that("only the listed pairs contribute, and min_pvalue comes from them alone", {
    r <- ss_select(pairs = setdiff(ALL_PAIRS, "bio_1::EMMAX"))
    expect_identical(r$snps$SNPID, c("1:100", "1:200", "2:300"))
    # 1:100's 1e-6 came from the dropped EMMAX bio_1 cell; LFMM bio_1 gives 1e-5.
    expect_equal(r$snps[SNPID == "1:100", min_pvalue], 1e-5)
})

test_that("filter_to_pairs matches the method exactly, not as a suffix", {
    sig <- list(MLM  = data.table::data.table(SNPID = "a", trait = "bio_1"),
                CMLM = data.table::data.table(SNPID = "b", trait = "bio_1"))
    got <- filter_to_pairs(sig, "bio_1::MLM")
    expect_named(got, "MLM")
})

test_that("filter_to_pairs ignores a malformed key instead of stopping", {
    sig <- list(EMMAX = data.table::data.table(SNPID = "a", trait = "bio_1"))
    expect_named(filter_to_pairs(sig, c("garbage", "::EMMAX", "bio_1::EMMAX")), "EMMAX")
    expect_length(filter_to_pairs(sig, "garbage"), 0L)
})

test_that("split_snp_set_pairs splits at the first ::", {
    sp <- split_snp_set_pairs(c("bio_1::EMMAX", "climate_multivariate::RDA"))
    expect_identical(sp$trait,  c("bio_1", "climate_multivariate"))
    expect_identical(sp$method, c("EMMAX", "RDA"))
    expect_error(split_snp_set_pairs("bio_1"), "malformed")
})

# ---------------------------------------------------------------- strategies

test_that("Union keeps every selected SNP with its minimum p", {
    r <- ss_select()
    expect_equal(r$snps$min_pvalue, c(1e-6, 1e-7, 1e-5))
})

test_that("Cross-method needs a second method nearby, any trait", {
    r <- ss_select(strategy = "Cross-method")
    # 2:300 pairs EMMAX bio_1 with LFMM bio_2: different traits, still kept.
    expect_identical(r$snps$SNPID, c("1:100", "1:200", "2:300"))
    # A single-method SNP is dropped.
    r1 <- ss_select(pairs = c("bio_1::EMMAX", "bio_2::EMMAX", "bio_1::LFMM"),
                    strategy = "Cross-method")
    expect_identical(r1$snps$SNPID, "1:100")
})

test_that("Cross-method per-trait needs the second method on the SAME trait", {
    r <- ss_select(strategy = "Cross-method per-trait")
    expect_identical(r$snps$SNPID, c("1:100", "1:200"))
})

test_that("legacy strategy aliases reach the same branch", {
    expect_identical(ss_select(strategy = "Sum")$snps, ss_select(strategy = "Union")$snps)
    expect_identical(ss_select(strategy = "Overlap")$snps,
                     ss_select(strategy = "Cross-method")$snps)
    expect_identical(ss_select(strategy = "PairOverlap")$snps,
                     ss_select(strategy = "Cross-method per-trait")$snps)
})

test_that("an unknown strategy stops", {
    expect_error(ss_select(strategy = "Nonsense"), "Unknown combine strategy")
    expect_error(ss_select(strategy = "Cross-method >= 3"), "Unknown combine strategy")
})

# ---------------------------------------------------------------- k-of-n

# Calls (long rows) of several methods; one row per (method, SNP, trait).
kn_list <- function(method, chr, pos, trait = "bio_1") {
    dt <- data.table::data.table(SNPID = paste0(chr, ":", pos), chr = chr,
                                 pos = as.integer(pos), pvalue = 1e-6,
                                 method = method, trait = trait)
    split(dt, by = "method", keep.by = TRUE)
}
kn_ids <- function(lst, strategy, gap) {
    sort(unique(combine_sigsnps_long(lst, strategy = strategy, gap = gap)$SNPID))
}
# chr 1: A, B, C inside one 1 kb window. chr 2: a chain, A-B and B-C 800 bp apart but
# A-C 1,600 bp. chr 3: A on two traits plus B: three cells, two methods.
KN <- kn_list(method = c("A", "B", "C", "A", "B", "C", "A", "A", "B"),
              chr    = c("1", "1", "1", "2", "2", "2", "3", "3", "3"),
              pos    = c(1000, 1500, 2000, 1000, 1800, 2600, 500, 500, 600),
              trait  = c(rep("bio_1", 7), "bio_2", "bio_1"))

test_that("Cross-method >=3 counts methods within the gap, the SNP's own call included", {
    expect_identical(kn_ids(KN, "Cross-method >=3", 1000L),
                     c("1:1000", "1:1500", "1:2000", "2:1800"))
    # the gap is inclusive: at 999 bp the chr-1 ends lose the call 1,000 bp away
    expect_identical(kn_ids(KN, "Cross-method >=3", 999L), c("1:1500", "2:1800"))
})

test_that("k-of-n counts methods, not cells, and never chains", {
    ids <- kn_ids(KN, "Cross-method >=3", 1000L)
    expect_false(any(c("3:500", "3:600") %in% ids))      # 3 cells, 2 methods
    expect_false(any(c("2:1000", "2:2600") %in% ids))    # A-B-C chain, A-C beyond the gap
})

test_that("Cross-method >=2 selects exactly what Cross-method does", {
    for (gap in c(0L, 500L, 1000L)) {
        expect_identical(combine_sigsnps_long(KN, "Cross-method >=2", gap),
                         combine_sigsnps_long(KN, "Cross-method", gap), info = gap)
    }
    expect_identical(ss_select(strategy = "Cross-method >=2"),
                     ss_select(strategy = "Cross-method"))
})

test_that("Cross-method >=K keeps every row of a qualifying SNP", {
    long <- combine_sigsnps_long(KN, "Cross-method >=2", 1000L)
    expect_identical(nrow(long[SNPID == "3:500"]), 2L)    # A's bio_1 and bio_2 cells
    expect_identical(unique(long[SNPID == "3:500"]$min_pvalue), 1e-6)
})

# The benchmark's own builder (benchmarks/mvp_build_snp_sets.R) built its sets from
# benchmarks/lib_detection.R: call_by_threshold() per method, then combine_support().
bench_lib <- function() {
    env <- new.env()
    sys.source(file.path(.clinego_root, "benchmarks", "lib_detection.R"), envir = env)
    env
}

test_that("k-of-n gives combine_support()'s n_methods >= k on random calls", {
    bench <- bench_lib()
    for (seed in 1:5) {
        set.seed(seed)
        parts <- lapply(c(A = 40L, B = 30L, C = 50L), function(n)
            unique(paste0(sample(c("1", "2"), n, TRUE), ":", sample.int(20000L, n, TRUE))))
        lst <- data.table::rbindlist(lapply(names(parts), function(m) {
            k <- data.table::tstrsplit(parts[[m]], ":", fixed = TRUE)
            data.table::data.table(SNPID = parts[[m]], chr = k[[1]], pos = as.integer(k[[2]]),
                                   pvalue = 1e-6, method = m, trait = "bio_1")
        }))
        lst <- split(lst, by = "method", keep.by = TRUE)
        for (gap in c(0L, 300L, 2000L)) {
            sup <- bench$combine_support(parts, gap)
            for (k in 2:3) {
                expect_identical(kn_ids(lst, paste0("Cross-method >=", k), gap),
                                 sort(sup[n_methods >= k, snp]),
                                 info = sprintf("seed %d gap %d k %d", seed, gap, k))
            }
        }
    }
})

# Three methods' wide p-value tables on one SNP set: LFMM and EMMAX test bio_1 and
# bio_2, RDA one multivariate column. SNPs within 3 kb of 6 shared loci get small
# p-values in a random half of each (method, trait) cell, so the top calls agree in
# part; one private locus per method gives calls no other method makes.
bench_tables <- function(seed, n = 4000L) {
    set.seed(seed)
    snp <- data.table::data.table(chr = sample(c("1", "2", "3"), n, TRUE))
    snp[, pos := sample.int(300000L, .N), by = chr]
    data.table::setorder(snp, chr, pos)
    snp[, SNPID := paste0(chr, ":", pos)]
    shared  <- data.table::data.table(chr = rep(c("1", "2", "3"), each = 2),
                                      pos = rep(c(50000L, 100000L), 3))
    private <- data.table::data.table(chr = c("1", "2", "3"), pos = 250000L,
                                      method = c("LFMM", "EMMAX", "RDA"))
    within_3kb <- function(loci) snp[, vapply(seq_len(.N), function(i)
        any(loci$chr == chr[i] & abs(loci$pos - pos[i]) <= 3000L), logical(1))]
    near_shared <- within_3kb(shared)
    col <- function(m) {
        p <- runif(n)
        hit <- near_shared & runif(n) < 0.5
        p[hit] <- runif(sum(hit), 1e-9, 1e-5)
        own <- within_3kb(private[method == m]) & runif(n) < 0.2
        p[own] <- runif(sum(own), 1e-12, 1e-10)
        p
    }
    list(LFMM  = cbind(snp[, .(SNPID, chr, pos)], bio_1 = col("LFMM"), bio_2 = col("LFMM")),
         EMMAX = cbind(snp[, .(SNPID, chr, pos)], bio_1 = col("EMMAX"), bio_2 = col("EMMAX")),
         RDA   = cbind(snp[, .(SNPID, chr, pos)], climate_multivariate = col("RDA")))
}

test_that("toppct cells and k-of-n reproduce the benchmark builder's six sets", {
    # mvp_build_snp_sets.R: top N = ceiling(0.25 % x rows) per predictor (RDA 0.5 %),
    # unioned over predictors, then combine_support() at the agreement window.
    bench <- bench_lib()
    old_sets <- function(pv, gap) {
        mult <- c(LFMM = 1, EMMAX = 1, RDA = 2)
        parts <- lapply(stats::setNames(names(mult), names(mult)), function(m) {
            l <- data.table::copy(pv[[m]])
            l[, key := paste(chr, pos, sep = ":")]
            tc <- setdiff(names(l), c("SNPID", "chr", "pos", "key"))
            bench$call_by_threshold(l, tc, "top", ceiling(0.0025 * mult[[m]] * nrow(l)))$called
        })
        sup <- bench$combine_support(parts, gap)
        list(solo_lfmm = parts$LFMM, solo_emmax = parts$EMMAX, solo_rda = parts$RDA,
             union = sup$snp, best = sup[n_methods >= 2L, snp],
             intersect3 = sup[n_methods >= 3L, snp])
    }
    new_sets <- function(pv, gap) {
        sel <- function(methods, strategy) {
            pairs <- unlist(lapply(methods, function(m)
                paste0(setdiff(names(pv[[m]]), c("SNPID", "chr", "pos")), "::", m)))
            cells <- stats::setNames(lapply(pairs, function(p)
                list(type = "toppct", value = if (endsWith(p, "::RDA")) 0.5 else 0.25)), pairs)
            quiet(select_snp_set(pv, cells, pairs, strategy = strategy, gap = gap))$snps$SNPID
        }
        all3 <- c("LFMM", "EMMAX", "RDA")
        list(solo_lfmm = sel("LFMM", "Union"), solo_emmax = sel("EMMAX", "Union"),
             solo_rda = sel("RDA", "Union"), union = sel(all3, "Union"),
             best = sel(all3, "Cross-method"), intersect3 = sel(all3, "Cross-method >=3"))
    }
    for (seed in 1:3) {
        pv <- bench_tables(seed)
        for (gap in c(0L, 5000L)) {
            old <- old_sets(pv, gap)
            new <- new_sets(pv, gap)
            for (nm in names(old)) {
                expect_identical(sort(new[[nm]]), sort(unique(old[[nm]])),
                                 info = sprintf("seed %d gap %d set %s", seed, gap, nm))
            }
            # the fixture must exercise every rule: 3/3 calls exist and are a strict subset
            if (gap == 5000L) {
                expect_gt(length(new$intersect3), 0L)
                expect_lt(length(new$intersect3), length(new$best))
                expect_lt(length(new$best), length(new$union))
            }
        }
    }
})

# ---------------------------------------------------------------- collapse

test_that("zero selected SNPs give a typed 4-column empty table", {
    r <- ss_select(cells = ss_cells(value = 1e-30))
    expect_identical(names(r$snps), c("SNPID", "chr", "pos", "min_pvalue"))
    expect_identical(nrow(r$snps), 0L)
    expect_identical(sum(r$cells$n_sig), 0L)
})

test_that("the collapse order does not depend on the method list's order", {
    pv <- ss_pvals()
    a <- quiet(select_snp_set(pv, ss_cells(), ALL_PAIRS))$snps
    b <- quiet(select_snp_set(rev(pv), ss_cells(), rev(ALL_PAIRS)))$snps
    expect_identical(a, b)
})

test_that("collapse_snp_set breaks a chr/pos tie by SNPID", {
    long <- data.table::data.table(SNPID = c("1:100_B", "1:100_A"), chr = "1",
                                   pos = 100L, min_pvalue = c(1e-5, 1e-6))
    expect_identical(collapse_snp_set(long)$SNPID, c("1:100_A", "1:100_B"))
})

# ---------------------------------------------------------------- loud failures

test_that("select_snp_set stops on a pair it cannot resolve", {
    pv <- ss_pvals()
    expect_error(select_snp_set(pv, ss_cells(), character(0)), "no trait::method pairs")
    expect_error(select_snp_set(pv, ss_cells("bio_1::EMMAX"), "bio_2::EMMAX"),
                 "no rule for pair")
    expect_error(select_snp_set(pv, ss_cells("bio_1::RDA"), "bio_1::RDA"),
                 "no p-value table for method")
    expect_error(select_snp_set(pv, ss_cells("bio_9::EMMAX"), "bio_9::EMMAX"),
                 "trait column missing")
    expect_error(select_snp_set(pv, ss_cells(), c("bio_1::EMMAX", "bio_1::EMMAX")),
                 "duplicate pairs")
})

# ---------------------------------------------------------------- declaration / record

ss_decl <- function(pairs = c("bio_1::EMMAX", "bio_2::LFMM")) {
    list(name = "gea.strict", origin = "config", source_module = "GEA",
         regime = "snp", strategy = "Cross-method", clumping_distance = 100000L,
         pairs = pairs, cells = ss_cells(pairs, type = "bonf", value = 0.0123456789012))
}

test_that("a record round-trips through YAML with full numeric precision", {
    d <- tempfile("ssrec_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    rec <- c(ss_decl(), list(k_best = 3L, n_snps = 12L, created = "2026-10-14T12:00:00"))
    path <- file.path(d, "snp_set.yaml")
    write_snp_set_record(rec, path)
    back <- read_snp_set_decl(path)
    expect_identical(back$pairs, rec$pairs)
    expect_identical(back$cells, rec$cells)       # 0.0123456789012 survives
    expect_identical(back$clumping_distance, 100000L)
    expect_identical(back$k_best, 3L)
    expect_identical(back$created, "2026-10-14T12:00:00")
    # Atomic write: nothing but the record is left in the directory.
    expect_identical(list.files(d), "snp_set.yaml")
})

test_that("a one-pair record keeps pairs as a sequence", {
    d <- tempfile("ssrec_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    path <- file.path(d, "snp_set.yaml")
    write_snp_set_record(ss_decl("bio_1::EMMAX"), path)
    expect_true(any(grepl("^- bio_1::EMMAX$", readLines(path))))
    expect_identical(read_snp_set_decl(path)$pairs, "bio_1::EMMAX")
})

test_that("the declaration fields lead the record, in a fixed order", {
    d <- tempfile("ssrec_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    path <- file.path(d, "snp_set.yaml")
    write_snp_set_record(c(list(n_snps = 1L), ss_decl()), path)
    top <- sub(":.*$", "", grep("^[a-z_]+:", readLines(path), value = TRUE))
    expect_identical(top, c(SNP_SET_DECL_FIELDS, "n_snps"))
})

test_that("a JSON declaration reads into the same normalized form", {
    d <- tempfile("ssdecl_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    path <- file.path(d, "declaration.json")
    # What python's json.dump writes: a one-element list stays a list, numbers
    # are plain floats.
    writeLines('{"cells": {"bio_1::EMMAX": {"type": "custom", "value": 1e-05}},
                 "clumping_distance": 50000, "name": "gea.one", "origin": "config",
                 "pairs": ["bio_1::EMMAX"], "regime": "snp",
                 "source_module": "GEA", "strategy": "Union"}', path)
    got <- read_snp_set_decl(path)
    expect_identical(got$pairs, "bio_1::EMMAX")
    expect_identical(got$cells, list(`bio_1::EMMAX` = list(type = "custom", value = 1e-05)))
    expect_identical(got$clumping_distance, 50000L)
})

test_that("a declaration missing a field or a cell value stops", {
    bad <- ss_decl(); bad$cells <- NULL
    expect_error(normalize_snp_set_decl(bad), "lacks field\\(s\\): cells")
    bad <- ss_decl(); bad$cells[[1]]$value <- NULL
    expect_error(normalize_snp_set_decl(bad), "needs a scalar type and a numeric value")
    bad <- ss_decl(); bad$pairs <- list()
    expect_error(normalize_snp_set_decl(bad), "has no pairs")
    bad <- ss_decl(); bad$pairs <- c("bio_1")
    expect_error(normalize_snp_set_decl(bad), "malformed")
})

test_that("a declaration read from disk drives select_snp_set end to end", {
    d <- tempfile("ssdecl_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
    path <- file.path(d, "snp_set.yaml")
    decl <- ss_decl(ALL_PAIRS)
    decl$cells <- ss_cells()
    decl$strategy <- "Cross-method per-trait"
    write_snp_set_record(decl, path)
    back <- read_snp_set_decl(path)
    r <- quiet(select_snp_set(ss_pvals(), back$cells, back$pairs,
                              back$strategy, back$clumping_distance))
    expect_identical(r$snps$SNPID, c("1:100", "1:200"))
})
