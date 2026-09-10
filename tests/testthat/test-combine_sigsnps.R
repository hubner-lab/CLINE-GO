# scripts/R/lib/combine_sigsnps.R
#
# Merges the per-method significant-SNP tables into the single selected_snps
# table every downstream region/gene step reads.

# --- fixture ---------------------------------------------------------------

# BUILDER: combine_sigsnps() does `dt[, chr := as.character(chr)]` on the list
# elements, which mutates them by reference — a shared fixture would carry that
# between tests.
#
# Hand-countable 3-method layout (positions chosen so that "same locus" means
# "identical position", i.e. clumping_distance = 0 is decisive):
#
#   SNPID    chr:pos    EMMAX    LFMM     RDA      in N methods
#   A        1:1000     bio_1    bio_1    -        2, same trait
#   B        1:2000     bio_2    -        -        1
#   C        1:3000     -        bio_1    bio_2    2, different traits
#   D        2:1000     -        -        bio_1    1
#   E        1:5000     bio_99   -        -        1, trait not in predictors
three_methods <- function() {
    row <- function(SNPID, chr, pos, trait, method, pvalue) {
        data.table::data.table(SNPID = SNPID, chr = as.character(chr),
                               pos = as.integer(pos), trait = trait,
                               method = method, pvalue = pvalue)
    }
    list(
        EMMAX = rbind(row("A", 1, 1000, "bio_1",  "EMMAX", 1e-8),
                      row("B", 1, 2000, "bio_2",  "EMMAX", 1e-7),
                      row("E", 1, 5000, "bio_99", "EMMAX", 1e-9)),
        LFMM  = rbind(row("A", 1, 1000, "bio_1",  "LFMM",  1e-6),
                      row("C", 1, 3000, "bio_1",  "LFMM",  1e-5)),
        RDA   = rbind(row("C", 1, 3000, "bio_2",  "RDA",   1e-4),
                      row("D", 2, 1000, "bio_1",  "RDA",   1e-3))
    )
}

PREDICTORS <- c("bio_1", "bio_2")

# .empty_combined_snps_dt() builds its method columns with `dt[[m]] <- character()`
# and then adds min_pvalue with `:=`, which makes data.table emit its "shallow copy
# was taken" warning. Harmless and internal to the library, but it fires on every
# empty-path call. Asserted once below, suppressed elsewhere so it does not bury
# the suite output.
quiet_empty <- function(expr) suppressWarnings(suppressMessages(expr))

# --- .normalise_strategy ---------------------------------------------------

test_that(".normalise_strategy maps every legacy alias to its canonical name", {
    expect_identical(.normalise_strategy("All"), "Union")
    expect_identical(.normalise_strategy("Sum"), "Union")
    expect_identical(.normalise_strategy("Overlap"), "Cross-method")
    expect_identical(.normalise_strategy("MethodOverlap"), "Cross-method per-trait")
    expect_identical(.normalise_strategy("PairOverlap"), "Cross-method per-trait")
})

test_that(".normalise_strategy passes canonical and single-method names through", {
    expect_identical(.normalise_strategy("Union"), "Union")
    expect_identical(.normalise_strategy("Cross-method"), "Cross-method")
    expect_identical(.normalise_strategy("EMMAX"), "EMMAX")
})

# --- strategies ------------------------------------------------------------

test_that("Union keeps every SNP significant in any method", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    expect_setequal(r$SNPID, c("A", "B", "C", "D"))   # E filtered out by predictors
})

test_that("the Sum and All aliases behave identically to Union", {
    u   <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    for (alias in c("Sum", "All")) {
        a <- quiet(combine_sigsnps(three_methods(), alias, 0L, PREDICTORS))
        expect_identical(a$SNPID, u$SNPID)
    }
})

test_that("Cross-method keeps only SNPs supported by >= 2 methods", {
    # A is EMMAX + LFMM; C is LFMM + RDA. B and D are single-method.
    r <- quiet(combine_sigsnps(three_methods(), "Cross-method", 0L, PREDICTORS))
    expect_setequal(r$SNPID, c("A", "C"))
})

test_that("Cross-method per-trait additionally requires the methods to share a trait", {
    # C is 2-method but LFMM calls it bio_1 and RDA calls it bio_2, so it drops
    # out; A is bio_1 in both methods and survives. This is the one assertion
    # that separates the two overlap strategies.
    r <- quiet(combine_sigsnps(three_methods(), "Cross-method per-trait", 0L, PREDICTORS))
    expect_setequal(r$SNPID, "A")
})

test_that("a single method name is a passthrough of that method's SNPs", {
    r <- quiet(combine_sigsnps(three_methods(), "EMMAX", 0L, PREDICTORS))
    expect_setequal(r$SNPID, c("A", "B"))             # E dropped by predictors
})

test_that("clumping_distance widens Cross-method from identical positions to nearby ones", {
    # At distance 0 only A and C pair up. At 1500 bp, EMMAX's B (1:2000) is within
    # reach of LFMM's A (1:1000) and C (1:3000), so B joins the set.
    near <- quiet(combine_sigsnps(three_methods(), "Cross-method", 1500L, PREDICTORS))
    expect_true("B" %in% near$SNPID)
    expect_false("D" %in% near$SNPID)                 # different chromosome
})

test_that("an unknown strategy errors rather than silently returning something", {
    expect_error(quiet(combine_sigsnps(three_methods(), "Nonsense", 0L, PREDICTORS)),
                 "Unknown strategy")
})

# --- output shape ----------------------------------------------------------

test_that("columns are SNPID, chr, pos, the methods in list order, then min_pvalue", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    expect_identical(colnames(r),
                     c("SNPID", "chr", "pos", "EMMAX", "LFMM", "RDA", "min_pvalue"))
})

test_that("each method cell is a comma-separated sorted trait list, empty when absent", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    expect_identical(r[SNPID == "A"]$EMMAX, "bio_1")
    expect_identical(r[SNPID == "A"]$LFMM,  "bio_1")
    expect_identical(r[SNPID == "A"]$RDA,   "")       # not NA
    expect_identical(r[SNPID == "C"]$RDA,   "bio_2")
    expect_false(any(is.na(unlist(r[, .(EMMAX, LFMM, RDA)]))))
})

test_that("multiple traits in one method collapse into one sorted comma string", {
    lst <- three_methods()
    lst$EMMAX <- rbind(lst$EMMAX,
                       data.table::data.table(SNPID = "A", chr = "1", pos = 1000L,
                                              trait = "bio_2", method = "EMMAX",
                                              pvalue = 1e-9))
    r <- quiet(combine_sigsnps(lst, "Union", 0L, PREDICTORS))
    expect_identical(r[SNPID == "A"]$EMMAX, "bio_1,bio_2")
})

test_that("one row per SNPID, and min_pvalue is the smallest across methods", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    expect_identical(anyDuplicated(r$SNPID), 0L)
    expect_equal(r[SNPID == "A"]$min_pvalue, 1e-8)    # EMMAX 1e-8 beats LFMM 1e-6
})

test_that("rows are sorted by chr then pos", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, PREDICTORS))
    expect_identical(r$SNPID, c("A", "B", "C", "D"))
})

test_that("predictors filters out traits that were not requested", {
    r <- quiet(combine_sigsnps(three_methods(), "Union", 0L, "bio_1"))
    expect_setequal(r$SNPID, c("A", "C", "D"))        # B is bio_2 only
    expect_identical(r[SNPID == "C"]$RDA, "")         # RDA called C bio_2
})

# --- empty paths -----------------------------------------------------------

test_that("the empty-table constructor warns about a data.table shallow copy", {
    lst <- lapply(three_methods(), function(dt) dt[0])
    expect_warning(quiet(combine_sigsnps(lst, "Union", 0L, PREDICTORS)),
                   "shallow copy")
})

test_that("no SNPs in any method returns the empty schema keyed by method names", {
    lst <- lapply(three_methods(), function(dt) dt[0])
    r <- quiet_empty(combine_sigsnps(lst, "Union", 0L, PREDICTORS))
    expect_identical(nrow(r), 0L)
    expect_identical(colnames(r),
                     c("SNPID", "chr", "pos", "EMMAX", "LFMM", "RDA", "min_pvalue"))
})

test_that("no SNP surviving the predictors filter returns the empty schema", {
    r <- quiet_empty(combine_sigsnps(three_methods(), "Union", 0L, "no_such_trait"))
    expect_identical(nrow(r), 0L)
})

test_that("a strategy that matches nothing returns the empty schema, not NULL", {
    # Every SNP is single-method, so Cross-method finds no pairs.
    lst <- three_methods()
    lst$LFMM <- lst$LFMM[0]
    lst$RDA  <- lst$RDA[0]
    r <- quiet_empty(combine_sigsnps(lst, "Cross-method", 0L, PREDICTORS))
    expect_identical(nrow(r), 0L)
    expect_true(data.table::is.data.table(r))
})

test_that("NULL entries in the list are tolerated", {
    lst <- c(three_methods()[c("EMMAX", "LFMM")], list(RDA = NULL))
    r <- quiet(combine_sigsnps(lst, "Union", 0L, PREDICTORS))
    expect_setequal(r$SNPID, c("A", "B", "C"))
})
