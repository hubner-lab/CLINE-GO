# R/fct_data_loading.R — the four pure functions in an otherwise I/O-heavy file.
# Everything else there reads a path through load_cached(); these four take data
# and return data, so they are the part that can be pinned without a fixture tree.
#
# cluster_pop_summary() is the one doing real science: there is no discrete
# cluster column on disk, so population membership is an ARGMAX over the Q-matrix
# columns. Every number the Structure tab shows about population sizes comes out
# of that one max.col() call.
#
# The other three each have a documented return that is not what a caller would
# guess — noted at each block.

# ---------------------------------------------------------------- parse_gff_attributes

test_that("parse_gff_attributes splits one row into key columns", {
    got <- parse_gff_attributes("ID=gene1;Name=ABC;biotype=protein_coding")
    expect_identical(nrow(got), 1L)
    expect_identical(got$ID, "gene1")
    expect_identical(got$Name, "ABC")
    expect_identical(got$biotype, "protein_coding")
})

test_that("parse_gff_attributes fills missing keys with NA across rows", {
    # fill = TRUE in the rbindlist is what makes variable key sets work; without
    # it this errors rather than producing NA.
    got <- parse_gff_attributes(c("ID=g1;Name=A", "ID=g2"))
    expect_identical(got$ID, c("g1", "g2"))
    expect_identical(got$Name, c("A", NA_character_))
})

test_that("parse_gff_attributes tolerates a different key ORDER between rows", {
    got <- parse_gff_attributes(c("ID=g1;Name=A", "Name=B;ID=g2"))
    expect_identical(got$ID, c("g1", "g2"))
    expect_identical(got$Name, c("A", "B"))
})

test_that("parse_gff_attributes keeps everything after the FIRST = in a value", {
    # sub("^[^=]+=", "") is non-greedy by construction, so an = inside the value
    # survives. GFF notes and URLs rely on this.
    got <- parse_gff_attributes("ID=g1;Note=a=b=c")
    expect_identical(got$Note, "a=b=c")
})

test_that("parse_gff_attributes gives a fragment with no = the same key and value", {
    # sub("=.*", "") leaves the token untouched when there is no =, and so does
    # sub("^[^=]+=", ""), so the column is named after its own content.
    got <- parse_gff_attributes("ID=g1;orphan")
    expect_identical(got$orphan, "orphan")
})

test_that("parse_gff_attributes SPLITS a value containing a semicolon", {
    # No URL-decoding: the GFF spec escapes ; as %3B, and this parser does not
    # decode it. A literal ; inside a value is therefore split into a new pair,
    # which here produces a column named after the value fragment.
    got <- parse_gff_attributes("ID=g1;Note=first;second")
    expect_identical(got$Note, "first")
    expect_true("second" %in% names(got))
})

test_that("parse_gff_attributes ignores empty fragments from doubled separators", {
    got <- parse_gff_attributes("ID=g1;;Name=A")
    expect_identical(nrow(got), 1L)
    expect_identical(got$ID, "g1")
    expect_identical(got$Name, "A")
})

test_that("parse_gff_attributes returns a zero-column table for empty input", {
    # Every row parses to list(), so there are no columns to bind. A caller
    # indexing a column by name gets NULL, not an error.
    got <- parse_gff_attributes(c(NA_character_, ""))
    expect_identical(nrow(got), 0L)
    expect_identical(ncol(got), 0L)
})

test_that("parse_gff_attributes DROPS empty rows instead of keeping them aligned", {
    # Characterisation of a live defect, filed 2026-09-12. An NA or empty
    # attribute string parses to list(), and rbindlist drops empty lists rather
    # than emitting an all-NA row — so the result has FEWER rows than the input
    # and no longer corresponds to it positionally.
    #
    # load_gff_genes() then does cbind(dt[, .(gene_id, chr, start, end)],
    # attr_dt) (:550), a POSITIONAL bind. Measured consequence: the attributes
    # recycle, so every gene after the attribute-less row is annotated with
    # another gene's attributes — right row count, plausible values, and for an
    # exact multiple not even a warning.
    got <- parse_gff_attributes(c("ID=g1", NA_character_, "ID=g3"))
    expect_identical(nrow(got), 2L)
    expect_identical(got$ID, c("g1", "g3"))
})

# ---------------------------------------------------------------- cluster_pop_summary

mk_q <- function(...) data.table::data.table(...)

test_that("cluster_pop_summary assigns each sample by argmax over C1..Ck", {
    dt <- mk_q(C1 = c(0.8, 0.1, 0.2), C2 = c(0.1, 0.8, 0.3), C3 = c(0.1, 0.1, 0.5))
    got <- cluster_pop_summary(dt, 3L)
    expect_identical(got$n_samples, 3L)
    expect_identical(got$n_pops, 3L)
    expect_identical(got$min_n, 1L)
    expect_identical(got$max_n, 1L)
})

test_that("cluster_pop_summary counts only OCCUPIED clusters, so n_pops can be < k", {
    # No sample has C3 as its maximum. n_pops is length(table(assign)), which
    # counts observed levels — the empty cluster silently vanishes, and the
    # Structure tab reports 2 populations for a K=3 run.
    dt <- mk_q(C1 = c(0.8, 0.1), C2 = c(0.1, 0.8), C3 = c(0.1, 0.1))
    got <- cluster_pop_summary(dt, 3L)
    expect_identical(got$n_pops, 2L)
    expect_identical(got$n_samples, 2L)
})

test_that("cluster_pop_summary breaks ties toward the leftmost column", {
    dt <- mk_q(C1 = 0.5, C2 = 0.5)
    got <- cluster_pop_summary(dt, 2L)
    expect_identical(got$n_pops, 1L)
    expect_identical(got$max_n, 1L)
})

test_that("cluster_pop_summary reports uneven population sizes", {
    dt <- mk_q(C1 = c(0.9, 0.9, 0.9, 0.1), C2 = c(0.1, 0.1, 0.1, 0.9))
    got <- cluster_pop_summary(dt, 2L)
    expect_identical(got$n_pops, 2L)
    expect_identical(got$min_n, 1L)
    expect_identical(got$max_n, 3L)
})

test_that("cluster_pop_summary summarises over the columns PRESENT when k overshoots", {
    # intersect() with names(dt) means a missing C3 is dropped rather than an
    # error — the summary is over C1/C2 alone.
    dt <- mk_q(C1 = c(0.8, 0.2), C2 = c(0.2, 0.8))
    got <- cluster_pop_summary(dt, 3L)
    expect_identical(got$n_pops, 2L)
    expect_identical(got$n_samples, 2L)
})

test_that("cluster_pop_summary returns NULL for an empty table or no Q columns", {
    expect_null(cluster_pop_summary(mk_q(C1 = numeric(0)), 1L))
    expect_null(cluster_pop_summary(mk_q(sample = "a", site = "b"), 3L))
})

test_that("cluster_pop_summary ignores non-Q columns sitting beside the matrix", {
    dt <- mk_q(sample = c("a", "b"), site = c("s1", "s2"),
               C1 = c(0.9, 0.1), C2 = c(0.1, 0.9))
    got <- cluster_pop_summary(dt, 2L)
    expect_identical(got$n_pops, 2L)
    expect_identical(got$n_samples, 2L)
})

# ---------------------------------------------------------------- wza_collapse_stats

test_that("wza_collapse_stats derives the window size from the SNPID range, inclusive", {
    dt <- data.table::data.table(SNPID = c("1:1-10000", "1:10001-20000"),
                                 n_snps = c(5L, 7L))
    got <- wza_collapse_stats(list(EMMAX = dt))
    expect_identical(got$n_windows, 2L)
    expect_identical(got$n_snps, 12L)
    expect_identical(got$window_bp, 10000L)   # end - start + 1
})

test_that("wza_collapse_stats inspects only the FIRST method's table", {
    # Window positions are identical across methods (positional fixed tiles), so
    # a second method with a different window count is ignored by design.
    a <- data.table::data.table(SNPID = "1:1-1000", n_snps = 3L)
    b <- data.table::data.table(SNPID = c("1:1-99", "1:100-199"), n_snps = c(1L, 1L))
    got <- wza_collapse_stats(list(EMMAX = a, LFMM = b))
    expect_identical(got$n_windows, 1L)
    expect_identical(got$n_snps, 3L)
})

test_that("wza_collapse_stats returns NA_integer_ for n_snps when the column is absent", {
    dt <- data.table::data.table(SNPID = "1:1-1000")
    got <- wza_collapse_stats(list(EMMAX = dt))
    expect_identical(got$n_snps, NA_integer_)
    expect_identical(got$window_bp, 1000L)
})

test_that("wza_collapse_stats leaves window_bp NA when the SNPID is not a range", {
    dt <- data.table::data.table(SNPID = "1:12345", n_snps = 1L)
    got <- wza_collapse_stats(list(EMMAX = dt))
    expect_identical(got$window_bp, NA_integer_)
    expect_identical(got$n_windows, 1L)
})

test_that("wza_collapse_stats returns NULL for no methods, a NULL table or zero rows", {
    expect_null(wza_collapse_stats(list()))
    expect_null(wza_collapse_stats(list(EMMAX = NULL)))
    expect_null(wza_collapse_stats(list(EMMAX = data.table::data.table(SNPID = character(0)))))
})

test_that("wza_collapse_stats ignores NA n_snps rather than propagating them", {
    dt <- data.table::data.table(SNPID = c("1:1-100", "1:101-200"),
                                 n_snps = c(4L, NA_integer_))
    expect_identical(wza_collapse_stats(list(EMMAX = dt))$n_snps, 4L)
})

# ---------------------------------------------------------------- design_metric

test_that("design_metric converts the stored character value to a number", {
    design <- list(n_sites = list(value = "9", flag = "", note = ""))
    expect_identical(design_metric(design, "n_sites"), 9)
})

test_that("design_metric returns NA_real_ for a metric that is absent", {
    expect_identical(design_metric(list(), "n_sites"), NA_real_)
    expect_identical(design_metric(list(other = list(value = "1")), "n_sites"), NA_real_)
})

test_that("design_metric returns NA_real_ for a non-numeric value, warning suppressed", {
    # The adequacy table deliberately mixes integers, rounded reals and the
    # literal "inf", so a non-numeric cell must degrade quietly.
    design <- list(vif = list(value = "not a number", flag = "FAIL", note = ""))
    expect_warning(got <- design_metric(design, "vif"), NA)
    expect_identical(got, NA_real_)
})

test_that("design_metric parses the literal \"inf\" the table can contain", {
    design <- list(vif = list(value = "inf", flag = "FAIL", note = ""))
    expect_identical(design_metric(design, "vif"), Inf)
})

test_that("design_metric returns numeric(0), NOT NA, for an entry with no value field", {
    # entry is non-NULL so the guard passes, then entry$value is NULL and
    # as.numeric(NULL) is numeric(0). A caller writing `if (design_metric(...) >
    # 5)` gets "argument is of length zero" rather than the NA it guards for.
    # Characterisation of a real length-0 escape, not endorsement.
    design <- list(n_sites = list(flag = "", note = "no value column"))
    got <- design_metric(design, "n_sites")
    expect_identical(got, numeric(0))
    expect_length(got, 0L)
    expect_false(isTRUE(is.na(got)))
})
