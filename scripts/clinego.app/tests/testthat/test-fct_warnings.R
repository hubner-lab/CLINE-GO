# Data-problem warnings (fct_warnings.R, release A2). Each rule is pure: these fixtures
# are SIMDATA-like (silent) or carry the problem (firing). The real trees were checked
# by hand, read-only: Trifolium d150 / K5ERpc5 EMMAX bio_15 has a 135-hit tie on six
# chromosomes; SIMDATA stays silent up to the top 100 SNPs per trait.

# The badge text of a rendered warning (the trigger span, not the tooltip body).
badge_text <- function(tag) {
    html <- as.character(tag)
    m <- regmatches(html, regexec(
        '<span class="badge note-badge[^"]*">(.*?)</span>', html))[[1]]
    if (length(m) < 2) return(NA_character_)
    trimws(gsub("\\s+", " ", gsub("<[^>]+>", "", m[2])))
}

sig_cell <- function(method, trait, snp, chr, p) {
    data.table::data.table(SNPID = snp, chr = chr, pos = seq_along(snp), pvalue = p,
                           method = method, trait = trait)
}

test_that("every warning id has its help note, and every label fits in 8 words", {
    for (id in WARNING_NOTE_IDS) {
        expect_false(is.null(HELP_NOTES[[id]]), info = id)
        expect_false(is.null(warning_note(id, "x")), info = id)
        # information, not teaching (A3): items, at most one short clause, the Config line
        html <- as.character(warning_note(id, "x", items = c("a: 1", "b: 2")))
        expect_match(html, "<li>a: 1</li>", fixed = TRUE, info = id)
        expect_false(grepl("What it (means|shows)", html), info = id)
        expect_lte(nchar(HELP_NOTES[[id]]$desc %||% ""), 100L)
    }
    expect_null(warning_note("warn_no_such_id", "x"))
    flags <- data.table::data.table(predictor = c("bio_8", "bio_11", "bio_17"),
                                    flag = c("identical", "identical", "binary"),
                                    detail = c("bio_11", "bio_8", "0 and 1"))
    stripe2 <- data.table::data.table(method = "EMMAX", traits = "bio_15", n_hits = 210L,
                                      n_stripes = 8L, n_stripe_hits = 203L, largest = 135L,
                                      largest_p = 2.5e-9, n_chr = 6L)
    notes <- list(
        stripe_note(stripe2),
        stripe_note(data.table::copy(stripe2)[, n_stripes := 1L]),
        ld_decay_note(list(status = "unreliable", reason = "r")),
        ld_decay_note(list(status = "not assessed", reason = "r")),
        climate_predictors_note(flags, c("bio_8", "bio_11")),
        climate_predictors_note(flags, "bio_17"),
        climate_predictors_note(rbind(flags, data.table::data.table(
            predictor = "bio_18", flag = "binary", detail = "1 and 2")),
            c("bio_8", "bio_11", "bio_17", "bio_18")),
        climate_dropped_note(list(no_coords = "GAL05",
                                  no_climate = data.table::data.table(sample = "X1",
                                                                      reason = "sea"))),
        chip_legend(flags, "bio_14", "bio_8", "leg", "warn_predictor_chips", "predictor"),
        cross_entropy_note(list(k_end = 7L, k_values = 2:7, ce = 7:2 / 10)),
        dbmem_selection_note(list(state = "fallback", mems = "MEM2")),
        dbmem_selection_note(list(state = "unknown")))
    for (n in notes) {
        expect_false(is.null(n))
        lab <- badge_text(n)
        expect_false(is.na(lab))
        expect_lte(label_words(lab), WARNING_LABEL_MAX_WORDS)
    }
})

test_that("stripes: SNPs of one scan tied on p across chromosomes", {
    p_rand <- c(1.2e-8, 3.4e-8, 5.6e-8, 7.8e-8, 9.1e-8)
    calm <- list(EMMAX = sig_cell("EMMAX", "bio_1", paste0("s", 1:5),
                                  c("1", "1", "2", "3", "4"), p_rand))
    expect_identical(nrow(manhattan_stripes(calm)), 0L)

    # an LD block: three SNPs with one p on one chromosome is no stripe
    ld <- list(EMMAX = sig_cell("EMMAX", "bio_1", paste0("s", 1:3), "6", rep(1.854419e-07, 3)))
    expect_identical(nrow(manhattan_stripes(ld)), 0L)

    # five SNPs sharing a p on three chromosomes (plus one ordinary hit)
    st <- manhattan_stripes(list(EMMAX = sig_cell("EMMAX", "bio_15", paste0("s", 1:6),
        c("1", "1", "2", "3", "3", "5"), c(rep(2.5e-9, 5), 4e-9))))
    expect_identical(nrow(st), 1L)
    expect_identical(c(st$n_hits, st$n_stripes, st$n_stripe_hits, st$largest, st$n_chr),
                     c(6L, 1L, 5L, 5L, 3L))
    expect_identical(stripe_label(st), "Stripe: 5 hits share one genotype")
    expect_match(badge_text(stripe_note(st)), "^Stripe: 5 hits")
    # one item per scan: the scan, its ties, the largest tie's p and hits; no teaching block
    html <- as.character(stripe_note(st))
    expect_match(html, "<li>EMMAX, bio_15: 5 of 6 hits in 1 tie; largest p = 2.5e-09, 5 hits on 3 chromosomes</li>",
                 fixed = TRUE)
    expect_false(grepl("What it means", html))
    expect_null(stripe_note(st, wza = TRUE))

    # ties read at 10 significant digits
    near <- sig_cell("LFMM", "bio_1", c("a", "b"), c("1", "2"), c(1.23456789012e-7, 1.23456789049e-7))
    expect_identical(nrow(manhattan_stripes(list(LFMM = near))), 1L)
    apart <- sig_cell("LFMM", "bio_1", c("a", "b"), c("1", "2"), c(1.234567891e-7, 1.234567892e-7))
    expect_identical(nrow(manhattan_stripes(list(LFMM = apart))), 0L)

    # p = 0 (underflow) and p = 1 ties are not read
    edge <- sig_cell("RDA", "climate_multivariate", paste0("s", 1:4), c("1", "2", "3", "4"),
                     c(0, 0, 1, 1))
    expect_identical(nrow(manhattan_stripes(list(RDA = edge))), 0L)
})

test_that("stripes: twin traits count once, and a card reads its own scan", {
    snp <- paste0("s", 1:4); chr <- c("1", "2", "3", "4"); p <- rep(3e-8, 4)
    sig <- list(EMMAX = rbind(sig_cell("EMMAX", "bio_16", snp, chr, p),
                              sig_cell("EMMAX", "bio_19", snp, chr, p),
                              sig_cell("EMMAX", "bio_1", snp[1:2], chr[1:2], c(1e-8, 2e-8))))
    st <- manhattan_stripes(sig)
    expect_identical(st$traits, "bio_16 = bio_19")
    expect_identical(nrow(stripes_for(st, "EMMAX", "bio_19")), 1L)
    expect_identical(nrow(stripes_for(st, "EMMAX", "bio_1")), 0L)
    expect_identical(nrow(stripes_for(st, "LFMM", "bio_16")), 0L)
    expect_null(stripes_for(NULL, "EMMAX", "bio_16"))
    expect_null(stripe_note(stripes_for(st, "EMMAX", "bio_1")))
})

# Per-sample climate: three sites of two samples each.
clim_fixture <- function(n_sites = 3) {
    site <- rep(paste0("S", seq_len(n_sites)), each = 2)
    v <- rep(c(10.5, 12.1, 14.9, 16.0)[seq_len(n_sites)], each = 2)
    list(
        clim = data.table::data.table(
            sample = paste0("x", seq_along(site)),
            bio_1 = v, bio_8 = v * 2 + 1, bio_11 = v * 2 + 1,   # bio_8 = bio_11; bio_1 linear in both
            bio_12 = rep(c(300, 410, 520, 600)[seq_len(n_sites)], each = 2) + c(0, 3),
            bio_14 = 0,                                              # invariant: not read here
            bio_17 = rep(c(0, 1, 1, 0)[seq_len(n_sites)], each = 2)),
        meta = data.table::data.table(sample = paste0("x", seq_along(site)), site = site))
}

test_that("identical and two-valued predictors, from three sites on", {
    f <- clim_fixture(3)
    fl <- predictor_flags(f$clim, f$meta)
    idn <- fl[fl$flag == "identical", ]
    expect_setequal(idn$predictor, c("bio_1", "bio_8", "bio_11"))
    expect_identical(idn$detail[idn$predictor == "bio_8"], "bio_1, bio_11")
    expect_identical(fl$predictor[fl$flag == "binary"], "bio_17")
    expect_identical(fl$detail[fl$flag == "binary"], "0 and 1")
    expect_false("bio_14" %in% fl$predictor)
    expect_false("bio_12" %in% fl$predictor)

    # two sites: every predictor is two-valued and perfectly correlated, so nothing is read
    f2 <- clim_fixture(2)
    expect_identical(nrow(predictor_flags(f2$clim, f2$meta)), 0L)
    expect_identical(nrow(predictor_flags(NULL, f$meta)), 0L)
})

test_that("predictor counters read the selected predictors only (e83e54)", {
    f <- clim_fixture(3)
    fl <- predictor_flags(f$clim, f$meta)
    # SIMDATA-like: the flagged layers are not selected
    expect_null(climate_predictors_note(fl, c("bio_12")))
    # one twin selected alone is no duplicate
    expect_null(climate_predictors_note(fl, c("bio_8", "bio_12")))
    expect_identical(selected_twin_sets(fl, c("bio_8", "bio_11", "bio_12")),
                     list(c("bio_8", "bio_11")))
    expect_identical(badge_text(climate_predictors_note(fl, c("bio_8", "bio_11"))),
                     "2 selected predictors: identical")
    expect_identical(badge_text(climate_predictors_note(fl, c("bio_1", "bio_8", "bio_11"))),
                     "3 selected predictors: identical")
    expect_identical(badge_text(climate_predictors_note(fl, c("bio_17", "bio_12"))),
                     "1 selected predictor: two-valued")
    # one badge for both problems (W4), one item per twin set and per two-valued predictor
    both <- climate_predictors_note(fl, c("bio_8", "bio_11", "bio_17"))
    expect_identical(badge_text(both), "3 selected predictors: identical or two-valued")
    expect_match(as.character(both), "<li>bio_8 = bio_11: identical at every sample</li>", fixed = TRUE)
    expect_match(as.character(both), "<li>bio_17: two values across the sites, 0 and 1</li>",
                 fixed = TRUE)
    expect_identical(selected_invariant(c("bio_14", "bio_15"), c("bio_1", "bio_15")), "bio_15")
    expect_identical(selected_invariant("bio_14", c("bio_1", "bio_2")), character(0))
})

test_that("the chips mark identical and two-valued predictors, selectable; invariant struck", {
    f <- clim_fixture(3)
    fl <- predictor_flags(f$clim, f$meta)
    html <- as.character(render_bio_chips("cfg-Climate_predictors", "bio_8,bio_12",
                                          invariant = "bio_14", flags = fl))
    chip <- function(b) regmatches(html, regexpr(
        paste0('<button[^>]*data-bio="', b, '"[^>]*>'), html))
    expect_match(chip("bio_8"), "bc-flagged")
    expect_match(chip("bio_8"), "bc-active")
    expect_match(chip("bio_8"), 'data-flag="identical"')
    expect_match(chip("bio_8"), "identical to bio_1, bio_11 at every sample")
    expect_match(chip("bio_17"), 'data-flag="binary"')
    expect_match(chip("bio_17"), "two values across the sites, 0 and 1")
    # no corner icons: a marked chip is orange, the legend names each problem (C1)
    expect_false(grepl("bc-flag-icons", html))
    expect_match(html, "<li>bio_8: identical to bio_1, bio_11 at every sample</li>", fixed = TRUE)
    expect_match(html, "<li>bio_14: no variance, cannot be selected</li>", fixed = TRUE)
    expect_false(grepl("bc-active", chip("bio_17")))
    expect_match(chip("bio_14"), "bc-invariant")
    expect_false(grepl("bc-flagged", chip("bio_14")))
    expect_false(grepl("bc-flagged", chip("bio_12")))
    expect_match(html, 'id="cfg-Climate_predictors_legend"')
    # a marked chip is selected (bio_8): the legend starts amber
    expect_match(html, 'class="badge note-badge bg-warning text-dark"')
    calm <- as.character(render_bio_chips("cfg-x", "bio_12", invariant = "bio_14", flags = fl))
    expect_false(grepl('class="badge note-badge bg-warning', calm))
    expect_match(calm, 'class="badge note-badge bg-secondary"')
    expect_match(calm, "marked predictors")
    # nothing marked: no legend
    none <- as.character(render_bio_chips("cfg-y", "bio_1", choices = c("bio_1", "bio_2")))
    expect_false(grepl('id="cfg-y_legend"', none))
})

test_that("identical traits need three distinct value pairs", {
    meta <- data.table::data.table(site = "A", sample = paste0("s", 1:6),
                                   latitude = 1, longitude = 1,
                                   height = c(1, 2, 3, 4, NA, 6),
                                   height_cm = c(10, 20, 30, 40, 50, 60),
                                   score = c(5, 3, 4, 1, 2, 2),
                                   site_a = c(1, 1, 1, 2, 2, 2), site_b = c(7, 7, 7, 9, 9, 9))
    fl <- trait_flags(meta, c("height", "height_cm", "score", "site_a", "site_b"))
    expect_setequal(fl$predictor, c("height", "height_cm"))
    expect_identical(fl$detail[fl$predictor == "height"], "height_cm")
    expect_identical(nrow(trait_flags(meta, "height")), 0L)
})

test_that("samples left out of the climate analyses", {
    none <- climate_dropped_samples(
        data.table::data.table(n_total = 57L, n_available = 57L, n_missing = 0L,
                               missing_samples = ""),
        data.table::data.table(sample = character(), site = character(), reason = character()))
    expect_null(none)
    expect_null(climate_dropped_note(none))
    d <- climate_dropped_samples(
        data.table::data.table(n_total = 47L, n_available = 45L, n_missing = 2L,
                               missing_samples = "GAL05,GAL06"),
        data.table::data.table(sample = "NEG03", site = "Negev",
                               reason = "falls on ocean/NoData pixel"))
    expect_identical(d$no_coords, c("GAL05", "GAL06"))
    expect_identical(badge_text(climate_dropped_note(d)), "3 samples left out of climate analyses")
    expect_match(as.character(climate_dropped_note(d)), "<li>no coordinates: GAL05, GAL06</li>",
                 fixed = TRUE)
    expect_match(as.character(climate_dropped_note(d)), "<li>falls on ocean/NoData pixel: NEG03</li>",
                 fixed = TRUE)
    expect_null(climate_dropped_samples(NULL, NULL))
})

test_that("cross-entropy without a minimum: lowest at the largest K, three K or more", {
    ce <- function(K, v) data.table::data.table(K = rep(K, each = 2), run = rep(1:2, length(K)),
                                                cross_entropy = rep(v, each = 2) + c(0, 0.01),
                                                best_run = rep(c(TRUE, FALSE), length(K)))
    # SIMDATA: lowest at K = 2, the first K tested
    expect_null(cross_entropy_no_minimum(ce(2:7, c(1.09, 1.15, 1.19, 1.25, 1.26, 1.34))))
    # a minimum inside the range
    expect_null(cross_entropy_no_minimum(ce(2:5, c(0.9, 0.8, 0.85, 0.88))))
    nm <- cross_entropy_no_minimum(ce(2:5, c(0.9, 0.8, 0.75, 0.7)))
    expect_identical(nm$k_end, 5L)
    expect_identical(badge_text(cross_entropy_note(nm)), "No minimum up to K = 5")
    # two K values, or one (the d150 tree's K5-5): nothing to read
    expect_null(cross_entropy_no_minimum(ce(2:3, c(0.9, 0.8))))
    expect_null(cross_entropy_no_minimum(ce(5L, 0.9)))
    expect_null(cross_entropy_no_minimum(NULL))
})

test_that("dbMEM: the fallback axis is amber, a table without the column unknown", {
    expect_null(dbmem_selection_state(data.table::data.table(
        mem = c("MEM1", "MEM2"), selected = c(TRUE, FALSE),
        selection = c("forward_selection", "not_selected"))))
    fb <- dbmem_selection_state(data.table::data.table(
        mem = c("MEM1", "MEM2", "MEM3"), selected = c(FALSE, TRUE, FALSE),
        selection = c("not_selected", "fallback_best_marginal_r2", "not_selected")))
    expect_identical(fb, list(state = "fallback", mems = "MEM2"))
    expect_match(as.character(dbmem_selection_note(fb)), "bg-warning")
    old <- dbmem_selection_state(data.table::data.table(mem = "MEM1", selected = TRUE))
    expect_identical(old$state, "unknown")
    expect_match(as.character(dbmem_selection_note(old)), 'class="badge note-badge bg-secondary"')
    expect_null(dbmem_selection_state(NULL))
})

test_that("the LD-decay KPI is green only for a reliable curve (8bbdee)", {
    expect_identical(ld_decay_kpi_theme(list(status = "reliable", reason = "")), "success")
    expect_identical(ld_decay_kpi_theme(list(status = "unreliable", reason = "x")), "warning")
    expect_identical(ld_decay_kpi_theme(list(status = "not assessed", reason = "x")), "warning")
    expect_identical(ld_decay_kpi_theme(NULL), "secondary")
    expect_null(ld_decay_note(list(status = "reliable", reason = "")))
    expect_null(ld_decay_note(NULL))
})
