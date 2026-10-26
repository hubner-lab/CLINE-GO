# R/utils_helpers.R — three small functions whose failure mode is silence.
#
# keep_merged_method_files() is the sharp one. The discovery layer globs
# methods/<M>/*_pvalues_K*.tsv, which matches BOTH the merged per-method table and
# the per-trait tables written beside it. It separates them by a path-token rule:
# the basename must start with the name of the `methods/` child directory. Get that
# wrong and the app reads a single trait's p-values as if they were the method's
# whole table — no error, a smaller SNP set, and nothing says which file was used.

mk_path <- function(...) paste(c(...), collapse = .Platform$file.sep)

test_that("keep_merged_method_files keeps the merged per-method table", {
    files <- c(
        mk_path("SIMDATA_results", "GEA", "tables", "methods", "EMMAX",
                "EMMAX_pvalues_K3.tsv"),
        mk_path("SIMDATA_results", "GEA", "tables", "methods", "LFMM",
                "LFMM_pvalues_K3.tsv"))
    expect_identical(keep_merged_method_files(files), files)
})

test_that("keep_merged_method_files drops the per-trait siblings", {
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    merged  <- mk_path(dir, "EMMAX_pvalues_K3.tsv")
    pertrait <- c(mk_path(dir, "bio_1_pvalues_K3.tsv"),
                  mk_path(dir, "height_pvalues_K3.tsv"))
    got <- keep_merged_method_files(c(pertrait[1], merged, pertrait[2]))
    expect_identical(got, merged)
})

test_that("keep_merged_method_files rejects a path with no methods/ segment", {
    # Without the segment there is no method name to compare the basename against,
    # so the only safe answer is "not a merged table".
    files <- mk_path("SIMDATA_results", "GEA", "tables", "EMMAX_pvalues_K3.tsv")
    expect_length(keep_merged_method_files(files), 0L)
})

test_that("keep_merged_method_files uses the WZA token in wza regime", {
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    snp <- mk_path(dir, "EMMAX_pvalues_K3.tsv")
    wza <- mk_path(dir, "EMMAX_wza_K3.tsv")

    expect_identical(keep_merged_method_files(c(snp, wza), regime = "snp"), snp)
    expect_identical(keep_merged_method_files(c(snp, wza), regime = "wza"), wza)
})

test_that("keep_merged_method_files handles an empty input", {
    expect_length(keep_merged_method_files(character(0)), 0L)
})

test_that("keep_merged_method_files is not fooled by a method-name PREFIX", {
    # A directory named EMMAX must not claim EMMAX2_pvalues_K3.tsv. startsWith on
    # "EMMAX" plus the token is what prevents it: the token has to follow the name
    # immediately.
    dir <- c("SIMDATA_results", "GEA", "tables", "methods", "EMMAX")
    expect_length(
        keep_merged_method_files(mk_path(dir, "EMMAX2_pvalues_K3.tsv")), 0L)
})

test_that("resolve_adjust finds the method's own rule", {
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"),
        list(method = "LFMM",  adjust = "qval", threshold = "0.1"))))
    expect_identical(resolve_adjust(cfg, "EMMAX"), "bonf_0.05")
    expect_identical(resolve_adjust(cfg, "LFMM"), "qval_0.1")
})

test_that("resolve_adjust returns NULL for a method that is not configured", {
    # NULL, not a default: the caller uses it to decide whether the method ran at
    # all, so inventing "bonf_0.05" here would fabricate an output path.
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))))
    expect_null(resolve_adjust(cfg, "RDA"))
    expect_null(resolve_adjust(list(), "EMMAX"))
})

test_that("resolve_adjust prefers the module's own configs block", {
    cfg <- list(
        GEA  = list(configs = list(list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))),
        GWAS = list(configs = list(list(method = "EMMAX", adjust = "top",  threshold = "500"))))
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GEA"),  "bonf_0.05")
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GWAS"), "top_500")
})

test_that("resolve_adjust falls back to GEA.configs when the module has none", {
    # The documented inheritance: GWAS inherits GEA's rules unless overridden.
    cfg <- list(GEA = list(configs = list(
        list(method = "EMMAX", adjust = "bonf", threshold = "0.05"))))
    expect_identical(resolve_adjust(cfg, "EMMAX", module = "GWAS"), "bonf_0.05")
})

test_that("tab_to_mode maps each tab to its Snakemake mode string", {
    expect_identical(tab_to_mode("processing"), "processing")
    expect_identical(tab_to_mode("prestructure"), "prestructure")
    expect_identical(tab_to_mode("climate"), "climate")
    expect_identical(tab_to_mode("traits"), "traits")
    expect_identical(tab_to_mode("pregea"), "pregea")
    expect_identical(tab_to_mode("structure"), "structure")
    expect_identical(tab_to_mode("gea"), "gea")
    expect_identical(tab_to_mode("gwas"), "gwas")
})

test_that("tab_to_mode returns NULL for a tab with no single mode", {
    # Haplotype is driven by two separate runner instances, and an unknown tab must
    # not silently become a mode string that Snakemake would reject at parse time.
    expect_null(tab_to_mode("haplotype"))
    expect_null(tab_to_mode("not_a_tab"))
})

test_that("file_ok requires a single existing non-empty path", {
    d <- withr::local_tempdir()
    full  <- file.path(d, "full.tsv");  writeLines("x", full)
    empty <- file.path(d, "empty.tsv"); file.create(empty)

    expect_true(file_ok(full))
    expect_false(file_ok(empty))
    expect_false(file_ok(file.path(d, "missing.tsv")))
    expect_false(file_ok(NULL))
    expect_false(file_ok(NA_character_))
    expect_false(file_ok(c(full, full)))
})

# --- default_clumping_distance(): the filter bar's starting distance --------------
#
# The project's own distance, read as the pipeline reads it: bp as given; an auto_*
# value or no key at all (the pipeline's default since 2026-10-26) is the genome-wide
# LD-decay distance under the project's clumping_ld_rule, via the pipeline's own
# resolve_clumping_distance(). 100 kb stands in only when that cannot be read.

ld_table_file <- function(rows = list(list(group = "All", critical = 3100, half = 1240, r2thr = 2600)),
                          env = parent.frame()) {
    dt <- data.table::rbindlist(lapply(rows, function(r) data.table::data.table(
        group = r$group, scope = "genome_wide", n_samples = 50L, n_pairs = 1e5,
        max_dist_bp = 5e5, short_r2 = 0.45, background_r2 = 0.06, background_mean_r2 = 0.02,
        decay_critical_bp = r$critical, decay_half_bp = r$half, decay_r2thr_bp = r$r2thr,
        r2_threshold = 0.2, decay_method = "loess_log")))
    path <- withr::local_tempfile(fileext = ".tsv", .local_envir = env)
    data.table::fwrite(dt, path, sep = "\t")
    path
}

test_that("default_clumping_distance uses a bp spec as given, table or not", {
    cfg <- list(GEA = list(snp_clumping_distance = 25000L))
    expect_identical(default_clumping_distance(cfg, MOD_GEA, "NULL"), 25000L)
    expect_identical(default_clumping_distance(cfg, MOD_GEA, ld_table_file()), 25000L)
})

test_that("an auto spec, or none, is the table's distance under the project's rule", {
    tbl <- ld_table_file()
    for (cfg in list(list(GEA = list(snp_clumping_distance = "auto_genome_wide")),
                     list(GEA = list(snp_clumping_distance = "auto_per_chromosome")),
                     list(GEA = list(configs = list())),
                     list())) {
        expect_identical(default_clumping_distance(cfg, MOD_GEA, tbl), 3100L)
    }
    half <- list(GEA = list(clumping_ld_rule = "half_decay"))
    expect_identical(default_clumping_distance(half, MOD_GEA, tbl), 1240L)
})

test_that("GWAS inherits GEA's distance and rule, and may set its own", {
    tbl <- ld_table_file()
    expect_identical(default_clumping_distance(list(GEA = list(snp_clumping_distance = 7000L)),
                                               MOD_GWAS, tbl), 7000L)
    expect_identical(default_clumping_distance(list(GEA = list(clumping_ld_rule = "half_decay")),
                                               MOD_GWAS, tbl), 1240L)
    own <- list(GEA = list(snp_clumping_distance = 7000L),
                GWAS = list(snp_clumping_distance = "auto_genome_wide",
                            clumping_ld_rule = "r2_threshold"))
    expect_identical(default_clumping_distance(own, MOD_GWAS, tbl), 2600L)
    expect_identical(default_clumping_distance(own, MOD_GEA, tbl), 7000L)
})

test_that("the legacy region_distance / combine_gap names count, region_distance first", {
    tbl <- ld_table_file()
    expect_identical(default_clumping_distance(
        list(GEA = list(combine_gap = 100000L, region_distance = 1000000L)), MOD_GEA, tbl), 1000000L)
    expect_identical(default_clumping_distance(
        list(GEA = list(combine_gap = 50000L)), MOD_GEA, tbl), 50000L)
    expect_identical(default_clumping_distance(
        list(GEA = list(snp_clumping_distance = 9000L, region_distance = 1000000L)), MOD_GEA, tbl), 9000L)
})

test_that("the ld_decay_group is read, its legacy name too", {
    tbl <- ld_table_file(list(list(group = "All", critical = 3100, half = 1240, r2thr = 2600),
                              list(group = "Tel Aviv", critical = 4200, half = 900, r2thr = 3000)))
    expect_identical(default_clumping_distance(list(GEA = list(ld_decay_group = "Tel Aviv")),
                                               MOD_GEA, tbl), 4200L)
    expect_identical(default_clumping_distance(list(GEA = list(region_ld_decay_group = "Tel Aviv")),
                                               MOD_GWAS, tbl), 4200L)
})

test_that("100 kb stands in when the LD distance cannot be read", {
    old <- withr::local_tempfile(fileext = ".tsv")
    data.table::fwrite(data.table::data.table(group = "All", scope = "genome_wide",
                                              half_decay_bp = 1500, r2_02_bp = 130811), old, sep = "\t")
    unreached <- ld_table_file(list(list(group = "All", critical = NA, half = 1240, r2thr = 2600)))
    cases <- list(missing   = list(list(), "NULL"),
                  absent    = list(list(), file.path(tempdir(), "no_such_ld_table.tsv")),
                  old_table = list(list(), old),          # predates the PLINK estimator
                  unreached = list(list(), unreached),    # critical r2 never reached
                  bad_group = list(list(GEA = list(ld_decay_group = "nowhere")), ld_table_file()),
                  bad_spec  = list(list(GEA = list(snp_clumping_distance = "wide")), ld_table_file()))
    for (nm in names(cases)) {
        expect_identical(default_clumping_distance(cases[[nm]][[1]], MOD_GEA, cases[[nm]][[2]]),
                         DEFAULT_CLUMPING_DISTANCE, info = nm)
    }
})

# --- clumping_distance_state(): one value per project, saved only when moved -------
#
# The filter bar used to render from its own input. Once each project started at its
# own distance, a project switch made the bar re-render in a loop and saved one
# project's value into another's region_params.json (cold app, 2026-10-26). The state
# helper keeps one value per project and saves only a move of the box.

distance_server <- function(id, project_data) {
    shiny::moduleServer(id, function(input, output, session) {
        dist <- clumping_distance_state(input, project_data, "dist", MOD_GEA,
                                        "snp_clumping_distance", MOD_GEA)
    })
}

saved_distance <- function(project) {
    get_global_param(read_region_params(project), MOD_GEA, "snp_clumping_distance")
}

test_that("the box starts at each project's distance and saves only a move, under its project", {
    a <- paste0("cds_a_", as.integer(Sys.time()))
    b <- paste0("cds_b_", as.integer(Sys.time()))
    withr::defer(unlink(file.path(CLINEGO_TEST_ROOT, paste0(c(a, b), "_results")), recursive = TRUE))
    proj <- shiny::reactiveVal(list(name = a, config = list(GEA = list(snp_clumping_distance = 5000L))))
    shiny::testServer(distance_server, args = list(project_data = proj), {
        expect_identical(dist(), 5000L)
        session$setInputs(dist = 5000)          # the bar binding with the app's own value
        expect_null(saved_distance(a))          # opening a tab saves nothing
        session$setInputs(dist = 12345)         # a move
        expect_identical(dist(), 12345L)
        expect_identical(as.integer(saved_distance(a)), 12345L)

        proj(list(name = b, config = list(GEA = list(snp_clumping_distance = 7000L))))
        session$flushReact()
        expect_identical(dist(), 7000L)         # the new project's own start, not a's move
        session$setInputs(dist = 7000)          # the re-rendered bar binding
        expect_null(saved_distance(b))
        expect_identical(as.integer(saved_distance(a)), 12345L)

        proj(list(name = a, config = list(GEA = list(snp_clumping_distance = 5000L))))
        session$flushReact()
        expect_identical(dist(), 12345L)        # the move survives a switch away and back
    })
    # a new session reads the saved move
    proj2 <- shiny::reactiveVal(list(name = a, config = list(GEA = list(snp_clumping_distance = 5000L))))
    shiny::testServer(distance_server, args = list(project_data = proj2), {
        expect_identical(dist(), 12345L)
    })
})

test_that("an empty or zero box is not a move", {
    a <- paste0("cds_c_", as.integer(Sys.time()))
    withr::defer(unlink(file.path(CLINEGO_TEST_ROOT, paste0(a, "_results")), recursive = TRUE))
    proj <- shiny::reactiveVal(list(name = a, config = list(GEA = list(snp_clumping_distance = 5000L))))
    shiny::testServer(distance_server, args = list(project_data = proj), {
        session$setInputs(dist = NA)
        session$setInputs(dist = 0)
        expect_identical(dist(), 5000L)
        expect_null(saved_distance(a))
    })
})

test_that("a late report of the previous project's box is an echo, not a move", {
    # the cold-app race: the box bound for project a reports 12345 only after the
    # switch to b, which used to save a's value into b's file and loop the bar
    a <- paste0("cds_d_", as.integer(Sys.time()))
    b <- paste0("cds_e_", as.integer(Sys.time()))
    withr::defer(unlink(file.path(CLINEGO_TEST_ROOT, paste0(c(a, b), "_results")), recursive = TRUE))
    proj <- shiny::reactiveVal(list(name = a, config = list(GEA = list(snp_clumping_distance = 12345L))))
    shiny::testServer(distance_server, args = list(project_data = proj), {
        proj(list(name = b, config = list(GEA = list(snp_clumping_distance = 7000L))))
        session$flushReact()
        session$setInputs(dist = 12345)        # a's box, bound late
        expect_identical(dist(), 7000L)
        expect_null(saved_distance(b))
        session$setInputs(dist = 7000)         # b's box binds
        session$setInputs(dist = 8000)         # now a move
        expect_identical(dist(), 8000L)
        expect_identical(as.integer(saved_distance(b)), 8000L)
        expect_null(saved_distance(a))
    })
})

test_that("a box that already holds the new project's value takes a move at once", {
    a <- paste0("cds_f_", as.integer(Sys.time()))
    b <- paste0("cds_g_", as.integer(Sys.time()))
    withr::defer(unlink(file.path(CLINEGO_TEST_ROOT, paste0(c(a, b), "_results")), recursive = TRUE))
    proj <- shiny::reactiveVal(list(name = a, config = list(GEA = list(snp_clumping_distance = 9000L))))
    shiny::testServer(distance_server, args = list(project_data = proj), {
        session$setInputs(dist = 9000)
        proj(list(name = b, config = list(GEA = list(snp_clumping_distance = 9000L))))
        session$flushReact()                   # same start: the re-rendered box reports nothing
        session$setInputs(dist = 4000)
        expect_identical(as.integer(saved_distance(b)), 4000L)
        expect_null(saved_distance(a))
    })
})
