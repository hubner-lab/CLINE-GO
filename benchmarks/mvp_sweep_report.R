#!/usr/bin/env Rscript
# =============================================================================
# mvp_sweep_report.R -- aggregate the per-seed sweep into the cross-seed tables and
# figures the lab journal reports.
#
#   Rscript benchmarks/mvp_sweep_report.R [--seeds=...] [--outdir=...]
#
# Inputs   benchmarks/mvp_eval/sweep/MVP{seed}/{tag}_threshold_sweep.tsv
#          benchmarks/mvp_eval/sweep/MVP{seed}/rank.tsv
#          benchmarks/mvp_sweep_cells.tsv, mvp_seeds.tsv, mvp_published_baseline.tsv
# Outputs  benchmarks/mvp_eval/report/*.tsv, *.png, *.svg
#
# THE THREE CONFIGURATIONS, and why the distinction is the whole report
# ---------------------------------------------------------------------
#   default  cell c3 -- every method at sNMF k_best, i.e. what the pipeline does if the
#            user changes nothing.
#   pregea   each method at the rung mode=pregea recommended. No truth table consulted.
#            THIS is the achievable configuration and the one the claim rests on.
#   oracle   each method at its AUC-PR-best rung, chosen against the truth table. An UPPER
#            BOUND that no user can reach; it exists to bound what `pregea` leaves behind.
# Reporting `oracle` as the pipeline's performance would reduce the claim to "knowing the
# answer helps you find the answer". Every table here keeps the three separate.
#
# WHY at_least_2 AND rank_sum LEAD THE COMBINE COLUMNS
# ----------------------------------------------------
# A 4-way `intersection` that includes RDA is degenerate on this data by construction, not
# by accident: rda.R combines its two fits with pmax (rda.R:478), and the max of two roughly
# uniform p-values is a valid ORDERING but not a valid p-value -- median 0.706 vs the
# expected 0.5, so no multiple-testing rule fires (MVP_README.md:283-306). RDA therefore
# contributes ranks but almost never calls, and anything requiring unanimity inherits that.
# `rank_sum` is threshold-free and sidesteps the calibration problem entirely; `at_least_2`
# tolerates one abstaining method. Both are reported, and the intersection column is kept
# so the emptying is visible rather than hidden.
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

ROOT <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
source(file.path(ROOT, "scripts/R/utils/theme_clinego.R"))

A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    a <- sub("^--", "", a)
    A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}
SWEEP  <- file.path(ROOT, "benchmarks/mvp_eval/sweep")
OUTDIR <- if (is.null(A$outdir)) file.path(ROOT, "benchmarks/mvp_eval/report") else A$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

METHODS <- c("EMMAX", "LFMM", "RDA", "BLINK")
seeds_man <- fread(file.path(ROOT, "benchmarks/mvp_seeds.tsv"), colClasses = c(seed = "character"))
cells_man <- fread(file.path(ROOT, "benchmarks/mvp_sweep_cells.tsv"), colClasses = c(seed = "character"))
if (!is.null(A$seeds)) {
    keep <- strsplit(A$seeds, ",")[[1]]
    seeds_man <- seeds_man[seed %in% keep]
}

# tags are MVP{seed}_{c<i>|oracle|pregea}_{temp|sal|any}
TAG_RE <- "^MVP([0-9]+)_(c[0-9]+|oracle|pregea)_(temp|sal|any)$"

# ------------------------------------------------------------------ 1. threshold surface
files <- unlist(lapply(seeds_man$seed, function(s)
    list.files(file.path(SWEEP, paste0("MVP", s)), pattern = "_threshold_sweep\\.tsv$", full.names = TRUE)))
files <- files[grepl(TAG_RE, sub("_threshold_sweep\\.tsv$", "", basename(files)))]
if (!length(files)) stop("No *_threshold_sweep.tsv found under ", SWEEP)

sw <- rbindlist(lapply(files, function(f) {
    d <- fread(f)
    t <- sub("_threshold_sweep\\.tsv$", "", basename(f))
    d[, `:=`(seed = sub(TAG_RE, "\\1", t), arm = sub(TAG_RE, "\\2", t), axis = sub(TAG_RE, "\\3", t))]
    d
}), fill = TRUE)
sw[, config_group := fifelse(arm == "c3", "default", fifelse(arm %in% c("oracle", "pregea"), arm, "ladder"))]
sw <- merge(sw, seeds_man[, .(seed, architecture, arch_level, arm_class = arm, meanFst,
                              r2_pc1_temp, r2_pc1_sal, n_snps, n_causal_temp, n_causal_sal, k_best)],
            by = "seed", all.x = TRUE)
fwrite(sw, file.path(OUTDIR, "mvp_sweep_all.tsv"), sep = "\t")
message("INFO: mvp_sweep_all.tsv — ", nrow(sw), " rows, ", uniqueN(sw$seed), " seeds")

# ------------------------------------------------------------------ 2. ranking ladders
rk <- rbindlist(lapply(seeds_man$seed, function(s) {
    f <- file.path(SWEEP, paste0("MVP", s), "rank.tsv")
    if (!file.exists(f)) return(NULL)
    fread(f)
}), fill = TRUE)
if (nrow(rk)) {
    rk <- rk[metric %in% c("auc_pr", "best_f1", "n_causal_testable")]
    rk[, c("proj", "method", "cell", "ax_parsed") := tstrsplit(config, "|", fixed = TRUE)]
    rk[, `:=`(seed = sub("^MVP", "", proj), value = as.numeric(value))]
    ladder <- dcast(rk, seed + method + cell + ax_parsed ~ metric, value.var = "value")
    setnames(ladder, "ax_parsed", "axis")
    # attach the actual rung value (n_pcs / K), which no output path records
    ladder <- merge(ladder, cells_man[, .(seed, cell, method, param, rung = value)],
                    by = c("seed", "method", "cell"), all.x = TRUE)
    ladder <- merge(ladder, seeds_man[, .(seed, architecture, arch_level, meanFst,
                                          r2_pc1_temp, k_best)], by = "seed", all.x = TRUE)
    fwrite(ladder, file.path(OUTDIR, "mvp_ladders.tsv"), sep = "\t")
    message("INFO: mvp_ladders.tsv — ", nrow(ladder), " (seed, method, rung, axis) rows")

    p <- ggplot(ladder[!is.na(rung)], aes(x = rung, y = auc_pr, colour = method)) +
        geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
        facet_grid(axis ~ arch_level, scales = "free_x") +
        scale_color_clinego() +
        labs(title = "Structure-correction ladder: detection ranking vs rung",
             subtitle = "rung = EMMAX/RDA/BLINK PC count, or LFMM latent factors",
             x = "Rung (n_pcs / condition_pcs / LFMM K)", y = "AUC-PR", colour = "Method") +
        theme_clinego()
    ggsave(file.path(OUTDIR, "ladder_aucpr.png"), p, width = 11, height = 7, dpi = 150)
    ggsave(file.path(OUTDIR, "ladder_aucpr.svg"), p, width = 11, height = 7)

    # The two Est-Clines controls (R2(PC1~temp) ~ 0.92-0.94) exist to collapse on purpose
    # once structure correction switches on. If they do not, the benchmark is not measuring
    # what the seed selection claims (MVP_README.md:413).
    ctrl <- ladder[seed %in% c("1232568", "1231578") & axis == "temp" & !is.na(rung)]
    if (nrow(ctrl)) {
        chk <- ctrl[order(rung), .(auc_at_min = auc_pr[1], auc_at_max = auc_pr[.N],
                                   collapsed = auc_pr[.N] < auc_pr[1]), by = .(seed, method)]
        fwrite(chk, file.path(OUTDIR, "control_collapse_check.tsv"), sep = "\t")
        message("INFO: control collapse — ", sum(chk$collapsed), "/", nrow(chk),
                " (seed, method) pairs decline from lowest to highest rung")
    }
} else message("WARN: no rank.tsv found — ladders and control check skipped")

# ---------------------------------------------------- 3. multi-method advantage
# For each (seed, axis, configuration): the best single method versus the best combine
# rule, at two operating points -- max F1, and max recall while precision stays >= 0.9.
best_of <- function(d, kinds) {
    d <- d[kind %in% kinds & n_called > 0]
    if (!nrow(d)) return(data.table(f1 = NA_real_, f1_kind = NA_character_, f1_at = NA_character_,
                                    hp_recall = NA_real_, hp_kind = NA_character_, hp_at = NA_character_))
    b  <- d[!is.na(f1)][order(-f1)][1]
    hp <- d[!is.na(precision) & precision >= 0.9][order(-recall_testable)][1]
    data.table(f1        = if (nrow(b))  b$f1 else NA_real_,
               f1_kind   = if (nrow(b))  b$kind else NA_character_,
               f1_at     = if (nrow(b))  paste0(b$adjust, " ", b$value) else NA_character_,
               hp_recall = if (nrow(hp) && !is.na(hp$recall_testable)) hp$recall_testable else NA_real_,
               hp_kind   = if (nrow(hp) && !is.na(hp$recall_testable)) hp$kind else NA_character_,
               hp_at     = if (nrow(hp) && !is.na(hp$recall_testable)) paste0(hp$adjust, " ", hp$value) else NA_character_)
}
COMBINES <- c("combine_union", "combine_at_least_2", "combine_rank_sum", "combine_intersection")
adv <- sw[config_group %in% c("default", "pregea", "oracle"),
          {
              s <- best_of(.SD, METHODS)
              c <- best_of(.SD, COMBINES)
              cbind(data.table(single_f1 = s$f1, single_f1_kind = s$f1_kind, single_f1_at = s$f1_at,
                               single_hp_recall = s$hp_recall, single_hp_kind = s$hp_kind,
                               combine_f1 = c$f1, combine_f1_kind = c$f1_kind, combine_f1_at = c$f1_at,
                               combine_hp_recall = c$hp_recall, combine_hp_kind = c$hp_kind))
          }, by = .(seed, axis, config_group, architecture, arch_level, n_causal_temp, n_causal_sal)]
adv[, `:=`(f1_gain     = combine_f1 - single_f1,
           recall_gain = combine_hp_recall - single_hp_recall)]
setorder(adv, seed, axis, config_group)
fwrite(adv, file.path(OUTDIR, "multimethod_advantage.tsv"), sep = "\t")
message("INFO: multimethod_advantage.tsv — combine beats best single on F1 in ",
        adv[!is.na(f1_gain) & f1_gain > 0, .N], "/", adv[!is.na(f1_gain), .N], " (seed, axis, config) cells")

# per-combine-rule detail, so the intersection emptying is visible rather than inferred
rule_detail <- sw[config_group %in% c("default", "pregea", "oracle") & kind %in% c(METHODS, COMBINES),
                  .(best_f1 = suppressWarnings(max(f1, na.rm = TRUE)),
                    max_called = max(n_called),
                    n_operating_points_with_calls = sum(n_called > 0)),
                  by = .(seed, axis, config_group, kind)]
rule_detail[!is.finite(best_f1), best_f1 := NA_real_]
fwrite(rule_detail, file.path(OUTDIR, "combine_rule_detail.tsv"), sep = "\t")

pd <- adv[!is.na(f1_gain)]
if (nrow(pd)) {
    p2 <- ggplot(pd, aes(x = single_f1, y = combine_f1, colour = axis, shape = config_group)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = CLINEGO_THRESHOLD) +
        geom_point(size = 2.6, alpha = 0.9) +
        scale_color_clinego() +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(title = "Best cross-method rule vs best single method",
             subtitle = "Above the diagonal: combining methods wins at that seed / axis / configuration",
             x = "Best single-method F1", y = "Best combine-rule F1",
             colour = "Axis", shape = "Configuration") +
        theme_clinego()
    ggsave(file.path(OUTDIR, "multimethod_advantage.png"), p2, width = 8, height = 7, dpi = 150)
    ggsave(file.path(OUTDIR, "multimethod_advantage.svg"), p2, width = 8, height = 7)
}

# --------------------------------------------------------- 4. false-positive controls
# The 1-trait seeds carry zero causal loci on bio_2. Precision and recall are undefined
# there by construction; the number that means something is how many calls a rule makes
# when every one of them is false.
fp_seeds <- seeds_man[n_causal_sal == 0, seed]
if (length(fp_seeds)) {
    fp <- sw[seed %in% fp_seeds & axis == "sal",
             .(seed, arm, config_group, kind, adjust, value, n_called, fp_background,
               expected_linked, architecture)]
    setorder(fp, seed, kind, adjust, value)
    fwrite(fp, file.path(OUTDIR, "fp_control_sal.tsv"), sep = "\t")
    message("INFO: fp_control_sal.tsv — ", length(fp_seeds), " zero-causal-sal seeds, ",
            nrow(fp), " operating points, all calls false by construction")
}

# ------------------------------------------------------- 5. published-baseline comparison
pub_f <- file.path(ROOT, "benchmarks/mvp_published_baseline.tsv")
if (file.exists(pub_f) && nrow(rk)) {
    pub <- fread(pub_f, colClasses = c(seed = "character"))
    # Report BOTH arms against the published number. Lotterhos ran one fixed setting; our
    # default rung (c3) is measurably the wrong one for EMMAX and RDA on this landscape, so
    # citing c3 alone would understate the pipeline, while citing the oracle would overstate
    # it. The achievable comparison is the PreGEA-recommended rung -- no truth table used.
    pregea_cells <- rbindlist(lapply(seeds_man$seed, function(s) {
        f <- file.path(SWEEP, paste0("MVP", s), "pregea_cells.tsv")
        if (file.exists(f)) fread(f, colClasses = c(seed = "character")) else NULL
    }), fill = TRUE)
    ours <- ladder[cell == "c3" & method %in% c("LFMM", "RDA"),
                   .(seed, method, axis, our_aucpr = auc_pr)]
    if (nrow(pregea_cells)) {
        pg <- merge(ladder[, .(seed, method, cell, axis, pregea_aucpr = auc_pr)],
                    pregea_cells[, .(seed, method, cell)], by = c("seed", "method", "cell"))
        ours <- merge(ours, pg[, .(seed, method, axis, pregea_aucpr)],
                      by = c("seed", "method", "axis"), all.x = TRUE)
    } else ours[, pregea_aucpr := NA_real_]
    cmp <- merge(ours, pub[, .(seed, pub_lfmm_aucpr_temp_neut, pub_lfmm_aucpr_sal_neut,
                               pub_rda_aucpr_neut, pub_rda_aucpr_neut_corr,
                               lfmm_temp_aucpr_suspect)], by = "seed", all.x = TRUE)
    # Always the *_neutSNPs columns: their neutral definition matches lib_detection.R's
    # background_neutral-only false-positive rule (MVP_README.md:203).
    cmp[, published := fifelse(method == "LFMM" & axis == "temp", pub_lfmm_aucpr_temp_neut,
                        fifelse(method == "LFMM" & axis == "sal",  pub_lfmm_aucpr_sal_neut,
                        fifelse(method == "RDA",                   pub_rda_aucpr_neut_corr, NA_real_)))]
    cmp[, caveat := fifelse(method == "LFMM" & axis == "temp" & lfmm_temp_aucpr_suspect,
                            "published temp AUC-PR not reproducible from the deposit's own temp p-value column — column swap suspected (MVP_README.md:220); compare against the re-scored value instead",
                            "")]
    fwrite(cmp[, .(seed, method, axis, published,
                   ours_default_c3 = our_aucpr, ours_pregea = pregea_aucpr, caveat)],
           file.path(OUTDIR, "published_comparison.tsv"), sep = "\t")
    message("INFO: published_comparison.tsv — ", nrow(cmp),
            " rows (published vs default c3 vs PreGEA-recommended rung)")
}

message("INFO: report written to ", OUTDIR)
