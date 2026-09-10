#!/usr/bin/env Rscript
# pregea_transfer_guard.R — C6/C.2:338: re-estimate calibration ONCE on the
# FULL marker set at the hyperparameters preGEA recommends, before the user
# commits to the full GEA run. Only the TREND across the ladder transfers;
# lambda_GC is not scale-free (Yang et al. 2011), so the pruned-set absolute
# value must not be reused as-is on the full set.
#
# OPT-IN (PreGEA.TransferGuard.enabled) — the ONE preGEA rule that pulls the
# full-set imputation chain into this otherwise LD-pruned-only mode.
# Scope: LFMM + EMMAX only. RDA's full-vs-pruned question is answered
# separately by A6 ("fit the final scan on the full set", not "transfer
# pruned hyperparameters") — RDA is not folded in here.
#
# 'auto' hyperparameters are resolved INSIDE this script from the
# recommendations file (passed as an input, not a Snakemake param — Snakemake
# cannot read a file that does not exist yet at DAG-build time).

suppressPackageStartupMessages({
    library(data.table)
    library(LEA)
    library(ggplot2)
})

source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/emmax_core.R")

args <- commandArgs(trailingOnly = TRUE)
################################################################################
LFMM_PRUNED_IMP <- args[1]
LFMM_FULL_IMP   <- args[2]
VCFSNP_FULL     <- args[3]
CLIMATE         <- args[4]
PREDICTORS      <- args[5]
LFMM_K_ARG      <- args[6]                   # "auto" or int
RECOMMENDATIONS <- args[7]
EMMAX_VCF_PRUNED       <- args[8]            # always a real path now (PreGEA runs LFMM+EMMAX+RDA
                                              # together, no per-block enabled switch) — the "NULL"
                                              # skip branch below is a vestigial safety net
EMMAX_TPED_PREFIX_PRUNED <- args[9]
EMMAX_KINSHIP_PRUNED   <- args[10]
EMMAX_VCF_FULL         <- args[11]           # "NULL" when EMMAX not in GEA.configs
EMMAX_TPED_PREFIX_FULL <- args[12]
EMMAX_KINSHIP_FULL     <- args[13]
EMMAX_NPCS_ARG  <- args[14]                  # "auto" or int
PCA             <- args[15]
SAMPLES_ORDER   <- args[16]
METADATA        <- args[17]
SEED            <- as.numeric(args[18])
OUT_TSV         <- args[19]
PLOT_DIR        <- args[20]
INTER_DIR       <- args[21]
################################################################################

for (d in c(dirname(OUT_TSV), PLOT_DIR, INTER_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
predictor_names <- strsplit(PREDICTORS, ",")[[1]] |> trimws()

recs <- if (file.exists(RECOMMENDATIONS)) fread(RECOMMENDATIONS, sep = "\t", header = TRUE) else data.table()
# NOTE: the data.table `[` scoping resolves bare `method`/`param` to RECS'S OWN
# COLUMNS first — a function arg of the same name is invisible inside the `[`
# call, so `recs[method == method & param == param]` silently compares each
# column to itself (always TRUE) and just returns row 1 regardless of what was
# asked. Args below are named target_method/target_param specifically so no
# name collides with a column in `recs`.
resolve_auto <- function(arg, target_method, target_param, fallback) {
    if (arg != "auto") return(as.numeric(arg))
    r <- recs[method == target_method & param == target_param]
    if (nrow(r) > 0) return(as.numeric(r$recommended_value[1]))
    message("WARNING: no recommendation for ", target_method, ".", target_param, " — falling back to ", fallback)
    fallback
}

# lambda_GC — same formula pregea_ladder_stats.R uses (median chi-sq(1)
# statistic over the median of the null chi-sq(1)).
lambda_gc <- function(pvals) {
    pv <- pvals[!is.na(pvals)]
    if (length(pv) < 3) return(NA_real_)
    chisq1 <- suppressWarnings(qchisq(pv, df = 1, lower.tail = FALSE))
    stats::median(chisq1, na.rm = TRUE) / qchisq(0.5, df = 1)
}

rows <- list()
add_row <- function(method, param, value, scope, metric, val, delta, verdict)
    rows[[length(rows) + 1]] <<- data.table(method = method, param = param, value = value, scope = scope,
                                            metric = metric, value_measured = val, delta_vs_pruned = delta, verdict = verdict)

climate_dt <- fread(CLIMATE, sep = "\t", header = TRUE)

################################################################################
# LFMM — pruned vs full lambda at the recommended K
################################################################################
lfmm_k <- resolve_auto(LFMM_K_ARG, "LFMM", "K", 3)
message("INFO: Transfer guard — LFMM K=", lfmm_k)

lfmm_lambda_for <- function(geno_path) {
    Y <- fread(geno_path, sep = " ", header = FALSE) |> as.matrix()
    lambdas <- numeric(0)
    for (tr in predictor_names) {
        env <- climate_dt[[tr]]
        if (length(env) != nrow(Y)) next
        env <- scale(env)
        mod <- tryCatch(LEA::lfmm2(Y, env = env, K = lfmm_k), error = function(e) NULL)
        if (is.null(mod)) next
        res <- tryCatch(LEA::lfmm2.test(mod, input = Y, env = env, genomic.control = FALSE),
                        error = function(e) NULL)
        if (!is.null(res)) lambdas <- c(lambdas, lambda_gc(res$pvalues))
    }
    if (length(lambdas) == 0) NA_real_ else mean(lambdas, na.rm = TRUE)
}

lambda_pruned <- tryCatch(lfmm_lambda_for(LFMM_PRUNED_IMP), error = function(e) { message("WARNING: LFMM pruned fit failed: ", e$message); NA_real_ })
lambda_full   <- tryCatch(lfmm_lambda_for(LFMM_FULL_IMP),   error = function(e) { message("WARNING: LFMM full fit failed: ", e$message); NA_real_ })
delta_lfmm <- if (!is.na(lambda_pruned) && !is.na(lambda_full)) lambda_full - lambda_pruned else NA_real_
verdict_lfmm <- if (is.na(delta_lfmm)) "unavailable" else if (abs(delta_lfmm) < 0.3) "transfers" else "diverges"
add_row("LFMM", "K", lfmm_k, "pruned", "lambda_gc", lambda_pruned, NA_real_, "reference")
add_row("LFMM", "K", lfmm_k, "full", "lambda_gc", lambda_full, delta_lfmm, verdict_lfmm)
message("INFO: LFMM lambda pruned=", round(lambda_pruned, 3), " full=", round(lambda_full, 3), " verdict=", verdict_lfmm)

################################################################################
# EMMAX — pruned vs full lambda at the recommended #PCs. The pruned arm
# (Block 0's own LD-pruned TPED/kinship) is always available now — PreGEA
# runs LFMM+EMMAX+RDA together, no per-block enabled switch. The full arm
# needs 'EMMAX' configured in GEA.configs (gea.smk's full-set TPED/kinship
# chain), so only that one may still be unavailable.
################################################################################
emmax_npcs <- resolve_auto(EMMAX_NPCS_ARG, "EMMAX", "n_pcs", 3)
message("INFO: Transfer guard — EMMAX n_pcs=", emmax_npcs)

emmax_lambda_for <- function(vcf_path, tped_prefix, kinship_path, inter_dir) {
    dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
    lambdas <- numeric(0)
    for (tr in predictor_names) {
        out_file <- file.path(inter_dir, paste0(tr, "_pvalues.tsv"))
        cmd <- sprintf(
            "Rscript /pipeline/scripts/emmax.R %s %d %s %s %s %s %s %s %s %s %s",
            shQuote(vcf_path), as.integer(emmax_npcs), shQuote(CLIMATE), shQuote(PCA), tr,
            shQuote(inter_dir), shQuote(METADATA), shQuote(out_file), shQuote(tped_prefix),
            shQuote(kinship_path), shQuote(SAMPLES_ORDER))
        status <- system(cmd)
        if (status == 0 && file.exists(out_file)) {
            pv <- fread(out_file, sep = "\t", header = TRUE)
            if (tr %in% names(pv)) lambdas <- c(lambdas, lambda_gc(pv[[tr]]))
        }
    }
    if (length(lambdas) == 0) NA_real_ else mean(lambdas, na.rm = TRUE)
}
emmax_arm_available <- function(vcf) !identical(vcf, "NULL") && file.exists(vcf)

lambda_pruned_e <- NA_real_
if (!emmax_arm_available(EMMAX_VCF_PRUNED)) {
    # Vestigial safety net — the pruned arm is always available now
    # (PreGEA runs LFMM+EMMAX+RDA together), so this branch should never fire.
    add_row("EMMAX", "n_pcs", emmax_npcs, "pruned", "status", NA_real_, NA_real_, "unavailable")
    message("INFO: EMMAX pruned arm skipped — input path missing")
} else {
    lambda_pruned_e <- tryCatch(
        emmax_lambda_for(EMMAX_VCF_PRUNED, EMMAX_TPED_PREFIX_PRUNED, EMMAX_KINSHIP_PRUNED,
                         file.path(INTER_DIR, "transfer_emmax_pruned")),
        error = function(e) { message("WARNING: EMMAX pruned transfer-guard fit failed: ", e$message); NA_real_ })
    add_row("EMMAX", "n_pcs", emmax_npcs, "pruned", "lambda_gc", lambda_pruned_e, NA_real_, "reference")
}

if (!emmax_arm_available(EMMAX_VCF_FULL)) {
    add_row("EMMAX", "n_pcs", emmax_npcs, "full", "status", NA_real_, NA_real_, "unavailable_no_gea_emmax_config")
    message("INFO: EMMAX full arm skipped — no GEA EMMAX config available")
} else {
    lambda_full_e <- tryCatch(
        emmax_lambda_for(EMMAX_VCF_FULL, EMMAX_TPED_PREFIX_FULL, EMMAX_KINSHIP_FULL,
                         file.path(INTER_DIR, "transfer_emmax_full")),
        error = function(e) { message("WARNING: EMMAX full transfer-guard fit failed: ", e$message); NA_real_ })
    delta_emmax <- if (!is.na(lambda_pruned_e) && !is.na(lambda_full_e)) lambda_full_e - lambda_pruned_e else NA_real_
    verdict_emmax <- if (is.na(delta_emmax)) "unavailable" else if (abs(delta_emmax) < 0.3) "transfers" else "diverges"
    add_row("EMMAX", "n_pcs", emmax_npcs, "full", "lambda_gc", lambda_full_e, delta_emmax, verdict_emmax)
    message("INFO: EMMAX lambda pruned=", round(lambda_pruned_e, 3), " full=", round(lambda_full_e, 3), " verdict=", verdict_emmax)
}

out <- rbindlist(rows, fill = TRUE)
fwrite(out, OUT_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote transfer guard table (", nrow(out), " rows) to ", OUT_TSV)

################################################################################
# Plot: pruned vs full lambda, LFMM only (the arm this rule fully computes)
################################################################################
plot_dt <- out[method == "LFMM"]
g <- ggplot(plot_dt, aes(x = scope, y = value_measured)) +
    geom_hline(yintercept = 1, color = CLINEGO_THRESHOLD, linetype = "dashed") +
    geom_col(fill = CLINEGO_NEUTRAL) +
    labs(x = NULL, y = expression(lambda[GC]),
        title = "Transfer guard: LFMM lambda, pruned vs full marker set",
        subtitle = sprintf("K=%s | verdict=%s", lfmm_k, verdict_lfmm)) +
    theme_clinego()
ggsave(file.path(PLOT_DIR, "transfer_guard_lambda.png"), g, width = 6, height = 4.5, dpi = 300)
ggsave(file.path(PLOT_DIR, "transfer_guard_lambda.svg"), g, width = 6, height = 4.5,
      device = svglite::svglite, bg = "transparent")

message("INFO: preGEA transfer guard complete.")
