#!/usr/bin/env Rscript
# pregea_varpart.R — shared climate/structure/geography variance-partitioning
# panel (preGEA Block 4b). Ports rda_varpart_lasky.Rmd's forward-selection,
# Px, varpart, and anova.cca blocks (see that file for the reference
# implementation this was verified against).
#
# CORRECTIONS applied on port (docs/rda_research.md C.0/C5):
#   - Px: the Rmd's CODE already divides by tot.chi (sigma^2, total SNP
#     variance) correctly — only a COMMENT in the Rmd was wrong ("/ r2").
#     Fixed here from the start; see compute_Px() below.
#   - ordiR2step gains Pin/R2permutations/R2scope explicitly (the Rmd omitted
#     the first two) — Blanchet's double stopping rule (A17/C4).
#   - A17: the FULL forward-selection path is written (cumulative adjusted R2
#     per step + the full-model ceiling), never just the final variable set.
#
# RESPONSE MATRIX defaults to genomic PCs (LEA .projections truncated at a
# cumulative-variance target) — a documented SCALABILITY ADAPTATION, not
# literal Lasky et al. 2012 methodology (which used raw SNPs). The
# partitioning FRAMEWORK (varpart/dbMEM/Px) is what is being ported, not the
# exact response matrix — see docs/rda_research.md C.0.
#
# STRUCTURE TABLE is the sNMF Q-matrix (C1..CK, LAST column dropped for
# compositional singularity) — explicitly NOT the LEA PCs used as the
# response, which would make X_struct = Y tautological.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)      # load_pca_covariates() (emmax_core.R) uses %>% internally
    library(vegan)
    library(ggplot2)
    library(ggrepel)    # dbMEM selection-path labels only (see below)
    library(VennDiagram)  # schematic (non-proportional) 2/3-set Venn for the
                          # variance-partition plot — bioinformatics-standard
                          # tool for exactly this diagram, computes its own
                          # region label positions (draw.triple.venn()'s
                          # region layout is fixed/schematic by default —
                          # "General scaling for three-set Venn diagrams are
                          # disabled due to potentially misleading visual
                          # representation", its own docs). Two hand-rolled
                          # attempts before this (donutsk nested donut, then
                          # ggforce::geom_arc_bar + ggrepel schematic circles)
                          # both needed manually-guessed anchor points; this
                          # doesn't.
})

source("/pipeline/scripts/R/utils/theme_clinego.R")
source("/pipeline/scripts/R/utils/emmax_core.R")   # load_pca_covariates()

args <- commandArgs(trailingOnly = TRUE)
################################################################################
PROJECTIONS      <- args[1]
EIGENVALUES      <- args[2]
DBMEM_VECTORS    <- args[3]
CLUSTERS         <- args[4]                  # "NULL" if structure_table=none
CLIMATE          <- args[5]
CLIMATE_VALID    <- args[6]
SAMPLES_ORDER    <- args[7]
PREDICTORS       <- args[8]
RESPONSE         <- args[9]                  # "pcs" | "snps"
RESPONSE_VAR     <- as.numeric(args[10])
RESPONSE_MAX_PCS <- as.integer(args[11])
RESPONSE_MIN_PCS <- as.integer(args[12])
LFMM_PRUNED      <- args[13]                 # "NULL" when RESPONSE="pcs"
ORDIR2STEP_PIN   <- as.numeric(args[14])
R2_PERMUTATIONS  <- as.numeric(args[15])     # ordiR2step's internal R2permutations (fixed constant, Climate.Varpart layer)
VARPART_PERMS    <- as.numeric(args[16])     # the one surviving Varpart permutations knob (varpart anova test)
SEED             <- as.numeric(args[17])
OUT_SELECTION_TSV<- args[18]
OUT_SELECTED_TSV <- args[19]
OUT_VARPART_TSV  <- args[20]                 # human-readable variance_partition.tsv (was OUT_FRACTIONS_TSV/OUT_ANOVA_TSV, merged)
OUT_CONFOUND_TSV <- args[21]                 # climate_confounding.tsv (was OUT_ANOVA_TSV slot)
OUT_PX_TSV       <- args[22]
PLOT_DIR         <- args[23]
INTER_DIR        <- args[24]
################################################################################
# NOTE: SEL_PERMUTATIONS (old arg 16, "ordir2step_permutations" config key) was
# dropped here — it was assigned but never referenced (ordiR2step's own
# R2permutations control comes from R2_PERMUTATIONS above, arg 15). Removed
# together with the config key (module split cleanup) rather than left dead.
# CONFOUND_FLAG (old arg 17) was dropped too — the confounding-flag diagnostic
# is now always computed, not switchable (Climate config simplification).

for (d in c(dirname(OUT_SELECTION_TSV), PLOT_DIR, INTER_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

predictor_names <- strsplit(PREDICTORS, ",")[[1]] |> trimws()
climate_valid_ids <- fread(CLIMATE_VALID, header = FALSE, colClasses = "character")$V2

################################################################################
# 1. Response matrix Y — genomic PCs (default) or raw pruned SNPs
################################################################################
if (RESPONSE == "snps") {
    Y <- fread(LFMM_PRUNED, sep = " ", header = FALSE) |> as.matrix()
    response_desc <- sprintf("raw pruned SNPs (m=%d)", ncol(Y))
} else {
    eig <- as.numeric(readLines(EIGENVALUES))
    cum_var <- cumsum(eig / sum(eig))
    k_target <- which(cum_var >= RESPONSE_VAR)[1]
    if (is.na(k_target)) k_target <- length(eig)
    k_use <- max(RESPONSE_MIN_PCS, min(k_target, RESPONSE_MAX_PCS, length(eig)))
    Y <- load_pca_covariates(PROJECTIONS, k_use, SAMPLES_ORDER, climate_valid_ids)
    if (is.null(Y)) stop("preGEA varpart: could not build genomic-PC response matrix (k_use=", k_use, ")")
    response_desc <- sprintf("genomic PCs (k=%d, %.1f%% cumulative variance target)", k_use, RESPONSE_VAR * 100)
}
message("INFO: Response matrix Y: ", nrow(Y), " x ", ncol(Y), " [", response_desc, "]")
n_samples <- nrow(Y)

################################################################################
# 2. Climate predictors X_clim (climate-valid, same row order as Y)
################################################################################
climate_dt <- fread(CLIMATE, sep = "\t", header = TRUE)
X_clim <- as.data.frame(climate_dt[, ..predictor_names])
if (nrow(X_clim) != n_samples) stop("preGEA varpart: climate/response row mismatch (", nrow(X_clim), " vs ", n_samples, ")")
rownames(X_clim) <- climate_dt$sample %||% seq_len(n_samples)

################################################################################
# 3. Structure X_struct — sNMF Q-matrix, last column dropped (compositional)
################################################################################
X_struct <- NULL
if (!identical(CLUSTERS, "NULL") && file.exists(CLUSTERS)) {
    clusters_dt <- fread(CLUSTERS, colClasses = c(sample = "character"))
    q_cols <- grep("^C\\d+$", names(clusters_dt), value = TRUE)
    if (length(q_cols) >= 2) {
        q_cols_use <- q_cols[-length(q_cols)]   # drop last column (compositional singularity)
        clusters_dt <- clusters_dt[match(climate_valid_ids, sample)]
        X_struct <- as.data.frame(clusters_dt[, ..q_cols_use])
        message("INFO: Structure table: ", length(q_cols_use), "/", length(q_cols), " Q-matrix columns (last dropped)")
    }
}

################################################################################
# 4. Geography X_geo — dbMEM vectors, forward-selected (A17: full path kept)
################################################################################
dbmem_dt <- fread(DBMEM_VECTORS, colClasses = c(sample = "character"))
mem_cols <- grep("^MEM\\d+$", names(dbmem_dt), value = TRUE)
dbmem_status <- if (length(mem_cols) == 0) "no_spatial_vectors" else "ok"

X_geo_sel <- NULL
selection_dt <- data.table(step = integer(), variable = character(),
                           r2_adj_cumulative = numeric(), r2_adj_full_ceiling = numeric())
selected_dt <- data.table(mem = character(), selected = logical())

if (dbmem_status == "ok") {
    dbmem_dt <- dbmem_dt[match(climate_valid_ids, sample)]
    X_geo_full <- as.data.frame(dbmem_dt[, ..mem_cols])
    full_fit <- tryCatch(rda(as.formula(paste("Y ~", paste(mem_cols, collapse = " + "))), data = X_geo_full),
                         error = function(e) NULL)
    full_ceiling <- if (!is.null(full_fit)) suppressWarnings(RsquareAdj(full_fit)$adj.r.squared) else NA_real_

    selected_names <- character(0)
    if (!is.null(full_fit) && !is.na(full_ceiling) && full_ceiling > 0) {
        null_fit <- rda(Y ~ 1, data = X_geo_full)
        set.seed(SEED)
        step_res <- tryCatch(
            ordiR2step(null_fit, scope = formula(full_fit), Pin = ORDIR2STEP_PIN,
                      R2scope = TRUE, R2permutations = R2_PERMUTATIONS, trace = FALSE),
            error = function(e) { message("WARNING: ordiR2step (dbMEM) failed: ", e$message); NULL })
        if (!is.null(step_res)) selected_names <- attr(terms(step_res), "term.labels")
    }
    if (length(selected_names) == 0) {
        message("WARNING: forward selection chose 0 MEMs; falling back to all ", length(mem_cols), " MEM(s).")
        selected_names <- mem_cols
    }

    cum_terms <- character(0)
    for (i in seq_along(selected_names)) {
        cum_terms <- c(cum_terms, selected_names[i])
        step_fit <- tryCatch(rda(as.formula(paste("Y ~", paste(cum_terms, collapse = " + "))),
                                 data = X_geo_full), error = function(e) NULL)
        if (is.null(step_fit)) next
        selection_dt <- rbind(selection_dt, data.table(
            step = i, variable = selected_names[i],
            r2_adj_cumulative = suppressWarnings(RsquareAdj(step_fit)$adj.r.squared),
            r2_adj_full_ceiling = full_ceiling))
    }
    selected_dt <- data.table(mem = mem_cols, selected = mem_cols %in% selected_names)
    X_geo_sel <- X_geo_full[, selected_names, drop = FALSE]
    message("INFO: dbMEM forward selection: ", length(selected_names), "/", length(mem_cols), " MEM(s) selected")
} else {
    message("INFO: dbMEM status=", dbmem_status, " — geography fraction unavailable, running climate-vs-structure only")
}
fwrite(selection_dt, OUT_SELECTION_TSV, sep = "\t", quote = FALSE)
fwrite(selected_dt,  OUT_SELECTED_TSV,  sep = "\t", quote = FALSE)

################################################################################
# 5. Lasky Px per climate variable — compute_Px(), verbatim port (denominator
#    is tot.chi / (tot.chi - pCCA$tot.chi for partial models), NOT r2 — the
#    Rmd's comment was wrong, its code was always correct).
################################################################################
compute_Px <- function(rda_obj, X) {
    eig <- rda_obj$CCA$eig
    if (is.null(eig) || length(eig) == 0) return(NULL)
    denom <- if (!is.null(rda_obj$pCCA)) rda_obj$tot.chi - rda_obj$pCCA$tot.chi else rda_obj$tot.chi
    if (is.null(denom) || denom <= 0) return(NULL)
    lc <- scores(rda_obj, display = "lc", choices = seq_along(eig), scaling = 0)
    if (is.vector(lc)) lc <- matrix(lc, ncol = 1)
    Xs <- scale(as.matrix(X))
    R  <- suppressWarnings(cor(Xs, lc))
    Px <- as.numeric(abs(R) %*% eig) / denom
    data.table(variable = colnames(Xs), Px = Px, Px_pct = 100 * Px, rank = rank(-Px, ties.method = "min"))
}

rda_B <- tryCatch(rda(as.formula(paste("Y ~", paste(predictor_names, collapse = " + "))), data = X_clim),
                  error = function(e) NULL)
px_rows <- list()
if (!is.null(rda_B)) {
    px_b <- compute_Px(rda_B, X_clim)
    if (!is.null(px_b)) px_rows[[length(px_rows) + 1]] <- px_b[, model := "unconditional"]
}
if (!is.null(X_geo_sel) && ncol(X_geo_sel) > 0) {
    combined_df <- cbind(X_clim, X_geo_sel)
    mem_names <- paste(colnames(X_geo_sel), collapse = " + ")
    rda_C <- tryCatch(
        rda(as.formula(paste("Y ~", paste(predictor_names, collapse = " + "),
                             "+ Condition(", mem_names, ")")), data = combined_df),
        error = function(e) { message("WARNING: partial RDA (Px model C) failed: ", e$message); NULL })
    if (!is.null(rda_C)) {
        px_c <- compute_Px(rda_C, X_clim)
        if (!is.null(px_c)) px_rows[[length(px_rows) + 1]] <- px_c[, model := "partial_geo"]
    }
}
px_dt <- if (length(px_rows) > 0) rbindlist(px_rows) else data.table(
    variable = character(), Px = numeric(), Px_pct = numeric(), rank = integer(), model = character())
fwrite(px_dt, OUT_PX_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote Px table (", nrow(px_dt), " rows)")

################################################################################
# 6. Variance partitioning — ONE unified tree (Unique-to-one-factor / Shared-
#    between-factors / Unexplained, each with its own leaves) per run,
#    populated from whichever varpart fit is available:
#      - 3-table (climate+structure+geography) when both auxiliary tables exist
#      - 2-table climate-vs-geography when structure_table="none"
#      - 2-table climate-vs-structure fallback when geography is unavailable
#      - climate alone, last resort, when neither auxiliary table exists
#    The climate-vs-geography CONFOUNDING CHECK is a separate, dedicated
#    2-table fit, ALWAYS computed when geography is available regardless of
#    which tree above ends up displayed (docs/rda_research.md A.2 point 4) —
#    written to its own file (climate_confounding.tsv), not mixed into the
#    variance-partition rows (that mixing was the old fraction_type=
#    "diagnostic" hack).
################################################################################
# vegan::varpart()'s indfract table names its adj-R2 column "Adj.R.squared" for a
# 2-table Venn but "Adj.R.square" (no trailing "d") for 3+ tables — confirmed on
# this Docker's vegan 2.6-8. data.frame[idx, "wrong_name"] returns NULL silently
# (no error), so hardcoding either spelling makes the OTHER table count go to NA
# with zero warning. Resolve the real column name at each call site instead.
adjr2_col <- function(df) {
    hit <- grep("^Adj\\.R\\.square", colnames(df), value = TRUE)
    if (!length(hit)) stop("varpart indfract table has no Adj.R.square(d) column — vegan API changed?")
    hit[1]
}

# p-value for one testable (unique) fraction via anova.cca on its partial-RDA
# object. Shared/unexplained fractions are never directly testable this way —
# callers simply don't call this for those, leaving p_value NA.
compute_p <- function(rda_obj) {
    if (is.null(rda_obj)) return(NA_real_)
    set.seed(SEED)
    av <- tryCatch(anova.cca(rda_obj, permutations = VARPART_PERMS), error = function(e) NULL)
    if (is.null(av)) return(NA_real_)
    as.data.frame(av)[["Pr(>F)"]][1]
}

# Component/group naming — SEMANTIC keys, independent of which vegan lettering
# produced the number (2-table [a]/[b]/[c] and 3-table [a]-[g] mean different
# things by letter — conflating them was the source of the bug below). Single
# source of truth for both variance_partition.tsv and the Venn plot labels
# (CLAUDE.md Rule 4). No "(unique)"/"(shared)" suffixes — the Venn's own
# geometry (own-circle-only vs. overlap region) already makes that obvious.
COMPONENT <- c(
    climate_u   = "Climate",
    structure_u = "Structure",
    geo_u       = "Geography",
    clim_struct = "Climate ∩ Structure",
    struct_geo  = "Structure ∩ Geography",
    clim_geo    = "Climate ∩ Geography",
    all_three   = "All three",
    unexplained = "Unexplained"
)
GROUP_OF <- c(
    climate_u = "Unique to one factor", structure_u = "Unique to one factor", geo_u = "Unique to one factor",
    clim_struct = "Shared between factors", struct_geo = "Shared between factors",
    clim_geo = "Shared between factors", all_three = "Shared between factors",
    unexplained = "Unexplained"
)
tree_rows <- list()
add_tree_row <- function(region_key, value, p = NA_real_) {
    # NOTE: the column is named region_key, NOT "key" — data.table()'s
    # constructor reserves `key` as a formal argument (sets the table's
    # index), so data.table(key = ..., ...) silently does NOT create a "key"
    # column; it errors downstream instead ("some columns are not in the
    # data.table"). Caught by testing before this shipped.
    tree_rows[[length(tree_rows) + 1]] <<- data.table(
        region_key = region_key, group = GROUP_OF[[region_key]], component = COMPONENT[[region_key]],
        variance_pct = 100 * value, p_value = p)
}

have_struct <- !is.null(X_struct) && ncol(X_struct) > 0
have_geo    <- !is.null(X_geo_sel) && ncol(X_geo_sel) > 0
status      <- if (!have_geo) "no_spatial_vectors" else if (!have_struct) "no_structure_table" else "ok"
model_label <- NA_character_
confound_dt <- data.table(confounded = logical(), shared_pct = numeric(), max_unique_pct = numeric())

# 2-table climate-vs-geography — the dedicated confounding check, ALWAYS
# computed when geography is available. Also doubles as the MAIN tree when
# structure is unavailable (structure_table="none") — same fit, no wasted work.
if (have_geo) {
    vp2 <- tryCatch(varpart(Y, X_clim, X_geo_sel), error = function(e) { message("WARNING: 2-table varpart failed: ", e$message); NULL })
    if (!is.null(vp2)) {
        fr <- vp2$part$indfract
        get_fr <- function(pat) { idx <- grep(pat, rownames(fr), fixed = TRUE); if (length(idx)) fr[idx, adjr2_col(fr)] else NA_real_ }
        # vegan 2-table convention (empirically verified with designed test
        # data): [a]=table1 unique, [b]=table2 unique, [c]=shared, [d]=residual.
        # BUG FIX: this script previously read [b] as "shared" and [c] as
        # "table2 unique" (swapped) — which fed the WRONG two numbers into the
        # confounding-flag comparison below. Fixed here.
        clim_unique <- get_fr("[a]"); geo_unique <- get_fr("[b]"); shared_2 <- get_fr("[c]")
        resid_2     <- get_fr("[d]")
        rda_clim_u <- tryCatch(rda(Y, X_clim, X_geo_sel), error = function(e) NULL)
        rda_geo_u  <- tryCatch(rda(Y, X_geo_sel, X_clim), error = function(e) NULL)
        p_clim <- compute_p(rda_clim_u); p_geo <- compute_p(rda_geo_u)

        if (is.finite(shared_2) && is.finite(clim_unique) && is.finite(geo_unique)) {
            max_unique <- max(clim_unique, geo_unique, na.rm = TRUE)
            confounded <- shared_2 > max_unique
            confound_dt <- data.table(confounded = confounded, shared_pct = 100 * shared_2, max_unique_pct = 100 * max_unique)
            if (confounded) message("WARNING: climate-geography CONFOUNDED — shared (", round(shared_2, 4),
                                    ") exceeds the largest unique fraction (", round(max_unique, 4), ")")
        }

        if (!have_struct) {
            add_tree_row("climate_u", clim_unique, p_clim)
            add_tree_row("geo_u",     geo_unique,  p_geo)
            add_tree_row("clim_geo",  shared_2)
            add_tree_row("unexplained", resid_2)
            model_label <- "2-way (climate + geography)"
        }
    }
}

# 3-table: climate / structure / geography, when both auxiliary tables exist
if (have_struct && have_geo) {
    vp3 <- tryCatch(varpart(Y, X_clim, X_struct, X_geo_sel), error = function(e) { message("WARNING: 3-table varpart failed: ", e$message); NULL })
    if (!is.null(vp3)) {
        fr3 <- vp3$part$indfract
        get_fr3 <- function(lt) { idx <- grep(paste0("\\[", lt, "\\]"), rownames(fr3)); if (length(idx)) fr3[idx, adjr2_col(fr3)] else NA_real_ }
        # Verified empirically (designed tests, isolating each pairwise overlap):
        # [a]=climate unique, [b]=structure unique, [c]=geography unique,
        # [d]=climate∩structure, [e]=structure∩geography,
        # [f]=climate∩geography, [g]=all three, [h]=residual.
        rda_clim3 <- tryCatch(rda(Y, X_clim, cbind(X_struct, X_geo_sel)), error = function(e) NULL)
        rda_str3  <- tryCatch(rda(Y, X_struct, cbind(X_clim, X_geo_sel)), error = function(e) NULL)
        rda_geo3  <- tryCatch(rda(Y, X_geo_sel, cbind(X_clim, X_struct)), error = function(e) NULL)
        add_tree_row("climate_u",   get_fr3("a"), compute_p(rda_clim3))
        add_tree_row("structure_u", get_fr3("b"), compute_p(rda_str3))
        add_tree_row("geo_u",       get_fr3("c"), compute_p(rda_geo3))
        add_tree_row("clim_struct", get_fr3("d"))
        add_tree_row("struct_geo",  get_fr3("e"))
        add_tree_row("clim_geo",    get_fr3("f"))
        add_tree_row("all_three",   get_fr3("g"))
        add_tree_row("unexplained", get_fr3("h"))
        model_label <- "3-way (climate + structure + geography)"
    }
} else if (have_struct && !have_geo) {
    # Degenerate-input fallback (dbmem status != ok): climate-vs-structure 2-table instead.
    vp2s <- tryCatch(varpart(Y, X_clim, X_struct), error = function(e) NULL)
    if (!is.null(vp2s)) {
        frs <- vp2s$part$indfract
        get_frs <- function(pat) { idx <- grep(pat, rownames(frs), fixed = TRUE); if (length(idx)) frs[idx, adjr2_col(frs)] else NA_real_ }
        # Same [a]/[b]/[c]/[d] convention as the 2-table block above (same fix applies).
        clim_unique_s   <- get_frs("[a]"); struct_unique_s <- get_frs("[b]")
        shared_s        <- get_frs("[c]"); resid_s <- get_frs("[d]")
        rda_clim_s <- tryCatch(rda(Y, X_clim, X_struct), error = function(e) NULL)
        rda_str_s  <- tryCatch(rda(Y, X_struct, X_clim), error = function(e) NULL)
        add_tree_row("climate_u",   clim_unique_s,   compute_p(rda_clim_s))
        add_tree_row("structure_u", struct_unique_s, compute_p(rda_str_s))
        add_tree_row("clim_struct", shared_s)
        add_tree_row("unexplained", resid_s)
        model_label <- "2-way (climate + structure)"
    }
} else if (!have_struct && !have_geo) {
    # Last resort: neither structure nor geography available — climate alone.
    rda_clim_only <- tryCatch(rda(Y, X_clim), error = function(e) NULL)
    if (!is.null(rda_clim_only)) {
        r2_only <- suppressWarnings(RsquareAdj(rda_clim_only)$adj.r.squared)
        if (is.finite(r2_only)) {
            add_tree_row("climate_u", r2_only)
            add_tree_row("unexplained", 1 - r2_only)
            model_label <- "climate only (no structure or geography available)"
        }
    }
}

tree_dt <- if (length(tree_rows) > 0) rbindlist(tree_rows) else data.table(
    region_key = character(), group = character(), component = character(), variance_pct = numeric(), p_value = numeric())
tree_dt[, `:=`(model = model_label, status = status)]

# `region_key` is an internal plotting aid (which predictor-table region a
# row is — see the Venn plot below); the on-disk table stays the clean,
# human-readable schema only.
fwrite(tree_dt[, .(group, component, variance_pct, p_value, model, status)], OUT_VARPART_TSV, sep = "\t", quote = FALSE)
fwrite(confound_dt, OUT_CONFOUND_TSV, sep = "\t", quote = FALSE)
message("INFO: Wrote variance partition (", nrow(tree_dt), " rows, status=", status,
       ", model=", model_label %||% "none", ") and confounding check (", nrow(confound_dt), " row)")

################################################################################
# 7. Plots. dbMEM selection path stays a line plot (low-N, Rule 6 doesn't
#    apply). The variance partition is a SCHEMATIC Venn (fixed circle
#    positions, NOT area-proportional) via VennDiagram::draw.{single,pairwise,
#    triple}.venn() — a purpose-built, widely-used bioinformatics package for
#    exactly this diagram, chosen after two hand-rolled attempts (donutsk
#    nested donut, then ggforce::geom_arc_bar + ggrepel schematic circles)
#    both needed manually-guessed label anchor points. VennDiagram computes
#    its own region-label positions; drawn non-proportionally by design
#    (adjusted-R2 fractions are routinely negative — a real Euler/proportional
#    diagram can't represent that as area).
#
#    The package's own auto-generated region labels are just numbers, with no
#    way to pass custom text directly — so each region is drawn with a
#    UNIQUE PLACEHOLDER value first, then that grob's $label is swapped for
#    our real "Component\nXX.X%" string post-hoc (a standard, documented
#    customization pattern for this package). draw.triple.venn()'s
#    area.vector order (verified from the package's own source, NOT the
#    argument order): [area1-only, n12-only, area2-only, n13-only, all-three,
#    n23-only, area3-only], with area1=Climate/area2=Structure/area3=Geography
#    (verified via the category-label coordinates it returns).
#
#    Unexplained itself is NOT drawn (it's everything outside the circles,
#    the standard Venn convention) — its magnitude is instead surfaced as the
#    "X% variance explained" badge next to the plot in Shiny (mod_climate.R).
#
#    Every LABEL always shows the real, unclamped, signed value (CLAUDE.md
#    Rule 8 — a value label is data, not baked-in commentary). The "read a
#    negative fraction as ~0 (Peres-Neto et al. 2006)" interpretation lives in
#    the Shiny help note, not on the plot.
################################################################################
save_both <- function(name, g, w = 7, h = 5) clinego_save_both(file.path(PLOT_DIR, name), g, w, h)
empty_plot <- clinego_empty_plot

if (nrow(selection_dt) > 0) {
    g_path <- ggplot(selection_dt, aes(x = step, y = r2_adj_cumulative)) +
        geom_hline(aes(yintercept = r2_adj_full_ceiling), color = CLINEGO_THRESHOLD, linetype = "dashed") +
        geom_line(color = CLINEGO_NEUTRAL) + geom_point(color = CLINEGO_NEUTRAL, size = 2) +
        ggrepel::geom_text_repel(aes(label = variable), size = 3, color = CLINEGO_COL$fg) +
        labs(x = "Forward-selection step", y = "Cumulative adjusted R2",
            title = "dbMEM forward-selection path (A17)") +
        theme_clinego()
    save_both("dbmem_selection_path", g_path)
} else save_both("dbmem_selection_path", empty_plot("No dbMEM selection path (geography unavailable)"))

# Swap placeholder numeric grob labels (as drawn by VennDiagram) for our real
# text, matched by exact placeholder string — order-independent, robust to
# whatever internal grob order the package uses.
replace_grob_labels <- function(g, mapping) {
    for (i in seq_along(g)) {
        if (inherits(g[[i]], "text")) {
            lbl <- as.character(g[[i]]$label)
            if (lbl %in% names(mapping)) g[[i]]$label <- mapping[[lbl]]
        }
    }
    g
}
region_label <- function(row) sprintf("%s\n%.1f%%", row$component, row$variance_pct)

venn_ok <- tryCatch({
    region_dt <- tree_dt[region_key != "unexplained"]
    if (nrow(region_dt) == 0) stop("no variance-partition data")

    has_clim     <- "climate_u"   %in% region_dt$region_key
    has_struct_r <- "structure_u" %in% region_dt$region_key
    has_geo_r    <- "geo_u"       %in% region_dt$region_key
    n_sets <- sum(has_clim, has_struct_r, has_geo_r)
    SET_COLOR <- c(Climate = CLINEGO_CATEGORICAL[1], Structure = CLINEGO_CATEGORICAL[2], Geography = CLINEGO_CATEGORICAL[3])
    by_key <- function(k) region_dt[region_key == k]

    v <- if (n_sets == 3) {
        order3 <- c("climate_u", "clim_struct", "structure_u", "clim_geo", "all_three", "struct_geo", "geo_u")
        mapping <- setNames(vapply(order3, function(k) region_label(by_key(k)), character(1)), as.character(1:7))
        draw.triple.venn(area1 = 100, area2 = 100, area3 = 100, n12 = 40, n23 = 40, n13 = 40, n123 = 20,
            category = c("Climate", "Structure", "Geography"),
            fill = unname(SET_COLOR[c("Climate", "Structure", "Geography")]), col = "white", alpha = 0.45,
            cat.col = CLINEGO_COL$fg, cat.cex = 1.15, cat.fontface = "bold",
            label.col = CLINEGO_COL$fg, cex = 0.9,
            direct.area = TRUE, area.vector = 1:7, ind = FALSE)
    } else if (n_sets == 2) {
        other <- if (has_struct_r) "Structure" else "Geography"
        other_key <- if (has_struct_r) "structure_u" else "geo_u"
        shared_key <- if (has_struct_r) "clim_struct" else "clim_geo"
        # No direct.area for pairwise — pick area1/area2/cross.area so the
        # three AUTO-COMPUTED displayed numbers (area1-cross, cross,
        # area2-cross) are distinct known placeholders to match on.
        v0 <- draw.pairwise.venn(area1 = 1000, area2 = 2000, cross.area = 100,
            category = c("Climate", other), euler.d = FALSE, scaled = FALSE,
            fill = unname(SET_COLOR[c("Climate", other)]), col = "white", alpha = 0.45,
            cat.col = CLINEGO_COL$fg, cat.cex = 1.15, cat.fontface = "bold",
            label.col = CLINEGO_COL$fg, cex = 0.9, ind = FALSE)
        mapping <- c(`900` = region_label(by_key("climate_u")), `1900` = region_label(by_key(other_key)),
                    `100` = region_label(by_key(shared_key)))
        v0
    } else {
        # Last resort: climate alone, no overlap geometry.
        v0 <- draw.single.venn(area = 1, category = "Climate",
            fill = unname(SET_COLOR["Climate"]), col = "white", alpha = 0.45,
            cat.col = CLINEGO_COL$fg, cat.cex = 1.15, cat.fontface = "bold",
            label.col = CLINEGO_COL$fg, cex = 0.9, ind = FALSE)
        mapping <- c(`1` = region_label(by_key("climate_u")))
        v0
    }
    v <- replace_grob_labels(v, mapping)
    v <- add.title(gList = v, x = "Variance partition", pos = c(0.5, 1.05),
                   cex = 1.3, fontface = "bold", col = CLINEGO_COL$fg)
    v <- add.title(gList = v, x = model_label, pos = c(0.5, 1.00), cex = 0.9, col = CLINEGO_COL$muted)

    grDevices::png(file.path(PLOT_DIR, "varpart_venn.png"), width = 7, height = 7.5, units = "in", res = 300, bg = "white")
    grid::grid.newpage(); grid::grid.draw(v); grDevices::dev.off()
    svglite::svglite(file.path(PLOT_DIR, "varpart_venn.svg"), width = 7, height = 7.5, bg = "white")
    grid::grid.newpage(); grid::grid.draw(v); grDevices::dev.off()
    TRUE
}, error = function(e) {
    message("WARNING: variance-partition Venn failed: ", e$message)
    FALSE
})
if (!venn_ok) save_both("varpart_venn", empty_plot(paste0("Variance partition unavailable (", status, ")")))

if (nrow(px_dt) > 0) {
    g_px <- ggplot(px_dt, aes(x = reorder(variable, Px), y = Px_pct, fill = model)) +
        geom_col(position = "dodge") + coord_flip() +
        scale_fill_clinego(name = NULL) +
        labs(x = NULL, y = "Px (% of available variance)", title = "Lasky Px per climate variable") +
        theme_clinego()
    save_both("px_barplot", g_px)
} else save_both("px_barplot", empty_plot("No Px values computed"))

message("INFO: preGEA varpart complete.")
