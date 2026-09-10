# Module path constants — match workflow/rules/common.smk MOD_* values
MOD_PROC      <- "Processing"
MOD_PRESTRUCT <- "PreStructure"
MOD_PREGEA    <- "PreGEA"
MOD_STRUCT    <- "Structure"
MOD_CLIMATE   <- "climate"
MOD_TRAITS    <- "Traits"
MOD_GEA       <- "GEA"
MOD_GWAS      <- "GWAS"
MOD_GEAXGWAS  <- "GEAxGWAS"
MOD_MALAD     <- "Maladaptation"
MOD_INTER     <- "_intermediate"

#' Base results directory for a project
#' @noRd
project_base <- function(project, pipeline_path = get_pipeline_path()) {
    file.path(pipeline_path, paste0(project, "_results"))
}

#' Path within a module's output directory
#' @noRd
mod_path <- function(project, module, ...) {
    file.path(project_base(project), module, ...)
}

# ─── Manhattan plots ─────────────────────────────────────────────────────────

#' Per-method Manhattan background PNG path
#' Background shows all SNPs (threshold-independent); threshold line drawn live in Shiny.
#' @noRd
manhattan_bg_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_background.png"))
}

#' Per-method coordinate mapping JSON path
#' @noRd
manhattan_coords_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_", trait, "_K", k, "_", adjust, "_coords.json"))
}

#' Combined Manhattan background PNG path
#' @noRd
combined_manhattan_bg_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_background.png"))
}

#' Combined Manhattan coordinate mapping JSON path
#' @noRd
combined_manhattan_coords_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_combined_K", k, "_coords.json"))
}

# ─── WZA Manhattan plots ──────────────────────────────────────────────────────

#' WZA per-method Manhattan background PNG path
#' @noRd
manhattan_wza_bg_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_wza_", trait, "_K", k, "_", adjust, "_background.png"))
}

#' WZA per-method coordinate mapping JSON path
#' @noRd
manhattan_wza_coords_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("manhattan_wza_", trait, "_K", k, "_", adjust, "_coords.json"))
}

#' WZA combined Manhattan background PNG path
#' @noRd
combined_manhattan_wza_bg_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_wza_combined_K", k, "_background.png"))
}

#' WZA combined Manhattan coordinate mapping JSON path
#' @noRd
combined_manhattan_wza_coords_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("manhattan_wza_combined_K", k, "_coords.json"))
}

#' WZA sig-windows table (exact K — for loading in Shiny overlay)
#' @noRd
method_wza_sigwindows_path <- function(project, module = MOD_GEA, method, k, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_wza_K", k, "_sig_windows_", adjust, ".tsv"))
}

#' WZA p-values TSV (windows, all traits)
#' @noRd
method_wza_pvalues_path <- function(project, module = MOD_GEA, method, k) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_wza_K", k, ".tsv"))
}

#' Directory for disk-cached interactive sig SNPs (keyed by MD5 hash of threshold params).
#' @noRd
interactive_sigsnps_dir <- function(project, module = MOD_GEA) {
    mod_path(project, MOD_INTER, "interactive_sigsnps", module)
}

# ─── RDA (glob patterns — K isn't user-selectable on this tab, resolved via
#     Sys.glob() in fct_data_loading.R, same convention as method_sigsnps_path) ─

#' RDA candidates side table (glob pattern for discovery)
#' @noRd
rda_candidates_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "methods", "RDA", "RDA_candidates_K*.tsv")
}

#' RDA diagnostics table (glob pattern for discovery)
#' @noRd
rda_diagnostics_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "methods", "RDA", "RDA_diagnostics_K*.tsv")
}

#' RDA anova table (glob pattern for discovery)
#' @noRd
rda_anova_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "methods", "RDA", "RDA_anova_K*.tsv")
}

#' RDA diagnostic plot path (screeplot | pvalue_histogram | loadings_biplot)
#' @noRd
rda_plot_path <- function(project, module = MOD_GEA, which) {
    glob <- mod_path(project, module, "plots", "rda", paste0("rda_", which, "_K*.png"))
    hits <- Sys.glob(glob)
    if (length(hits) > 0) hits[[1]] else glob
}

#' QQ plot path (per-method)
#' @noRd
qq_plot_path <- function(project, module = MOD_GEA, method, trait, k, adjust) {
    mod_path(project, module, "plots", "manhattan", method,
             paste0("qq_", trait, "_K", k, "_", adjust, ".png"))
}

#' Combined QQ plot path
#' @noRd
qq_combined_path <- function(project, module = MOD_GEA, k) {
    mod_path(project, module, "plots", "manhattan", "combined",
             paste0("qq_combined_K", k, ".png"))
}

# ─── Miami plot ───────────────────────────────────────────────────────────────

#' Miami background PNG path
#' @noRd
miami_bg_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_combined_K", k, "_background.png"))
}

#' Miami coordinate mapping JSON path
#' @noRd
miami_coords_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_combined_K", k, "_coords.json"))
}

#' WZA Miami background PNG path
#' @noRd
miami_wza_bg_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_wza_combined_K", k, "_background.png"))
}

#' WZA Miami coordinate mapping JSON path
#' @noRd
miami_wza_coords_path <- function(project, k) {
    mod_path(project, MOD_GEAXGWAS, "plots",
             paste0("miami_wza_combined_K", k, "_coords.json"))
}

#' WZA pairwise overlap table (pipeline-computed)
#' @noRd
pairwise_wza_table_path <- function(project) {
    mod_path(project, MOD_GEAXGWAS, "tables", "pairwise_overlap_wza.tsv")
}

# ─── Association tables ───────────────────────────────────────────────────────

#' Glob all per-method sig SNPs files for a module (across all methods and adjustments)
#' @noRd
find_method_sigsnps_files <- function(project, module = MOD_GEA) {
    Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_sig_snps_*.tsv"))
}

#' Pairwise trait overlap table (pipeline-computed)
#' @noRd
pairwise_table_path <- function(project) {
    mod_path(project, MOD_GEAXGWAS, "tables", "pairwise_overlap_table.tsv")
}

#' Pairwise collapsed sig SNPs (long format, one row per SNP per trait)
#' @noRd
pairwise_collapsed_path <- function(project) {
    mod_path(project, MOD_GEAXGWAS, "tables", "pairwise_collapsed_snps.tsv")
}

#' Selected SNPs table (combined all methods)
#' @noRd
selected_snps_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "selected_snps.tsv")
}

#' Per-method sig SNPs table (glob pattern for discovery)
#' @noRd
method_sigsnps_path <- function(project, module = MOD_GEA, method, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_pvalues_K*_sig_snps_", adjust, ".tsv"))
}

#' Per-method sig SNPs table (exact K — for loading in Shiny overlay)
#' @noRd
method_sigsnps_direct_path <- function(project, module = MOD_GEA, method, k, adjust) {
    mod_path(project, module, "tables", "methods", method,
             paste0(method, "_pvalues_K", k, "_sig_snps_", adjust, ".tsv"))
}

#' Per-trait regions table
#' @noRd
regions_per_trait_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "regions_per_trait.tsv")
}

#' Combined regions table
#' @noRd
regions_combined_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "regions_combined.tsv")
}

#' Genes per region table
#' @noRd
genes_per_region_path <- function(project, module = MOD_GEA) {
    mod_path(project, module, "tables", "genes_per_region_collapsed.tsv")
}

#' GO enrichment table for a specific region and trait
#' @noRd
enrichment_table_path <- function(project, module = MOD_GEA, trait, region_id) {
    mod_path(project, module, "tables", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_enrichment.tsv"))
}

# ─── Enrichment plots ─────────────────────────────────────────────────────────

#' Enrichment plot path
#' @noRd
enrichment_plot_path <- function(project, module = MOD_GEA, trait, region_id, plot_type) {
    mod_path(project, module, "plots", "enrichment", trait,
             paste0("Region_", region_id, "_", trait, "_", plot_type, ".png"))
}

# ─── PreStructure ─────────────────────────────────────────────────────────────

#' PCA plot path
#' @noRd
pca_path <- function(project) {
    mod_path(project, MOD_PRESTRUCT, "plots", "pca.png")
}

#' Tracy-Widom path
#' @noRd
tracy_widom_path <- function(project) {
    mod_path(project, MOD_PRESTRUCT, "plots", "tracy_widom.png")
}

#' Cross-entropy plot path
#' @noRd
cross_entropy_path <- function(project, k_start, k_end) {
    mod_path(project, MOD_PRESTRUCT, "plots",
             paste0("cross_entropy_K", k_start, "-", k_end, ".png"))
}

#' Structure barplot path for a given K
#' @noRd
structure_k_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("structure_K", k, ".png"))
}

#' PCA structure plot path for a given K
#' @noRd
pca_structure_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("pca_structure_K", k, ".png"))
}

#' Pop differentiation path for a given K
#' @noRd
pop_diff_path <- function(project, k) {
    mod_path(project, MOD_PRESTRUCT, "plots", paste0("K", k),
             paste0("pop_diff_K", k, ".png"))
}

# ─── Structure ────────────────────────────────────────────────────────────────

#' Piemap path
#' @param points if TRUE, resolve the points ("clear map") companion instead of
#'   the pie chart. Points are trait/metric-independent, so `metric` is ignored.
#' @noRd
piemap_path <- function(project, bio, metric = NULL, zoom = NULL, points = FALSE) {
    dir <- mod_path(project, MOD_STRUCT, "plots", "piemap")
    metric <- if (points) NULL else metric
    suffix <- if (points) "_points" else if (!is.null(metric) && nzchar(metric) && metric != "none") paste0("_", metric) else ""
    if (!is.null(zoom) && nzchar(zoom) && zoom != "none") {
        dir <- file.path(dir, "zoom", zoom)
    }
    fname <- paste0("piemap_bio_", bio, suffix, ".png")
    file.path(dir, fname)
}

#' Climate correlation heatmap path
#' @noRd
climate_heatmap_path <- function(project) {
    mod_path(project, MOD_CLIMATE, "plots", "correlation_heatmap.png")
}

#' Climate density plot path
#' @noRd
climate_density_path <- function(project, bio = NULL) {
    if (is.null(bio)) {
        mod_path(project, MOD_CLIMATE, "plots", "density_plot_present.png")
    } else {
        mod_path(project, MOD_CLIMATE, "plots", paste0("density_plot_present_bio", bio, ".png"))
    }
}

#' Phenotype density plot path
#' Produced by mode=traits since Phase 3 (was mode=climate) — it describes
#' traits, and had to survive Climate.enabled: false.
#' @noRd
phenotype_density_path <- function(project) {
    mod_path(project, MOD_TRAITS, "plots", "density_plot_phenotypes.png")
}

# ─── Traits (mode=traits) ─────────────────────────────────────────────────────

#' Generic Traits plot path: Traits/plots/{stem}.png
#' @noRd
trait_plot_path <- function(project, stem) {
    mod_path(project, MOD_TRAITS, "plots", paste0(stem, ".png"))
}

#' Generic Traits table path: Traits/tables/{stem}.tsv
#' @noRd
trait_table_path <- function(project, stem) {
    mod_path(project, MOD_TRAITS, "tables", paste0(stem, ".tsv"))
}

#' Trait correlogram path
#' @param with_climate when TRUE, resolve the joint traits x climate correlogram
#'   (standard regime only — absent when the project has no climate)
#' @noRd
trait_corr_path <- function(project, with_climate = FALSE) {
    trait_plot_path(project, if (with_climate) "correlation_heatmap_traits_climate"
                             else "correlation_heatmap_traits")
}

#' Climate invariant predictors file path
#' @noRd
climate_invariant_path <- function(project) {
    mod_path(project, MOD_CLIMATE, "tables", "present", "climate_invariant_predictors.tsv")
}

#' Phenotype missing summary path (GWAS mode)
#' @noRd
pheno_missing_summary_path <- function(project) {
    mod_path(project, MOD_GWAS, "tables", "phenotype_missing_summary.tsv")
}

#' LEA .removed file path — glob from _work/ (LD-pruned variant)
#' Returns empty string if not found.
#' @noRd
removed_snps_path <- function(project) {
    base  <- project_base(project)
    files <- Sys.glob(file.path(base, "_work", "maf*", paste0(project, ".removed")))
    if (length(files) == 0) "" else files[1]
}

#' Glob first per-SNP pvalue TSV for a module/method (any K, any adjust).
#' Used to extract trait column names. Excludes snp-mode `_sig_snps_` and
#' @noRd
find_pvalue_tsv <- function(project, module = MOD_GEA) {
    files <- Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_pvalues_K*.tsv"))
    files <- grep("_sig_snps_|_wza_|_block_pvalues_", files, value = TRUE, invert = TRUE)
    files <- keep_merged_method_files(files, regime = "snp")
    if (length(files) == 0) "" else files[1]
}

#' All per-method pvalue TSV paths for a module (one per method).
#' Used by regionplot subprocess discovery instead of sig_snps files.
#' @noRd
find_method_pvalue_files <- function(project, module = MOD_GEA) {
    files <- Sys.glob(mod_path(project, module, "tables", "methods", "*", "*_pvalues_K*.tsv"))
    files <- grep("_sig_snps_|_wza_|_block_pvalues_", files, value = TRUE, invert = TRUE)
    keep_merged_method_files(files, regime = "snp")
}

#' LD decay plot path
#' @noRd
ld_decay_path <- function(project, per_chr = FALSE) {
    fname <- if (per_chr) "ld_decay_per_chromosome.png" else "ld_decay_genome_wide.png"
    mod_path(project, MOD_STRUCT, "plots", "ld_decay", fname)
}

#' LD decay half-distances TSV path
#' @noRd
ld_decay_table_path <- function(project) {
    mod_path(project, MOD_STRUCT, "tables", "ld_decay_half_distances.tsv")
}

# ─── PreGEA (hyperparameter-exploration mode) ─────────────────────────────────
# LFMM/EMMAX artifacts are grid-level only — plot_pregea_ladder.R emits one
# histogram grid / QQ grid / lambda-vs-K/#PCs / hits-vs-K/#PCs PNG per engine,
# never one file per K or #PC rung — so there is no per-rung selector for
# those two ladders; the ladder tables (below) carry the per-rung numbers.
# EXCEPTION: RDA IS per-model (one biplot/screeplot/axis-anova set per swept
# condition_pcs value, see rda_model_*_path() above) — the Shiny RDA panel
# selects between them (module-split follow-up: predictor characterization
# moved to the Climate tab, and RDA's model comparison gained per-rung detail).
# Two generic helpers below (not one function per artifact, unlike the rest of
# this file) because every other PreGEA output is a uniform {subdir}/{stem}.{ext}
# — see mod_pregea.R for the concrete subdir/stem values actually used.

#' Generic PreGEA plot path: PreGEA/plots/{subdir}/{stem}.png
#' @noRd
pregea_plot_path <- function(project, subdir, stem) {
    mod_path(project, MOD_PREGEA, "plots", subdir, paste0(stem, ".png"))
}

#' Generic PreGEA table path: PreGEA/tables/{subdir}/{stem}.tsv
#' subdir = "" resolves to PreGEA/tables/{stem}.tsv (recommendations, transfer guard)
#' @noRd
pregea_table_path <- function(project, subdir, stem) {
    if (nzchar(subdir)) mod_path(project, MOD_PREGEA, "tables", subdir, paste0(stem, ".tsv"))
    else mod_path(project, MOD_PREGEA, "tables", paste0(stem, ".tsv"))
}

#' RDA per-model artifact paths — one set per condition_pcs value the PreGEA
#' RDA ladder swept (pregea_rda_setup.R writes these inside one rule's
#' internal loop, not a Snakemake fan-out — see that rule's header). Mirrors
#' common.smk's rda_model_biplot()/rda_model_screeplot()/rda_model_axis_anova()
#' Python helpers exactly: PreGEA/{plots,tables}/rda/models/pc{n}/...
#' @noRd
rda_model_biplot_path <- function(project, n) {
    mod_path(project, MOD_PREGEA, "plots", "rda", "models", paste0("pc", n), "biplot.png")
}
rda_model_screeplot_path <- function(project, n) {
    mod_path(project, MOD_PREGEA, "plots", "rda", "models", paste0("pc", n), "axis_screeplot.png")
}
rda_model_axis_anova_path <- function(project, n) {
    mod_path(project, MOD_PREGEA, "tables", "rda", "models", paste0("pc", n), "axis_anova.tsv")
}

# ─── Climate (predictor characterization + spatial/varpart — mode=climate) ───
# Two generic helpers, same convention as pregea_plot_path()/pregea_table_path()
# above — every Climate output is a uniform {subdir}/{stem}.{ext}.

#' Generic Climate plot path: climate/plots/{subdir}/{stem}.png
#' subdir = "" resolves to climate/plots/{stem}.png (correlation_heatmap, density plots)
#' @noRd
climate_plot_path <- function(project, subdir, stem) {
    if (nzchar(subdir)) mod_path(project, MOD_CLIMATE, "plots", subdir, paste0(stem, ".png"))
    else mod_path(project, MOD_CLIMATE, "plots", paste0(stem, ".png"))
}

#' Generic Climate table path: climate/tables/{subdir}/{stem}.tsv
#' @noRd
climate_table_path <- function(project, subdir, stem) {
    mod_path(project, MOD_CLIMATE, "tables", subdir, paste0(stem, ".tsv"))
}

# ─── Maladaptation ───────────────────────────────────────────────────────────

#' Directory for a curated SNP set (produced by Shiny GEA Save)
#' @noRd
snp_sets_dir <- function(project) {
    mod_path(project, MOD_INTER, "snp_sets")
}

#' Path to a specific SNP set's selected_snps.tsv
#' @noRd
snp_set_path <- function(project, name) {
    mod_path(project, MOD_INTER, "snp_sets", name, "selected_snps.tsv")
}

#' Path to the SNP sets manifest JSON
#' @noRd
snp_sets_manifest_path <- function(project) {
    mod_path(project, MOD_INTER, "snp_sets", "manifest.json")
}

#' GF selected SNPs path (intermediate file under _intermediate/{method}/{suffix}/)
#' @noRd
gf_selected_snps_path <- function(project, suffix, method = "gradient_forest") {
    mod_path(project, MOD_INTER, method, suffix, "selected_snps.tsv")
}

#' GF importance plot path
#' @noRd
gf_importance_path <- function(project, suffix, type = "overall", method = "gradient_forest") {
    mod_path(project, MOD_MALAD, "plots", method, suffix,
             paste0(type, "_importance.png"))
}

#' Geometric-offset conditioning diagnostics table (long key/value)
#' @noRd
geometric_offset_diagnostics_path <- function(project, suffix) {
    mod_path(project, MOD_MALAD, "tables", "geometric_offset", suffix,
             "geometric_offset_diagnostics.tsv")
}

#' Gradient-Forest imputation-sensitivity diagnostics table (long key/value)
#'
#' Written by scripts/gf_sensitivity_check.R. Same two-column (quantity, value) shape as
#' geometric_offset_diagnostics.tsv.
#' @noRd
gf_sensitivity_path <- function(project, suffix, method = "gradient_forest") {
    mod_path(project, MOD_MALAD, "tables", method, suffix, "imputation_sensitivity.tsv")
}

#' Scenario names available for one method/SNP-set, newest layout only
#'
#' Offset products gained a {scenario} level when the pipeline learned to project
#' one model onto several futures. Returns character(0) for a pre-scenario tree.
#' @noRd
find_offset_scenarios <- function(project, suffix, method = "gradient_forest") {
    base <- mod_path(project, MOD_MALAD, "tables", method, suffix)
    if (!dir.exists(base)) return(character(0))
    d <- list.dirs(base, full.names = FALSE, recursive = FALSE)
    sort(d[nzchar(d)])
}

#' Resolve one offset artefact, tolerating the pre-scenario layout
#'
#' Offsets moved from `{suffix}/{file}` to `{suffix}/{scenario}/{file}`. Results
#' trees written before that change have no scenario directory, so this is the
#' single place that knows about both shapes:
#'   - an explicit `scenario` wins when that file exists;
#'   - otherwise the flat legacy path is used when it exists;
#'   - otherwise the first available scenario is used, so a scenario-era tree
#'     renders without the caller having to pick one.
#' When nothing exists it returns the modern shape, so "missing" reads as the
#' location the pipeline would write today.
#' @noRd
mala_offset_file <- function(project, kind, method, suffix, fname, scenario = NULL) {
    base <- mod_path(project, MOD_MALAD, kind, method, suffix)
    if (!is.null(scenario) && nzchar(scenario)) {
        p <- file.path(base, scenario, fname)
        if (file.exists(p)) return(p)
    }
    legacy <- file.path(base, fname)
    if (file.exists(legacy)) return(legacy)
    scen <- find_offset_scenarios(project, suffix, method)
    if (length(scen)) return(file.path(base, scen[1], fname))
    if (!is.null(scenario) && nzchar(scenario)) file.path(base, scenario, fname) else legacy
}

#' GF genetic offset piemap path
#' @noRd
gf_offset_piemap_path <- function(project, suffix, variant = "base", method = "gradient_forest",
                                  scenario = NULL) {
    # "base" maps to "notrait" suffix; other variants used as-is
    piemap_variant <- if (is.null(variant) || variant == "base" || !nzchar(variant)) "notrait" else variant
    fname <- paste0("genetic_offset_piemap_", piemap_variant, ".png")
    mala_offset_file(project, "plots", method, suffix, fname, scenario)
}

#' GF site offset table path
#' @noRd
gf_site_table_path <- function(project, suffix, method = "gradient_forest", scenario = NULL) {
    mala_offset_file(project, "tables", method, suffix, "genetic_offset_site.tsv", scenario)
}

#' GF map offset table path (landscape cells)
#' @noRd
gf_map_table_path <- function(project, suffix, method = "gradient_forest", scenario = NULL) {
    mala_offset_file(project, "tables", method, suffix, "genetic_offset_map.tsv", scenario)
}

# ─── Model comparison paths ───────────────────────────────────────────────────

#' Stable comparison key for two models (sorted so A/B order doesn't matter for cache)
#' @noRd
model_compare_key <- function(key_a, key_b) {
    sorted <- sort(c(key_a, key_b))
    paste0(
        gsub(":::", "_", sorted[1], fixed = TRUE), "__vs__",
        gsub(":::", "_", sorted[2], fixed = TRUE)
    )
}

#' Cache directory for a two-model comparison
#' @noRd
model_compare_dir <- function(project, key_a, key_b) {
    mod_path(project, MOD_INTER, "model_compare", model_compare_key(key_a, key_b))
}

#' Novelty cache directory (per scenario — shared across comparisons)
#' @noRd
novelty_cache_dir <- function(project, scenario_label) {
    mod_path(project, MOD_INTER, "novelty", scenario_label)
}

#' stats.json path for a model comparison
#' @noRd
model_compare_stats_path <- function(project, key_a, key_b) {
    file.path(model_compare_dir(project, key_a, key_b), "stats.json")
}

#' Disagreement raster path
#' @noRd
model_compare_disagree_raster <- function(project, key_a, key_b) {
    file.path(model_compare_dir(project, key_a, key_b), "disagree_rank.tif")
}

#' Site offsets table path (merged, two columns offset_a / offset_b)
#' @noRd
model_compare_site_offsets <- function(project, key_a, key_b) {
    file.path(model_compare_dir(project, key_a, key_b), "site_offsets.tsv")
}

#' Rank stability table path
#' @noRd
model_compare_rank_stability <- function(project, key_a, key_b) {
    file.path(model_compare_dir(project, key_a, key_b), "rank_stability.tsv")
}

#' ExDet novelty raster path
#' @noRd
novelty_raster_path <- function(project, scenario_label) {
    file.path(novelty_cache_dir(project, scenario_label), "exdet_novelty.tif")
}

# ─── Phenotype piemap ────────────────────────────────────────────────────────

#' Phenotype association piemap (phenomap) path
#' @param points if TRUE, resolve the points ("clear map") companion instead of
#'   the pie chart. Points are trait-independent, so `trait` is ignored.
#' @noRd
pheno_piemap_path <- function(project, trait, points = FALSE) {
    fname <- if (points) "phenomap_points.png" else paste0("phenomap_", trait, ".png")
    mod_path(project, MOD_GWAS, "plots", "piemap", fname)
}

# ─── Input data paths ─────────────────────────────────────────────────────────

#' Resolve full path to the project's GFF annotation file
#' @noRd
gff_path <- function(project, pipeline_path = get_pipeline_path()) {
    config  <- read_project_config(project, pipeline_path)
    inp_dir <- config_get(config, "Input", "dir", default = "data")
    gff_rel <- config_get(config, "Input", "gff", default = "")
    if (is.null(gff_rel) || gff_rel == "") return("")
    file.path(pipeline_path, inp_dir, gff_rel)
}

# ─── Regional Manhattan (topr) ────────────────────────────────────────────────

#' Directory for on-demand regionplots persisted to disk.
#' Subdirectory per combo_hash keeps files for different filter combinations separate.
#' @noRd
ondemand_regionplot_dir <- function(project, module = MOD_GEA, combo_hash) {
    mod_path(project, "regionplot", module, combo_hash)
}

# ─── On-demand pairwise trait overlap ─────────────────────────────────────────

#' Directory for on-demand pairwise trait overlap results.
#' Keyed by a parameter hash so different threshold/clumping combinations are cached separately.
#' Lives under GEAxGWAS/pairwise_ondemand/{hash}/ in the project results dir.
#' @noRd
ondemand_pairwise_dir <- function(project, hash) {
    mod_path(project, MOD_GEAXGWAS, "pairwise_ondemand", hash)
}

#' Regionplot image path
#' @noRd
regionplot_path <- function(project, region_id, trait = NULL) {
    # Sanitize region_id: pipeline replaces colons and hyphens with underscores
    safe_id <- gsub("[:\\-]", "_", region_id)
    fname <- if (is.null(trait)) {
        paste0("regionplot_", safe_id, ".png")
    } else {
        paste0("regionplot_", safe_id, "_", trait, ".png")
    }
    mod_path(project, "regionplot", fname)
}

# ─── Haplotype ────────────────────────────────────────────────────────────────
# All haplotype outputs live under _intermediate/haplotype/{tag}/ (Shiny on-demand only;
# no pipeline haplotype modes since Phase 3).

#' Haplotype clustree path (per-region)
#' @noRd
hap_clustree_path <- function(project, tag, region_id, type = "MG") {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots", "clustree",
             paste0("Region_", region_id, "_clustree_", type, ".png"))
}

#' Haplotype crosshap_viz path (per-trait)
#' @noRd
crosshap_viz_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
             paste0("Region_", region_id, "_crosshap_viz_", trait, ".png"))
}

#' Haplotype boxplot path (per-trait)
#' @noRd
hap_boxplot_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
             paste0("Region_", region_id, "_boxplot_", trait, ".png"))
}

#' Haplotype piemap path (per-trait)
#' @noRd
hap_piemap_path <- function(project, tag, region_id, trait) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots",
             paste0("haplotype_piemap_region_", region_id, "_", trait, ".png"))
}

#' Haplotype scan status path
#' @noRd
hap_scan_status_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "scan_status.tsv")
}

#' Haplotype selected regions path
#' @noRd
hap_selected_regions_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "selected_regions.tsv")
}

#' Haplotype assignments table path
#' @noRd
hap_assignments_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "Haplotype_assignments.tsv")
}

#' Haplotype frequencies table path
#' @noRd
hap_frequencies_path <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables", "Haplotype_frequencies.tsv")
}

#' Haplotype intermediate directory (HapObject .qs files)
#' @noRd
hap_intermediate_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag)
}

#' Haplotype scan clustree plots output directory
#' @noRd
hap_scan_plots_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots", "clustree")
}

#' Haplotype viz plots output directory
#' @noRd
hap_viz_plots_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "plots")
}

#' Haplotype viz tables output directory
#' @noRd
hap_viz_tables_dir <- function(project, tag) {
    mod_path(project, MOD_INTER, "haplotype", tag, "tables")
}

#' Derive VCF basename from config (mirrors pipeline's VCF_BASE).
#' Falls back to `project` if config is missing an Input.vcf entry.
#' @noRd
.vcf_base <- function(project, config = NULL) {
    vcf <- if (!is.null(config)) config_get(config, "Input", "vcf", default = NULL) else NULL
    if (is.null(vcf) || !nzchar(vcf)) return(project)
    name <- basename(vcf)
    name <- sub("\\.vcf\\.gz$", "", name)
    name <- sub("\\.vcf$",     "", name)
    name
}

#' Filtered VCF path (glob from _work, returns first match)
#' @noRd
hap_filtered_vcf_path <- function(project, config = NULL) {
    base <- project_base(project)
    vb <- .vcf_base(project, config)
    vcfs <- Sys.glob(file.path(base, "_work", "maf*", paste0(vb, ".vcf")))
    if (length(vcfs) == 0) "" else vcfs[1]
}

#' Imputed VCF path for LD computation (glob from _work, returns first match)
#' @noRd
hap_imputed_vcf_path <- function(project, k_best, config = NULL) {
    base <- project_base(project)
    vb <- .vcf_base(project, config)
    vcfs <- Sys.glob(file.path(base, "_work", "maf*",
                               paste0(vb, "_K", k_best, "imp.vcf")))
    if (length(vcfs) == 0) "" else vcfs[1]
}

#' Aligned metadata path
#' @noRd
hap_metadata_path <- function(project) {
    mod_path(project, MOD_PROC, "tables", "metadata.tsv")
}

#' Clusters table path for k_best
#' @noRd
hap_clusters_path <- function(project, k_best) {
    mod_path(project, MOD_PRESTRUCT, "tables", paste0("K", k_best),
             paste0("clusters_K", k_best, ".tsv"))
}

# ── Processing QC paths ──────────────────────────────────────────────────────

#' QC plot path (Processing/plots/<filename>)
#' @noRd
qc_plot_path <- function(project, filename) {
    mod_path(project, MOD_PROC, "plots", filename)
}

#' QC table path (Processing/tables/<filename>)
#' @noRd
qc_table_path <- function(project, filename) {
    mod_path(project, MOD_PROC, "tables", filename)
}

#' Pipeline summary TSV path
#' @noRd
pipeline_summary_path <- function(project) {
    file.path(project_base(project), "pipeline_summary.tsv")
}
