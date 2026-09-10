#' Central registry of plot help-notes: what a plot shows + which config controls it.
#'
#' One entry per note-ready plot card, feeding help_note() below. Centralizing this
#' prose here (instead of inlining a description + config path at ~24 module call
#' sites) keeps it a single source of truth, per the CLAUDE.md "compute once, reuse
#' downstream" / "single source of truth" rules. Each entry:
#'   - desc:   one-sentence description of what the plot shows
#'   - config: config dot-path(s) that control it, or NULL if none applies
#' @noRd
HELP_NOTES <- list(
    snp_density = list(
        desc   = "SNPs retained per chromosome after all filters.",
        config = NULL
    ),
    pca = list(
        desc   = "Sample scatter on the first two genotype PCs, labeled by sampling site — unsupervised structure overview before any clustering model is fit.",
        config = NULL
    ),
    tracy_widom = list(
        desc   = "Tracy-Widom test statistic per PC — the eigenvalue variance each PC explains, in decreasing order. PCs before the curve flattens capture real structure; flat PCs are noise.",
        config = NULL
    ),
    cross_entropy = list(
        desc   = "sNMF cross-entropy per tested K (lower = better fit) — the standard criterion for choosing the number of ancestral populations.",
        config = "sNMF.k_start, sNMF.k_end, sNMF.repeats"
    ),
    structure_bar = list(
        desc   = "sNMF ancestry proportions per sample at the selected K, one bar per sample.",
        config = "sNMF.k_best"
    ),
    pca_structure = list(
        desc   = "PCA scatter with each sample's genotype PCs pie-chart colored by its sNMF ancestry proportions at the selected K.",
        config = "sNMF.k_best"
    ),
    pop_diff = list(
        desc   = "sNMF population-differentiation p-values per SNP among K clusters (histogram + genome-wide -log10(p) track).",
        config = "sNMF.k_best"
    ),
    pca_structure_k = list(
        desc   = "PCA colored by cluster at the best K.",
        config = "sNMF.k_best"
    ),
    climate_heatmap = list(
        desc   = "Pairwise correlation among climate predictors — collinearity check before GEA.",
        config = "Climate.predictors"
    ),
    climate_density = list(
        desc   = "Distribution of present-day climate values across sampling sites.",
        config = "Climate.predictors, Map.climate_extent"
    ),
    phenotype_density = list(
        desc   = "Distribution of phenotype/trait values across samples.",
        config = NULL
    ),
    trait_correlogram = list(
        desc   = "Pairwise correlation among phenotypic traits — which traits carry redundant information before they are taken to GWAS. Constant or all-missing traits are dropped automatically.",
        config = NULL
    ),
    trait_correlogram_climate = list(
        desc   = "The same correlogram with the climate predictors added, so a trait's environmental association is visible next to its trait-trait correlations. Built only when the project has climate; the two blocks are joined on sample, so climate-excluded samples drop out of both.",
        config = "Climate.predictors"
    ),
    trait_pairs = list(
        desc   = "Pairwise scatter (lower), density (diagonal) and correlation coefficient (upper) across traits. Above Traits.pairs_max_factors traits the grid is unreadable, so the pipeline writes a placeholder instead — raise the cap only if you intend to read a grid that large.",
        config = "Traits.pairs_max_factors"
    ),
    mantel = list(
        desc   = "Isolation by distance vs isolation by environment: Mantel and partial-Mantel correlations of the ancestry (Hellinger) distance against geographic and climate distance, squared into a geography/climate/shared/unexplained budget. Computed at SITE level — samples are collapsed to one row per site first, since coordinates and climate are site properties and per-sample rows would let sampling effort inflate both the distances and the permutation p-value. Descriptive rather than a test when there are few sites (the permutation p cannot go below 1/n_sites!). The variance-partition Venn on the Climate tab answers the same question with dbMEM + partial RDA and adjusted R2 — prefer it where both exist.",
        config = "Population.calc_stats, Population.window_size"
    ),
    amova = list(
        desc   = "Variance partition among vs within populations.",
        config = "Population.calc_stats"
    ),
    qq_plot_gea = list(
        desc   = "Observed vs expected p-value quantiles — calibration check for this method/trait.",
        config = "GEA.configs"
    ),
    qq_plot_gwas = list(
        desc   = "Observed vs expected p-value quantiles — calibration check for this method/trait.",
        config = "GWAS.configs"
    ),
    enrich_dotplot = list(
        desc   = "Top enriched GO terms for this region's genes.",
        config = "Enrichment.top_terms, GFF.go_field"
    ),
    enrich_emapplot = list(
        desc   = "Network of related enriched GO terms (needs ≥ 2 terms).",
        config = "Enrichment.top_terms, GFF.go_field"
    ),
    regionplot_img = list(
        desc   = "Zoomed regional Manhattan with gene track.",
        config = "Filter bar → Clumping distance (bp); fixed 100 kb default, tune interactively"
    ),
    hap_clustree_mg = list(
        desc   = "Marker-group cluster stability across scanned epsilon values.",
        config = "haplotype.scan.epsilon_range, .min_group_size, .min_snps"
    ),
    hap_clustree_hap = list(
        desc   = "Haplotype-group cluster stability across epsilon.",
        config = "haplotype.scan.epsilon_range, .min_group_size, .min_snps, .min_haplotype_size"
    ),
    hap_crosshap_viz = list(
        desc   = "Marker-group × haplotype × phenotype panel at the selected epsilon.",
        config = "haplotype.epsilon_selected"
    ),
    hap_boxplot = list(
        desc   = "Phenotype distribution per haplotype group.",
        config = "haplotype.epsilon_selected"
    ),
    hap_piemap = list(
        desc   = "Geographic haplotype-group frequency map.",
        config = "haplotype.epsilon_selected"
    ),
    overall_importance = list(
        desc   = paste0("Variable importance ranking. Gradient Forest: permutation importance. ",
                        "Geometric offset: diag(cov(B)) — the offset produced by a unit shift ",
                        "in each predictor, summed over every eigen-direction including the ",
                        "numerically null ones. When more predictors than distinct site ",
                        "environments are supplied, the ridge inverse inside LEA::lfmm2() ",
                        "amplifies the lowest-variance directions and these bars become an ",
                        "artefact of the predictor block rather than a marker result."),
        config = "Maladaptation.methods.gradient_forest.ntree, .cor_threshold; Climate.predictors"
    ),
    cumulative_importance = list(
        desc   = "Gradient Forest cumulative importance / turnover along each predictor.",
        config = "Maladaptation.methods.gradient_forest.ntree, .cor_threshold"
    ),
    gf_imputation_sensitivity = list(
        desc   = paste0("Whether the imputation decided the answer. Gradient Forest cannot ",
                        "read missing genotypes, so it is always fitted on an imputed matrix. ",
                        "This check refits the SAME SNP set a second way — as site allele ",
                        "frequencies computed from OBSERVED calls only, inventing nothing — ",
                        "and compares the two. Agreement means the offset is a property of ",
                        "the data; disagreement means it is a property of the imputer. Site ",
                        "frequencies are also the canonical Gradient Forest response unit ",
                        "(Fitzpatrick & Keller 2015), so this is the textbook parameterisation ",
                        "run as a control on the individual-level one."),
        states = list(
            pass = list(
                label = "Pass", class = "bg-success",
                desc  = paste0("Site ranking agrees (Spearman ≥ 0.80) and at least 2 of the ",
                               "top 3 climate drivers are shared. The imputed result is ",
                               "reproduced without imputing anything — safe to report.")),
            warning = list(
                label = "Warning", class = "bg-warning text-dark",
                desc  = paste0("Partial agreement (Spearman 0.60–0.80, or only 1 shared top ",
                               "driver). The broad pattern survives but individual site ranks ",
                               "move. Quote the ranking with a caveat, and do not lean on the ",
                               "exact order of adjacent sites.")),
            fail = list(
                label = "Fail", class = "bg-danger",
                desc  = paste0("The two fits disagree (Spearman < 0.60, below the p<0.05 ",
                               "critical value at typical site counts, or no shared top ",
                               "driver). The offset map is substantially an artefact of the ",
                               "imputation and should not be reported as a marker result.")),
            not_run = list(
                label = "Not run", class = "bg-secondary",
                desc  = paste0("The check could not be evaluated — too few sites (<8), too ",
                               "few SNPs surviving observed-only aggregation (<30), or the ",
                               "check is disabled. Absence of a warning here is NOT evidence ",
                               "that the imputation is harmless."))
        ),
        config = paste0("Maladaptation.methods.gradient_forest.sensitivity_check, ",
                        ".freq_min_calls")
    ),
    offset_piemap = list(
        desc   = "Predicted genetic offset under future climate.",
        config = "Future.ssp, .year, .models, Maladaptation.methods.gradient_forest.spatial_correction"
    ),
    structure_piemap = list(
        desc   = "Geographic map of per-population ancestry/climate pie charts.",
        config = "Piemap.alpha, .show_labels, .pie_scale, Map.zoom_extent"
    ),
    manhattan_gea_combined = list(
        desc   = "Genome-wide GEA associations from all methods combined, significant SNPs highlighted over the full SNP cloud.",
        config = "GEA.configs"
    ),
    manhattan_gea_method = list(
        desc   = "Genome-wide GEA associations for the selected method and trait, over the full SNP cloud.",
        config = "GEA.configs"
    ),
    manhattan_gwas_combined = list(
        desc   = "Genome-wide GWAS associations from all methods combined, significant SNPs over the full SNP cloud.",
        config = "GWAS.configs"
    ),
    manhattan_gwas_method = list(
        desc   = "Genome-wide GWAS associations for the selected method and trait, over the full SNP cloud.",
        config = "GWAS.configs"
    ),
    rda_screeplot = list(
        desc   = "Constrained RDA eigenvalues from the PARTIAL (structure-corrected) fit; retained axes (green) vs dropped (red), per-axis anova.cca(by=\"axis\") p annotated above each bar.",
        config = "GEA.configs (RDA method params: axes, axis_alpha, condition_pcs)"
    ),
    rda_pval_hist = list(
        desc   = "rdadapt p-value histogram for the B6-combined candidate rule — max(p_partial, p_unconstrained), the continuous-p-value analogue of intersecting two independent outlier lists (a SNP only stays significant if BOTH the structure-corrected and unconstrained fits call it). GIF lambda computed with genomic.control implicit in the robust-Mahalanobis calibration (not LFMM-style genomic control). Should be flat with a spike near 0; a U-shape or hump signals miscalibration (see docs/rda_research.md A4, B6).",
        config = "GEA.configs (RDA method params)"
    ),
    rda_biplot = list(
        desc   = "SNP loadings on the first two constrained axes of the PARTIAL fit (binned density, not raw points — WGS-scale marker counts), with predictor biplot vectors and candidate SNPs colored by their most-correlated predictor. Candidate status reflects the B6-combined rule (both fits must agree).",
        config = "GEA.configs (RDA method params: predictor_set)"
    ),
    miami_plot = list(
        desc   = "GEA associations (upward) mirrored against GWAS associations (downward) to reveal loci shared by climate and phenotype.",
        config = "GEAxGWAS.pairwise.window_size, .min_snps"
    ),
    pairwise_miami = list(
        desc   = "Miami plot for one selected trait pair — the two traits' associations mirrored to compare locus overlap.",
        config = "GEAxGWAS.pairwise.window_size, .min_snps"
    ),
    phenomap = list(
        desc   = "Geographic map of phenotype values across sampling sites as scaled pie charts.",
        config = "Piemap.alpha, .show_labels, .pie_scale, Map.zoom_extent"
    ),
    compare_novelty = list(
        desc   = "ExDet climate-novelty surface (NT1 univariate, NT2 multivariate) — where future climate leaves the training range and offset extrapolation is least reliable.",
        config = "Future.ssp, .year, .models"
    ),
    compare_disagree = list(
        desc   = "Rank disagreement between the two offset models mapped against climate novelty.",
        config = "Future.ssp, .year, .models"
    ),
    compare_concordance = list(
        desc   = "Per-site scatter of offset ranks from the two models with a linear fit — overall ranking agreement.",
        config = NULL
    ),
    compare_stability = list(
        desc   = "The 20 sites with the largest offset-rank difference between the two models.",
        config = NULL
    ),
    compare_nway = list(
        desc   = "Kendall's W concordance across all offset models — a single 0-1 agreement statistic over site rankings.",
        config = NULL
    ),

    # ── PreGEA (hyperparameter-exploration mode: "how many K / how many PCs") ──
    pregea_screeplot = list(
        desc   = "LD-pruned genotype PCA screeplot (eigenvalue bars + broken-stick null) with the swept LFMM K band marked. Anchors the K range — the structure-axis count is an UPPER BOUND on K, not a target.",
        config = "sNMF.k_best, PreGEA.k_offset"
    ),
    pregea_ladder_table = list(
        desc   = "One row per (K or #PC value, trait). Trait: which predictor, or __pooled__ for all predictors combined. Histogram shape: flat_spike = well-calibrated; u_shape/hump = value too low; depleted = value too high. Flatness (KS): the actual rule used to pick LFMM's K — lower is flatter. Hits (FDR 0.1): the over-correction detector — a steep drop as the value rises means real signal is being absorbed. lambda_GC: EMMAX's selection rule, target near 1.0.",
        config = NULL
    ),
    pregea_lfmm_hist = list(
        desc   = "P-value histogram per swept K — flat is well-calibrated, U-shaped/peaked at 0 signals residual structure or overcorrection.",
        config = "PreGEA.k_offset"
    ),
    pregea_lfmm_qq = list(
        desc   = "Observed vs expected -log10(p) quantiles per swept K — deviation from the diagonal at the same K as the histogram check.",
        config = "PreGEA.k_offset"
    ),
    pregea_lfmm_lambda = list(
        desc   = "Genomic inflation factor (lambda_GC) vs K, fit with genomic_control OFF (fixed, not configurable — see README 'Fixed constants') so lambda is actually informative. Recommended K picks the flattest p-value histogram, NOT the lambda closest to 1 — see the note on the histogram panel.",
        config = NULL
    ),
    pregea_lfmm_hits = list(
        desc   = "Read the TREND, not the height — this is an over-correction detector, not a power measure. A steep DROP in hits as K rises means the model is absorbing real signal (over-correction onset); a stable plateau means the correction is doing its job. On small test datasets (e.g. SIMDATA) this line is often flat at 1 hit across the whole range — the metric only separates values at WGS-scale SNP density, that is not a bug.",
        config = NULL
    ),
    pregea_emmax_hist = list(
        desc   = "P-value histogram per swept #PCs (BN kinship matrix held fixed, same one GAPIT uses) — flat is well-calibrated.",
        config = "PreGEA.n_pcs_max"
    ),
    pregea_emmax_qq = list(
        desc   = "Observed vs expected -log10(p) quantiles per swept #PCs.",
        config = "PreGEA.n_pcs_max"
    ),
    pregea_emmax_lambda = list(
        desc   = "Genomic inflation factor (lambda_GC) vs #PCs — the recommended value is the smallest #PCs landing inside the calibration band [0.90, 1+0.15] (fixed constants, not configurable), i.e. the cheapest covariate count that still calibrates.",
        config = "PreGEA.n_pcs_max"
    ),
    pregea_emmax_hits = list(
        desc   = "Read the TREND, not the height — same over-correction detector as the LFMM panel. A steep DROP in hits as #PCs rises means real signal is being absorbed. On small test datasets (e.g. SIMDATA) this line is often flat — the metric needs WGS-scale SNP density to separate values.",
        config = NULL
    ),
    pregea_rda_model_comparison = list(
        desc   = "R2adj, max VIF, #significant axes, and full-model p across every swept Condition()-PC value (dashed lines: vif_max, axis_alpha) — the same sweep-and-compare pattern as the LFMM-K and EMMAX-#PC ladders, applied to RDA's own structure-correction strength. Which predictors were pre-screened for collinearity is shown on the Climate tab's correlation heatmap badge.",
        config = "PreGEA.n_pcs_max, .Advanced.vif_max, .Advanced.axis_alpha"
    ),
    pregea_rda_model_screeplot = list(
        desc   = "Constrained-axis eigenvalue screeplot for the SELECTED model (pick a Condition()-PC value above) — retained axes (green) vs non-significant (plum).",
        config = "PreGEA.Advanced.axis_alpha"
    ),
    pregea_rda_model_biplot = list(
        desc   = "RDA site-scores biplot for the SELECTED model — sample scores + predictor arrows. Plain scatter is correct here (N=samples, not SNPs).",
        config = "PreGEA.predictors"
    ),
    pregea_rda_model_axis_anova = list(
        desc   = "Per-axis eigenvalue and anova.cca significance for the SELECTED model.",
        config = "PreGEA.Advanced.axis_alpha, .permutations"
    ),
    pregea_rda_ordir2step = list(
        desc   = "ordiR2step forward-selection path (Blanchet double stopping rule) — which predictors entered the model, in order, and the cumulative R2adj at each step vs the full-model ceiling. An empty path means no variable passed the entry criterion, not a missing plot.",
        config = NULL
    ),
    climate_dbmem_screeplot = list(
        desc   = "dbMEM spatial eigenvector screeplot (Moran's I per MEM) — positive-autocorrelation MEMs are kept as the spatial predictor set for varpart and for Gradient Forest's spatial correction.",
        config = "No config — always site-level (one point per site centroid); skipped with a warning below 3 sites."
    ),
    climate_dbmem_selection_path = list(
        desc   = "Forward-selected MEM path (Blanchet double stopping rule) — which spatial eigenvectors entered the varpart geography term, in order.",
        config = NULL
    ),
    climate_varpart_venn = list(
        desc   = "Schematic Venn (fixed circle positions — NOT area-proportional) of the genomic-PC response's variance: one circle per available factor (Climate, Structure, Geography), each region labeled with its adjusted-R2 share. 3 circles when both auxiliary tables are available, falling back to 2 (climate+geography, or climate+structure) or 1 (climate only) otherwise. \"Unexplained\" is deliberately NOT drawn — standard Venn convention, it's everything outside the circles — its size is instead reported by the \"% variance explained\" badge above the plot. A region's % label can be NEGATIVE — adjusted R2 corrects for the number of predictors used and overcorrects when a fraction explains ~nothing at low N; read a negative value as ≈0 (Peres-Neto et al. 2006), not an error. The separate climate/geography confounding badge is its own dedicated check, independent of which Venn variant is shown.",
        config = "Climate.Varpart.response, .structure_table (confounding flag is always computed; dbMEM geography has no config — always site-level)"
    ),
    climate_px_barplot = list(
        desc   = "Lasky et al. 2012 Px metric per climate predictor — each predictor's share of the RDA-explained variance, ranked, for the unconditional and geography-partialled models.",
        config = "Climate.Varpart.response"
    ),
    pregea_transfer_guard = list(
        desc   = "LFMM/EMMAX lambda re-estimated on the FULL marker set at the recommended hyperparameters, compared against the LD-pruned ladder value. Only the TREND is expected to transfer — lambda is not scale-free (Yang et al. 2011), so a large jump does not necessarily mean the recommended hyperparameter is wrong.",
        config = "PreGEA.TransferGuard.enabled, .lfmm_k, .emmax_n_pcs"
    )
)

#' Build a neutral hover-note badge for a plot: what it shows + results + config link.
#'
#' Icon-only trigger (no visible label) — this is always-present reference info,
#' not a status alert, so it stays visually quiet next to the card title. Reveals
#' a structured tooltip body on hover via filter_note() (R/utils_ui.R):
#'   - "What it shows": entry$desc (always present)
#'   - "Results": caller-supplied, data-derived from tables already on disk
#'     (e.g. best K, N samples, pops-per-K) — omitted when NULL, never invented
#'   - "Config": entry$config, the dot-path(s) that control the plot — omitted
#'     when NULL
#'
#' Returns NULL silently for an unknown id — a missing help-note should never
#' break a card's rendering.
#'
#' @param id key into HELP_NOTES
#' @param results optional data-derived results line (string or htmltools tag)
#' @param extra optional extra line appended to the tooltip body
#' @param label optional visible badge label (e.g. "K = 3") — same convention as
#'   the MAF/missingness filter badges, which show the live value on the badge
#'   itself, not just in the tooltip. Pass a value the caller already has —
#'   don't recompute it here.
#' @noRd
help_note <- function(id, results = NULL, extra = NULL, label = "",
                      class = "bg-secondary", states = FALSE) {
    entry <- HELP_NOTES[[id]]
    if (is.null(entry)) return(NULL)

    body <- htmltools::tagList(
        htmltools::p(htmltools::strong("What it shows: "), entry$desc),
        if (!is.null(results))
            htmltools::p(htmltools::strong("Results: "), results),
        # `states` renders an entry's per-status legend (green/amber/red), so a status badge
        # explains its own colour scale in the same tooltip that reports the result.
        if (isTRUE(states) && !is.null(entry$states))
            htmltools::div(
                class = "small mt-1",
                lapply(names(entry$states), function(s) {
                    st <- entry$states[[s]]
                    htmltools::p(
                        class = "mb-1",
                        htmltools::tags$span(class = paste("badge me-1", st$class), st$label),
                        st$desc
                    )
                })
            ),
        if (!is.null(entry$config))
            htmltools::p(class = "text-muted small",
                htmltools::strong("Config: "), htmltools::code(entry$config)),
        if (!is.null(extra)) htmltools::p(extra)
    )

    filter_note(label = label, body = body, class = class)
}
