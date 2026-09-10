# MODULE: CLIMATE — predictor characterization (mode=climate)
#
# Owns the full predictor story: correlation/density plots, invariant-
# predictor detection, and spatial structure (dbMEM + variance partitioning).
# Split out of the old preGEA module (which mixed predictor characterization
# with LFMM-K/EMMAX-#PC/RDA hyperparameter sweeps) so each module answers one
# question: Climate answers "what do my predictors look like and how much do
# they overlap with geography/structure", PreGEA answers "how many K / #PCs".
#
# Climate DATA ACQUISITION (WorldClim download / custom-table staging) stays
# in structure.smk — climate_site is an input to Structure's own piemaps, and
# the boundary here is "predictor characterization", not "data acquisition".
#
# Two blocks:
#   Block 1 — predictor characterization (moved from structure.smk unchanged)
#   Block 2 — dbMEM (shared spatial artifact) + varpart (moved from pregea.smk,
#             renamed climate_dbmem/climate_varpart)

# =============================================================================
# Block 1 — predictor characterization
# =============================================================================

rule check_climate_variance:
    """Detect invariant (zero-variance or all-NA) bioclimatic predictors at sample sites."""
    input:  site = O['climate_site']
    output: O['climate_invariant']
    params: columns = "bio"
    log:    f"{LOGDIR}climate/check_climate_variance.log"
    shell:
        """
        Rscript /pipeline/scripts/check_climate_variance.R \
            {input.site} {output} {params.columns} > {log} 2>&1
        """

rule design_adequacy:
    """What the SAMPLING DESIGN can support, computed before any GEA/offset fit.

    Reports the environmental degrees of freedom (n_sites - 1, NOT n_samples - 1),
    per-site balance, and the conditioning of the site-level predictor block
    (numerical rank, participation-ratio effective dimensionality, condition
    number, PC1 share). These bound varpart, RDA's constrained rank and every
    Mahalanobis-type genomic offset.

    The d.f./balance/pseudoreplication half is computed nowhere else. Rank and
    condition number ARE already reported downstream — geometric_offset.R writes
    env_cov_rank_numeric / env_cov_condition_number post-hoc, and
    pregea_rda_setup.R screens site-level max|r| per rung — but only after a fit,
    only when those optional modes run, and on the COVARIANCE matrix. This rule
    reports them in mode=climate, before anything is fitted, on the CORRELATION
    matrix; see scripts/design_adequacy.py for why the two do not agree.

    Warn-only by design — every flag here is a judgement call about the study,
    not a code error, so it is recorded and surfaced rather than used to abort."""
    input:
        meta    = O['metadata'],
        climate = O['climate_site'],
    output: O['climate_design']
    params: predictors = ",".join(get_predictors_list())
    log:    f"{LOGDIR}climate/design_adequacy.log"
    shell:
        """
        python3 /pipeline/scripts/design_adequacy.py \
            {input.meta} {input.climate} {params.predictors} {output} > {log} 2>&1
        """

rule density_plot:
    """Generate combined density plot for all climate predictors (all BIO columns)."""
    input:  climate = O['climate_site']
    output: DENSITY_PLOT_COMBINED
    params:
        predictors = "all",
        inter_dir = INTER
    log:    f"{LOGDIR}climate/density_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

# NOTE: the phenotype density plot used to live here (density_plot_phenotypes).
# It moved to traits.smk in Phase 3 — it describes traits, not climate, and had
# to survive Climate.enabled: false. Same for the trait columns that were
# cbind()ed into the correlogram below: this heatmap is climate-only now, and
# the traits x climate version is traits.smk's trait_correlogram_climate.

rule correlation_heatmap:
    """Correlogram of the climate predictors alone (Phase 3: traits dropped —
    they have their own correlogram, plus a joint traits x climate one, in
    mode=traits). Second table passed as NULL."""
    input:  climate = O['climate_site']
    output: O['corr_heatmap']
    params: inter_dir = INTER, second = 'NULL', title = "Correlogram of climate factors"
    log:    f"{LOGDIR}climate/correlation_heatmap.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_correlation_heatmap.R \
            {input.climate} {params.second} {output} {params.inter_dir} \
            '{params.title}' > {log} 2>&1
        """

# =============================================================================
# Block 2 — dbMEM (shared spatial artifact) + varpart. Split into two rules
# specifically so the shared artifact is cheap and standalone: climate_dbmem
# depends only on coordinates, so Gradient Forest's 'spatial' variant can
# request it in maladaptation mode without dragging any genotype work into
# the DAG.
# =============================================================================

rule climate_dbmem:
    """SHARED spatial-vector artifact — the single dbMEM computation reused
    by varpart AND Gradient Forest (A20/C4). Replaces gradient_forest_model.R's
    old ad hoc raw-Euclidean-dist + "keep half the positive eigenvalues"
    heuristic — see scripts/pregea_dbmem.R header."""
    input:
        metadata = W['metadata_climate_valid'],
        samples  = W['climate_valid_samples'],
    output:
        vectors = O['climate_dbmem_vectors'], diag = O['climate_dbmem_diag'],
        png = O['climate_dbmem_png'], svg = O['climate_dbmem_svg'],
    params: inter_dir = INTER
    log:    f"{LOGDIR}climate/dbmem.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_dbmem.R \
            {input.metadata} {input.samples} \
            {output.vectors} {output.diag} {output.png} {params.inter_dir} > {log} 2>&1
        """

rule climate_varpart:
    """Shared climate/structure/geography variance-partitioning panel.
    scripts/pregea_varpart.R — see its header for the full port-and-correction
    notes from rda_varpart_lasky.Rmd."""
    input:
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues'],
        dbmem       = O['climate_dbmem_vectors'],
        clusters    = clusters_table(K_BEST) if CLIMATE_VP_STRUCT == 'qmatrix' else [],
        climate     = O['climate_site_scaled'],
        climate_valid = W['climate_valid_samples'],
        samples_order = W['samples_order'],
        lfmm_pruned = W['lfmm_imp_climate'] if CLIMATE_VP_RESPONSE == 'snps' else [],
    output:
        selection  = O['climate_vp_selection'], selected = O['climate_vp_selected'],
        table      = O['climate_vp_table'], confound = O['climate_vp_confound'],
        px         = O['climate_vp_px'],
        p_path = O['climate_vp_path_png'], p_venn = O['climate_vp_venn_png'],
        p_px   = O['climate_vp_px_png'],
    params:
        clusters_arg = lambda wc, input: input.clusters if CLIMATE_VP_STRUCT == 'qmatrix' else 'NULL',
        lfmm_arg     = lambda wc, input: input.lfmm_pruned if CLIMATE_VP_RESPONSE == 'snps' else 'NULL',
        predictors = PREDICTORS_SELECTED,
        response = CLIMATE_VP_RESPONSE, response_var = CLIMATE_VP_VAR_CUTOFF,
        response_max = CLIMATE_VP_MAX_PCS, response_min = CLIMATE_VP_MIN_PCS,
        pin = CLIMATE_VP_PIN, r2perms = CLIMATE_VP_R2PERMS,
        varpartperms = CLIMATE_VP_PERMS,
        seed = PREGEA_SEED, plot_dir = f"{MOD_CLIMATE}plots/varpart/", inter_dir = INTER,
    log: f"{LOGDIR}climate/varpart.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_varpart.R \
            {input.projections} {input.eigenvalues} {input.dbmem} {params.clusters_arg} \
            {input.climate} {input.climate_valid} {input.samples_order} {params.predictors} \
            {params.response} {params.response_var} {params.response_max} {params.response_min} \
            {params.lfmm_arg} {params.pin} {params.r2perms} {params.varpartperms} \
            {params.seed} \
            {output.selection} {output.selected} {output.table} {output.confound} {output.px} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """
