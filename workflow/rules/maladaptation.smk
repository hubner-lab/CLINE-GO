#=============================================================================
# MODULE 6: MALADAPTATION
#=============================================================================

# Wildcard constraints for all Gradient Forest rules.
# run_label = set name (matches Shiny charset [A-Za-z0-9_.-]+, no slashes).
# spatial_tag = exactly 'spatial' or 'nospatial' — disambiguates trailing token
#   even when run_label itself contains underscores (e.g. 'EMMAX_bonf005').
# scenario = one future-climate projection (see SCENARIOS in common.smk). Names
#   come from future-table basenames, so the charset matches a filename stem.
wildcard_constraints:
    run_label   = r"[A-Za-z0-9_.-]+",
    spatial_tag = r"spatial|nospatial",
    scenario    = r"[A-Za-z0-9_.+-]+",
    method      = r"gradient_forest|geometric_offset|rda_offset|rda_offset_corrected"

if CLIMATE_SOURCE == 'worldclim':
    # Per-model future climate download (runs in parallel via Snakemake)
    rule download_climate_future_model:
        """Download CMIP6 future climate data for a single model.
        Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples)."""
        input: samples = W['metadata_climate']
        output: f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{{model}}.tif"
        wildcard_constraints: model = r"[A-Za-z0-9_-]+"
        params:
            crop = CLIMATE_EXTENT,
            gap = GAP,
            resolution = RESOLUTION,
            data_dir = INDIR
        log: f"{LOGDIR}maladaptation/download_climate_future_{{model}}.log"
        shell:
            """
            Rscript /pipeline/scripts/download_climate_future_model.R \
                {input.samples} {params.crop} {params.gap} {params.data_dir} \
                {SSP} {YEAR} {wildcards.model} {params.resolution} \
                {output} > {log} 2>&1
            """

    # Merge per-model rasters into averaged future climate
    rule merge_climate_future:
        """Average future climate across models and extract site values.
        Uses metadata_climate_valid.tsv (climate-valid samples -- already narrowed by
        filter_climate_valid_samples to exclude any present-climate NA/ocean-pixel samples;
        see download_climate_present). Still runs its own NA check as defense-in-depth against
        a residual future-only NA (e.g. a different NoData mask in a CMIP6 model raster)."""
        input:
            samples = W['metadata_climate_valid'],
            model_rasters = [f"{MOD_CLIMATE}rasters/future/climate_future_year{YEAR}_ssp{SSP}_{model}.tif" for model in MODELS_LIST],
            present_raster = W['climate_raster'],
            present_all = O['climate_all']
        output:
            raster = climate_future_raster(DEFAULT_SCENARIO),
            all_vals = climate_future_all(DEFAULT_SCENARIO),
            site_vals = climate_future_site(DEFAULT_SCENARIO),
            na_excluded = O['climate_future_na_excluded']
        params:
            raster_str = lambda wc, input: ','.join(input.model_rasters),
            n_models = len(MODELS_LIST)
        log: f"{LOGDIR}maladaptation/merge_climate_future.log"
        shell:
            """
            Rscript /pipeline/scripts/merge_climate_future.R \
                {input.samples} {params.raster_str} {params.n_models} \
                {input.present_raster} {input.present_all} \
                {output.raster} {output.all_vals} {output.site_vals} \
                {output.na_excluded} > {log} 2>&1
            """
else:
    rule stage_custom_climate_future:
        """Stage one user-supplied future climate table into pipeline-standard outputs.
        Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples).

        env_table is an INPUT, not a param: as a param Snakemake could not see the
        table change, so editing a future scenario left stale offsets in place and
        reported "Nothing to be done"."""
        input:
            samples        = W['metadata_climate'],
            present_all    = O['climate_all'],
            present_raster = W['climate_raster'],
            env_table      = lambda wc: SCENARIOS[wc.scenario]
        output:
            raster    = climate_future_raster('{scenario}'),
            all_vals  = climate_future_all('{scenario}'),
            site_vals = climate_future_site('{scenario}')
        params:
            columns   = CUSTOM_CLIMATE_COLUMNS,
            key       = CUSTOM_CLIMATE_KEY
        log: f"{LOGDIR}maladaptation/stage_custom_climate_future_{{scenario}}.log"
        shell:
            """
            Rscript /pipeline/scripts/stage_custom_climate.R future \
                {input.samples} {input.env_table} {params.columns} {params.key} \
                {input.present_raster} {input.present_all} \
                {output.raster} {output.all_vals} {output.site_vals} > {log} 2>&1
            """

rule density_plot_future:
    """Generate combined density plot for future climate predictors."""
    input: climate = O['climate_future_site']
    output: O['density_future']
    params:
        predictors = PREDICTORS_SELECTED,
        inter_dir = INTER
    log: f"{LOGDIR}maladaptation/density_plot_future.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.climate} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

# Gradient Forest - adaptive model
# sigsnps is an ancestor-less source file produced by the Shiny GEA tab.
# {run_label} = saved SNP-set name; {spatial_tag} = spatial|nospatial.
# pcnm param translates spatial_tag -> 'with'/'without' for the R script.
# dbmem/dbmem_selected: the climate module's forward-selected dbMEM vectors
# (climate.smk Block 2, adespatial::dbmem + ordiR2step) — only pulled in when
# spatial_tag == 'spatial'. Varpart is unconditional inside mode=climate (no
# more enabled switch to gate on), so these paths always resolve to real
# climate outputs whenever this branch of the DAG is actually reachable —
# the user just needs to have run mode=climate first.
rule gradient_forest_adaptive:
    """Build adaptive Gradient Forest model using a user-curated SNP set.

    lfmm is the IMPUTED matrix, matching geometric_offset/rda_offset and the GEA rules.
    [FIXED 2026-08-10] This rule and gradient_forest_random previously took
    W['lfmm_full_climate'] -- the raw matrix, in which a missing call is the literal value 9.
    gradientForest() takes the SNP columns as RESPONSE variables and has no na.action, and
    gradient_forest_model.R never converted 9 to NA, so 9 entered the regression forest as a
    genotype 3.5x further from the reference homozygote than the real alternate homozygote.
    On Trifolium86FinalPC3 (45.05% missing) the fit retained 18 of 185 SNPs, and those 18 were
    selected by the between-site variance of their MISSINGNESS (Wilcoxon p = 6.6e-11), not by
    allele frequency (p = 0.36). Refitting on the imputed matrix retains 184 of 185 and gives an
    unrelated map (rho = -0.064) and an unrelated importance profile (rho = -0.044).
    See docs/pipeline_improvement_requests.md bug O."""
    input:
        lfmm    = W['lfmm_imp_full_climate'],
        sigsnps = lambda wc: snp_set_file(wc.run_label),
        vcfsnp  = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = W['metadata_climate_valid'],
        climate = O['climate_site'],
        dbmem          = lambda wc: O['climate_dbmem_vectors'] if wc.spatial_tag == 'spatial' else [],
        dbmem_selected = lambda wc: O['climate_vp_selected']   if wc.spatial_tag == 'spatial' else []
    output: mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive')
    params:
        predictors    = PREDICTORS_SELECTED,
        ntree         = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm          = lambda wc: 'with' if wc.spatial_tag == 'spatial' else 'without',
        dbmem_arg     = lambda wc, input: input.dbmem          if wc.spatial_tag == 'spatial' else 'NULL',
        dbmem_sel_arg = lambda wc, input: input.dbmem_selected if wc.spatial_tag == 'spatial' else 'NULL'
    log: f"{LOGDIR}maladaptation/gradient_forest_adaptive_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            adaptive {output} {params.dbmem_arg} {params.dbmem_sel_arg} > {log} 2>&1
        """

# Gradient Forest - random/neutral model (optional)
rule gradient_forest_random:
    """Build neutral Gradient Forest model using random SNPs.

    Imputed matrix, same fix as gradient_forest_adaptive -- see the note there and bug O."""
    input:
        lfmm    = W['lfmm_imp_full_climate'],
        sigsnps = lambda wc: snp_set_file(wc.run_label),
        vcfsnp  = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = W['metadata_climate_valid'],
        climate = O['climate_site'],
        dbmem          = lambda wc: O['climate_dbmem_vectors'] if wc.spatial_tag == 'spatial' else [],
        dbmem_selected = lambda wc: O['climate_vp_selected']   if wc.spatial_tag == 'spatial' else []
    output: mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random')
    params:
        predictors    = PREDICTORS_SELECTED,
        ntree         = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm          = lambda wc: 'with' if wc.spatial_tag == 'spatial' else 'without',
        dbmem_arg     = lambda wc, input: input.dbmem          if wc.spatial_tag == 'spatial' else 'NULL',
        dbmem_sel_arg = lambda wc, input: input.dbmem_selected if wc.spatial_tag == 'spatial' else 'NULL'
    log: f"{LOGDIR}maladaptation/gradient_forest_random_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            random {output} {params.dbmem_arg} {params.dbmem_sel_arg} > {log} 2>&1
        """

# Gradient Forest - site-allele-frequency model (imputation-sensitivity control)
rule gradient_forest_frequency:
    """Refit the adaptive SNP set on site allele frequencies from OBSERVED calls only.

    The ONLY difference from gradient_forest_adaptive is the genotype matrix: this takes
    W['lfmm_full_climate'] (non-imputed, missing = LEA's sentinel 9) instead of the imputed
    one, and gradient_forest_model.R aggregates it to one row per site using
    mean(x[x != 9]) / 2 -- so no genotype is ever invented. The plumbing already exists:
    subset_lfmm_climate (processing.smk) is generic over its `variant` wildcard and already
    accepts 'lfmm_full'; it simply had no consumer until now.

    Site frequencies are the canonical GF response unit (Fitzpatrick & Keller 2015,
    Ecol Lett 18:1-16, doi:10.1111/ele.12376), so this is the textbook parameterisation run
    as a control on the imputed one -- not an exotic variant. Individual-level fitting is
    NOT wrong (Lind & Lotterhos 2025, Am Nat 207(3), doi:10.1086/739098, find genotype-based
    models can beat frequency-based ones off strict clines); the point of running both is
    that agreement between them is evidence the imputation did not decide the answer."""
    input:
        lfmm    = W['lfmm_full_climate'],
        sigsnps = lambda wc: snp_set_file(wc.run_label),
        vcfsnp  = W['vcfsnp_full'],
        removed = W['removed_full'],
        samples = W['metadata_climate_valid'],
        climate = O['climate_site'],
        dbmem          = lambda wc: O['climate_dbmem_vectors'] if wc.spatial_tag == 'spatial' else [],
        dbmem_selected = lambda wc: O['climate_vp_selected']   if wc.spatial_tag == 'spatial' else []
    output: mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'frequency')
    params:
        predictors    = PREDICTORS_SELECTED,
        ntree         = NTREE,
        cor_threshold = COR_THRESHOLD,
        pcnm          = lambda wc: 'with' if wc.spatial_tag == 'spatial' else 'without',
        dbmem_arg     = lambda wc, input: input.dbmem          if wc.spatial_tag == 'spatial' else 'NULL',
        dbmem_sel_arg = lambda wc, input: input.dbmem_selected if wc.spatial_tag == 'spatial' else 'NULL',
        min_calls     = GF_FREQ_MIN_CALLS
    log: f"{LOGDIR}maladaptation/gradient_forest_frequency_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_model.R \
            {input.lfmm} {input.sigsnps} {input.vcfsnp} {input.removed} \
            {input.samples} {input.climate} {params.predictors} \
            {params.ntree} {params.cor_threshold} {params.pcnm} \
            adaptive {output} {params.dbmem_arg} {params.dbmem_sel_arg} \
            site_frequency {params.min_calls} > {log} 2>&1
        """

# Imputation-sensitivity check — compares the two fits above.
rule gf_sensitivity_check:
    """Compare the imputed individual-level GF with the observed-only site-frequency GF.

    Writes a long key/value TSV (same shape as geometric_offset_diagnostics.tsv) carrying a
    status of pass / warning / fail / not_run, which the Shiny maladaptation tab renders as a
    traffic-light badge. This is a TABLE, not a plot, so it is never gated on
    Maladaptation.emit_plots -- the status is the deliverable."""
    input:
        imputed      = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        frequency    = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'frequency'),
        clim_present = O['climate_site'],
        clim_future  = climate_future_site(DEFAULT_SCENARIO)
    output: mala_sensitivity('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        predictors = PREDICTORS_SELECTED,
        samples    = W['metadata_climate_valid']
    log: f"{LOGDIR}maladaptation/gf_sensitivity_check_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gf_sensitivity_check.R \
            {input.imputed} {input.frequency} {input.clim_present} {input.clim_future} \
            {params.predictors} {params.samples} {output} > {log} 2>&1
        """

# Genetic offset calculation — all scenarios in one job.
# Multi-output by design: the model load and the present-climate prediction are
# scenario-invariant and dominate runtime on small SNP panels, so one job projects
# N futures instead of N jobs each re-loading the same model.
rule gradient_forest_offset:
    """Project one Gradient Forest model onto every future scenario."""
    input:
        gf             = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        future_all     = [climate_future_all(s) for s in SCENARIO_NAMES],
        present_all    = O['climate_all'],
        present_raster = W['climate_raster'],
        samples        = W['metadata_climate_valid']
    output:
        raster      = expand(mala_offset_raster('gradient_forest', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        map_values  = expand(mala_offset_map_values('gradient_forest', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        site_values = expand(mala_offset_site_values('gradient_forest', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES)
    params:
        predictors = PREDICTORS_SELECTED,
        scenarios  = ','.join(SCENARIO_NAMES),
        future_arg = lambda wc, input:  ','.join(input.future_all),
        raster_arg = lambda wc, output: ','.join(output.raster),
        map_arg    = lambda wc, output: ','.join(output.map_values),
        site_arg   = lambda wc, output: ','.join(output.site_values)
    log: f"{LOGDIR}maladaptation/gradient_forest_offset_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/gradient_forest_offset.R \
            {input.gf} {params.predictors} {params.future_arg} {input.present_all} \
            {input.present_raster} {input.samples} \
            {params.raster_arg} {params.map_arg} {params.site_arg} \
            {params.scenarios} > {log} 2>&1
        """

# Geometric Genetic Offset (Gain et al. 2023, MBE) — single-call rule.
# No separate model artifact: genetic.gap() fits LFMM2 + computes offset in one pass.
# Nospatial-only: genetic.gap() has no spatial-correction argument.
# candidate.loci = INTEGER INDEX vector into the full imputed matrix (never subset the matrix).
rule geometric_offset:
    """Compute geometric genetic offset using LEA::genetic.gap()."""
    input:
        lfmm_full   = W['lfmm_imp_full_climate'],
        vcfsnp      = W['vcfsnp_full'],
        removed     = W['removed_full'],
        sigsnps     = lambda wc: snp_set_file(wc.run_label),
        env_site_pres  = O['climate_site'],
        env_site_fut   = [climate_future_site(s) for s in SCENARIO_NAMES],
        env_all_pres   = O['climate_all'],
        env_all_fut    = [climate_future_all(s) for s in SCENARIO_NAMES],
        pres_raster    = W['climate_raster'],
        samples        = W['metadata_climate_valid']
    output:
        site_values = expand(mala_offset_site_values('geometric_offset', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        map_values  = expand(mala_offset_map_values('geometric_offset', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        raster      = expand(mala_offset_raster('geometric_offset', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        # Scenario-free: both describe the FIT, not the projection, so they are written once
        # from the first scenario and are NOT expanded over SCENARIO_NAMES.
        importance  = mala_importance('geometric_offset', '{run_label}', '{spatial_tag}'),
        diagnostics = f"{mala_table_dir('geometric_offset', '{run_label}', '{spatial_tag}')}geometric_offset_diagnostics.tsv"
    wildcard_constraints:
        spatial_tag = r"nospatial"   # geometric_offset is nospatial-only
    params:
        predictors = PREDICTORS_SELECTED,
        k          = lambda wc: GO_K if GO_K != '' else K_BEST,
        scale      = GO_SCALE,
        scenarios  = ','.join(SCENARIO_NAMES),
        fut_site_arg = lambda wc, input:  ','.join(input.env_site_fut),
        fut_all_arg  = lambda wc, input:  ','.join(input.env_all_fut),
        site_arg     = lambda wc, output: ','.join(output.site_values),
        map_arg      = lambda wc, output: ','.join(output.map_values),
        raster_arg   = lambda wc, output: ','.join(output.raster)
    log: f"{LOGDIR}maladaptation/geometric_offset_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/geometric_offset.R \
            {input.lfmm_full} {input.vcfsnp} {input.removed} {input.sigsnps} \
            {input.env_site_pres} {params.fut_site_arg} \
            {input.env_all_pres} {params.fut_all_arg} \
            {input.pres_raster} {input.samples} \
            {params.predictors} {params.k} {params.scale} \
            {params.site_arg} {params.map_arg} \
            {params.raster_arg} {output.importance} \
            {params.scenarios} {output.diagnostics} > {log} 2>&1
        """

# RDA Genetic Offset (Capblancq & Forester 2021) — single-call rule.
# Second, independent RDA fit vs. the GEA-scan RDA (scripts/rda.R): candidate
# SNPs only (the curated {run_label} set), no Condition() by default (B7).
# "Compute once, reuse downstream" does NOT apply here — see docs/rda_research.md
# Part B.0. Nospatial-only: this construction never carries structure covariates.
rule rda_offset:
    """Compute RDA genetic offset (second RDA fit on candidate SNPs, no structure
    conditioning by default) — see scripts/rda_offset.R and docs/rda_research.md Part B."""
    input:
        lfmm_full     = W['lfmm_imp_full_climate'],
        vcfsnp        = W['vcfsnp_full'],
        removed       = W['removed_full'],
        sigsnps       = lambda wc: snp_set_file(wc.run_label),
        env_site_pres = O['climate_site'],
        env_site_fut  = [climate_future_site(s) for s in SCENARIO_NAMES],
        env_all_pres  = O['climate_all'],
        env_all_fut   = [climate_future_all(s) for s in SCENARIO_NAMES],
        pres_raster   = W['climate_raster'],
        samples       = W['metadata_climate_valid'],
        pca           = W['pca_projections'],
        samples_order = W['samples_order'],
        climate_valid = W['climate_valid_samples']
    output:
        site_values = expand(mala_offset_site_values('{{rda_method}}', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        map_values  = expand(mala_offset_map_values('{{rda_method}}', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        raster      = expand(mala_offset_raster('{{rda_method}}', '{{run_label}}', '{{spatial_tag}}', '{scenario}'), scenario=SCENARIO_NAMES),
        # Scenario-free: the RDA is fitted on present climate; the future only
        # enters the projection, so importance/diagnostics/screeplot are written once.
        importance  = mala_importance('{rda_method}', '{run_label}', '{spatial_tag}'),
        diagnostics = f"{mala_table_dir('{rda_method}', '{run_label}', '{spatial_tag}')}rda_offset_diagnostics.tsv",
        screeplot     = f"{mala_plot_dir('{rda_method}', '{run_label}', '{spatial_tag}')}rda_axis_screeplot.png",
        screeplot_svg = f"{mala_plot_dir('{rda_method}', '{run_label}', '{spatial_tag}')}rda_axis_screeplot.svg"
    wildcard_constraints:
        spatial_tag = r"nospatial",  # rda_offset is nospatial-only (B7 — no Condition() by default)
        # One rule serves both registered RDA methods; they differ only in params.
        rda_method  = r"rda_offset|rda_offset_corrected"
    params:
        predictors   = PREDICTORS_SELECTED,
        axes         = lambda wc: rdo_cfg(wc.rda_method)['axes'],
        axis_alpha   = lambda wc: rdo_cfg(wc.rda_method)['axis_alpha'],
        permutations = lambda wc: rdo_cfg(wc.rda_method)['permutations'],
        condition_pcs = lambda wc: rdo_cfg(wc.rda_method)['condition_pcs'],
        seed         = lambda wc: rdo_cfg(wc.rda_method)['seed'],
        plot_dir     = lambda wc: mala_plot_dir(wc.rda_method, wc.run_label, wc.spatial_tag),
        scenarios    = ','.join(SCENARIO_NAMES),
        fut_site_arg = lambda wc, input:  ','.join(input.env_site_fut),
        fut_all_arg  = lambda wc, input:  ','.join(input.env_all_fut),
        site_arg     = lambda wc, output: ','.join(output.site_values),
        map_arg      = lambda wc, output: ','.join(output.map_values),
        raster_arg   = lambda wc, output: ','.join(output.raster)
    threads: CPU
    log: f"{LOGDIR}maladaptation/{{rda_method}}_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/rda_offset.R \
            {input.lfmm_full} {input.vcfsnp} {input.removed} {input.sigsnps} \
            {input.env_site_pres} {params.fut_site_arg} \
            {input.env_all_pres} {params.fut_all_arg} \
            {input.pres_raster} {input.samples} \
            {params.predictors} {params.axes} {params.axis_alpha} \
            {params.permutations} {params.condition_pcs} \
            {input.pca} {input.samples_order} {input.climate_valid} \
            {params.seed} {threads} \
            {params.site_arg} {params.map_arg} \
            {params.raster_arg} {output.importance} \
            {output.diagnostics} {params.plot_dir} \
            {params.scenarios} > {log} 2>&1
        """

# Cumulative importance plot
rule plot_gf_cumimp:
    """Plot cumulative importance curves for adaptive (and optionally neutral) GF model."""
    input:
        gf        = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        gf_random = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random') if GF_RANDOM_MODEL else []
    output: mala_cumimp('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        gf_random_path = lambda wc: mala_model('gradient_forest', wc.run_label, wc.spatial_tag, 'random') if GF_RANDOM_MODEL else 'NULL',
        predictors     = PREDICTORS_SELECTED,
        inter_dir      = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_cumimp_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_cumimp.R \
            {input.gf} {params.gf_random_path} {params.predictors} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# Overall importance plot
rule plot_gf_importance:
    """Plot R2-weighted importance for adaptive (and optionally neutral) GF model."""
    input:
        gf        = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'adaptive'),
        gf_random = mala_model('gradient_forest', '{run_label}', '{spatial_tag}', 'random') if GF_RANDOM_MODEL else []
    output: mala_importance('gradient_forest', '{run_label}', '{spatial_tag}')
    params:
        gf_random_path = lambda wc: mala_model('gradient_forest', wc.run_label, wc.spatial_tag, 'random') if GF_RANDOM_MODEL else 'NULL',
        inter_dir      = INTER
    log: f"{LOGDIR}maladaptation/plot_gf_importance_{{run_label}}_{{spatial_tag}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_gf_importance.R \
            {input.gf} {params.gf_random_path} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# Genetic offset PieMaps (collapsed: notrait / tajima_d / pi_diversity)
def _piemap_trait_input(wc):
    """Return trait file path for size-scaled piemaps, or [] for notrait."""
    if wc.size_trait == 'tajima_d':
        return O['tajima']
    elif wc.size_trait == 'pi_diversity':
        return O['pi_div']
    return []

def _piemap_trait_path(wc):
    """Return trait file arg for plot_piemap.R shell command."""
    if wc.size_trait == 'tajima_d':
        return f"{O['tajima']} \"Tajima's D\""
    elif wc.size_trait == 'pi_diversity':
        return f"{O['pi_div']} \"Pi Diversity\""
    return "NULL NULL"

def _piemap_output_prefix(wc):
    """Return OUTPUT_PREFIX arg for plot_piemap.R."""
    return f'genetic_offset_piemap_{wc.size_trait}'

rule plot_gf_offset_piemap:
    """Plot genetic offset piemap for any maladaptation method, optionally scaled by a population statistic."""
    wildcard_constraints: size_trait = "notrait|tajima_d|pi_diversity"
    input:
        offset_raster = lambda wc: mala_offset_raster(wc.method, wc.run_label, wc.spatial_tag, wc.scenario),
        samples       = W['metadata_climate_valid'],
        clusters      = clusters_table(K_BEST),
        trait_file    = _piemap_trait_input
    output: mala_offset_piemap('{method}', '{run_label}', '{spatial_tag}', '{size_trait}', '{scenario}')
    params:
        pie_alpha       = PIEMAP_ALPHA,
        pop_label       = PIEMAP_SHOW_LABELS,
        pop_label_size  = PIEMAP_LABEL_SIZE,
        pie_scale       = PIEMAP_PIE_SCALE,
        plot_dir        = lambda wc: mala_offset_plot_dir(wc.method, wc.run_label, wc.spatial_tag, wc.scenario),
        inter_dir       = INTER,
        regionmap_extent = REGIONMAP_EXTENT,
        trait_path      = _piemap_trait_path,
        output_prefix   = _piemap_output_prefix
    log: f"{LOGDIR}maladaptation/plot_offset_piemap_{{method}}_{{run_label}}_{{spatial_tag}}_{{scenario}}_{{size_trait}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            {params.trait_path} \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            {params.output_prefix} {params.regionmap_extent} {params.pie_scale} Clusters > {log} 2>&1
        """

rule plot_gf_offset_piemap_points:
    """Points ('clear map') companion for the genetic offset piemap — trait-independent
    (points only mark sample locations), one per method/run_label/spatial_tag."""
    input:
        offset_raster = lambda wc: mala_offset_raster(wc.method, wc.run_label, wc.spatial_tag, wc.scenario),
        samples       = W['metadata_climate_valid'],
        clusters      = clusters_table(K_BEST)
    output: mala_offset_piemap_points('{method}', '{run_label}', '{spatial_tag}', '{scenario}')
    params:
        pie_alpha       = PIEMAP_ALPHA,
        pop_label       = PIEMAP_SHOW_LABELS,
        pop_label_size  = PIEMAP_LABEL_SIZE,
        pie_scale       = PIEMAP_PIE_SCALE,
        plot_dir        = lambda wc: mala_offset_plot_dir(wc.method, wc.run_label, wc.spatial_tag, wc.scenario),
        inter_dir       = INTER,
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}maladaptation/plot_offset_piemap_points_{{method}}_{{run_label}}_{{spatial_tag}}_{{scenario}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.offset_raster} 1 "Genetic Offset" \
            {input.samples} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            genetic_offset_piemap_points {params.regionmap_extent} {params.pie_scale} Clusters T > {log} 2>&1
        """
