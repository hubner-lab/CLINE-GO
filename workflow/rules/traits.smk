# MODULE: TRAITS — phenotypic factor characterization (mode=traits)
#
# The phenotypic counterpart of the climate module: "what do my traits look
# like, how do they covary, and do any of them carry no signal at all".
#
# Runs in BOTH regimes. Climate characterization is unavailable in gwas_only
# (mode=climate raises when Climate.enabled is false) — which is exactly the
# regime where traits are the ONLY factors a user has, so nothing here may
# depend on climate or on coordinates. The one exception is the traits x
# climate correlogram, which is a separate rule gated on CLIMATE_ENABLED.
#
# Everything reads O['metadata'] (the full cohort, coordinates optional). The
# only rule that touches a climate-side file is trait_correlogram_climate.

rule trait_density:
    """Density plot for every phenotypic trait (metadata columns 5+).
    Moved here from climate.smk (Phase 3): the density plot describes traits,
    not climate, and must survive Climate.enabled: false."""
    input:  meta = O['metadata']
    output: DENSITY_PLOT_PHENOTYPES
    params:
        predictors = "all",
        inter_dir = INTER
    log:    f"{LOGDIR}traits/trait_density.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_density.R \
            {input.meta} {params.predictors} {output} {params.inter_dir} > {log} 2>&1
        """

rule trait_correlogram:
    """Correlogram of the traits alone. Always produced — this is the file the
    Shiny Phenotypic tab falls back to when there is no climate."""
    input:  meta = O['metadata']
    output: O['trait_corr']
    params: inter_dir = INTER, second = 'NULL', title = "Correlogram of traits"
    log:    f"{LOGDIR}traits/trait_correlogram.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_correlation_heatmap.R \
            {input.meta} {params.second} {output} {params.inter_dir} \
            '{params.title}' > {log} 2>&1
        """

rule trait_correlation_table:
    """The matrix trait_correlogram draws, as a pair table (bffcfc). Its own rule: the
    correlogram is an input of write_summary, so a new output there would rerun both."""
    input:  meta = O['metadata']
    output: O['trait_corr_table']
    params: second = 'NULL'
    log:    f"{LOGDIR}traits/trait_correlation_table.log"
    shell:
        """
        Rscript /pipeline/scripts/correlation_table.R \
            {input.meta} {params.second} {output} > {log} 2>&1
        """

if CLIMATE_ENABLED:
    rule trait_correlogram_climate:
        """Correlogram of traits AND climate predictors together — how much of
        a trait's variation tracks the environment. Standard regime only.
        Uses the climate-valid metadata subset so the two blocks describe the
        same samples; plot_correlation_heatmap.R joins them on `sample`. Every
        trait beside the selected predictors (Climate.predictors cuts the climate
        block only, c9edc9)."""
        input:
            meta    = W['metadata_climate_valid'],
            climate = O['climate_site'],
        output: O['trait_corr_climate']
        params: inter_dir = INTER, title = "Correlogram of traits and climate factors",
                factors = ",".join(get_predictors_list())
        log:    f"{LOGDIR}traits/trait_correlogram_climate.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_correlation_heatmap.R \
                {input.meta} {input.climate} {output} {params.inter_dir} \
                '{params.title}' {params.factors} > {log} 2>&1
            """

    rule trait_climate_correlation_table:
        """The matrix trait_correlogram_climate draws, as a pair table (bffcfc), from the
        same two tables joined on `sample`."""
        input:
            meta    = W['metadata_climate_valid'],
            climate = O['climate_site'],
        output: O['trait_corr_climate_table']
        log:    f"{LOGDIR}traits/trait_climate_correlation_table.log"
        shell:
            """
            Rscript /pipeline/scripts/correlation_table.R \
                {input.meta} {input.climate} {output} > {log} 2>&1
            """

rule trait_pairs:
    """Pairwise scatter/density/correlation grid across traits. Above
    Traits.pairs_max_factors the script writes a placeholder instead (both PNG
    and SVG) rather than an unreadable k x k grid."""
    input:  meta = O['metadata']
    output: O['trait_pairs']
    params: max_factors = TRAITS_PAIRS_MAX
    log:    f"{LOGDIR}traits/trait_pairs.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_trait_pairs.R \
            {input.meta} {output} {params.max_factors} > {log} 2>&1
        """

rule trait_summary:
    """Per-trait n / missingness / mean / sd / min / median / max."""
    input:  meta = O['metadata']
    output: O['trait_summary']
    log:    f"{LOGDIR}traits/trait_summary.log"
    shell:
        """
        Rscript /pipeline/scripts/trait_summary.R \
            {input.meta} {output} > {log} 2>&1
        """

rule trait_invariant:
    """Detect invariant (zero-variance or all-NA) traits — the phenotypic twin
    of check_climate_variance, same script with COLUMNS=traits."""
    input:  meta = O['metadata']
    output: O['trait_invariant']
    params: columns = "traits"
    log:    f"{LOGDIR}traits/trait_invariant.log"
    shell:
        """
        Rscript /pipeline/scripts/check_climate_variance.R \
            {input.meta} {output} {params.columns} > {log} 2>&1
        """
