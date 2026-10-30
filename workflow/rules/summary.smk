#=============================================================================
# MODULE 9: PIPELINE SUMMARY
#=============================================================================

# Summary rule uses ruleorder to resolve ambiguity - only the matching mode runs
if MODE == 'processing':
    rule write_summary:
        """Write processing mode summary to Pipeline_summary.tsv."""
        input:
            vcf_filt           = W['vcf_filt'],
            vcf_ld             = W['vcf_ld'],
            samples_list       = W['samples_list'],
            # The actual final sample list -- same resolution filter_vcf uses -- so the
            # "Samples" summary box reports the true post-ALL-filtering count (including
            # relatedness/het removal), not just post-missingness.
            samples_filtered   = W['samples_rel_filtered'] if RELATEDNESS_REMOVE else (W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered']),
            samples_removed    = W['samples_removed'],
            qc_raw_summary     = O['qc_raw_summary'],
            filtering_summary  = O['qc_filtering_summary'],
            plot_done          = O['qc_plot_attrition'],
            depth_summary      = O['qc_depth_site'] if HAS_FORMAT_DP else [],
            relatedness_removed= O['qc_relatedness_removed'] if PI_HAT is not None else [],
            coord_missing_summary = O['coord_missing_summary']
        output: W['summary_done']
        params:
            summary_tsv    = O['summary'],
            het_outlier_sd = HET_OUTLIER_SD if HET_OUTLIER_SD is not None else 'NULL',
            has_dp         = 'TRUE' if HAS_FORMAT_DP else 'FALSE',
            depth_summary  = f"{MOD_PROCESSING}tables/depth_summary.tsv" if HAS_FORMAT_DP else 'NULL',
            pihat_thresh   = PI_HAT if PI_HAT is not None else 'NULL',
            relatedness_removed = O['qc_relatedness_removed'] if PI_HAT is not None else 'NULL',
            relatedness_action  = RELATEDNESS_ACTION
        log: f"{LOGDIR}processing/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                processing {params.summary_tsv} \
                {input.vcf_filt} {input.vcf_ld} \
                {input.samples_list} {input.samples_filtered} {input.samples_removed} \
                {input.qc_raw_summary} {input.filtering_summary} \
                {params.het_outlier_sd} {params.has_dp} {params.depth_summary} \
                {params.pihat_thresh} {params.relatedness_removed} {params.relatedness_action} \
                {input.coord_missing_summary} > {log} 2>&1
            touch {output}
            """

elif MODE == 'prestructure':
    rule write_summary:
        """Write prestructure mode summary to Pipeline_summary.tsv."""
        input:
            cross_entropy = O['cross_entropy']
        output: W['summary_done']
        params: ks = K_START, ke = K_END, summary_tsv = O['summary']
        log: f"{LOGDIR}prestructure/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                prestructure {params.summary_tsv} \
                {params.ks} {params.ke} > {log} 2>&1
            touch {output}
            """

elif MODE == 'pregea':
    rule write_summary:
        """Write pregea mode summary. A mode cannot complete without
        W['summary_done'], so this is what makes 'pregea' a real mode.
        Varpart/dbMEM metrics moved to the 'climate' mode summary branch
        below (module split — PreGEA no longer produces those artifacts)."""
        input:
            recs   = O['pregea_recommendations'],
            lfmm   = O['pregea_lfmm_ladder'],
            emmax  = O['pregea_emmax_ladder'],
            rda    = O['pregea_rda_ladder'],
            collin = O['pregea_rda_collin'],
            guard  = O['pregea_transfer_guard'] if PREGEA_TRANSFER_GUARD else [],
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            guard_p     = O['pregea_transfer_guard'] if PREGEA_TRANSFER_GUARD else 'NULL',
        log: f"{LOGDIR}pregea/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                pregea {params.summary_tsv} \
                {input.recs} {input.lfmm} {input.emmax} {input.rda} \
                {input.collin} {params.guard_p} > {log} 2>&1
            touch {output}
            """

elif MODE == 'climate':
    rule write_summary:
        """Write climate mode summary. A mode cannot complete without
        W['summary_done'], so this is what makes 'climate' a real mode."""
        input:
            invariant = O['climate_invariant'],
            dbmem_diag = O['climate_dbmem_diag'],
            varpart    = O['climate_vp_table'],
            confound   = O['climate_vp_confound'],
            px         = O['climate_vp_px'],
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            predictors  = ",".join(get_predictors_list())
        log: f"{LOGDIR}climate/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                climate {params.summary_tsv} \
                {input.invariant} {input.dbmem_diag} {input.varpart} {input.confound} {input.px} \
                {params.predictors} > {log} 2>&1
            touch {output}
            """

elif MODE == 'traits':
    rule write_summary:
        """Write traits mode summary. Climate correlogram path is passed as
        'NULL' in gwas_only — its presence is itself a reported metric."""
        input:
            summary   = O['trait_summary'],
            invariant = O['trait_invariant'],
            pairs     = O['trait_pairs'],
            corr      = O['trait_corr'],
            corr_clim = O['trait_corr_climate'] if CLIMATE_ENABLED else [],
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            corr_clim_p = lambda wc, input: input.corr_clim if CLIMATE_ENABLED else 'NULL',
        log: f"{LOGDIR}traits/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                traits {params.summary_tsv} \
                {input.summary} {input.invariant} {input.pairs} {params.corr_clim_p} > {log} 2>&1
            touch {output}
            """

elif MODE == 'structure':
    rule write_summary:
        """Write structure mode summary to Pipeline_summary.tsv."""
        input:
            climate_site = O['climate_site'] if CLIMATE_ENABLED else [],
            climate_na_excluded = O['climate_na_excluded'] if CLIMATE_ENABLED and CLIMATE_SOURCE == 'worldclim' else [],
            ld_decay = O['ld_decay_table'],
            imputed_maf_report = W['imputed_maf_report']
        output: W['summary_done']
        params:
            k = K_BEST,
            predictors = ALL_BIO_STR if CLIMATE_ENABLED else 'NULL',
            climate_site_path = O['climate_site'] if CLIMATE_ENABLED else 'NULL',
            climate_na_excluded_path = O['climate_na_excluded'] if CLIMATE_ENABLED and CLIMATE_SOURCE == 'worldclim' else 'NULL',
            summary_tsv = O['summary'],
            ld_decay_path = O['ld_decay_table'],
            ld_decay_group_by = LD_DECAY_GROUP_BY,
            ld_decay_scope = LD_DECAY_SCOPE
        log: f"{LOGDIR}structure/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                structure {params.summary_tsv} \
                {params.k} {params.climate_site_path} {params.predictors} \
                {params.ld_decay_path} {params.ld_decay_group_by} {params.ld_decay_scope} \
                {params.climate_na_excluded_path} {input.imputed_maf_report} > {log} 2>&1
            touch {output}
            """

elif MODE == 'gea':
    rule write_summary:
        """Write gea mode summary to Pipeline_summary.tsv."""
        input:
            selected_snps = O['selected_snps'],
            regions_per_trait = O['regions_per_trait'],
            regions_combined = O['regions_combined'],
            genes = O['genes_per_region'],
            ld_decay = ld_decay_input(REGION_DISTANCE_MODE),
            imputed_maf_report = W['imputed_maf_report_full']
        output: W['summary_done']
        params:
            summary_tsv      = O['summary'],
            region_dist_spec = REGION_DISTANCE,
            ld_rule          = REGION_LD_RULE,
            ld_decay_group   = REGION_LD_DECAY_GROUP,
            ld_decay_path    = O.get('ld_decay_table', 'NULL') if REGION_DISTANCE_MODE != 'fixed' else 'NULL'
        log: f"{LOGDIR}GEA/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                gea {params.summary_tsv} \
                {input.selected_snps} {input.regions_per_trait} {input.regions_combined} \
                {input.genes} {params.region_dist_spec} {params.ld_rule} \
                "{params.ld_decay_group}" {params.ld_decay_path} {input.imputed_maf_report} > {log} 2>&1
            touch {output}
            """

elif MODE == 'maladaptation':
    # Build fan-out input lists over all (method, set_name, spatial_tag) combinations.
    # resolve_active_snp_sets() was already called in get_targets; call again here
    # (cheap — just globs/validates) so summary.smk stays self-contained. The same
    # function, not a copy of its glob: a copy would keep listing the empty and
    # not-yet-built declared sets the resolver drops or refuses. It never returns
    # an empty list (it raises instead), so no placeholder set is needed below.
    _summary_sets = resolve_active_snp_sets()

    # Collect model artifacts only for methods that build a separate model file
    _summary_adaptive = [
        mala_model(m, s, t, 'adaptive')
        for m in ACTIVE_MALA_METHODS
        if MALADAPTATION_METHODS[m]['builds_model']
        for s in _summary_sets
        for t in mala_spatial_tags(m)
    ]
    # Offset-site TSVs for all methods × scenarios
    _summary_offset_sites = [
        mala_offset_site_values(m, s, t, sc)
        for m in ACTIVE_MALA_METHODS
        for s in _summary_sets
        for t in mala_spatial_tags(m)
        for sc in SCENARIO_NAMES
    ]

    rule write_summary:
        """Write maladaptation mode summary to Pipeline_summary.tsv (one row per method × SNP set × spatial tag × scenario).

        The offset paths go through a MANIFEST FILE, not the command line. Once
        offsets gained the scenario dimension this list became methods × sets ×
        scenarios — 1176 paths (~130 KB) on a 42-scenario sweep, which overflows
        the shell and fails the rule with exit 126."""
        input:
            mala_adaptive = _summary_adaptive,
            offset_site   = _summary_offset_sites
        output: W['summary_done']
        params:
            summary_tsv = O['summary'],
            manifest    = f"{INTER}flags/summary_offset_sites.txt"
        log: f"{LOGDIR}maladaptation/write_summary.log"
        run:
            import os
            os.makedirs(os.path.dirname(params.manifest), exist_ok=True)
            with open(params.manifest, 'w') as fh:
                fh.write("\n".join(input.offset_site) + "\n")
            shell(
                "Rscript /pipeline/scripts/write_summary.R "
                "maladaptation {params.summary_tsv} {params.manifest} > {log} 2>&1"
            )
            shell("touch {output}")

elif MODE == 'gwas':
    rule write_summary:
        """Write gwas mode summary to Pipeline_summary.tsv."""
        input:
            missing_summary = O['pheno_missing_summary'],
            selected_snps = O['pheno_selected_snps'],
            regions_per_trait = O['pheno_regions_per_trait'],
            regions_combined = O['pheno_regions_combined'],
            genes = O['pheno_genes_per_region'],
            ld_decay = ld_decay_input(PHENO_REGION_DISTANCE_MODE),
            imputed_maf_report = W['imputed_maf_report_full']
        output: W['summary_done']
        params:
            summary_tsv      = O['summary'],
            region_dist_spec = PHENO_REGION_DISTANCE,
            ld_rule          = PHENO_REGION_LD_RULE,
            ld_decay_group   = PHENO_REGION_LD_DECAY_GROUP,
            ld_decay_path    = O.get('ld_decay_table', 'NULL') if PHENO_REGION_DISTANCE_MODE != 'fixed' else 'NULL'
        log: f"{LOGDIR}GWAS/write_summary.log"
        shell:
            """
            Rscript /pipeline/scripts/write_summary.R \
                gwas {params.summary_tsv} \
                {input.missing_summary} {input.selected_snps} {input.regions_per_trait} \
                {input.regions_combined} {input.genes} {params.region_dist_spec} \
                {params.ld_rule} "{params.ld_decay_group}" {params.ld_decay_path} \
                {input.imputed_maf_report} > {log} 2>&1
            touch {output}
            """

elif MODE == 'gea_x_gwas':
    rule write_summary:
        """Overlapping mode: region analysis is fully interactive (Shiny). Just mark done."""
        input:
            [O['overlap_miami']] if (GEA_RUNS and GWAS_CONFIGS) else [],
            [O['pairwise_overlap_table']] if (GEA_RUNS or GWAS_CONFIGS) else [],
        output: W['summary_done']
        shell: "touch {output}"

