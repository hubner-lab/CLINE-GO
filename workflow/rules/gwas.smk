#=============================================================================
# MODULE 4B: PHENOTYPE ASSOCIATION (GWAS)
#=============================================================================

# --- PATH A: MEAN/MEDIAN mode (single sample set, all traits together) ---
if GWAS_CONFIGS and PHENO_MISSING != 'DROP':

    rule prepare_phenotypes:
        """Extract phenotype traits, impute missing values."""
        input: metadata = O['metadata']
        output:
            summary = O['pheno_missing_summary'],
            all_pheno = W['pheno_all_phenotypes'],
            site_means = [f"{WORK_FILT}phenotypes/{trait}_site_means.tsv" for trait in PHENO_TRAITS]
        params:
            strategy = PHENO_MISSING,
            work_dir = f"{WORK_FILT}phenotypes/"
        log: f"{LOGDIR}gwas/prepare_phenotypes.log"
        shell:
            """
            Rscript /pipeline/scripts/prepare_phenotypes.R \
                {input.metadata} {params.strategy} {params.work_dir} {output.summary} > {log} 2>&1
            """

    rule tped_gwas:
        """Convert filtered VCF to TPED/TFAM for phenotype EMMAX."""
        input: vcf = W['vcf_filt']
        output: tped = W['pheno_tped'], tfam = W['pheno_tfam']
        params: prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gwas/tped_pheno.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_gwas:
        """Compute BN kinship matrix for phenotype EMMAX."""
        input: tped = W['pheno_tped']
        output: W['pheno_kinship']
        params: prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gwas/kinship_pheno.log"
        shell:
            "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

    # ==========================================================================
    # GWAS Path A — per-factor caching (symmetric with GEA Phase 2)
    #
    # Each phenotype trait is computed independently and cached under
    #   _intermediate/gwas_per_trait/{method}/{trait}_pvalues_K{k}.tsv
    # The wide table is assembled from the currently-selected traits using
    # combine_pheno_pvalues.R (which also computes q-values).
    # ==========================================================================

    # --- EMMAX per-trait (Path A) ---
    if "EMMAX" in GWAS_OTHER_CONFIGS:

        rule gwas_a_emmax_trait:
            """Run EMMAX for a single phenotype trait in MEAN/MEDIAN mode (per-factor caching)."""
            input:
                vcf        = W['vcf_filt'],
                tped       = W['pheno_tped'],
                kinship    = W['pheno_kinship'],
                pca        = W['pca_projections'],
                phenotypes = W['pheno_all_phenotypes'],
            output: f"{INTER}gwas_per_trait/EMMAX/{{pheno_trait}}_pvalues_K{K_BEST}.tsv"
            wildcard_constraints:
                pheno_trait = r"[a-zA-Z]\w*"
            params:
                tped_prefix = f"{WORK_FILT}phenotypes/emmax/{VCF_BASE}",
                k           = K_BEST,
                tables_dir  = f"{INTER}gwas_per_trait/",
            log: f"{LOGDIR}gwas/gwas_a_emmax_{{pheno_trait}}.log"
            shell:
                "Rscript /pipeline/scripts/emmax_phenotypes.R "
                "{input.vcf} {params.tped_prefix} {input.kinship} {input.pca} "
                "{params.k} {input.phenotypes} {params.tables_dir} NULL "
                "{output} {wildcards.pheno_trait} > {log} 2>&1"

        rule gwas_a_emmax_assemble:
            """Assemble per-trait EMMAX pvalue files + compute q-values (Path A)."""
            input: [pheno_pvalues_trait("EMMAX", t) for t in PHENO_TRAITS]
            output:
                pvals = pheno_pvalues("EMMAX"),
                qvals = pheno_qvalues("EMMAX"),
            params:
                files_str = lambda wc, input: " ".join(input)
            log: f"{LOGDIR}gwas/gwas_a_emmax_assemble.log"
            shell:
                'Rscript /pipeline/scripts/combine_pheno_pvalues.R '
                '"{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1'

    # --- GAPIT per-trait (Path A) ---
    if GWAS_GAPIT_CONFIGS:

        rule gapit_gwas_trait_a:
            """Run GAPIT for a single phenotype trait in MEAN/MEDIAN mode (per-factor caching)."""
            input:
                gd         = W['gapit_gd'],
                gm         = W['gapit_gm'],
                pca        = W['pca_projections'],
                kinship    = W['pheno_kinship'],
                phenotypes = W['pheno_all_phenotypes'],
                metadata   = O['metadata'],
            output:
                [f"{INTER}gwas_per_trait/{model}/{{pheno_trait}}_pvalues_K{K_BEST}.tsv"
                 for model in GWAS_GAPIT_CONFIGS]
            wildcard_constraints:
                pheno_trait = r"[a-zA-Z]\w*"
            params:
                k             = K_BEST,
                models        = ','.join(GWAS_GAPIT_CONFIGS.keys()),
                workdir       = lambda wc: f"{INTER}gapit/gwas_a/{wc.pheno_trait}/",
                tables_dir    = f"{INTER}gwas_per_trait/",
                native_outdir = f"{MOD_GWAS}GAPIT_native_output/",
            log: f"{LOGDIR}gwas/gapit_gwas_a_{{pheno_trait}}.log"
            shell:
                """
                Rscript /pipeline/scripts/gapit.R \
                    {input.gd} {input.gm} {input.phenotypes} {input.pca} \
                    {input.kinship} {params.k} {params.models} \
                    {params.workdir} {params.tables_dir} {wildcards.pheno_trait} \
                    {input.metadata} {params.native_outdir} NULL {wildcards.pheno_trait} > {log} 2>&1
                """

        for _gwas_model in GWAS_GAPIT_CONFIGS:
            _gwas_trait_files = [pheno_pvalues_trait(_gwas_model, t) for t in PHENO_TRAITS]
            rule:
                name:   f"gwas_a_{_gwas_model.lower()}_assemble"
                input:  _gwas_trait_files
                output:
                    pvals = pheno_pvalues(_gwas_model),
                    qvals = pheno_qvalues(_gwas_model),
                params: files_str = lambda wc, input: " ".join(input)
                log:    f"{LOGDIR}gwas/gwas_a_{_gwas_model.lower()}_assemble.log"
                shell:
                    'Rscript /pipeline/scripts/combine_pheno_pvalues.R '
                    '"{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1'

# --- PATH B: DROP mode (per-trait sample sets) ---
if GWAS_CONFIGS and PHENO_MISSING == 'DROP':

    rule prepare_phenotypes:
        """Extract phenotype traits, create per-trait sample lists for DROP mode."""
        input: metadata = O['metadata']
        output:
            summary = O['pheno_missing_summary'],
            samples = [f"{WORK_FILT}phenotypes/{trait}_samples.list" for trait in PHENO_TRAITS],
            phenotypes = [f"{WORK_FILT}phenotypes/{trait}_phenotype.tsv" for trait in PHENO_TRAITS],
            site_means = [f"{WORK_FILT}phenotypes/{trait}_site_means.tsv" for trait in PHENO_TRAITS]
        params:
            strategy = PHENO_MISSING,
            work_dir = f"{WORK_FILT}phenotypes/"
        log: f"{LOGDIR}gwas/prepare_phenotypes.log"
        shell:
            """
            Rscript /pipeline/scripts/prepare_phenotypes.R \
                {input.metadata} {params.strategy} {params.work_dir} {output.summary} > {log} 2>&1
            """

    rule subset_vcf_gwas:
        """Subset filtered VCF to per-trait samples."""
        input:
            vcf = W['vcf_filt'],
            samples = f"{WORK_FILT}phenotypes/{{pheno_trait}}_samples.list"
        output: f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}"
        log: f"{LOGDIR}gwas/subset_vcf_pheno_{{pheno_trait}}.log"
        shell:
            """
            plink --vcf {input.vcf} --keep {input.samples} --const-fid 0 --allow-extra-chr \
                --recode vcf --out {params.prefix} > {log} 2>&1
            sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}
            """

    rule tped_gwas_trait:
        """Convert per-trait VCF to TPED/TFAM."""
        input: vcf = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf"
        output:
            tped = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tped",
            tfam = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tfam"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gwas/tped_pheno_trait_{{pheno_trait}}.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_gwas_trait:
        """Compute BN kinship matrix for per-trait sample set."""
        input: tped = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.tped"
        output: f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aBN.kinf"
        wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
        params: prefix = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gwas/kinship_pheno_trait_{{pheno_trait}}.log"
        shell:
            "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

    # Non-GAPIT phenotype Path B: per-trait rule + combine rule, one pair per method.
    # Shell strings are inlined per engine to avoid deferred-evaluation closure capture.
    for _method in GWAS_OTHER_CONFIGS:
        _engine           = GWAS_METHODS[_method]["engine"]
        _out_path         = f"{MOD_GWAS}tables/methods/{_method}/{{pheno_trait}}_pvalues_K{K_BEST}.tsv"
        _logpath          = f"{LOGDIR}gwas/gwas_b_{_method.lower()}_{{pheno_trait}}.log"
        _tped_pfx         = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}"
        _per_trait_inputs = expand(
            f"{MOD_GWAS}tables/methods/{_method}/{{trait}}_pvalues_K{K_BEST}.tsv",
            trait=PHENO_TRAITS,
        )

        if _engine == "emmax":
            rule:
                name:   f"gwas_b_{_method.lower()}_trait"
                input:
                    vcf       = f"{WORK_FILT}phenotypes/{{pheno_trait}}/{VCF_BASE}.vcf",
                    kinship   = f"{WORK_FILT}phenotypes/{{pheno_trait}}/emmax/{VCF_BASE}.aBN.kinf",
                    pca       = W['pca_projections'],
                    phenotype = f"{WORK_FILT}phenotypes/{{pheno_trait}}_phenotype.tsv",
                output: _out_path
                wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
                params:
                    tped_prefix   = _tped_pfx,
                    k             = K_BEST,
                    tables_dir    = f"{MOD_GWAS}tables/methods/{_method}/",
                    samples_order = W['samples_order'],
                log: _logpath
                shell:
                    "Rscript /pipeline/scripts/emmax_phenotypes.R "
                    "{input.vcf} {params.tped_prefix} {input.kinship} {input.pca} "
                    "{params.k} {input.phenotype} {params.tables_dir} "
                    "{params.samples_order} {output} > {log} 2>&1"

        rule:
            name:  f"combine_gwas_pvalues_{_method.lower()}"
            input: _per_trait_inputs
            output:
                pvals = pheno_pvalues(_method),
                qvals = pheno_qvalues(_method),
            params: files_str = lambda wc, input: ' '.join(input)
            log: f"{LOGDIR}gwas/combine_gwas_pvalues_{_method.lower()}.log"
            shell:
                "Rscript /pipeline/scripts/combine_pheno_pvalues.R "
                '"{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1'

    if GWAS_GAPIT_CONFIGS:
        rule gapit_gwas_trait:
            """Run GAPIT for a single phenotype trait (DROP mode).
            Uses full-dataset GD and kinship; gapit.R subsets internally."""
            input:
                gd = W['gapit_gd'],
                gm = W['gapit_gm'],
                pca = W['pca_projections'],
                kinship = W['assoc_kinship'],
                phenotype = f"{WORK_FILT}phenotypes/{{pheno_trait}}_phenotype.tsv",
                samples = f"{WORK_FILT}phenotypes/{{pheno_trait}}_samples.list",
                metadata = O['metadata']
            output: [f"{MOD_GWAS}tables/methods/{model}/{{pheno_trait}}_pvalues_K{K_BEST}.tsv" for model in GWAS_GAPIT_CONFIGS]
            wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
            params:
                k = K_BEST,
                models = ','.join(GWAS_GAPIT_CONFIGS.keys()),
                workdir = lambda wc: f"{INTER}gapit/phenotype_association/{wc.pheno_trait}/",
                tables_dir = f"{MOD_GWAS}tables/methods/",
                native_outdir = f"{MOD_GWAS}GAPIT_native_output/"
            log: f"{LOGDIR}gwas/gapit_pheno_trait_{{pheno_trait}}.log"
            shell:
                """
                Rscript /pipeline/scripts/gapit.R \
                    {input.gd} {input.gm} {input.phenotype} {input.pca} \
                    {input.kinship} {params.k} {params.models} \
                    {params.workdir} {params.tables_dir} {wildcards.pheno_trait} \
                    {input.metadata} {params.native_outdir} {input.samples} \
                    {wildcards.pheno_trait} > {log} 2>&1
                """

        rule combine_gapit_gwas_pvalues:
            """Merge per-trait GAPIT p-value files into combined table (DROP mode)."""
            input: lambda wc: expand(f"{MOD_GWAS}tables/methods/{wc.method}/{{trait}}_pvalues_K{K_BEST}.tsv", trait=PHENO_TRAITS)
            output:
                pvals = f"{MOD_GWAS}tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv",
                qvals = f"{MOD_GWAS}tables/methods/{{method}}/{{method}}_qvalues_K{K_BEST}.tsv"
            wildcard_constraints:
                method = '|'.join(GWAS_GAPIT_CONFIGS.keys())
            params: files_str = lambda wc, input: ' '.join(input)
            log: f"{LOGDIR}gwas/combine_gwas_pvalues_{{method}}.log"
            shell:
                """
                Rscript /pipeline/scripts/combine_pheno_pvalues.R \
                    "{params.files_str}" {output.pvals} {output.qvals} > {log} 2>&1
                """

# piemap_pheno is GWAS-specific (no GEA equivalent); downstream rules moved to _assoc_downstream.smk
if GWAS_CONFIGS:

  rule piemap_gwas:
      """Generate pie map for phenotype trait with trait-controlled pie sizes.
      Uses metadata_climate.tsv for coordinates (samples with missing lat/lon can't be
      plotted spatially, but remain in the phenotype association analysis itself)."""
      input:
          raster = W['climate_raster'],
          metadata = W['metadata_climate'],
          clusters = clusters_table(K_BEST),
          trait = f"{WORK_FILT}phenotypes/{{pheno_trait}}_site_means.tsv"
      output:
          png = f"{MOD_GWAS}plots/piemap/phenomap_{{pheno_trait}}.png",
          svg = f"{MOD_GWAS}plots/piemap/phenomap_{{pheno_trait}}.svg"
      wildcard_constraints: pheno_trait = r"[a-zA-Z]\w*"
      params:
          raster_layer = get_predictors_list()[0] if get_predictors_list() else "1",
          pie_alpha = PIEMAP_ALPHA,
          pop_label = PIEMAP_SHOW_LABELS,
          pop_label_size = PIEMAP_LABEL_SIZE,
          plot_dir = f"{MOD_GWAS}plots/piemap/",
          inter_dir = INTER,
          regionmap_extent = REGIONMAP_EXTENT,
          pie_scale = PIEMAP_PIE_SCALE
      log: f"{LOGDIR}gwas/piemap_pheno_{{pheno_trait}}.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_piemap.R \
              {input.raster} {params.raster_layer} {params.raster_layer} \
              {input.metadata} {input.clusters} \
              {input.trait} "{wildcards.pheno_trait}" \
              {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
              {params.plot_dir} {params.inter_dir} \
              phenomap_{wildcards.pheno_trait} {params.regionmap_extent} {params.pie_scale} Clusters > {log} 2>&1
          """

  rule piemap_gwas_points:
      """Points ('clear map') companion for the phenotype piemap — trait-independent
      (points only mark sample locations), one per project."""
      input:
          raster = W['climate_raster'],
          metadata = W['metadata_climate'],
          clusters = clusters_table(K_BEST)
      output:
          png = f"{MOD_GWAS}plots/piemap/phenomap_points.png",
          svg = f"{MOD_GWAS}plots/piemap/phenomap_points.svg"
      params:
          raster_layer = get_predictors_list()[0] if get_predictors_list() else "1",
          pie_alpha = PIEMAP_ALPHA,
          pop_label = PIEMAP_SHOW_LABELS,
          pop_label_size = PIEMAP_LABEL_SIZE,
          plot_dir = f"{MOD_GWAS}plots/piemap/",
          inter_dir = INTER,
          regionmap_extent = REGIONMAP_EXTENT,
          pie_scale = PIEMAP_PIE_SCALE
      log: f"{LOGDIR}gwas/piemap_pheno_points.log"
      shell:
          """
          Rscript /pipeline/scripts/plot_piemap.R \
              {input.raster} {params.raster_layer} {params.raster_layer} \
              {input.metadata} {input.clusters} \
              NULL NULL \
              {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
              {params.plot_dir} {params.inter_dir} \
              phenomap_points {params.regionmap_extent} {params.pie_scale} Clusters T > {log} 2>&1
          """
