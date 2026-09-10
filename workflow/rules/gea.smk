#=============================================================================
# MODULE 4: ASSOCIATION (GEA)
# Non-GAPIT methods (EMMAX, LFMM, RDA) are each a hand-written rule block
# gated on `if "<METHOD>" in GEA_OTHER_CONFIGS:` — GEA_METHODS registry drives
# config parsing/validation/params (see common.smk), not rule generation for
# these three (contrast the registry-driven loop idiom in gwas.smk:207-240,
# which fits N methods sharing one shell template; these three don't share one).
# GAPIT models are declared as a single multi-output rule (gapit_gea_trait).
#=============================================================================

# Full dataset processing (non-LD pruned) - needed for LFMM
rule vcf_to_lfmm_full:
    """Convert full filtered VCF to LEA formats for association analysis."""
    input:  vcf = W['vcf_filt']
    output:
        geno = W['geno_full'],
        lfmm = W['lfmm_full'],
        vcfsnp = W['vcfsnp_full'],
        removed = W['removed_full'],
        nmissing = W['lfmm_full_nmissing']
    log: f"{LOGDIR}gea/vcf_to_lfmm_full.log"
    shell: "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

rule snmf_full:
    """Run sNMF on full (non-LD pruned) dataset for imputation."""
    input:  geno = W['geno_full']
    output: W['snmf_full']
    params: ks = K_BEST, ke = K_BEST, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}gea/snmf_full.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/snmf.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule impute_full:
    """Impute missing genotypes in full dataset using SNMF."""
    input:  snmf = W['snmf_full'], lfmm = W['lfmm_full'], nmissing = W['lfmm_full_nmissing']
    output: W['lfmm_imp_full']
    params: k = K_BEST
    log:    f"{LOGDIR}gea/impute_full.log"
    shell:
        """
        Rscript /pipeline/scripts/impute.R {input.snmf} {input.lfmm} {params.k} {output} {input.nmissing} > {log} 2>&1
        """

rule lfmm2vcf_full:
    """Convert imputed full LFMM to VCF format."""
    input:  vcf = W['vcf_filt'], lfmm_imp = W['lfmm_imp_full']
    output: W['vcf_imp_full']
    params: blocksize = 20000
    log:    f"{LOGDIR}gea/lfmm2vcf_full.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """

# TPED/kinship for association (shared by EMMAX and GAPIT)
rule tped_gea:
    """Convert filtered VCF to TPED/TFAM for association EMMAX."""
    input: vcf = W['vcf_filt']
    output: tped = W['assoc_tped'], tfam = W['assoc_tfam']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}gea/tped_assoc.log"
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
            --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
        awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
            {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
        """

rule kinship_gea:
    """Compute BN kinship matrix for association analysis."""
    input: tped = W['assoc_tped']
    output: W['assoc_kinship']
    params: prefix = f"{WORK_FILT}emmax/{VCF_BASE}"
    log: f"{LOGDIR}gea/kinship_assoc.log"
    shell:
        "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

# --- EMMAX climate-coord subset (climate coordinate NA handling) ---
# EMMAX binds climate values to genotypes positionally (cbind on TFAM row order),
# so when samples with missing lat/lon or NA-climate raster extraction are dropped from
# the climate table (filter_coord_samples / filter_climate_valid_samples), the genotype-side
# TPED/kinship must be rebuilt on the same climate-valid sample set. Mirrors
# subset_vcf_gwas/tped_gwas_trait exactly.
if "EMMAX" in GEA_OTHER_CONFIGS:

    rule subset_vcf_gea_climate:
        """Subset filtered VCF to climate-valid samples for GEA EMMAX."""
        input:
            vcf = W['vcf_filt'],
            samples = W['climate_valid_samples']
        output: W['vcf_filt_climate']
        params: prefix = f"{WORK_FILT}climate/{VCF_BASE}"
        log: f"{LOGDIR}gea/subset_vcf_gea_climate.log"
        shell:
            """
            plink --vcf {input.vcf} --keep {input.samples} --const-fid 0 --allow-extra-chr \
                --recode vcf --out {params.prefix} > {log} 2>&1
            sed -i '/^#CHROM/s/\\t0_/\\t/g' {output}
            """

    rule tped_gea_climate:
        """Convert coord-valid-subset VCF to TPED/TFAM for association EMMAX."""
        input: vcf = W['vcf_filt_climate']
        output: tped = W['assoc_tped_climate'], tfam = W['assoc_tfam_climate']
        params: prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gea/tped_gea_climate.log"
        shell:
            """
            plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
                --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
            awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
                {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
            """

    rule kinship_gea_climate:
        """Compute BN kinship matrix for the coord-valid subset. Always built
        (GAPIT consumes this one unconditionally, regardless of EMMAX.kinship)."""
        input: tped = W['assoc_tped_climate']
        output: W['assoc_kinship_climate']
        params: prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}"
        log: f"{LOGDIR}gea/kinship_gea_climate.log"
        shell:
            "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"

    # EMMAX.kinship='IBS' flavor — a SEPARATE file, never overwriting the BN one
    # above. Only built when configured (GEA_PARAMS gate below), so a config
    # that never sets kinship never adds this rule to the DAG.
    if GEA_PARAMS.get("EMMAX", {}).get("kinship") == "IBS":

        rule kinship_gea_climate_ibs:
            """Compute IBS kinship matrix for the coord-valid subset (EMMAX.kinship='IBS').
            emmax-kin-intel64 -s switches the estimator to IBS (default, no flag: BN) and
            names its own output file [tped].[aBN|aIBS].kinf accordingly — no -o needed."""
            input: tped = W['assoc_tped_climate']
            output: W['assoc_kinship_climate_ibs']
            params: prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}"
            log: f"{LOGDIR}gea/kinship_gea_climate_ibs.log"
            shell:
                "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -s -d 10 -x {params.prefix} > {log} 2>&1"

# =============================================================================
# GEA Association rules — per-factor caching (Phase 2)
#
# Each bioclimatic factor is computed independently and persisted as:
#   _intermediate/gea_per_trait/{method}/{trait}_pvalues_K{k}.tsv
#
# The wide table consumed by downstream rules is assembled from the selected
# factors. Changing the selection only reruns the assembly + any NEW factors;
# factors already computed are reused. Changing K or filter params changes the
# imputation input path, cascading invalidation across all factors.
# =============================================================================

# --- EMMAX per-trait ---
if "EMMAX" in GEA_OTHER_CONFIGS:

    rule assoc_emmax_gea_trait:
        """Run EMMAX for a single GEA bioclimatic trait (per-factor caching).
        Uses the climate-valid-subset VCF/TPED/kinship (see subset_vcf_gea_climate) since
        climate values (traits) exclude samples with missing lat/lon or NA-climate raster
        extraction; PCA covariates stay on the full cohort and are row-subset internally by
        load_pca_covariates() via samples_order (see emmax.R). Kinship file and #PCs are
        both driven by EMMAX.kinship / EMMAX.n_pcs (default: BN / sNMF.k_best — byte-
        identical to pre-hyperparameter-system behavior)."""
        input:
            vcf        = W['vcf_filt_climate'],
            tped       = W['assoc_tped_climate'],
            kinship    = emmax_kinship_climate_path(),
            traits     = O['climate_site_scaled'],
            covariates = W['pca_projections'],
            metadata   = W['metadata_climate_valid'],
            samples_order = W['samples_order'],
        output: f"{INTER}gea_per_trait/EMMAX/{{trait}}_pvalues_K{K_BEST}.tsv"
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            k           = GEA_PARAMS.get("EMMAX", {}).get("n_pcs", K_BEST),
            inter_dir   = INTER,
            tped_prefix = f"{WORK_FILT}climate/emmax/{VCF_BASE}",
        log: f"{LOGDIR}gea/assoc_emmax_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/emmax.R \
                {input.vcf} {params.k} {input.traits} {input.covariates} \
                {wildcards.trait} {params.inter_dir} {input.metadata} \
                {output} {params.tped_prefix} {input.kinship} {input.samples_order} > {log} 2>&1
            """

    rule assoc_emmax_gea_assemble:
        """Assemble per-trait EMMAX files into the wide table consumed by downstream."""
        input: [assoc_pvalues_trait("EMMAX", t) for t in get_predictors_list()]
        output: assoc_pvalues("EMMAX")
        params:
            traits_str = lambda wc, input: " ".join(input)
        log: f"{LOGDIR}gea/assoc_emmax_assemble.log"
        shell:
            "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"

# --- LFMM per-trait ---
if "LFMM" in GEA_OTHER_CONFIGS:

    rule assoc_lfmm_gea_trait:
        """Run LFMM for a single GEA bioclimatic trait (per-factor caching).
        Uses coord-valid-subset genotype matrices (see subset_lfmm_climate in
        processing.smk) since climate (traits) excludes samples with missing lat/lon
        and LFMM binds climate rows to genotype-matrix rows positionally."""
        input:
            lfmm_ld  = W['lfmm_imp_climate'],
            lfmm_full = W['lfmm_imp_full_climate'],
            climate  = O['climate_site_scaled'],
            vcfsnp   = W['vcfsnp_full'],
        output: f"{INTER}gea_per_trait/LFMM/{{trait}}_pvalues_K{K_BEST}.tsv"
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            # GEA.configs[LFMM].params.K, NOT sNMF.k_best. The registry declares K with
            # the @k_best sentinel as its default (gea.py), so an unset K still resolves
            # to k_best here — but an explicitly configured K now actually reaches
            # lfmm2(). Before this, the param was parsed, range-validated and recommended
            # by pregea_recommend.R, yet silently discarded at the rule.
            k = GEA_PARAMS.get('LFMM', {}).get('K', K_BEST),
            gc = "TRUE" if GEA_PARAMS.get('LFMM', {}).get('genomic_control', True) else "FALSE",
        log: f"{LOGDIR}gea/assoc_lfmm_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/lfmm.R \
                {input.lfmm_ld} {input.lfmm_full} {input.climate} \
                {params.k} {wildcards.trait} {input.vcfsnp} \
                {output} {params.gc} > {log} 2>&1
            """

    rule assoc_lfmm_gea_assemble:
        """Assemble per-trait LFMM files into the wide table consumed by downstream."""
        input: [assoc_pvalues_trait("LFMM", t) for t in get_predictors_list()]
        output: assoc_pvalues("LFMM")
        params:
            traits_str = lambda wc, input: " ".join(input)
        log: f"{LOGDIR}gea/assoc_lfmm_assemble.log"
        shell:
            "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"

# --- RDA (multivariate genome scan, Capblancq & Forester 2021 / rdadapt) ---
# RDA fits ONE multivariate model across all predictors jointly and emits ONE
# p-value column (the "climate_multivariate" pseudo-trait) — there is no
# per-trait decomposition to cache, so unlike EMMAX/LFMM above this writes the
# wide p-value table directly and has no per-trait + assemble_pvalues.R pair.
if "RDA" in GEA_OTHER_CONFIGS:
    _rda_p = GEA_PARAMS["RDA"]

    rule assoc_rda_gea:
        """Multivariate RDA genome scan. Fits a partial RDA (Condition(PC1..PCk)
        from the LD-pruned LEA PCA) on the FULL marker set by default (literature
        default per Capblancq & Forester 2021 — see docs/rda_research.md A6),
        calls candidates via robust-Mahalanobis rdadapt (Capblancq et al. 2018),
        and writes the wide pvalue table plus a candidates/diagnostics/anova side
        table set and 3 diagnostic plots. See scripts/rda.R for the full algorithm."""
        input:
            lfmm_full     = W['lfmm_imp_full_climate'],
            lfmm_pruned   = W['lfmm_imp_climate'],
            climate       = O['climate_site_scaled'],
            vcfsnp_full   = W['vcfsnp_full'],
            vcfsnp_pruned = W['vcfsnp'],
            pca           = W['pca_projections'],
            samples_order = W['samples_order'],
            climate_valid = W['climate_valid_samples'],
        output:
            pvalues       = assoc_pvalues("RDA"),
            candidates    = O['rda_candidates'],
            diagnostics   = O['rda_diagnostics'],
            anova         = O['rda_anova'],
            screeplot     = O['rda_screeplot'],     screeplot_svg = O['rda_screeplot_svg'],
            pval_hist     = O['rda_pval_hist'],     pval_hist_svg = O['rda_pval_hist_svg'],
            biplot        = O['rda_biplot'],        biplot_svg    = O['rda_biplot_svg'],
        params:
            predictors    = _rda_p["predictor_set"] if _rda_p["predictor_set"] != "auto" else PREDICTORS_SELECTED,
            cond_pcs      = _rda_p["condition_pcs"],
            axes          = _rda_p["axes"],
            axis_alpha    = _rda_p["axis_alpha"],
            perms         = _rda_p["permutations"],
            fit_mode      = _rda_p["fit_mode"],
            max_snps      = _rda_p["full_fit_max_snps"],
            max_gb        = _rda_p["full_fit_max_gb"],
            seed          = _rda_p["seed"],
            plot_dir      = f"{MOD_GEA}plots/rda/",
            k_best        = K_BEST,
            # GEA_CONFIGS["RDA"] is the fused "{adjust}_{threshold}" token built
            # at common.smk's parse_association_configs(). Splitting it here
            # (rather than plumbing cfg['adjust']/cfg['threshold'] separately)
            # keeps that function's return shape untouched and guarantees these
            # two args reconstruct exactly the string that names the sig_snps
            # file downstream — unambiguous under the widened adjust grammar
            # since neither field may contain '_'.
            cand_adjust    = GEA_CONFIGS["RDA"].split("_", 1)[0],
            cand_threshold = GEA_CONFIGS["RDA"].split("_", 1)[1],
        threads: CPU
        log: f"{LOGDIR}gea/assoc_rda.log"
        shell:
            """
            Rscript /pipeline/scripts/rda.R \
                {input.lfmm_full} {input.lfmm_pruned} {input.climate} \
                {input.vcfsnp_full} {input.vcfsnp_pruned} \
                {input.pca} {input.samples_order} {input.climate_valid} \
                {params.predictors} {params.cond_pcs} {params.axes} {params.axis_alpha} \
                {params.perms} {params.fit_mode} {params.max_snps} {params.max_gb} \
                {params.seed} {threads} {params.k_best} \
                {output.pvalues} {output.candidates} {output.diagnostics} {output.anova} \
                {params.plot_dir} {params.cand_adjust} {params.cand_threshold} > {log} 2>&1
            """

# VCF to GAPIT numeric (shared by GEA and GWAS GAPIT models)
if GEA_GAPIT_CONFIGS or GWAS_GAPIT_CONFIGS:

    rule vcf_to_numeric:
        """Convert imputed VCF to GAPIT numeric format (GD + GM)."""
        input: vcf = W['vcf_imp_full']
        output:
            gd = W['gapit_gd'],
            gm = W['gapit_gm']
        log: f"{LOGDIR}gea/vcf_to_numeric.log"
        shell:
            "Rscript /pipeline/scripts/vcf_to_gapit_numeric.R {input.vcf} {output.gd} {output.gm} > {log} 2>&1"

# --- GAPIT GEA per-trait ---
if GEA_GAPIT_CONFIGS:

    def _gapit_shared_npcs():
        """The number of PCA covariates handed to GAPIT, from GEA.configs params.n_pcs.

        gapit_gea_trait runs every configured GAPIT model in ONE gapit.R call and
        gapit.R takes a single scalar (arg 6 -> pca_raw[, 1:K]), so per-model n_pcs is
        not expressible without splitting the rule. Rather than silently honouring one
        model's value and discarding the others, disagreement is a hard error: the user
        gets told to split the run instead of getting a number that is wrong for all but
        one model. Unset params resolve to the @k_best sentinel, so the all-defaults case
        agrees trivially and keeps the previous behaviour.
        """
        vals = {m: GEA_PARAMS.get(m, {}).get('n_pcs', K_BEST) for m in GEA_GAPIT_CONFIGS}
        distinct = set(vals.values())
        if len(distinct) > 1:
            raise ValueError(
                "GAPIT models in GEA.configs request different params.n_pcs "
                f"({vals}), but all configured GAPIT models share one gapit.R call and "
                "therefore one PCA covariate count. Give them the same n_pcs, or run "
                "them in separate pipeline invocations."
            )
        return distinct.pop()

    rule gapit_gea_trait:
        """Run GAPIT for a single GEA bioclimatic trait, all configured models in one call.
        Output files land under _intermediate/gea_per_trait/{model}/{trait}_pvalues_K{k}.tsv.
        GAPIT writes {TABLES_DIR}/{model}/{OUTPUT_PREFIX}_pvalues_K{k}.tsv (existing behaviour
        when OUTPUT_PREFIX is set), so we reuse gapit.R unchanged.
        Uses full-dataset GD/kinship; gapit.R subsets internally via coord_samples to the
        climate-valid sample set (coord-valid minus ocean/NoData raster NAs), matching the
        set already baked into O['climate_site_scaled'] — same mechanism as GWAS DROP mode."""
        input:
            gd       = W['gapit_gd'],
            gm       = W['gapit_gm'],
            traits   = O['climate_site_scaled'],
            pca      = W['pca_projections'],
            kinship  = W['assoc_kinship'],
            metadata = O['metadata'],
            coord_samples = W['climate_valid_samples'],
        output:
            [f"{INTER}gea_per_trait/{model}/{{trait}}_pvalues_K{K_BEST}.tsv"
             for model in GEA_GAPIT_CONFIGS]
        wildcard_constraints:
            trait = r"bio_\d+"
        params:
            # k is the FILENAME tag only (gapit.R writes "{prefix}_pvalues_K{k}.tsv", which
            # this rule declares as its output), so it must stay sNMF k_best. The PCA
            # covariate count travels separately as n_pcs.
            k             = K_BEST,
            n_pcs         = _gapit_shared_npcs(),
            models        = ','.join(GEA_GAPIT_CONFIGS.keys()),
            workdir       = lambda wc: f"{INTER}gapit/gea/{wc.trait}/",
            tables_dir    = f"{INTER}gea_per_trait/",
            native_outdir = f"{MOD_GEA}GAPIT_native_output/",
        log: f"{LOGDIR}gea/gapit_gea_{{trait}}.log"
        shell:
            """
            Rscript /pipeline/scripts/gapit.R \
                {input.gd} {input.gm} {input.traits} {input.pca} \
                {input.kinship} {params.k} {params.models} \
                {params.workdir} {params.tables_dir} {wildcards.trait} \
                {input.metadata} {params.native_outdir} {input.coord_samples} \
                {wildcards.trait} {params.n_pcs} > {log} 2>&1
            """

    for _gapit_model in GEA_GAPIT_CONFIGS:
        _gapit_trait_files = [assoc_pvalues_trait(_gapit_model, t) for t in get_predictors_list()]
        rule:
            name:   f"assoc_{_gapit_model.lower()}_gea_assemble"
            input:  _gapit_trait_files
            output: assoc_pvalues(_gapit_model)
            params: traits_str = lambda wc, input: " ".join(input)
            log:    f"{LOGDIR}gea/assoc_{_gapit_model.lower()}_assemble.log"
            shell:
                "Rscript /pipeline/scripts/assemble_pvalues.R {params.traits_str} {output} > {log} 2>&1"
