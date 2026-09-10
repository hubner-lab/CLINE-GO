# MODULE 2b: preGEA — GEA hyperparameter exploration (LD-pruned dataset only)
#
# Every rule here runs on the LD-PRUNED marker set. This is the cheap-
# exploration stage (docs/rda_research.md Part C, C.2). It deliberately does
# NOT extend to the GEA module's final RDA scan, which fits on the FULL set
# per A6. Redundant by design: sweep the ladders once so the user picks GEA
# hyperparameters with full visibility before the expensive full-SNP GEA run.
#
# preGEA is now ONE decision — LFMM/EMMAX/RDA always run together (no
# per-block enabled switches; only the opt-in TransferGuard stays separate,
# since it is the one rule that drags the full-set imputation chain into an
# otherwise LD-pruned-only mode). Predictor characterization (correlation,
# density, dbMEM, varpart) moved to climate.smk (mode=climate) — see that
# file's header for the module-split rationale.
#
# Three blocks, joined by pregea_recommend:
#   Block 0 — shared prep (LD-pruned + climate-valid VCF -> TPED -> BN kinship)
#   Block 1 — LFMM-K ladder                  (reuses scripts/lfmm.R)
#   Block 2 — EMMAX-#PC ladder               (reuses scripts/emmax.R)
#   Block 3 — RDA setup                      (scripts/pregea_rda_setup.R)
#   pregea_recommend                          (scripts/pregea_recommend.R)
#   pregea_transfer_guard (opt-in)            (scripts/pregea_transfer_guard.R)

# =============================================================================
# Block 0 — shared prep (EMMAX-only: LFMM reuses W['lfmm_imp_climate'] as-is,
# already LD-pruned + climate-valid + K_BEST-imputed — no new prep needed).
# Kinship is always BN — GAPIT uses the same BN matrix (scripts/gapit.R), and
# an IBS alternative was never exercised; not a user-facing choice.
# =============================================================================

rule pregea_subset_vcf_pruned_climate:
    """LD-pruned VCF subset to climate-valid samples. Mirrors
    gea.smk:91-104 (subset_vcf_gea_climate) but on WORK_LD, not WORK_FILT —
    the pruned+climate-valid VCF does not exist anywhere else in the
    pipeline (the EMMAX GEA chain is full-set only).

    W['vcf_ld'] (ld_prune's output) already carries a SINGLE '0_' FID
    prefix on every sample name in its #CHROM header — ld_prune reads
    already-'0_'-prefixed W['vcf_filt'] with its own --const-fid pass,
    and since plink's --const-fid only overrides FID (never IID), the
    written IID is the literal, unmodified VCF sample string it read —
    i.e. the pre-existing '0_' prefix survives into the IID column
    untouched (ld_prune's own 0_0_ sed only fires when TWO const-fid
    passes stack; here there's only one, so nothing to strip). Reading
    vcf_ld AGAIN with --const-fid here is exactly that second stacked
    pass: the --keep file's IID column must be pre-matched to the
    existing single prefix (awk below), and the resulting OUTPUT then
    carries a genuine double '0_0_' prefix — which the sed below strips,
    mirroring the exact convention documented at ld_prune (processing.smk)."""
    input:  vcf = W['vcf_ld'], samples = W['climate_valid_samples']
    output: W['vcf_ld_climate']
    params: prefix = f"{WORK_LD}climate/{VCF_BASE}"
    log:    f"{LOGDIR}pregea/subset_vcf_pruned_climate.log"
    shell:
        """
        awk '{{print $1, "0_"$2}}' {input.samples} > {params.prefix}_keep.tmp
        plink --vcf {input.vcf} --keep {params.prefix}_keep.tmp --const-fid 0 --allow-extra-chr \
            --recode vcf --out {params.prefix} > {log} 2>&1
        rm -f {params.prefix}_keep.tmp
        sed -i '/^#CHROM/s/\\t0_0_/\\t/g' {output}
        """

rule pregea_tped_pruned:
    """TPED/TFAM for the pruned+climate-valid VCF. Mirrors gea.smk:106-118."""
    input:  vcf = W['vcf_ld_climate']
    output: tped = W['pregea_tped'], tfam = W['pregea_tfam']
    params: prefix = f"{WORK_LD}climate/emmax/{VCF_BASE}"
    log:    f"{LOGDIR}pregea/tped_pruned.log"
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --recode12 transpose \
            --output-missing-genotype 0 --out {params.prefix} > {log} 2>&1
        awk '{{split($1,a,"_"); split($2,b,"_"); if(a[1]==b[1]){{$1=a[1];$2=a[1]}} print}}' \
            {output.tfam} > {params.prefix}_tmp.tfam && mv {params.prefix}_tmp.tfam {output.tfam}
        """

rule pregea_kinship_pruned:
    """BN kinship on the pruned+climate-valid TPED. Mirrors gea.smk:120-130.
    Kinship is held FIXED across the #PC ladder — only the PCA covariate
    count varies (C.2: 'with a kinship matrix already in the model')."""
    input:  tped = W['pregea_tped']
    output: W['pregea_kinship']
    params: prefix = f"{WORK_LD}climate/emmax/{VCF_BASE}"
    log:    f"{LOGDIR}pregea/kinship_pruned.log"
    shell:  "/pipeline/scripts/emmax_run.sh /pipeline/scripts/emmax-kin-intel64 -v -d 10 -x {params.prefix} > {log} 2>&1"


rule pregea_screeplot:
    """LD-pruned genotype PCA screeplot: eigenvalue bars + broken-stick null,
    with the swept LFMM K band marked. Anchors the K range (C.1) and labels
    the structure-axis count as an UPPER BOUND, not a target (C.2) — the
    caveat text itself lives in the Shiny note (Code Structure Rule 8)."""
    input:  eigenvalues = W['pca_eigenvalues']
    output: png = O['pregea_screeplot'], svg = O['pregea_screeplot_svg'],
            tsv = O['pregea_screeplot_tsv']
    params: k_best = K_BEST, k_range = PREGEA_LFMM_K_STR,
            n_pcs_max = PREGEA_NPC_MAX, inter_dir = INTER
    log:    f"{LOGDIR}pregea/screeplot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pregea_screeplot.R \
            {input.eigenvalues} {params.k_best} {params.k_range} \
            {params.n_pcs_max} {output.png} {output.tsv} {params.inter_dir} > {log} 2>&1
        """

# =============================================================================
# Block 1 — LFMM K ladder (reuses scripts/lfmm.R + scripts/assemble_pvalues.R)
# genomic_control is ALWAYS FALSE here (fixed, not configurable) so lambda
# stays informative — docs/rda_research.md C.0. This is independent of the
# GEA module's own LFMM genomic_control switch (workflow/methods/gea.py),
# which defaults TRUE for the production run.
# =============================================================================

rule pregea_lfmm_rung_trait:
    """One LFMM rung: lfmm2(K={k}) for one trait, trained AND tested on the
    LD-pruned imputed matrix (W['lfmm_imp_climate'] — already LD-pruned +
    climate-valid + K_BEST-imputed). Imputation stays FIXED at K_BEST —
    sweeping imputation too would cascade into impute/snmf/lfmm2vcf, out
    of scope for this module. REUSES scripts/lfmm.R verbatim except for
    the optional arg 8 (GENOMIC_CONTROL), always FALSE here."""
    input:
        lfmm_ld   = W['lfmm_imp_climate'],
        lfmm_test = W['lfmm_imp_climate'],
        climate   = O['climate_site_scaled'],
        vcfsnp    = W['vcfsnp'],
    output: f"{INTER}pregea/lfmm/K{{k}}/{{trait}}_pvalues.tsv"
    wildcard_constraints: k = r"\d+", trait = r"[A-Za-z0-9_.]+"
    log:    f"{LOGDIR}pregea/lfmm_K{{k}}_{{trait}}.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm.R \
            {input.lfmm_ld} {input.lfmm_test} {input.climate} \
            {wildcards.k} {wildcards.trait} {input.vcfsnp} \
            {output} FALSE > {log} 2>&1
        """

rule pregea_lfmm_rung_assemble:
    """Wide per-rung p-value table. REUSES scripts/assemble_pvalues.R unchanged."""
    input:  lambda wc: [f"{INTER}pregea/lfmm/K{wc.k}/{t}_pvalues.tsv" for t in PREGEA_PREDICTORS]
    output: f"{INTER}pregea/lfmm/K{{k}}/pvalues.tsv"
    wildcard_constraints: k = r"\d+"
    params: files = lambda wc, input: " ".join(input)
    log:    f"{LOGDIR}pregea/lfmm_assemble_K{{k}}.log"
    shell:  "Rscript /pipeline/scripts/assemble_pvalues.R {params.files} {output} > {log} 2>&1"

# =============================================================================
# Generic rung-statistics rule — serves BOTH ladder engines (one p-value
# table shape, SNPID/chr/pos/<traits>, regardless of engine). FDR is fixed
# at 0.1 for both engines (hit-count series is a trend read, not a threshold
# decision — no longer user-configurable; Bonferroni is computed but never
# plotted downstream, kept for the TSV's completeness only).
# =============================================================================
rule pregea_rung_stats:
    """Per-rung calibration diagnostics for ANY ladder engine.
    scripts/pregea_ladder_stats.R — see its header for the full metric list
    (lambda_gc, hits_bonf/qval, histogram-shape classification)."""
    input:  f"{INTER}pregea/{{engine}}/{{rung}}/pvalues.tsv"
    output: f"{INTER}pregea/{{engine}}/{{rung}}/rung_stats.tsv"
    wildcard_constraints: engine = "lfmm|emmax", rung = r"K\d+|npc\d+"
    params:
        rung_param = lambda wc: "K" if wc.engine == "lfmm" else "n_pcs",
        rung_value = lambda wc: wc.rung.lstrip("Knpc"),
        fdr        = PREGEA_FDR,
        bonf       = 0.05,
        # apply_gc is keyed on ENGINE NAME, not on a config flag: LFMM's rung
        # table is fit with genomic.control=FALSE (Block 1 above) so the raw
        # p-values are informative for lambda/histogram diagnosis; GC
        # recalibration is applied HERE so hit counts reflect what the
        # production genomic.control=TRUE GEA run would produce. EMMAX has no
        # equivalent recalibration step in the production run, so FALSE.
        apply_gc   = lambda wc: "TRUE" if wc.engine == "lfmm" else "FALSE",
    log:    f"{LOGDIR}pregea/rung_stats_{{engine}}_{{rung}}.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_ladder_stats.R \
            {input} {wildcards.engine} {params.rung_param} {params.rung_value} \
            {params.fdr} {params.bonf} {params.apply_gc} {output} > {log} 2>&1
        """

# =============================================================================
# Ladder assembly + plots — one script (plot_pregea_ladder.R), both engines.
# =============================================================================

rule pregea_lfmm_ladder:
    """Assemble the K ladder: histogram grid (THE selection signal), QQ
    grid, lambda-vs-K (informative only because genomic.control=FALSE,
    NOT a selection criterion), hit-count-vs-K at fixed FDR=0.1 (C.2:315)."""
    input:
        stats = [f"{INTER}pregea/lfmm/K{k}/rung_stats.tsv" for k in PREGEA_LFMM_KS],
        pvals = [f"{INTER}pregea/lfmm/K{k}/pvalues.tsv"    for k in PREGEA_LFMM_KS],
    output:
        ladder = O['pregea_lfmm_ladder'],
        hist   = O['pregea_lfmm_hist'],   qq     = O['pregea_lfmm_qq'],
        lam    = O['pregea_lfmm_lambda'], hits   = O['pregea_lfmm_hits'],
    params:
        stats_str = lambda wc, input: " ".join(input.stats),
        pvals_str = lambda wc, input: " ".join(input.pvals),
        values    = PREGEA_LFMM_K_STR, fdr = PREGEA_FDR,
        plot_dir  = f"{MOD_PREGEA}plots/lfmm/", inter_dir = INTER,
    log: f"{LOGDIR}pregea/lfmm_ladder.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pregea_ladder.R \
            lfmm K "{params.stats_str}" "{params.pvals_str}" {params.values} \
            NULL {params.fdr} {output.ladder} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

# =============================================================================
# Block 2 — EMMAX #PC ladder (reuses scripts/emmax.R + scripts/assemble_pvalues.R)
# Sweeps the SAME PC range as RDA's Condition()-PC ladder (PREGEA_NPC_MAX) —
# user decision: "most of the time the number will be equal", one shared range.
# =============================================================================

rule pregea_emmax_rung_trait:
    """One EMMAX rung: n_pcs={n} PCA covariates on top of a FIXED BN kinship,
    one trait, LD-pruned + climate-valid genotypes. n=0 is the kinship-only
    reference panel — load_pca_covariates() already returns NULL for k<=0
    (emmax_core.R, whose comment already names this preGEA ladder).

    INTER_DIR is PER RUNG (not shared): emmax.R writes
    'EMMAX_phenotype_<trait>.tsv' under INTER_DIR — a shared dir would
    make rungs clobber each other's intermediates (silent wrong output,
    not a crash).

    REUSES scripts/emmax.R unchanged (11 positional args, all satisfied)."""
    input:
        vcf        = W['vcf_ld_climate'],
        tped       = W['pregea_tped'],
        kinship    = W['pregea_kinship'],
        traits     = O['climate_site_scaled'],
        covariates = W['pca_projections'],
        metadata   = W['metadata_climate_valid'],
        samples_order = W['samples_order'],
    output: f"{INTER}pregea/emmax/npc{{n}}/{{trait}}_pvalues.tsv"
    wildcard_constraints: n = r"\d+", trait = r"[A-Za-z0-9_.]+"
    params:
        inter_dir   = lambda wc: f"{INTER}pregea/emmax/npc{wc.n}/",
        tped_prefix = f"{WORK_LD}climate/emmax/{VCF_BASE}",
    log: f"{LOGDIR}pregea/emmax_npc{{n}}_{{trait}}.log"
    shell:
        """
        mkdir -p {params.inter_dir}
        Rscript /pipeline/scripts/emmax.R \
            {input.vcf} {wildcards.n} {input.traits} {input.covariates} \
            {wildcards.trait} {params.inter_dir} {input.metadata} \
            {output} {params.tped_prefix} {input.kinship} {input.samples_order} > {log} 2>&1
        """

rule pregea_emmax_rung_assemble:
    """Identical shape to pregea_lfmm_rung_assemble. REUSES
    scripts/assemble_pvalues.R unchanged."""
    input:  lambda wc: [f"{INTER}pregea/emmax/npc{wc.n}/{t}_pvalues.tsv" for t in PREGEA_PREDICTORS]
    output: f"{INTER}pregea/emmax/npc{{n}}/pvalues.tsv"
    wildcard_constraints: n = r"\d+"
    params: files = lambda wc, input: " ".join(input)
    log:    f"{LOGDIR}pregea/emmax_assemble_npc{{n}}.log"
    shell:  "Rscript /pipeline/scripts/assemble_pvalues.R {params.files} {output} > {log} 2>&1"

rule pregea_emmax_ladder:
    """Same script as pregea_lfmm_ladder, engine='emmax'. lambda IS a
    valid criterion here (no forced recalibration, C.0) — the panel gets
    both the 1.0 line and the deflation floor (C.2/C3)."""
    input:
        stats = [f"{INTER}pregea/emmax/npc{n}/rung_stats.tsv" for n in PREGEA_EMMAX_NPCS],
        pvals = [f"{INTER}pregea/emmax/npc{n}/pvalues.tsv"    for n in PREGEA_EMMAX_NPCS],
    output:
        ladder = O['pregea_emmax_ladder'],
        hist   = O['pregea_emmax_hist'],   qq   = O['pregea_emmax_qq'],
        lam    = O['pregea_emmax_lambda'], hits = O['pregea_emmax_hits'],
    params:
        stats_str = lambda wc, input: " ".join(input.stats),
        pvals_str = lambda wc, input: " ".join(input.pvals),
        values    = PREGEA_EMMAX_NPC_STR, floor = PREGEA_DEFLATION_FLOOR,
        fdr       = PREGEA_FDR,
        plot_dir  = f"{MOD_PREGEA}plots/emmax/", inter_dir = INTER,
    log: f"{LOGDIR}pregea/emmax_ladder.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pregea_ladder.R \
            emmax n_pcs "{params.stats_str}" "{params.pvals_str}" {params.values} \
            {params.floor} {params.fdr} {output.ladder} {params.plot_dir} \
            {params.inter_dir} > {log} 2>&1
        """

# =============================================================================
# Block 3 — RDA setup (ONE rule, internal loop over the Condition()-PC ladder)
# =============================================================================

rule pregea_rda_setup:
    """RDA exploratory setup on the LD-pruned set: predictor collinearity
    pre-screen + post-fit VIF (A16, now a REAL gate — see script header);
    a ladder of partial-RDA candidate models varying the Condition() PC
    count (the SAME range as the EMMAX #PC ladder), each producing its own
    biplot + axis screeplot + by-axis/by-term anova.cca — plus a combined
    model-comparison panel across the whole ladder; plain-vs-partial
    candidate-count comparison; and the ordiR2step predictor-selection path
    drawn as cumulative adjusted R2 per step with the full-model ceiling (A17).

    scripts/pregea_rda_setup.R — REUSES scripts/R/lib/rdadapt.R (extracted
    from rda.R), emmax_core.R::load_pca_covariates(), theme_clinego.R.
    Diagnostics only (Rule 5) — does not call find_sig_snps.R or touch any
    GEA output.

    Deliberately NOT a wildcard fan-out: every rung refits the same Y (the
    LD-pruned matrix) with a different Condition() PC count — a fan-out
    would re-read the genotype matrix N times for no caching benefit; the
    per-model output FILES below are just this one rule writing more of
    them, not a Snakemake-level fan-out."""
    input:
        lfmm_pruned   = W['lfmm_imp_climate'],
        vcfsnp_pruned = W['vcfsnp'],
        climate       = O['climate_site_scaled'],
        pca           = W['pca_projections'],
        samples_order = W['samples_order'],
        climate_valid = W['climate_valid_samples'],
    output:
        collin = O['pregea_rda_collin'],  ladder = O['pregea_rda_ladder'],
        axis   = O['pregea_rda_axis'],    fwd    = O['pregea_rda_fwd'],
        p_fwd  = O['pregea_rda_fwd_png'], p_comparison = O['pregea_rda_comparison_png'],
        p_biplots     = [rda_model_biplot(n)      for n in PREGEA_RDA_COND_PCS],
        p_biplots_svg = [rda_model_biplot_svg(n)  for n in PREGEA_RDA_COND_PCS],
        p_screes      = [rda_model_screeplot(n)     for n in PREGEA_RDA_COND_PCS],
        p_screes_svg  = [rda_model_screeplot_svg(n) for n in PREGEA_RDA_COND_PCS],
        t_axis_models = [rda_model_axis_anova(n)    for n in PREGEA_RDA_COND_PCS],
    params:
        predictors = ",".join(PREGEA_PREDICTORS),
        cmin = min(PREGEA_RDA_COND_PCS), cmax = max(PREGEA_RDA_COND_PCS),
        collin_r = PREGEA_RDA_COLLIN_R, vif_max = PREGEA_RDA_VIF_MAX,
        min_preds = PREGEA_RDA_MIN_PREDS, axis_alpha = PREGEA_RDA_AXIS_ALPHA,
        perms = PREGEA_RDA_PERMS, pin = CLIMATE_VP_PIN,
        r2perms = CLIMATE_VP_R2PERMS, seed = PREGEA_SEED,
        plot_dir = f"{MOD_PREGEA}plots/rda/", table_dir = f"{MOD_PREGEA}tables/rda/",
        inter_dir = INTER,
    threads: CPU
    log: f"{LOGDIR}pregea/rda_setup.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_rda_setup.R \
            {input.lfmm_pruned} {input.vcfsnp_pruned} {input.climate} {input.pca} \
            {input.samples_order} {input.climate_valid} {params.predictors} \
            {params.cmin} {params.cmax} {params.collin_r} {params.vif_max} \
            {params.min_preds} {params.axis_alpha} {params.perms} \
            {params.pin} {params.r2perms} {params.seed} {threads} \
            {output.collin} {output.ladder} {output.axis} {output.fwd} \
            {params.plot_dir} {params.table_dir} {params.inter_dir} > {log} 2>&1
        """

# =============================================================================
# Recommendation handoff — the GEA method picker's pre-fill source. See
# scripts/pregea_recommend.R header for why this is NOT wired as a @pregea
# sentinel in gea.py's _resolve_param_default() (parse-time coupling risk).
# =============================================================================
rule pregea_recommend:
    """Emit the per-method hyperparameter recommendation artifact the GEA
    method picker pre-fills from (docs/rda_research.md C.2:340)."""
    input:
        lfmm    = O['pregea_lfmm_ladder'],
        emmax   = O['pregea_emmax_ladder'],
        rda     = O['pregea_rda_ladder'],
        collin  = O['pregea_rda_collin'],
        axis    = O['pregea_rda_axis'],
    output: O['pregea_recommendations']
    params:
        k_best    = K_BEST,
        lambda_tol = PREGEA_LAMBDA_TOL,
        deflation_floor = PREGEA_DEFLATION_FLOOR,
        vif_max = PREGEA_RDA_VIF_MAX,
        axis_alpha = PREGEA_RDA_AXIS_ALPHA,
    log: f"{LOGDIR}pregea/recommend.log"
    shell:
        """
        Rscript /pipeline/scripts/pregea_recommend.R \
            {input.lfmm} {input.emmax} {input.rda} {input.collin} \
            {input.axis} {params.k_best} \
            {params.lambda_tol} {params.deflation_floor} \
            {params.vif_max} {params.axis_alpha} {output} > {log} 2>&1
        """

# =============================================================================
# Transfer guard (opt-in, PreGEA.TransferGuard.enabled) — the ONE preGEA rule
# that pulls the full-set imputation chain (gea.smk) into this otherwise
# LD-pruned-only mode. See scripts/pregea_transfer_guard.R header.
# =============================================================================

if PREGEA_TRANSFER_GUARD:

    rule pregea_transfer_guard:
        """C6/C.2:338 — re-estimate calibration ONCE on the FULL marker set at
        the hyperparameters preGEA recommends, before the user commits to the
        full GEA run. Only the TREND across the ladder transfers; lambda_GC is
        not scale-free (Yang et al. 2011).

        Scope: LFMM + EMMAX only. RDA's full-vs-pruned question is answered
        separately by A6 — its answer is "fit the final scan on the full set",
        not "transfer pruned hyperparameters"; RDA is not folded in here."""
        input:
            lfmm_pruned   = W['lfmm_imp_climate'],
            lfmm_full     = W['lfmm_imp_full_climate'],
            vcfsnp_full   = W['vcfsnp_full'],
            climate       = O['climate_site_scaled'],
            recs          = O['pregea_recommendations'],
            pca           = W['pca_projections'],
            samples_order = W['samples_order'],
            metadata      = W['metadata_climate_valid'],
        output: tsv = O['pregea_transfer_guard'], png = O['pregea_transfer_png']
        params:
            predictors = ",".join(PREGEA_PREDICTORS),
            lfmm_k = PREGEA_TG_LFMM_K,
            # Pruned arm: Block 0's own LD-pruned+climate-valid TPED/kinship
            # (same artifacts the EMMAX ladder itself used).
            emmax_vcf_pruned         = W['vcf_ld_climate'],
            emmax_tped_prefix_pruned = f"{WORK_LD}climate/emmax/{VCF_BASE}",
            emmax_kinship_pruned     = W['pregea_kinship'],
            # Full arm: gea.smk's own full-set climate-valid TPED/kinship chain
            # (assoc_emmax_gea_trait's exact inputs) — independent of whether
            # EMMAX is configured in GEA.configs proper.
            emmax_vcf_full         = W['vcf_filt_climate'] if 'EMMAX' in GEA_OTHER_CONFIGS else 'NULL',
            emmax_tped_prefix_full = f"{WORK_FILT}climate/emmax/{VCF_BASE}" if 'EMMAX' in GEA_OTHER_CONFIGS else 'NULL',
            emmax_kinship_full     = (emmax_kinship_climate_path() if 'EMMAX' in GEA_OTHER_CONFIGS else 'NULL'),
            emmax_npcs = PREGEA_TG_EMMAX_NPCS,
            seed = PREGEA_SEED,
            plot_dir = f"{MOD_PREGEA}plots/transfer/", inter_dir = INTER,
        log: f"{LOGDIR}pregea/transfer_guard.log"
        shell:
            """
            Rscript /pipeline/scripts/pregea_transfer_guard.R \
                {input.lfmm_pruned} {input.lfmm_full} {input.vcfsnp_full} {input.climate} \
                {params.predictors} {params.lfmm_k} {input.recs} \
                {params.emmax_vcf_pruned} {params.emmax_tped_prefix_pruned} {params.emmax_kinship_pruned} \
                {params.emmax_vcf_full} {params.emmax_tped_prefix_full} {params.emmax_kinship_full} \
                {params.emmax_npcs} \
                {input.pca} {input.samples_order} {input.metadata} {params.seed} \
                {output.tsv} {params.plot_dir} {params.inter_dir} > {log} 2>&1
            """
