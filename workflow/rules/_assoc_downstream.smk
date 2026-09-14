# _assoc_downstream.smk
# Phase 4: unified GEA/GWAS downstream rules.
#
# Each rule is parameterized by the {source} wildcard whose value is the
# output directory name ("GEA" or "GWAS"). Source-specific config (module dir,
# predictors, distances, etc.) is looked up via _src(wc.source, key) defined
# in common.smk.
#
# These seven rules replace the following per-source pairs that were deleted:
#   GEA (association.smk)           GWAS (phenotype_assoc.smk)
#   find_sig_snps                   find_sig_snps_pheno
#   combine_selected_snps           combine_selected_snps_pheno
#   create_regions                  create_regions_pheno
#   find_genes_around_regions       find_genes_pheno
#   find_genes_combined_regions     find_genes_combined_pheno
#   manhattan_plot                  manhattan_pheno
#   manhattan_combined              manhattan_combined_pheno
#
# log: uses static wildcard paths (f"{LOGDIR}{{source}}/...") rather than
# lambdas — Snakemake 9 does not support lambda in log: directives.

if ASSOC_SOURCES:

    # Union of all method names across GEA + GWAS, for wildcard constraints.
    _ALL_ASSOC_METHODS_REGEX = "|".join(sorted({
        m for src in ASSOC_SOURCES.values() for m in src["configs"].keys()
    }))

    rule assoc_find_sig_snps:
        """Find significant SNPs for a (source, method, adjust) combination."""
        input:
            assoc = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv"
        output:
            f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}_sig_snps_{{adjust}}.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX,
            method = _ALL_ASSOC_METHODS_REGEX,
            # Alphabetic rule name + numeric value admitting scientific notation
            # (e.g. custom_1e-05 — R's as.character(1e-5)). \w+ was widened to
            # [A-Za-z]+ so digits/underscores can never ambiguously bleed into
            # the name half; find_sig_snps.R still splits on the single '_' into
            # exactly two fields since -/e/E/+ are not '_'.
            adjust = r"[A-Za-z]+_[0-9.eE+-]+"
        params:
            clumping_dist = lambda wc: _src(wc.source, "clumping_distance")
        log:
            f"{LOGDIR}{{source}}/find_sig_snps_{{method}}_{{adjust}}.log"
        threads: CPU
        shell:
            """
            Rscript /pipeline/scripts/find_sig_snps.R \
                {input.assoc} {wildcards.adjust} {params.clumping_dist} \
                {wildcards.method} {threads} {output} > {log} 2>&1
            """

    rule assoc_combine_selected_snps:
        """Combine significant SNPs from all methods for a given source.

        Every argument the script reads positionally is quoted: it checks for 5-6
        argv entries, so an unquoted value that is empty or carries a space (a
        trait column named "disease score" reaches {params.predictors} verbatim)
        would shift the positionals and be refused at exit 1."""
        input:
            sigsnps = lambda wc: [
                f"{OUTDIR}{wc.source}/tables/methods/{m}/{m}_pvalues_K{K_BEST}_sig_snps_{a}.tsv"
                for m, a in _src(wc.source, "configs").items()
            ]
        output:
            snps      = f"{OUTDIR}{{source}}/tables/selected_snps.tsv",
            # per-(SNP, trait) minimum p — the trait-scoped companion of selected_snps.tsv's
            # trait-agnostic min_pvalue; create_regions.R reads it for the per-trait table
            per_trait = f"{OUTDIR}{{source}}/tables/selected_snps_per_trait.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            sigsnps_str   = lambda wc, input: " ".join(input.sigsnps),
            method        = lambda wc: _src(wc.source, "combine_method"),
            clumping_dist = lambda wc: _src(wc.source, "clumping_distance"),
            # combine_predictors (not the plain "predictors" key) so a multivariate
            # method's pseudo-trait (e.g. RDA's climate_multivariate) isn't filtered
            # out of selected_snps.tsv/regions_*/genes_* by combine_sigsnps.R's
            # `trait %in% predictors` check.
            predictors    = lambda wc: _src(wc.source, "combine_predictors")
        log:
            f"{LOGDIR}{{source}}/combine_selected_snps.log"
        shell:
            """
            Rscript /pipeline/scripts/combine_selected_snps.R \
                "{params.sigsnps_str}" "{params.method}" {params.clumping_dist} \
                "{params.predictors}" {output.snps} {output.per_trait} > {log} 2>&1
            """

    rule assoc_create_regions:
        """Merge nearby significant SNPs into per-trait and combined regions."""
        input:
            selected_snps = f"{OUTDIR}{{source}}/tables/selected_snps.tsv",
            trait_pvalues = f"{OUTDIR}{{source}}/tables/selected_snps_per_trait.tsv",
            ld_decay      = lambda wc: ld_decay_input(_src(wc.source, "clumping_distance_mode"))
        output:
            per_trait = f"{OUTDIR}{{source}}/tables/regions_per_trait.tsv",
            combined  = f"{OUTDIR}{{source}}/tables/regions_combined.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            clumping_dist  = lambda wc: _src(wc.source, "clumping_distance"),
            ld_decay_path  = lambda wc: (
                O.get("ld_decay_table", "NULL")
                if _src(wc.source, "clumping_distance_mode") != "fixed" else "NULL"
            ),
            r2_threshold   = lambda wc: _src(wc.source, "clumping_r2_threshold"),
            ld_decay_group = lambda wc: _src(wc.source, "ld_decay_group")
        log:
            f"{LOGDIR}{{source}}/create_regions.log"
        shell:
            """
            Rscript /pipeline/scripts/create_regions.R \
                {input.selected_snps} {params.clumping_dist} \
                {output.per_trait} {output.combined} \
                {params.ld_decay_path} {params.r2_threshold} \
                {params.ld_decay_group} {input.trait_pvalues} > {log} 2>&1
            """

    rule assoc_wza:
        """Compute WZA p-values per window for a (source, method)."""
        input:
            pvalues  = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv",
            maf      = O["qc_maf_pos"],
            ld_decay = lambda wc: wza_ld_decay_input(_src(wc.source, "wza"))
        output:
            f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_wza_K{K_BEST}.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX,
            method = _ALL_ASSOC_METHODS_REGEX
        params:
            window_size      = lambda wc: _src(wc.source, "wza")["window_size"],
            fallback_bp      = lambda wc: _src(wc.source, "wza")["fallback_window_bp"],
            ld_decay_path    = lambda wc: (
                O.get("ld_decay_table", "NULL")
                if _src(wc.source, "wza")["window_size_mode"] != "fixed" else "NULL"
            ),
            ld_decay_group   = lambda wc: _src(wc.source, "ld_decay_group")
        log:
            f"{LOGDIR}{{source}}/wza_{{method}}.log"
        shell:
            """
            Rscript /pipeline/scripts/compute_wza.R \
                {input.pvalues} {input.maf} {params.ld_decay_path} \
                {params.window_size} {params.fallback_bp} {params.ld_decay_group} \
                {output} > {log} 2>&1
            """

    rule assoc_wza_sig_windows:
        """Find significant WZA windows for a (source, method, adjust) combination."""
        input:
            wza = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_wza_K{K_BEST}.tsv"
        output:
            f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_wza_K{K_BEST}_sig_windows_{{adjust}}.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX,
            method = _ALL_ASSOC_METHODS_REGEX,
            adjust = r"[A-Za-z]+_[0-9.eE+-]+"
        params:
            clumping_dist = lambda wc: _src(wc.source, "clumping_distance")
        log:
            f"{LOGDIR}{{source}}/wza_sig_windows_{{method}}_{{adjust}}.log"
        shell:
            """
            Rscript /pipeline/scripts/find_sig_snps.R \
                {input.wza} {wildcards.adjust} {params.clumping_dist} \
                {wildcards.method} 1 {output} --wza > {log} 2>&1
            """

    rule assoc_wza_manhattan:
        """WZA Manhattan + QQ plots for a (source, method, trait, adjust)."""
        input:
            wza = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_wza_K{K_BEST}.tsv"
        output:
            png         = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_wza_{{trait}}_K{K_BEST}_{{adjust}}.png",
            svg         = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_wza_{{trait}}_K{K_BEST}_{{adjust}}.svg",
            qq_png      = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/qq_wza_{{trait}}_K{K_BEST}_{{adjust}}.png",
            qq_svg      = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/qq_wza_{{trait}}_K{K_BEST}_{{adjust}}.svg",
            background  = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_wza_{{trait}}_K{K_BEST}_{{adjust}}_background.png",
            coords_json = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_wza_{{trait}}_K{K_BEST}_{{adjust}}_coords.json"
        wildcard_constraints:
            source = SOURCE_REGEX,
            method = _ALL_ASSOC_METHODS_REGEX,
            trait  = TRAIT_REGEX_ANY,
            adjust = r"[A-Za-z]+_[0-9.eE+-]+"
        params:
            k          = K_BEST,
            plot_dir   = lambda wc: f"{_src(wc.source, 'mod')}plots/manhattan/{wc.method}/",
            # combine_predictors so a multivariate method's pseudo-trait gets a
            # stable, distinct color (positional in get_trait_colors()) instead
            # of falling back to the same color as the first real predictor.
            predictors = lambda wc: _src(wc.source, "combine_predictors")
        log:
            f"{LOGDIR}{{source}}/wza_manhattan_{{method}}_{{trait}}_{{adjust}}.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_manhattan.R \
                {input.wza} {wildcards.adjust} {params.k} {wildcards.method} \
                {wildcards.trait} {params.plot_dir} {params.predictors} wza > {log} 2>&1
            """

    rule assoc_wza_manhattan_combined:
        """Combined WZA Manhattan for all traits/methods of a source.
        wza_methods(source) is the opt-out mechanism (registry supports_wza=False) —
        a method excluded there never had its WZA table built, so it must also be
        excluded here or Snakemake requests a target no rule produces."""
        input:
            wza_tables = lambda wc: [
                f"{OUTDIR}{wc.source}/tables/methods/{m}/{m}_wza_K{K_BEST}.tsv"
                for m in wza_methods(wc.source)
            ]
        output:
            simple_png  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}.png",
            simple_svg  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}.svg",
            qq_png      = f"{OUTDIR}{{source}}/plots/manhattan/combined/qq_wza_combined_K{K_BEST}.png",
            qq_svg      = f"{OUTDIR}{{source}}/plots/manhattan/combined/qq_wza_combined_K{K_BEST}.svg",
            background  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}_background.png",
            coords_json = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_wza_combined_K{K_BEST}_coords.json"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            assoc_str  = lambda wc: ",".join([
                f"{m}:{_src(wc.source, 'configs')[m]}:{OUTDIR}{wc.source}/tables/methods/{m}/{m}_wza_K{K_BEST}.tsv"
                for m in wza_methods(wc.source)
            ]),
            predictors = lambda wc: _src(wc.source, "combine_predictors"),
            k          = K_BEST,
            plot_dir   = lambda wc: f"{_src(wc.source, 'mod')}plots/manhattan/combined/"
        log:
            f"{LOGDIR}{{source}}/wza_manhattan_combined.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_manhattan_combined.R \
                "{params.assoc_str}" {params.predictors} {params.k} \
                {params.plot_dir} wza > {log} 2>&1
            touch {output.simple_png} {output.simple_svg} {output.qq_png} {output.qq_svg} \
                  {output.background} {output.coords_json}
            """

    rule assoc_find_genes_per_region:
        """Find genes overlapping per-trait regions for a given source."""
        input:
            regions = f"{OUTDIR}{{source}}/tables/regions_per_trait.tsv",
            gff     = W["gff_normalized"],
            vcfsnp  = W["vcfsnp_full"]
        output:
            genes     = f"{OUTDIR}{{source}}/tables/genes_per_region.tsv",
            collapsed = f"{OUTDIR}{{source}}/tables/genes_per_region_collapsed.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            feature      = GFF_FEATURE,
            promoter_len = lambda wc: _src(wc.source, "promoter_length")
        log:
            f"{LOGDIR}{{source}}/find_genes_per_region.log"
        threads: CPU
        shell:
            """
            Rscript /pipeline/scripts/find_genes_around_regions.R \
                {input.gff} {input.regions} {params.feature} \
                {params.promoter_len} {input.vcfsnp} {threads} \
                {output.genes} {output.collapsed} > {log} 2>&1
            """

    rule assoc_find_genes_combined:
        """Find genes overlapping combined regions for a given source."""
        input:
            regions = f"{OUTDIR}{{source}}/tables/regions_combined.tsv",
            gff     = W["gff_normalized"],
            vcfsnp  = W["vcfsnp_full"]
        output:
            genes = f"{OUTDIR}{{source}}/tables/genes_combined.tsv"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            feature      = GFF_FEATURE,
            promoter_len = lambda wc: _src(wc.source, "promoter_length")
        log:
            f"{LOGDIR}{{source}}/find_genes_combined.log"
        threads: CPU
        shell:
            """
            TEMP_COLLAPSED=$(mktemp)
            Rscript /pipeline/scripts/find_genes_around_regions.R \
                {input.gff} {input.regions} {params.feature} \
                {params.promoter_len} {input.vcfsnp} {threads} \
                {output.genes} $TEMP_COLLAPSED > {log} 2>&1
            rm -f $TEMP_COLLAPSED
            """

    rule assoc_manhattan_plot:
        """Per-(source, method, trait, adjust) Manhattan + QQ plots.
        Background PNG and coords include {adjust} in filename (Snakemake wildcard
        consistency), but the background shows ALL SNPs (threshold-independent).
        Shiny draws the threshold line live via the interactive threshold_y param."""
        input:
            assoc = f"{OUTDIR}{{source}}/tables/methods/{{method}}/{{method}}_pvalues_K{K_BEST}.tsv"
        output:
            png         = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.png",
            svg         = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}.svg",
            qq_png      = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/qq_{{trait}}_K{K_BEST}_{{adjust}}.png",
            qq_svg      = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/qq_{{trait}}_K{K_BEST}_{{adjust}}.svg",
            background  = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_background.png",
            coords_json = f"{OUTDIR}{{source}}/plots/manhattan/{{method}}/manhattan_{{trait}}_K{K_BEST}_{{adjust}}_coords.json"
        wildcard_constraints:
            source = SOURCE_REGEX,
            method = _ALL_ASSOC_METHODS_REGEX,
            trait  = TRAIT_REGEX_ANY,
            adjust = r"[A-Za-z]+_[0-9.eE+-]+"
        params:
            k          = K_BEST,
            plot_dir   = lambda wc: f"{_src(wc.source, 'mod')}plots/manhattan/{wc.method}/",
            predictors = lambda wc: _src(wc.source, "combine_predictors")
        log:
            f"{LOGDIR}{{source}}/manhattan_{{method}}_{{trait}}_{{adjust}}.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_manhattan.R \
                {input.assoc} {wildcards.adjust} {params.k} {wildcards.method} \
                {wildcards.trait} {params.plot_dir} {params.predictors} > {log} 2>&1
            """

    rule assoc_manhattan_combined:
        """Combined Manhattan + QQ plots for all traits/methods of a source."""
        input:
            assoc_tables = lambda wc: [
                f"{OUTDIR}{wc.source}/tables/methods/{m}/{m}_pvalues_K{K_BEST}.tsv"
                for m in _src(wc.source, "configs")
            ]
        output:
            simple_png  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_combined_K{K_BEST}.png",
            simple_svg  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_combined_K{K_BEST}.svg",
            qq_png      = f"{OUTDIR}{{source}}/plots/manhattan/combined/qq_combined_K{K_BEST}.png",
            qq_svg      = f"{OUTDIR}{{source}}/plots/manhattan/combined/qq_combined_K{K_BEST}.svg",
            background  = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_combined_K{K_BEST}_background.png",
            coords_json = f"{OUTDIR}{{source}}/plots/manhattan/combined/manhattan_combined_K{K_BEST}_coords.json"
        wildcard_constraints:
            source = SOURCE_REGEX
        params:
            assoc_str  = lambda wc: ",".join([
                f"{m}:{a}:{OUTDIR}{wc.source}/tables/methods/{m}/{m}_pvalues_K{K_BEST}.tsv"
                for m, a in _src(wc.source, "configs").items()
            ]),
            # combine_predictors — without it, a multivariate method's pseudo-trait
            # is absent from the traits list load_assoc_data() requests, so its
            # signal never appears in the "all methods" combined view at all.
            predictors = lambda wc: _src(wc.source, "combine_predictors"),
            k          = K_BEST,
            plot_dir   = lambda wc: f"{_src(wc.source, 'mod')}plots/manhattan/combined/"
        log:
            f"{LOGDIR}{{source}}/manhattan_combined.log"
        shell:
            """
            Rscript /pipeline/scripts/plot_manhattan_combined.R \
                "{params.assoc_str}" {params.predictors} {params.k} \
                {params.plot_dir} > {log} 2>&1
            touch {output.simple_png} {output.simple_svg} {output.qq_png} {output.qq_svg} \
                  {output.background} {output.coords_json}
            """
