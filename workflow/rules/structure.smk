#=============================================================================
# MODULE 3: STRUCTURE_K (requires K_BEST selection)
#=============================================================================

rule impute_ld:
    """Impute missing genotypes using SNMF results for K_BEST."""
    input:  snmf = W['snmf'], lfmm = W['lfmm'], nmissing = W['lfmm_nmissing']
    output: W['lfmm_imp']
    params: k = K_BEST
    log:    f"{LOGDIR}structure/impute_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/impute.R {input.snmf} {input.lfmm} {params.k} {output} {input.nmissing} > {log} 2>&1
        """

rule lfmm2vcf_ld:
    """Convert imputed LFMM back to VCF format."""
    input:  vcf = W['vcf_ld'], lfmm_imp = W['lfmm_imp']
    output: W['vcf_imp']
    params: blocksize = 20000
    log:    f"{LOGDIR}structure/lfmm2vcf_ld.log"
    shell:
        """
        Rscript /pipeline/scripts/lfmm2vcf.R {input.lfmm_imp} {input.vcf} {params.blocksize} > {log} 2>&1
        """
    # NOTE: lfmm2vcf.R already exists, outputs .vcf with same base as .lfmm input

if CLIMATE_SOURCE == 'worldclim':
    rule download_climate_present:
        """Download and process present climate data for sampling locations.
        Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples).
        Always emits climate_na_excluded.tsv (empty when clean) diagnosing any sample whose
        raster extraction returned NA (e.g. ocean/NoData pixel) and always excludes those
        samples from climate-VALUE-dependent steps only (see filter_climate_valid_samples
        below) -- never a hard stop; the user can fix coordinates or ignore it."""
        input:  meta = W['metadata_climate']
        output:
            site = O['climate_site'],
            site_scaled = O['climate_site_scaled'],
            all_vals = O['climate_all'],
            raster = W['climate_raster'],
            excluded = O['climate_na_excluded']
        params:
            crop = CLIMATE_EXTENT,
            gap = GAP,
            resolution = RESOLUTION,
            data_dir = f"{INDIR}",
            raster_dir = f"{MOD_CLIMATE}rasters/present/",
            tables_dir = f"{MOD_CLIMATE}tables/present/"
        log: f"{LOGDIR}structure/download_climate_present.log"
        shell:
            """
            Rscript /pipeline/scripts/download_climate_present.R \
                {input.meta} {params.crop} {params.gap} {params.data_dir} {params.resolution} \
                {params.raster_dir} {params.tables_dir} > {log} 2>&1
            """

    rule filter_climate_valid_samples:
        """Further narrow coord-valid samples to exclude any whose raster climate extraction
        returned NA (see download_climate_present). No-op (climate-valid == coord-valid) when
        climate_na_excluded.tsv is empty -- the common case, no samples excluded."""
        input:
            coord_samples = W['coord_valid_samples'],
            meta          = W['metadata_climate'],
            excluded      = O['climate_na_excluded']
        output:
            samples_list = W['climate_valid_samples'],
            metadata     = W['metadata_climate_valid']
        log: f"{LOGDIR}structure/filter_climate_valid_samples.log"
        shell:
            """
            Rscript /pipeline/scripts/filter_climate_valid_samples.R \
                {input.coord_samples} {input.meta} {input.excluded} \
                {output.samples_list} {output.metadata} > {log} 2>&1
            """
else:
    rule stage_custom_climate_present:
        """Stage user-supplied present climate table into pipeline-standard outputs.
        Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples)."""
        input:  meta = W['metadata_climate']
        output:
            site        = O['climate_site'],
            site_scaled = O['climate_site_scaled'],
            all_vals    = O['climate_all'],
            raster      = W['climate_raster']
        params:
            env_table  = CUSTOM_PRESENT_TABLE,
            columns    = CUSTOM_CLIMATE_COLUMNS,
            key        = CUSTOM_CLIMATE_KEY,
            grid_res   = CUSTOM_GRID_RES,
            raster_dir = f"{MOD_CLIMATE}rasters/present/",
            tables_dir = f"{MOD_CLIMATE}tables/present/"
        log: f"{LOGDIR}structure/stage_custom_climate_present.log"
        shell:
            """
            Rscript /pipeline/scripts/stage_custom_climate.R present \
                {input.meta} {params.env_table} {params.columns} {params.key} {params.grid_res} \
                {output.raster} {output.all_vals} {output.site} {output.site_scaled} > {log} 2>&1
            """

rule tajima_d:
    """Calculate Tajima's D per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['tajima']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}structure/tajima_d.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/tajima_d.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule pi_diversity:
    """Calculate nucleotide diversity (Pi) per population."""
    input:  vcf = W['vcf_filt'], meta = O['metadata']
    output: O['pi_div']
    params: winsize = METRICS_WINSIZE, inter_dir = INTER
    log:    f"{LOGDIR}structure/pi_diversity.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/pi_diversity.R \
            {input.vcf} {input.meta} {params.winsize} {threads} \
            {output} {params.inter_dir} > {log} 2>&1
        """

# `rule ibd` (scripts/ibd.R, ibd_raw.tsv / ibd_pairs.tsv) was removed 2026-09-13: a
# per-site-pair Mantel test with a two-valued geographic matrix is a differentiation
# contrast, not an isolation-by-distance test, and it silently dropped any site with a
# coordinate-missing sample (audit B3/SC5). IBD/IBE is mantel_test below, at site level.

rule mantel_test:
    """Perform Mantel test for IBD/IBE, at SITE level (mantel_test.R collapses
    samples to one row per site before scaling and computing distances -- geography
    and climate are site properties, and per-sample rows both weight scale()'s
    centre/sd by sampling effort and feed Mantel's permutation test pseudoreplicated
    units). Needs >= 4 sites; mantel_test.R stops with a clear message below that.
    Uses metadata_climate_valid.tsv (climate-valid samples, see filter_climate_valid_samples) --
    mantel_test.R's own complete.cases(env) guard drops geo/clust in lockstep with env, which
    requires meta and climate to already have the same row count/order."""
    input:  meta = W['metadata_climate_valid'], climate = O['climate_site'], clusters = clusters_table(K_BEST)
    output: plot = O['mantel'], stats = O['mantel_stats']
    params:
        predictors = ALL_BIO_STR,
        plot_dir = f"{MOD_STRUCT}plots/pop_stats/",
        inter_dir = INTER,
        tables_dir = f"{MOD_STRUCT}tables/pop_stats/"
    log: f"{LOGDIR}structure/mantel_test.log"
    shell:
        """
        Rscript /pipeline/scripts/mantel_test.R \
            {input.meta} {input.clusters} {input.climate} {params.predictors} \
            {params.plot_dir} {params.inter_dir} {params.tables_dir} > {log} 2>&1
        """

rule amova:
    """Perform AMOVA analysis."""
    input:  vcf = W['vcf_imp'], meta = O['metadata']
    output: table = O['amova'], plot = O['amova_plot']
    params: plot_dir = f"{MOD_STRUCT}plots/pop_stats/", tables_dir = f"{MOD_STRUCT}tables/pop_stats/", inter_dir = INTER
    log:    f"{LOGDIR}structure/amova.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/amova.R \
            {input.vcf} {input.meta} {threads} \
            {params.plot_dir} {params.tables_dir} {params.inter_dir} > {log} 2>&1
        """

rule piemap_plot:
    """Generate PieMap plots for a climate predictor with trait-scaled pie sizes.
    Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples)."""
    input:
        meta = W['metadata_climate'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster'],
        tajima = O['tajima'],
        diversity = O['pi_div']
    output:
        tajima_plot = piemap_tajima("{bio}"),
        diversity_plot = piemap_diversity("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        plot_dir = f"{MOD_STRUCT}plots/piemap/",
        inter_dir = INTER,
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}structure/piemap_{{bio}}.log"
    shell:
        """
        # TajimaD piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.tajima} "Tajima's D" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio}_tajima_d {params.regionmap_extent} {params.pie_scale} Clusters > {log} 2>&1

        # PiDiversity piemap
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            {input.diversity} "Pi Diversity" \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio}_pi_diversity {params.regionmap_extent} {params.pie_scale} Clusters >> {log} 2>&1
        """

rule piemap_simple:
    """Generate basic PieMap plots with autoscaled uniform pie size (no trait overlay).
    Uses metadata_climate.tsv (coord-valid samples only, see filter_coord_samples)."""
    input:
        meta = W['metadata_climate'],
        clusters = clusters_table(K_BEST),
        raster = W['climate_raster']
    output:
        pie = piemap_notrait("{bio}"),
        points = piemap_notrait_points("{bio}")
    wildcard_constraints: bio = r"bio_\d+"
    params:
        bio = lambda wc: wc.bio,
        pie_alpha = PIEMAP_ALPHA,
        pop_label = PIEMAP_SHOW_LABELS,
        pop_label_size = PIEMAP_LABEL_SIZE,
        pie_scale = PIEMAP_PIE_SCALE,
        plot_dir = f"{MOD_STRUCT}plots/piemap/",
        inter_dir = INTER,
        regionmap_extent = REGIONMAP_EXTENT
    log: f"{LOGDIR}structure/piemap_simple_{{bio}}.log"
    shell:
        """
        # Pie chart (main render)
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio} {params.regionmap_extent} {params.pie_scale} Clusters > {log} 2>&1

        # Points ("clear map") companion — same background, tiny dots instead of pies
        Rscript /pipeline/scripts/plot_piemap.R \
            {input.raster} {params.bio} {params.bio} \
            {input.meta} {input.clusters} \
            NULL NULL \
            {params.pie_alpha} {params.pop_label} {params.pop_label_size} \
            {params.plot_dir} {params.inter_dir} \
            piemap_{params.bio}_points {params.regionmap_extent} {params.pie_scale} Clusters T >> {log} 2>&1
        """

#=============================================================================
# LD DECAY ANALYSIS
#=============================================================================

rule ld_decay_prepare:
    """Create sample lists per group and split VCF by chromosome if needed."""
    input:
        vcf = W['vcf_filt'],
        meta = O['metadata'],
        clusters = clusters_table(K_BEST) if LD_DECAY_GROUP_BY == 'cluster' else []
    output:
        flag = touch(W['ld_decay_prep_done'])
    params:
        group_by = LD_DECAY_GROUP_BY,
        min_samples = LD_DECAY_MIN_SAMPLES,
        scope = LD_DECAY_SCOPE,
        sample_lists_dir = W['ld_decay_sample_lists'],
        chr_vcfs_dir = W['ld_decay_chr_vcfs'],
        clusters_path = clusters_table(K_BEST) if LD_DECAY_GROUP_BY == 'cluster' else 'NULL'
    log: f"{LOGDIR}structure/ld_decay_prepare.log"
    shell:
        """
        Rscript /pipeline/scripts/ld_decay_prepare.R \
            {input.vcf} {input.meta} {params.group_by} \
            {params.min_samples} {params.scope} \
            {params.sample_lists_dir} {params.chr_vcfs_dir} \
            {params.clusters_path} > {log} 2>&1
        """

rule ld_decay_run_gw:
    """Run PopLDdecay genome-wide for each group + All."""
    input:
        vcf = W['vcf_filt'],
        prep = W['ld_decay_prep_done']
    output:
        flag = touch(W['ld_decay_gw_done'])
    params:
        max_dist = LD_DECAY_MAX_DISTANCE,
        sample_lists_dir = W['ld_decay_sample_lists'],
        stat_dir = W['ld_decay_stat_gw']
    log: f"{LOGDIR}structure/ld_decay_run_gw.log"
    shell:
        """
        mkdir -p {params.stat_dir}
        LDFLAGS="-MAF 0 -Miss 1 -Het 1"

        # Run for "All" (no -SubPop)
        PopLDdecay -InVCF {input.vcf} \
            -OutStat {params.stat_dir}/All \
            -MaxDist {params.max_dist} $LDFLAGS \
            >> {log} 2>&1

        # Run per group
        for sample_list in {params.sample_lists_dir}/*.txt; do
            group=$(basename "$sample_list" .txt)
            [ "$group" = "manifest" ] && continue
            [ "$group" = "All" ] && continue
            PopLDdecay -InVCF {input.vcf} \
                -OutStat {params.stat_dir}/"$group" \
                -MaxDist {params.max_dist} $LDFLAGS \
                -SubPop "$sample_list" \
                >> {log} 2>&1
        done
        """

if LD_DECAY_SCOPE in ('per_chromosome', 'both'):
    rule ld_decay_run_chr:
        """Run PopLDdecay per chromosome for each group."""
        input:
            prep = W['ld_decay_prep_done']
        output:
            flag = touch(W['ld_decay_chr_done'])
        params:
            max_dist = LD_DECAY_MAX_DISTANCE,
            sample_lists_dir = W['ld_decay_sample_lists'],
            chr_vcfs_dir = W['ld_decay_chr_vcfs'],
            stat_dir = W['ld_decay_stat_chr']
        log: f"{LOGDIR}structure/ld_decay_run_chr.log"
        shell:
            """
            mkdir -p {params.stat_dir}
            LDFLAGS="-MAF 0 -Miss 1 -Het 1"

            while IFS= read -r chr; do
                # Run for "All" (no -SubPop)
                PopLDdecay -InVCF {params.chr_vcfs_dir}/"$chr".vcf \
                    -OutStat {params.stat_dir}/All_"$chr" \
                    -MaxDist {params.max_dist} $LDFLAGS \
                    >> {log} 2>&1

                # Run per group
                for sample_list in {params.sample_lists_dir}/*.txt; do
                    group=$(basename "$sample_list" .txt)
                    [ "$group" = "manifest" ] && continue
                    [ "$group" = "All" ] && continue
                    PopLDdecay -InVCF {params.chr_vcfs_dir}/"$chr".vcf \
                        -OutStat {params.stat_dir}/"$group"_"$chr" \
                        -MaxDist {params.max_dist} $LDFLAGS \
                        -SubPop "$sample_list" \
                        >> {log} 2>&1
                done
            done < {params.chr_vcfs_dir}/chromosomes.txt
            """

rule ld_decay_analyze:
    """Read PopLDdecay outputs, fit Hill-Weir curves, extract half-decay distances, plot."""
    input:
        gw = W['ld_decay_gw_done'] if LD_DECAY_SCOPE in ('genome_wide', 'both') else [],
        chr = W['ld_decay_chr_done'] if LD_DECAY_SCOPE in ('per_chromosome', 'both') else []
    output:
        table = O['ld_decay_table'],
        plot_gw = O['ld_decay_plot_gw'] if LD_DECAY_SCOPE in ('genome_wide', 'both') else [],
        plot_gw_svg = O['ld_decay_plot_gw_svg'] if LD_DECAY_SCOPE in ('genome_wide', 'both') else [],
        plot_chr = O['ld_decay_plot_chr'] if LD_DECAY_SCOPE in ('per_chromosome', 'both') else [],
        plot_chr_svg = O['ld_decay_plot_chr_svg'] if LD_DECAY_SCOPE in ('per_chromosome', 'both') else []
    params:
        scope = LD_DECAY_SCOPE,
        stat_gw_dir = W['ld_decay_stat_gw'],
        stat_chr_dir = W['ld_decay_stat_chr'],
        manifest = f"{W['ld_decay_sample_lists']}manifest.tsv",
        chr_list = f"{W['ld_decay_chr_vcfs']}chromosomes.txt",
        plot_gw_path = O.get('ld_decay_plot_gw', 'NULL'),
        plot_chr_path = O.get('ld_decay_plot_chr', 'NULL')
    log: f"{LOGDIR}structure/ld_decay_analyze.log"
    shell:
        """
        Rscript /pipeline/scripts/ld_decay_analyze.R \
            {params.scope} {params.stat_gw_dir} {params.stat_chr_dir} \
            {params.manifest} {params.chr_list} \
            {output.table} {params.plot_gw_path} {params.plot_chr_path} \
            > {log} 2>&1
        """
