#=============================================================================
# MODULE 1: VCF PROCESSING + QC
#=============================================================================

rule raw_vcf_stats:
    """Run bcftools stats on raw VCF; output raw stats file parsed later by R."""
    input:  vcf = f"{INDIR}{VCF_RAW}"
    output: O['qc_raw_summary']
    log:    f"{LOGDIR}processing/raw_vcf_stats.log"
    shell:  "bcftools stats {input.vcf} > {output} 2> {log}"

rule extract_samples:
    """Extract sample IDs from metadata for VCF subsetting."""
    input:  samples = f"{INDIR}{SAMPLES}"
    output: W['samples_list']
    log:    f"{LOGDIR}processing/extract_samples.log"
    # -F'\t' is load-bearing: the metadata is a TSV and awk's default FS splits on ANY run
    # of blanks, so a site name with a space ("Tel Aviv") made $2 the second WORD ("Aviv"),
    # plink --keep silently dropped every sample of that site, and the QC accounting still
    # reported 100 % retention (audit 2026-09-13 B41).
    shell:  "tail -n +2 {input.samples} | awk -F'\\t' '{{print 0, $2}}' > {output} 2> {log}"

rule calculate_sample_missing:
    """Calculate per-sample and per-SNP missing genotype rates. Filter samples by sample_miss threshold."""
    input:
        vcf = f"{INDIR}{VCF_RAW}",
        samples = W['samples_list']
    output:
        stats    = W['samples_missing_stats'],
        filtered = W['samples_filtered'],
        removed  = W['samples_removed'],
        lmiss    = O['qc_snp_miss_raw']
    params:
        threshold = SAMPLE_MISS,
        prefix    = f"{INTER}samples/sample_miss_tmp"
    log: f"{LOGDIR}processing/calculate_sample_missing.log"
    threads: CPU
    shell:
        """
        # Calculate per-sample and per-SNP missingness using plink
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --missing --out {params.prefix} > {log} 2>&1

        # Every ID in the keep-list must have survived --keep. plink drops IDs it cannot find
        # in the VCF with exit 0 ("--keep: N people remaining"), which is silent sample loss;
        # downstream counts are all taken from the .imiss, so nothing later can notice.
        n_keep=$(wc -l < {input.samples})
        n_imiss=$(($(wc -l < {params.prefix}.imiss) - 1))
        if [ "$n_keep" -ne "$n_imiss" ]; then
            echo "ERROR: keep-list has $n_keep samples but plink --keep retained $n_imiss." >> {log}
            echo "ERROR: metadata sample IDs not found in the VCF header:" >> {log}
            awk 'NR==FNR {{seen[$2]=1; next}} !($2 in seen) {{print "  " $2}}' {params.prefix}.imiss {input.samples} >> {log}
            grep -m1 -- '--keep:' {params.prefix}.log >> {log} || true
            exit 1
        fi

        # Create stats file with header (convert plink space-separated to TSV)
        echo -e "FID\tIID\tMISS_PHENO\tN_MISS\tN_GENO\tF_MISS" > {output.stats}
        tail -n +2 {params.prefix}.imiss | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' >> {output.stats}

        # Preserve per-SNP missingness for QC (convert plink space-sep to TSV)
        awk 'NR==1 {{print "CHR\tSNP\tN_MISS\tN_GENO\tF_MISS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5}}' \
            {params.prefix}.lmiss > {output.lmiss}

        # Filter samples: keep those with F_MISS <= threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 <= thresh {{print $1, $2}}' {params.prefix}.imiss > {output.filtered}

        # List removed samples: those with F_MISS > threshold
        awk -v thresh={params.threshold} 'NR>1 && $6 > thresh {{print $1, $2, $6}}' {params.prefix}.imiss > {output.removed}

        # Log summary
        n_total=$(wc -l < {input.samples})
        n_kept=$(wc -l < {output.filtered})
        n_removed=$(wc -l < {output.removed})
        echo "INFO: Sample missingness filtering (threshold: {params.threshold})" >> {log}
        echo "INFO: Total samples: $n_total" >> {log}
        echo "INFO: Samples passing: $n_kept" >> {log}
        echo "INFO: Samples removed: $n_removed" >> {log}

        # Cleanup temp files (keep .imiss and .lmiss — moved above)
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.imiss {params.prefix}.lmiss
        """

rule compute_het_raw:
    """Compute per-sample heterozygosity on the raw VCF (filtered samples).
    Used for het_vs_missingness plot and optional het outlier removal."""
    input:
        vcf     = f"{INDIR}{VCF_RAW}",
        samples = W['samples_filtered']
    output: W['het_stats_raw']
    params: prefix = f"{INTER}qc/het_raw_tmp"
    log:    f"{LOGDIR}processing/compute_het_raw.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --het --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "FID\tIID\tO_HOM\tE_HOM\tN_NM\tF"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.het > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.het
        """

if HET_OUTLIER_SD is not None:
    rule het_outlier_filter:
        """Remove samples with heterozygosity beyond het_outlier_sd SDs from mean.
        Produces updated samples list used by filter_vcf."""
        input:
            het     = W['het_stats_raw'],
            samples = W['samples_filtered']
        output: W['samples_het_filtered']
        params: sd_thresh = HET_OUTLIER_SD
        log:    f"{LOGDIR}processing/het_outlier_filter.log"
        shell:
            """
            # Compute obs_het per sample (N_NM - O_HOM) / N_NM, then filter beyond N SDs
            awk -v sd_thresh={params.sd_thresh} '
            BEGIN {{ OFS="\t" }}
            NR == 1 {{ next }}
            {{
                fid=$1; iid=$2; o_hom=$3; e_hom=$4; n_nm=$5
                het = (n_nm - o_hom) / n_nm
                samples[NR] = fid " " iid
                het_vals[NR] = het
                sum += het; count++
            }}
            END {{
                mean = sum / count
                for (i in het_vals) sq_sum += (het_vals[i] - mean)^2
                sd = sqrt(sq_sum / count)
                lo = mean - sd_thresh * sd
                hi = mean + sd_thresh * sd
                n_removed = 0
                for (i in samples) {{
                    if (het_vals[i] >= lo && het_vals[i] <= hi)
                        print samples[i]
                    else
                        n_removed++
                }}
                print "INFO: Het outlier filter: mean=" mean " sd=" sd " lo=" lo " hi=" hi > "/dev/stderr"
                print "INFO: Removed " n_removed " samples as het outliers" > "/dev/stderr"
            }}' {input.het} > {output} 2>> {log}
            """

rule compute_relatedness:
    """Compute pairwise IBS allele-sharing (plink --genome, DST) on LD-pruned genotypes
    (filtered/het-filtered samples). Always runs — feeds the relatedness QC histogram so the
    user can pick a sensible threshold. Sample removal only happens when Filter.relatedness is
    set AND Filter.relatedness_action is 'remove' (see relatedness_filter below; default action
    is 'keep' — visualize/count only).
    NOTE: IBS is model-free (no Hardy-Weinberg/outbred assumption), so it works for both
    high-selfing and outcrossing species — unlike pi-hat, which collapses under high selfing.
    Selfers show a higher IBS baseline overall; duplicates/clones cluster near 1.0."""
    input:
        vcf     = f"{INDIR}{VCF_RAW}",
        samples = W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered']
    output:
        genome = W['relatedness_genome'],
        pairs  = O['qc_relatedness_pairs']
    params:
        prefix = f"{INTER}qc/relatedness_tmp", win = LD_WIN, step = LD_STEP, r2 = LD_R2
    log:    f"{LOGDIR}processing/compute_relatedness.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --indep-pairwise {params.win} {params.step} {params.r2} \
            --out {params.prefix} > {log} 2>&1

        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --extract {params.prefix}.prune.in \
            --genome --out {params.prefix} >> {log} 2>&1

        echo -e "IID1\tIID2\tIBS" > {output.pairs}
        tail -n +2 {params.prefix}.genome | awk '{{print $2"\t"$4"\t"$12}}' >> {output.pairs}
        cp {params.prefix}.genome {output.genome}

        n_pairs=$(tail -n +2 {output.pairs} | wc -l)
        echo "WARNING: Relatedness (IBS allele-sharing) computed for $n_pairs sample pairs. Set Filter.relatedness from the histogram (Processing QC) — duplicates/clones cluster near 1.0." >> {log}

        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.prune.in {params.prefix}.prune.out {params.prefix}.genome
        """

if PI_HAT is not None:
    rule relatedness_filter:
        """Compute which related samples (IBS allele-sharing) would be dropped above the configured
        threshold. Always runs when Filter.relatedness is set (even in 'keep' mode) so the Shiny
        hover note can preview the removal count. For each over-threshold pair, drops the member
        with higher missingness (simple pairwise rule — does not attempt maximal-set optimization).
        The resulting samples list is only fed into filter_vcf when Filter.relatedness_action
        is 'remove' (see RELATEDNESS_REMOVE gate below) — in 'keep' mode nothing is excluded."""
        input:
            pairs   = O['qc_relatedness_pairs'],
            imiss   = W['samples_missing_stats'],
            samples = W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered']
        output:
            filtered = W['samples_rel_filtered'],
            removed  = O['qc_relatedness_removed']
        params:
            thresh    = PI_HAT,
            tmp_imiss = f"{INTER}samples/rel_tmp_imiss.tsv",
            tmp_pairs = f"{INTER}samples/rel_tmp_sorted_pairs.tsv"
        log:    f"{LOGDIR}processing/relatedness_filter.log"
        shell:
            """
            tail -n +2 {input.imiss} > {params.tmp_imiss}
            tail -n +2 {input.pairs} | awk -v thresh={params.thresh} '$3 > thresh' | sort -k3,3gr > {params.tmp_pairs}

            awk -v thresh={params.thresh} -v removed_out={output.removed} '
            BEGIN {{ OFS="\t"; file_num = 0; print "IID\tIBS\tF_MISS" > removed_out }}
            FNR==1 {{ file_num++ }}
            file_num==1 {{ fmiss[$2]=$6; next }}
            file_num==2 {{ keep[$2]=1; fid[$2]=$1; order[++n]=$2; next }}
            file_num==3 {{
                iid1=$1; iid2=$2; ibs=$3
                if (keep[iid1]==1 && keep[iid2]==1) {{
                    if (fmiss[iid1]+0 >= fmiss[iid2]+0) {{ drop=iid1 }} else {{ drop=iid2 }}
                    keep[drop]=0
                    print drop, ibs, fmiss[drop] >> removed_out
                }}
            }}
            END {{
                n_removed = 0
                for (i=1;i<=n;i++) {{
                    iid=order[i]
                    if (keep[iid]==1) print fid[iid], iid
                    else n_removed++
                }}
                print "WARNING: Relatedness filter (IBS>" thresh "): removed " n_removed " related samples" > "/dev/stderr"
            }}
            ' {params.tmp_imiss} {input.samples} {params.tmp_pairs} > {output.filtered} 2>> {log}

            rm -f {params.tmp_imiss} {params.tmp_pairs}
            """

rule compute_snp_freq_raw:
    """Compute allele frequencies on the raw VCF (filtered samples, before MAF filter)."""
    input:
        vcf     = f"{INDIR}{VCF_RAW}",
        samples = W['samples_rel_filtered'] if RELATEDNESS_REMOVE else (W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered'])
    output: O['qc_maf_raw']
    params: prefix = f"{INTER}qc/maf_raw_tmp"
    log:    f"{LOGDIR}processing/compute_snp_freq_raw.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --freq --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.frq > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.frq
        """

rule compute_snp_density_raw:
    """Compute SNP density in 1Mb windows on the raw VCF."""
    input: f"{INDIR}{VCF_RAW}"
    output: O['qc_snp_density_raw']
    params: prefix = f"{INTER}qc/snp_density_raw_tmp"
    log:    f"{LOGDIR}processing/compute_snp_density_raw.log"
    shell:
        """
        case "{input}" in
            *.gz) VCF_FLAG=--gzvcf ;;
            *)    VCF_FLAG=--vcf ;;
        esac
        /app/vcftools-0.1.16/bin/vcftools $VCF_FLAG {input} --SNPdensity 1000000 --out {params.prefix} > {log} 2>&1
        cp {params.prefix}.snpden {output}
        rm -f {params.prefix}.snpden {params.prefix}.log
        """

if (MIN_DEPTH is not None or MAX_DEPTH is not None) and HAS_FORMAT_DP:
    _depth_expr_parts = []
    if MIN_DEPTH is not None:
        _depth_expr_parts.append(f"FMT/DP < {MIN_DEPTH}")
    if MAX_DEPTH is not None:
        _depth_expr_parts.append(f"FMT/DP > {MAX_DEPTH}")
    _DEPTH_EXPR = " | ".join(_depth_expr_parts)

    rule depth_filter:
        """Mask genotypes with FORMAT/DP outside the configured min/max range.
        Sets out-of-range genotypes to missing (.) before MAF/missingness filtering."""
        input:  vcf = f"{INDIR}{VCF_RAW}"
        output: W['vcf_depth_filtered']
        params: expr = _DEPTH_EXPR
        log:    f"{LOGDIR}processing/depth_filter.log"
        shell:
            """
            bcftools filter -S . -e '{params.expr}' {input.vcf} -o {output} > {log} 2>&1
            echo "INFO: Depth filter applied: {params.expr}" >> {log}
            n_masked=$(bcftools stats {output} | awk '/^SN.*missing genotypes/ {{print $4}}')
            echo "INFO: Total missing genotypes after depth masking: $n_masked" >> {log}
            """

rule filter_vcf:
    """Filter VCF by MAF, missingness, and sample list (after sample missingness filter).
    Also normalizes chromosome names by removing 'chr' prefix to match LEA behavior."""
    input:
        vcf = W['vcf_depth_filtered'] if ((MIN_DEPTH is not None or MAX_DEPTH is not None) and HAS_FORMAT_DP) else f"{INDIR}{VCF_RAW}",
        samples = W['samples_rel_filtered'] if RELATEDNESS_REMOVE else (W['samples_het_filtered'] if HET_OUTLIER_SD is not None else W['samples_filtered'])
    output: W['vcf_filt']
    params: prefix = W['vcf_filt'].replace('.vcf', ''), maf = MAF, miss = MISS
    log:    f"{LOGDIR}processing/filter_vcf.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --keep {input.samples} \
            --maf {params.maf} --geno {params.miss} \
            --recode vcf --out {params.prefix} > {log} 2>&1
        sed -i '/^#CHROM/s/\t0_/\t/g' {output}

        # Normalize chromosome names — the pipeline's contract is: `chr` prefix stripped in
        # ANY case, X/Y/XY/MT kept as letters. plink already strips Chr/CHR/chr from the codes
        # it recognises and, with --output-chr MT above, writes X/Y/MT as letters instead of
        # 23/24/26; unrecognised contigs (chr2H, Chr5H) pass through verbatim, so the strip
        # here must be case-insensitive too. Body lines only, plus the ##contig header IDs so
        # header and body agree. scripts/normalize_gff.py applies the SAME rule to the GFF
        # and then compares the two contig sets (audit 2026-09-13 B1).
        sed -E -i '/^#/! s/^[Cc][Hh][Rr]//' {output}
        sed -E -i 's/^(##contig=<ID=)[Cc][Hh][Rr]/\1/' {output}

        echo "INFO: Normalized chromosome names (case-insensitive 'chr' strip; X/Y/MT kept as letters)" >> {log}
        """

rule compute_sample_het:
    """Compute per-sample heterozygosity on the filtered VCF (for het_vs_missingness plot)."""
    input:  vcf = W['vcf_filt']
    output: W['het_stats_filtered']
    params: prefix = f"{INTER}qc/het_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_sample_het.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# \
            --het --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "FID\tIID\tO_HOM\tE_HOM\tN_NM\tF"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.het > {output}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh {params.prefix}.het
        """

rule compute_snp_freq_filtered:
    """Compute allele frequencies on the filtered VCF (post MAF+missingness filter)."""
    input:  vcf = W['vcf_filt']
    output:
        frq = O['qc_maf_filtered'],
        pos = O['qc_maf_pos']
    params: prefix = f"{INTER}qc/maf_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_snp_freq_filtered.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# \
            --freq --make-bed --out {params.prefix} > {log} 2>&1
        awk 'NR==1 {{print "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS"; next}} {{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' \
            {params.prefix}.frq > {output.frq}
        paste <(awk '{{print $1"\t"$4}}' {params.prefix}.bim) \
              <(awk 'NR>1 {{print $5}}' {params.prefix}.frq) | \
            awk 'BEGIN{{print "CHR\tPOS\tMAF"}} {{print}}' > {output.pos}
        rm -f {params.prefix}.nosex {params.prefix}.log {params.prefix}.hh \
              {params.prefix}.frq {params.prefix}.bim {params.prefix}.bed {params.prefix}.fam
        """

rule compute_snp_density_filtered:
    """Compute SNP density in 1Mb windows on the filtered VCF."""
    input:  vcf = W['vcf_filt']
    output: O['qc_snp_density_filt']
    params: prefix = f"{INTER}qc/snp_density_filtered_tmp"
    log:    f"{LOGDIR}processing/compute_snp_density_filtered.log"
    shell:
        """
        /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --SNPdensity 1000000 --out {params.prefix} > {log} 2>&1
        cp {params.prefix}.snpden {output}
        rm -f {params.prefix}.snpden {params.prefix}.log
        """

if HAS_FORMAT_DP:
    rule compute_depth_stats:
        """Compute per-sample and per-site mean depth (only when FORMAT/DP field is present)."""
        input:  vcf = W['vcf_filt']
        output:
            sample = O['qc_depth_sample'],
            site   = O['qc_depth_site']
        params: prefix = f"{INTER}qc/depth_tmp"
        log:    f"{LOGDIR}processing/compute_depth_stats.log"
        shell:
            """
            /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --depth --out {params.prefix} > {log} 2>&1
            /app/vcftools-0.1.16/bin/vcftools --vcf {input.vcf} --site-mean-depth --out {params.prefix} >> {log} 2>&1
            cp {params.prefix}.idepth {output.sample}
            cp {params.prefix}.ldepth.mean {output.site}
            rm -f {params.prefix}.idepth {params.prefix}.ldepth.mean {params.prefix}.log
            """

rule ld_prune:
    """LD prune VCF (PCA is done separately via LEA)."""
    input:  vcf = W['vcf_filt']
    output: vcf = W['vcf_ld'], prune = W['prune_in']
    params: prefix = W['vcf_ld'].replace('.vcf', ''), win = LD_WIN, step = LD_STEP, r2 = LD_R2
    log:    f"{LOGDIR}processing/ld_prune.log"
    threads: CPU
    shell:
        """
        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# \
            --indep-pairwise {params.win} {params.step} {params.r2} \
            --out {params.prefix} > {log} 2>&1

        plink --vcf {input.vcf} --const-fid --allow-extra-chr --output-chr MT \
            --set-missing-var-ids @:# --extract {output.prune} \
            --make-bed --recode vcf \
            --out {params.prefix} >> {log} 2>&1

        sed -i '/^#CHROM/s/\t0_0_/\t/g' {output.vcf}
        """

rule extract_vcf_sample_order:
    """Get sample order from VCF header for metadata alignment."""
    input:  vcf = W['vcf_filt']
    output: W['samples_order']
    log:    f"{LOGDIR}processing/extract_vcf_sample_order.log"
    shell:  "grep -m1 CHROM {input.vcf} | cut -f10- | tr '\t' '\n' > {output} 2> {log}"

rule align_metadata:
    """Align metadata rows to match VCF sample order."""
    input:  meta = f"{INDIR}{SAMPLES}", order = W['samples_order']
    output: O['metadata']
    log:    f"{LOGDIR}processing/align_metadata.log"
    shell:  "Rscript /pipeline/scripts/filter_arrange_metadata.R {input.meta} {input.order} {output} > {log} 2>&1"

rule filter_coord_samples:
    """Filter samples to those with valid (non-NA) coordinates, for climate/GEA analysis.
    Unlike phenotypes, climate values cannot be meaningfully imputed, so missing
    coordinates are always dropped (no config option). Samples missing coordinates
    are retained in the full metadata.tsv (O['metadata']) for non-spatial analyses
    such as GWAS/phenotype association. Runs unconditionally (not gated on
    Climate.enabled) since piemaps/pop-stats also depend on valid coordinates."""
    input:  meta = O['metadata']
    output:
        samples_list = W['coord_valid_samples'],
        metadata     = W['metadata_climate'],
        summary      = O['coord_missing_summary']
    log:    f"{LOGDIR}processing/filter_coord_samples.log"
    shell:
        """
        Rscript /pipeline/scripts/filter_coord_samples.R \
            {input.meta} {output.samples_list} {output.metadata} {output.summary} > {log} 2>&1
        """

# Sources for LFMM-format matrix coord subsetting (climate coordinate NA handling).
# Each maps a {variant} wildcard to its full-cohort source matrix. Defined once here
# (rather than per-consumer in gea.smk/maladaptation.smk) to guarantee exactly one
# rule produces each output, since lfmm_imp_full is shared by GEA-LFMM and
# Maladaptation's geometric_offset.
_LFMM_CLIMATE_SOURCES = {
    'lfmm_imp':      W['lfmm_imp'],       # GEA-LFMM, LD-pruned imputed
    'lfmm_imp_full': W['lfmm_imp_full'],  # GEA-LFMM (full) + ALL Maladaptation methods
                                          # (gradient_forest, geometric_offset, rda_offset)
    'lfmm_full':     W['lfmm_full'],      # [UNUSED since 2026-08-10] was Maladaptation
                                          # gradient_forest. GF now takes the imputed matrix like
                                          # every other consumer -- feeding it the raw matrix meant
                                          # gradientForest() read LEA's missing code 9 as a
                                          # genotype (bug O, docs/pipeline_improvement_requests.md).
                                          # Plumbing kept so the variant stays buildable on demand;
                                          # no rule requests it, so it is no longer produced.
}

rule subset_lfmm_climate:
    """Row-subset an LFMM-format matrix to climate-valid samples (climate coordinate NA
    handling). EMMAX/LFMM/gradient_forest/geometric_offset all bind climate values to
    genotype-matrix rows positionally, so once filter_coord_samples/filter_climate_valid_samples
    drop samples with missing lat/lon or NA-climate raster extraction from the climate table,
    every genotype-side matrix consumed alongside it must be row-subset to match."""
    input:
        matrix        = lambda wc: _LFMM_CLIMATE_SOURCES[wc.variant],
        samples_order = W['samples_order'],
        coord_samples = W['climate_valid_samples']
    output: f"{INTER}climate_subset/{{variant}}_climate.lfmm"
    wildcard_constraints:
        variant = "lfmm_imp|lfmm_imp_full|lfmm_full"
    log: f"{LOGDIR}processing/subset_lfmm_climate_{{variant}}.log"
    shell:
        """
        Rscript /pipeline/scripts/subset_lfmm_matrix.R \
            {input.matrix} {input.samples_order} {input.coord_samples} {output} > {log} 2>&1
        """

rule normalize_gff:
    """Normalize GFF chromosome names to the VCF's contract and verify the two agree.

    scripts/normalize_gff.py strips the `chr` prefix case-insensitively (the same rule
    filter_vcf applies to the VCF body), then compares the GFF seqids with the contigs
    actually present in the FILTERED VCF: disjoint sets are a hard error, a partial
    overlap is logged. Until 2026-09-13 only a lowercase `chr` was stripped here while
    plink canonicalised `Chr1`/`CHR1`/`X`/`MT` on the VCF side, so a TAIR10-style or
    animal genome got an empty annotation with exit 0 (audit B1). The filtered VCF is an
    input purely for that check — every consumer of normalized.gff3 sits far downstream."""
    input:
        gff = f"{INDIR}{GFF}" if GFF else [],
        vcf = W['vcf_filt']
    output: W['gff_normalized']
    params: gff_arg = lambda wc, input: input.gff if GFF else 'NULL'
    log:    f"{LOGDIR}processing/normalize_gff.log"
    shell:
        """
        python3 /pipeline/scripts/normalize_gff.py {params.gff_arg} {input.vcf} {output} > {log} 2>&1
        """

rule vcf_to_lfmm:
    """Convert VCF to LEA formats (geno, lfmm)."""
    input:  vcf = W['vcf_ld']
    output:
        geno = W['geno'], lfmm = W['lfmm'], vcfsnp = W['vcfsnp'], removed = W['removed'],
        nmissing = W['lfmm_nmissing']
    log:    f"{LOGDIR}processing/vcf_to_lfmm.log"
    shell:  "Rscript /pipeline/scripts/vcf2lfmm.R {input.vcf} > {log} 2>&1"

rule plot_qc_processing:
    """Generate all processing QC plots and derived tables."""
    input:
        raw_summary      = O['qc_raw_summary'],
        sample_miss      = W['samples_missing_stats'],
        snp_miss_raw     = O['qc_snp_miss_raw'],
        maf_raw          = O['qc_maf_raw'],
        maf_filtered     = O['qc_maf_filtered'],
        het_raw          = W['het_stats_raw'],
        het_filtered     = W['het_stats_filtered'],
        snp_density_raw  = O['qc_snp_density_raw'],
        snp_density_filt = O['qc_snp_density_filt'],
        metadata         = O['metadata'],
        vcf_filt         = W['vcf_filt'],
        vcf_ld           = W['vcf_ld'],
        depth_sample     = O['qc_depth_sample'] if HAS_FORMAT_DP else [],
        depth_site       = O['qc_depth_site'] if HAS_FORMAT_DP else [],
        relatedness_pairs= O['qc_relatedness_pairs'],
        relatedness_removed = O['qc_relatedness_removed'] if PI_HAT is not None else []
    output:
        plot_sample_miss = O['qc_plot_sample_miss'],
        plot_het_miss    = O['qc_plot_het_miss'],
        plot_maf         = O['qc_plot_maf'],
        plot_snp_miss    = O['qc_plot_snp_miss'],
        plot_attrition   = O['qc_plot_attrition'],
        plot_snp_density = O['qc_plot_snp_density'],
        filtering_summary= O['qc_filtering_summary'],
        sample_het       = O['qc_sample_het'],
        plot_relatedness = O['qc_plot_relatedness'],
        plot_relatedness_mds = O['qc_plot_relatedness_mds'],
        **({"plot_depth": O['qc_plot_depth']} if HAS_FORMAT_DP else {})
    params:
        sample_miss_thresh = SAMPLE_MISS,
        snp_miss_thresh    = MISS,
        maf_thresh         = MAF,
        het_outlier_sd     = HET_OUTLIER_SD if HET_OUTLIER_SD is not None else 'NULL',
        has_dp             = 'TRUE' if HAS_FORMAT_DP else 'FALSE',
        depth_sample       = O['qc_depth_sample'] if HAS_FORMAT_DP else 'NULL',
        depth_site         = O['qc_depth_site'] if HAS_FORMAT_DP else 'NULL',
        plots_dir          = f"{MOD_PROCESSING}plots/",
        tables_dir         = f"{MOD_PROCESSING}tables/",
        pihat_thresh       = PI_HAT if PI_HAT is not None else 'NULL',
        relatedness_removed = O['qc_relatedness_removed'] if PI_HAT is not None else 'NULL',
        relatedness_action = RELATEDNESS_ACTION
    log: f"{LOGDIR}processing/plot_qc_processing.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_qc_processing.R \
            {input.raw_summary} \
            {input.sample_miss} \
            {input.snp_miss_raw} \
            {input.maf_raw} \
            {input.maf_filtered} \
            {input.het_raw} \
            {input.het_filtered} \
            {input.snp_density_raw} \
            {input.snp_density_filt} \
            {input.metadata} \
            {input.vcf_filt} \
            {input.vcf_ld} \
            {params.sample_miss_thresh} \
            {params.snp_miss_thresh} \
            {params.maf_thresh} \
            {params.het_outlier_sd} \
            {params.has_dp} \
            {params.depth_sample} \
            {params.depth_site} \
            {params.plots_dir} \
            {params.tables_dir} \
            {input.relatedness_pairs} \
            {params.pihat_thresh} \
            {params.relatedness_removed} \
            {params.relatedness_action} > {log} 2>&1
        """
