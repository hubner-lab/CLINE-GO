#!/usr/bin/env bash
# s_qc (supplement build 2026-10-25, B1): fig3/plot_qc_ms.R with the 25 arguments the pipeline's plot_qc_processing rule
# received for this run (Snakemake metadata shellcmd, 2026-10-14 20:53, written under the project's earlier name
# Trifolium86MainK5ER; the same files resolved under Trifolium86MainK5ERpc5_results/), then output stem and counts TSV.
# Run inside the image from /pipeline: bash ANALYSIS/Trifolium86MainK5ERpc5/fig3/run_qc_ms.sh
set -euo pipefail
R=/pipeline/Trifolium86MainK5ERpc5_results
W=$R/_work/maf0.1_miss0.3_smiss0.8_rel0.999
AN=/pipeline/ANALYSIS/Trifolium86MainK5ERpc5
Rscript $AN/fig3/plot_qc_ms.R \
    $R/Processing/tables/vcf_raw_summary.tsv \
    $R/_intermediate/samples/samples_missing_stats.tsv \
    $R/_intermediate/qc/snp_missingness_raw.lmiss \
    $R/_intermediate/qc/maf_raw.frq \
    $R/_intermediate/qc/maf_filtered.frq \
    $R/_intermediate/qc/het_raw.het \
    $R/_intermediate/qc/het_filtered.het \
    $R/_intermediate/qc/snp_density_raw.snpden \
    $R/_intermediate/qc/snp_density_filtered.snpden \
    $R/Processing/tables/metadata.tsv \
    $W/Trifolium86_chr1-6.vcf \
    $W/ld0.2_win100_step10/Trifolium86_chr1-6.vcf \
    0.8 0.3 0.1 NULL TRUE \
    $R/_intermediate/qc/depth_per_sample.idepth \
    $R/_intermediate/qc/depth_per_site.ldepth.mean \
    $R/Processing/plots/ \
    $R/Processing/tables/ \
    $R/Processing/tables/relatedness_pairs.tsv \
    0.999 \
    $R/Processing/tables/relatedness_removed.tsv \
    remove \
    $AN/manuscript_fig3/draft/s_qc \
    $AN/manuscript_fig3/draft/s_qc_counts.tsv
