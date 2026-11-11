#!/usr/bin/env bash
# j2_snp_missing.sh -- Job 2 (s_qc panel e): per-SNP missing-call fraction as filter_vcf sees it: the run's input VCF,
# the 57 plants that passed the sample-missingness and relatedness filters (the run's --keep list), before --geno and
# --maf. PLINK 1.9 --missing (.lmiss F_MISS) with the pipeline's flags (+ --chr-set 95, which the 2026-10-14 run's own
# call predates; contigs 1-6 are autosomes under either set). Runs inside cline-go:latest from /pipeline.
set -euo pipefail
R=/pipeline/Trifolium86MainK5ERpc5_results
O=/pipeline/work/final_redraw_2026-11/inputs
plink --vcf /pipeline/data/Trifolium86_chr1-6.vcf.gz --const-fid --chr-set 95 --allow-extra-chr --output-chr MT \
    --set-missing-var-ids @:# --keep $R/_intermediate/samples/samples_rel_filtered.list \
    --missing --memory 32000 --threads 8 --out /tmp/snp_missing_57
cp /tmp/snp_missing_57.log $O/snp_missing_57.log
cp /tmp/snp_missing_57.imiss $O/snp_missing_57.imiss
gzip -c /tmp/snp_missing_57.lmiss > $O/snp_missing_57.lmiss.gz
echo "j2 plink done"
