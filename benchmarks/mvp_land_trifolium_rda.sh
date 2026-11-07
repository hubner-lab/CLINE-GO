#!/usr/bin/env bash
# mvp_land_trifolium_rda.sh <control|uncorrected> [threads] -- Job C of the landscapes pack: RDA on Trifolium with
# condition_pcs 0, beside the pipeline. The GEA RDA block of config_Trifolium86MainK5d150.yaml, run as
# workflow/rules/gea.smk rule assoc_rda_gea runs it: scripts/rda.R on the paper's image with the paper's code.
#   control      condition_pcs 5, the published fit, to prove this path reproduces its p-values
#   uncorrected  condition_pcs 0, everything else equal
# Code: `git archive 3d4ccab scripts` at /mnt/data/eugene/TEMP/landscapes/trifolium_code, the last rda.R change before
# the fit of 2026-10-16 11:17 (equal to HEAD of that morning, db6b7ab, in rda.R and every file it sources); overlaid
# read-only because rda.R sources absolute /pipeline/scripts/... paths. The project tree is mounted read-only.
# Writes landscape_ms/trifolium_rda_uncorrected/<which>/ only. cwd = /tmp of the container, so a stray Rplots.pdf
# cannot land in the repo.
set -uo pipefail
WHICH="${1:?control|uncorrected}"; THREADS="${2:-8}"
case $WHICH in control) PCS=5 ;; uncorrected) PCS=0 ;; *) echo "control or uncorrected"; exit 2 ;; esac
REPO=/mnt/data/eugene/ADAPTOGENE; cd "$REPO"
CODE=/mnt/data/eugene/TEMP/landscapes/trifolium_code
IMG=cline-go:paper-d150
WANT=sha256:93bd6025e77b18e9b97ae1bea4a37e9c43a5b64c8acc067819632e09d550546f
DOCKER=(nix shell nixpkgs#docker-client -c docker)
[[ -s "$CODE/scripts/rda.R" ]] || { echo "code export missing at $CODE"; exit 2; }
id=$("${DOCKER[@]}" images --no-trunc --format '{{.ID}}' "$IMG" 2>/dev/null | head -1)
[[ "$id" == "$WANT" ]] || { echo "$IMG is $id, expected $WANT"; exit 2; }
P=Trifolium86MainK5d150_results; C=/pipeline/$P
W=$C/_work/maf0.1_miss0.3_smiss0.8_rel0.999; L=$W/ld0.2_win100_step10
OH=benchmarks/mvp_eval/landscape_ms/trifolium_rda_uncorrected/$WHICH; O=/pipeline/$OH
mkdir -p "$OH/plots"
T="/mnt/data/eugene/TEMP/landscapes/trifolium_$WHICH"; mkdir -p "$T"
"${DOCKER[@]}" run --rm --name "lpack_tri_$WHICH" -u "$(id -u):$(id -g)" -e USER=pipeline \
    --cpus="$THREADS" --memory=32g \
    -v "$T":/tmp -e TMPDIR=/tmp -e HOME=/tmp \
    -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e R_DATATABLE_NUM_THREADS=2 \
    -v "$REPO":/pipeline \
    -v "$CODE/scripts":/pipeline/scripts:ro \
    -v "$REPO/$P":"$C":ro \
    -w /tmp "$IMG" \
    Rscript /pipeline/scripts/rda.R \
    $C/_intermediate/climate_subset/lfmm_imp_full_climate.lfmm $C/_intermediate/climate_subset/lfmm_imp_climate.lfmm \
    $C/climate/tables/present/climate_present_site_scaled.tsv \
    $W/Trifolium86_chr1-6_K5imp_maf0.1.vcfsnp $L/Trifolium86_chr1-6_K5imp_maf0.1.vcfsnp \
    $L/Trifolium86_chr1-6.pca/Trifolium86_chr1-6.projections \
    $C/_intermediate/samples/samples_order.list $C/_intermediate/samples/climate_valid_samples.list \
    bio_1,bio_3,bio_5 "$PCS" auto 0.05 99 auto 2000000 64 42 "$THREADS" 5 \
    $O/RDA_pvalues_K5.tsv $O/RDA_candidates_K5.tsv $O/RDA_diagnostics_K5.tsv $O/RDA_anova_K5.tsv $O/plots/ qval 0.05
rc=$?
rm -rf "$T"
exit $rc
