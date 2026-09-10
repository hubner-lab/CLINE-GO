#!/usr/bin/env bash
# =============================================================================
# Prepare Arabidopsis 1001 Genomes benchmark dataset
# =============================================================================
# Source: The 1001 Genomes Consortium, 2016, Cell 166(2):481-491
#   - VCF: 1001genomes.org (v3.1, 1,135 accessions, ~11M SNPs)
#   - Phenotypes: AraPheno (FT10, FT16 — flowering time at 10C and 16C)
#   - GFF: Ensembl Plants (Arabidopsis_thaliana.TAIR10)
#
# Pipeline test: Full pipeline (GEA + GWAS + maladaptation)
# Ground truth: FRI (chr4:269k), FLC (chr5:3.17M) for flowering time
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CACHE_DIR="$SCRIPT_DIR/.cache/arabidopsis"
DATA_DIR="$PROJECT_DIR/data"

FORCE=false
[[ "${1:-}" == "--force" ]] && FORCE=true

# Final output files
VCF_OUT="$DATA_DIR/Arabidopsis_1001G.vcf"
META_OUT="$DATA_DIR/Arabidopsis_1001G_metadata.tsv"
GFF_OUT="$DATA_DIR/Arabidopsis_1001G.gff3"

# Target ~50K SNPs after thinning
TARGET_SNPS=50000

# AraPheno phenotype IDs (largest studies with ~1100 values each)
PHENO_FT10_ID=261
PHENO_FT16_ID=262

# Check if already done
if [[ "$FORCE" == false && -f "$VCF_OUT" && -f "$META_OUT" && -f "$GFF_OUT" ]]; then
    echo "All output files already exist. Use --force to re-run."
    exit 0
fi

mkdir -p "$CACHE_DIR" "$DATA_DIR"

# Check dependencies
for cmd in bcftools python3 wget; do
    if ! command -v "$cmd" &>/dev/null; then
        echo "ERROR: $cmd is required but not found."
        echo "  Install via: nix shell nixpkgs#$cmd"
        exit 1
    fi
done

echo "=== Arabidopsis 1001 Genomes Dataset Preparation ==="
echo ""

# =============================================================================
# Step 1: Download full VCF from 1001 Genomes
# =============================================================================
echo ">>> Step 1: Downloading Arabidopsis 1001 Genomes VCF..."

VCF_URL="https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
VCF_FULL="$CACHE_DIR/1001genomes_full.vcf.gz"

if [[ ! -f "$VCF_FULL" ]]; then
    echo "    Downloading ~18 GB VCF (this will take a while)..."
    echo "    URL: $VCF_URL"
    wget -c -q --show-progress -O "$VCF_FULL" "$VCF_URL"
else
    echo "    Full VCF already cached, skipping download."
fi

# Verify file is not truncated (should be >10GB)
FILE_SIZE=$(stat -c%s "$VCF_FULL" 2>/dev/null || stat -f%z "$VCF_FULL" 2>/dev/null)
if [[ "$FILE_SIZE" -lt 10000000000 ]]; then
    echo "WARNING: VCF file seems small ($FILE_SIZE bytes). Download may be incomplete."
    echo "    Delete $VCF_FULL and re-run, or use --force."
fi

# =============================================================================
# Step 2: Thin VCF to ~50K SNPs
# =============================================================================
echo ">>> Step 2: Thinning VCF to ~${TARGET_SNPS} SNPs..."

VCF_THINNED="$CACHE_DIR/1001genomes_thinned.vcf"

if [[ ! -f "$VCF_THINNED" || "$FORCE" == true ]]; then
    # 1001 Genomes v3.1 has ~10.7M SNPs. Calculate thinning step.
    # Single-pass: keep header + every Nth data line
    ESTIMATED_TOTAL=10700000
    STEP=$((ESTIMATED_TOTAL / TARGET_SNPS))
    [[ "$STEP" -lt 1 ]] && STEP=1
    echo "    Keeping every ${STEP}th variant (estimated total: ~${ESTIMATED_TOTAL})"
    echo "    This reads the full 18 GB file — may take 10-30 minutes..."

    zcat "$VCF_FULL" | awk -v step="$STEP" '
        /^#/ { print; next }
        { n++; if (n % step == 0) print }
    ' > "$VCF_THINNED"

    KEPT=$(grep -c -v '^#' "$VCF_THINNED")
    echo "    Kept $KEPT variants"
else
    echo "    Thinned VCF already cached, skipping."
fi

# Copy to final output
cp "$VCF_THINNED" "$VCF_OUT"

N_SAMPLES=$(grep '^#CHROM' "$VCF_OUT" | tr '\t' '\n' | tail -n+10 | wc -l)
N_VARIANTS=$(grep -c -v '^#' "$VCF_OUT")
echo "    Final VCF: $N_SAMPLES samples, $N_VARIANTS variants"

# =============================================================================
# Step 3: Build metadata with coordinates and phenotypes from AraPheno
# =============================================================================
echo ">>> Step 3: Building metadata from AraPheno API..."

# Extract VCF sample IDs
SAMPLES_FILE="$CACHE_DIR/vcf_samples.txt"
grep '^#CHROM' "$VCF_OUT" | tr '\t' '\n' | tail -n+10 > "$SAMPLES_FILE"

# Download accession list (coordinates) and phenotype values
ACCESSION_JSON="$CACHE_DIR/accessions.json"
PHENO_FT10_JSON="$CACHE_DIR/pheno_ft10.json"
PHENO_FT16_JSON="$CACHE_DIR/pheno_ft16.json"

echo "    Downloading accession metadata..."
if [[ ! -f "$ACCESSION_JSON" || "$FORCE" == true ]]; then
    wget -q -O "$ACCESSION_JSON" "https://arapheno.1001genomes.org/rest/accession/list.json"
fi

echo "    Downloading FT10 phenotype values (ID=$PHENO_FT10_ID)..."
if [[ ! -f "$PHENO_FT10_JSON" || "$FORCE" == true ]]; then
    wget -q -O "$PHENO_FT10_JSON" "https://arapheno.1001genomes.org/rest/phenotype/${PHENO_FT10_ID}/values.json"
fi

echo "    Downloading FT16 phenotype values (ID=$PHENO_FT16_ID)..."
if [[ ! -f "$PHENO_FT16_JSON" || "$FORCE" == true ]]; then
    wget -q -O "$PHENO_FT16_JSON" "https://arapheno.1001genomes.org/rest/phenotype/${PHENO_FT16_ID}/values.json"
fi

# Build metadata TSV using Python
echo "    Building metadata TSV..."
python3 << 'PYEOF'
import json
import os
import sys

cache_dir = os.environ.get('CACHE_DIR', 'benchmarks/.cache/arabidopsis')
meta_out = os.environ.get('META_OUT', 'data/Arabidopsis_1001G_metadata.tsv')

# Read VCF sample IDs
with open(os.path.join(cache_dir, 'vcf_samples.txt')) as f:
    vcf_samples = [line.strip() for line in f if line.strip()]

print(f"INFO: {len(vcf_samples)} samples from VCF", file=sys.stderr)

# Read accession metadata (coordinates)
with open(os.path.join(cache_dir, 'accessions.json')) as f:
    accessions = json.load(f)

# Build lookup: pk -> {country, lat, lon}
acc_lookup = {}
for acc in accessions:
    pk = str(acc.get('pk', ''))
    lat = acc.get('latitude')
    lon = acc.get('longitude')
    country = acc.get('country', 'unknown')
    if lat is not None and lon is not None:
        acc_lookup[pk] = {
            'country': country.replace('\t', ' ').strip() if country else 'unknown',
            'lat': float(lat),
            'lon': float(lon)
        }

print(f"INFO: {len(acc_lookup)} accessions with coordinates from AraPheno", file=sys.stderr)

# Read phenotype values
def load_phenotype(filepath):
    """Load phenotype JSON from AraPheno, return {accession_id: value}"""
    with open(filepath) as f:
        data = json.load(f)
    pheno = {}
    for entry in data:
        acc_id = str(entry.get('accession_id', entry.get('obs_unit_id', '')))
        value = entry.get('phenotype_value')
        if value is not None and acc_id:
            try:
                pheno[acc_id] = float(value)
            except (ValueError, TypeError):
                continue
    return pheno

ft10 = load_phenotype(os.path.join(cache_dir, 'pheno_ft10.json'))
ft16 = load_phenotype(os.path.join(cache_dir, 'pheno_ft16.json'))

print(f"INFO: FT10 values: {len(ft10)}, FT16 values: {len(ft16)}", file=sys.stderr)

# Write metadata TSV
# CLINE-GO format: site, sample, latitude, longitude, [trait1, trait2, ...]
matched_coords = 0
matched_ft10 = 0
matched_ft16 = 0

with open(meta_out, 'w') as f:
    f.write('site\tsample\tlatitude\tlongitude\tFT10\tFT16\n')
    for sample in vcf_samples:
        acc = acc_lookup.get(sample, {'country': 'unknown', 'lat': 0, 'lon': 0})
        country = acc['country']
        lat = acc['lat']
        lon = acc['lon']

        if lat != 0 or lon != 0:
            matched_coords += 1

        ft10_val = ft10.get(sample, '')
        ft16_val = ft16.get(sample, '')

        if ft10_val != '':
            matched_ft10 += 1
        if ft16_val != '':
            matched_ft16 += 1

        f.write(f"{country}\t{sample}\t{lat}\t{lon}\t{ft10_val}\t{ft16_val}\n")

print(f"INFO: Matched coordinates: {matched_coords}/{len(vcf_samples)}", file=sys.stderr)
print(f"INFO: Matched FT10: {matched_ft10}/{len(vcf_samples)}", file=sys.stderr)
print(f"INFO: Matched FT16: {matched_ft16}/{len(vcf_samples)}", file=sys.stderr)

# Warn about zero-coordinate samples
zero_coords = len(vcf_samples) - matched_coords
if zero_coords > 0:
    print(f"WARNING: {zero_coords} samples have zero coordinates", file=sys.stderr)
PYEOF

echo "    Metadata written to $META_OUT"

# =============================================================================
# Step 4: Download GFF3 from Ensembl Plants
# =============================================================================
echo ">>> Step 4: Downloading Arabidopsis thaliana GFF3 from Ensembl Plants..."

ENSEMBL_RELEASE="62"
GFF_URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${ENSEMBL_RELEASE}/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.${ENSEMBL_RELEASE}.gff3.gz"
GFF_GZ="$CACHE_DIR/Arabidopsis_thaliana.TAIR10.gff3.gz"

if [[ ! -f "$GFF_GZ" || "$FORCE" == true ]]; then
    echo "    Downloading from Ensembl Plants release $ENSEMBL_RELEASE..."
    wget -c -q --show-progress -O "$GFF_GZ" "$GFF_URL"
else
    echo "    GFF already cached, skipping download."
fi

echo "    Decompressing..."
gunzip -k -f "$GFF_GZ"
GFF_DECOMPRESSED="${GFF_GZ%.gz}"
cp "$GFF_DECOMPRESSED" "$GFF_OUT"
echo "    GFF copied to $GFF_OUT"

# =============================================================================
# Step 5: Validate
# =============================================================================
echo ">>> Step 5: Validating..."

N_META=$(tail -n+2 "$META_OUT" | wc -l)
echo "    VCF samples: $N_SAMPLES"
echo "    VCF variants: $N_VARIANTS"
echo "    Metadata rows: $N_META"

if [[ "$N_SAMPLES" -ne "$N_META" ]]; then
    echo "WARNING: Sample count mismatch! VCF=$N_SAMPLES, Metadata=$N_META"
fi

# Check coordinate coverage
N_ZERO=$(awk -F'\t' 'NR>1 && $3==0 && $4==0' "$META_OUT" | wc -l)
if [[ "$N_ZERO" -gt 0 ]]; then
    echo "WARNING: $N_ZERO samples have zero coordinates"
fi

# Check phenotype coverage
N_FT10=$(awk -F'\t' 'NR>1 && $5!=""' "$META_OUT" | wc -l)
N_FT16=$(awk -F'\t' 'NR>1 && $6!=""' "$META_OUT" | wc -l)
echo "    Phenotype FT10: $N_FT10 values"
echo "    Phenotype FT16: $N_FT16 values"

# GFF stats
N_GENES=$(grep -c 'mRNA' "$GFF_OUT" || echo "0")
echo "    GFF mRNA features: $N_GENES"

# Check chromosome names match
VCF_CHROMS=$(grep -v '^#' "$VCF_OUT" | cut -f1 | sort -u | tr '\n' ',' | sed 's/,$//')
GFF_CHROMS=$(grep -v '^#' "$GFF_OUT" | cut -f1 | sort -u | grep -v '^$' | tr '\n' ',' | sed 's/,$//')
echo "    VCF chromosomes: $VCF_CHROMS"
echo "    GFF chromosomes (first 10): $(echo "$GFF_CHROMS" | cut -d',' -f1-10)"

echo ""
echo "=== Arabidopsis 1001 Genomes dataset preparation complete ==="
echo ""
echo "Files:"
echo "  VCF:      $VCF_OUT ($N_SAMPLES samples, $N_VARIANTS variants)"
echo "  Metadata: $META_OUT (FT10: $N_FT10, FT16: $N_FT16 values)"
echo "  GFF:      $GFF_OUT ($N_GENES mRNA features)"
echo ""
echo "NOTE: Ensembl Plants GFF3 does NOT include GO terms."
echo "      Enrichment analysis will not produce results."
echo "      The key benchmark is detecting FRI (chr4:269k) and FLC (chr5:3.17M)"
echo "      in association with flowering time phenotypes."
echo ""
echo "Next: run pipeline with config_arabidopsis.yaml"
echo "  docker run ... snakemake --config mode=processing --configfile config_arabidopsis.yaml"
