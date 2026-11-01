#!/usr/bin/env bash
# =============================================================================
# Prepare Populus balsamifera benchmark dataset
# =============================================================================
# Source: Fitzpatrick et al. 2021 Molecular Ecology Resources
#   - VCF: github.com/stephenrkeller/Fitzpatrick_etal_MER_2021
#   - GFF: Ensembl Plants (Populus trichocarpa reference)
#
# Pipeline test: GEA + maladaptation (genetic offset validated by common garden)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CACHE_DIR="$SCRIPT_DIR/.cache/populus"
DATA_DIR="$PROJECT_DIR/data"

FORCE=false
[[ "${1:-}" == "--force" ]] && FORCE=true

# Final output files
VCF_OUT="$DATA_DIR/Populus_balsamifera.vcf"
META_OUT="$DATA_DIR/Populus_balsamifera_metadata.tsv"
GFF_OUT="$DATA_DIR/Populus_balsamifera.gff3"

# Check if already done
if [[ "$FORCE" == false && -f "$VCF_OUT" && -f "$META_OUT" && -f "$GFF_OUT" ]]; then
    echo "All output files already exist. Use --force to re-run."
    exit 0
fi

mkdir -p "$CACHE_DIR" "$DATA_DIR"

# =============================================================================
# Step 1: Clone GitHub repo
# =============================================================================
REPO_DIR="$CACHE_DIR/Fitzpatrick_etal_MER_2021"
REPO_URL="https://github.com/stephenrkeller/Fitzpatrick_etal_MER_2021.git"

echo ">>> Step 1: Cloning Fitzpatrick et al. 2021 data repository..."
if [[ ! -d "$REPO_DIR/.git" ]]; then
    git clone --depth 1 "$REPO_URL" "$REPO_DIR"
else
    echo "    Repository already cloned, skipping."
fi

# =============================================================================
# Step 2: Locate and prepare VCF
# =============================================================================
echo ">>> Step 2: Preparing VCF..."

# Find VCF file in repo (may be compressed)
VCF_SRC=""
for pattern in "*.vcf.gz" "*.vcf"; do
    found=$(find "$REPO_DIR" -name "$pattern" -type f 2>/dev/null | head -1)
    if [[ -n "$found" ]]; then
        VCF_SRC="$found"
        break
    fi
done

if [[ -z "$VCF_SRC" ]]; then
    echo "ERROR: No VCF file found in repository. Contents:"
    find "$REPO_DIR" -type f -name "*.vcf*" -o -name "*.gz" | head -20
    echo ""
    echo "Please check the repository structure and update this script."
    exit 1
fi

echo "    Found VCF: $VCF_SRC"

if [[ "$VCF_SRC" == *.gz ]]; then
    echo "    Decompressing..."
    gunzip -k -f "$VCF_SRC"
    VCF_SRC="${VCF_SRC%.gz}"
fi

cp "$VCF_SRC" "$VCF_OUT"
echo "    VCF copied to $VCF_OUT"

# Count samples and variants
N_SAMPLES=$(grep '^#CHROM' "$VCF_OUT" | tr '\t' '\n' | tail -n+10 | wc -l)
N_VARIANTS=$(grep -c -v '^#' "$VCF_OUT")
echo "    Samples: $N_SAMPLES, Variants: $N_VARIANTS"

# =============================================================================
# Step 3: Build metadata TSV
# =============================================================================
echo ">>> Step 3: Building metadata..."

# Extract sample names from VCF
SAMPLES_FILE="$CACHE_DIR/vcf_samples.txt"
grep '^#CHROM' "$VCF_OUT" | tr '\t' '\n' | tail -n+10 > "$SAMPLES_FILE"

# Look for population/coordinate data in the repo
POP_FILE=""
for pattern in "*pop*" "*coord*" "*meta*" "*sample*" "*loc*"; do
    found=$(find "$REPO_DIR" -iname "${pattern}.csv" -o -iname "${pattern}.tsv" -o -iname "${pattern}.txt" 2>/dev/null | head -1)
    if [[ -n "$found" ]]; then
        POP_FILE="$found"
        break
    fi
done

# Also check for RData files with population data
RDATA_FILES=$(find "$REPO_DIR" -name "*.RData" -o -name "*.rds" -o -name "*.Rdata" 2>/dev/null)

echo "    VCF samples file: $SAMPLES_FILE ($N_SAMPLES samples)"
if [[ -n "$POP_FILE" ]]; then
    echo "    Found population data: $POP_FILE"
    echo "    Preview:"
    head -3 "$POP_FILE" | sed 's/^/      /'
fi

# Build metadata using Python (handles CSV parsing, coordinate extraction)
python3 << 'PYEOF'
import csv
import os
import sys
import re

cache_dir = os.environ.get('CACHE_DIR', 'benchmarks/.cache/populus')
repo_dir = os.path.join(cache_dir, 'Fitzpatrick_etal_MER_2021')
samples_file = os.path.join(cache_dir, 'vcf_samples.txt')
meta_out = os.environ.get('META_OUT', 'data/Populus_balsamifera_metadata.tsv')

# Read VCF sample names
with open(samples_file) as f:
    vcf_samples = [line.strip() for line in f if line.strip()]

print(f"INFO: {len(vcf_samples)} samples from VCF", file=sys.stderr)

# Search for population/coordinate files
pop_data = {}  # sample -> {site, lat, lon}

# Try common file patterns
for root, dirs, files in os.walk(repo_dir):
    for fname in files:
        fpath = os.path.join(root, fname)
        if not fname.endswith(('.csv', '.tsv', '.txt')):
            continue

        try:
            # Try to read and detect columns
            with open(fpath) as f:
                first_line = f.readline()
            sep = ',' if fname.endswith('.csv') else '\t'
            header = first_line.strip().lower().split(sep)

            # Look for files with lat/lon and sample/pop columns
            has_coords = any(h in header for h in ['lat', 'latitude', 'y', 'lat_dd'])
            has_sample = any(h in header for h in ['sample', 'ind', 'individual', 'id', 'accession'])
            has_pop = any(h in header for h in ['pop', 'population', 'site', 'location', 'deme'])

            if has_coords and (has_sample or has_pop):
                print(f"INFO: Parsing coordinate file: {fname}", file=sys.stderr)
                print(f"  Header: {header}", file=sys.stderr)

                with open(fpath) as f:
                    reader = csv.DictReader(f, delimiter=sep)
                    fieldnames = [n.lower() for n in reader.fieldnames]
                    # Re-open with lowercase mapping
                    f.seek(0)
                    reader = csv.reader(f, delimiter=sep)
                    orig_header = next(reader)
                    header_map = {h.lower(): i for i, h in enumerate(orig_header)}

                    # Find column indices
                    lat_col = next((header_map[h] for h in ['latitude', 'lat', 'lat_dd', 'y'] if h in header_map), None)
                    lon_col = next((header_map[h] for h in ['longitude', 'lon', 'long', 'lon_dd', 'x'] if h in header_map), None)
                    pop_col = next((header_map[h] for h in ['pop', 'population', 'site', 'deme', 'location'] if h in header_map), None)
                    sample_col = next((header_map[h] for h in ['sample', 'ind', 'individual', 'id', 'accession'] if h in header_map), None)

                    if lat_col is not None and lon_col is not None:
                        for row in reader:
                            if len(row) <= max(lat_col, lon_col):
                                continue
                            try:
                                lat = float(row[lat_col])
                                lon = float(row[lon_col])
                            except (ValueError, IndexError):
                                continue

                            sample_id = row[sample_col] if sample_col is not None else None
                            pop_id = row[pop_col] if pop_col is not None else 'unknown'

                            if sample_id:
                                pop_data[sample_id] = {'site': pop_id, 'lat': lat, 'lon': lon}
                            else:
                                # Population-level coords (no per-sample IDs)
                                pop_data[f'_pop_{pop_id}'] = {'site': pop_id, 'lat': lat, 'lon': lon}

        except Exception as e:
            continue

# Check if we have per-sample or per-population data
per_sample = {k: v for k, v in pop_data.items() if not k.startswith('_pop_')}
per_pop = {k.replace('_pop_', ''): v for k, v in pop_data.items() if k.startswith('_pop_')}

print(f"INFO: Per-sample matches: {len(per_sample)}, Per-pop entries: {len(per_pop)}", file=sys.stderr)

# Write metadata
with open(meta_out, 'w') as f:
    f.write('site\tsample\tlatitude\tlongitude\n')

    matched = 0
    for sample in vcf_samples:
        if sample in per_sample:
            d = per_sample[sample]
            f.write(f"{d['site']}\t{sample}\t{d['lat']}\t{d['lon']}\n")
            matched += 1
        else:
            # Try to extract population from sample name (common pattern: POP_001 or POP001)
            pop_match = None
            for pop_id, d in per_pop.items():
                if sample.startswith(pop_id) or pop_id in sample:
                    pop_match = d
                    break

            if pop_match:
                f.write(f"{pop_match['site']}\t{sample}\t{pop_match['lat']}\t{pop_match['lon']}\n")
                matched += 1
            else:
                # Fallback: unknown site, no coords
                f.write(f"unknown\t{sample}\t0\t0\n")

print(f"INFO: Matched {matched}/{len(vcf_samples)} samples with coordinates", file=sys.stderr)

if matched == 0:
    print("WARNING: No samples matched with coordinate data!", file=sys.stderr)
    print("WARNING: You may need to manually create the metadata file.", file=sys.stderr)
    print("WARNING: Check available data files:", file=sys.stderr)
    for root, dirs, files in os.walk(repo_dir):
        for fname in sorted(files):
            if fname.endswith(('.csv', '.tsv', '.txt', '.RData', '.rds')):
                fpath = os.path.join(root, fname)
                size = os.path.getsize(fpath)
                print(f"  {os.path.relpath(fpath, repo_dir)} ({size} bytes)", file=sys.stderr)
PYEOF

echo "    Metadata written to $META_OUT"

# =============================================================================
# Step 4: Download GFF3 from Ensembl Plants
# =============================================================================
echo ">>> Step 4: Downloading Populus trichocarpa GFF3 from Ensembl Plants..."

# Pop_tri_v3, not v4: the VCF's positions are on v3 and run past v4's chromosome ends (a
# chromosome-1 SNP at 50,484,875; chromosome 1 is 49,788,581 bp in v4, 50,495,391 in v3).
# Ensembl Plants release 52 carries Pop_tri_v3: 39,514 genes, chromosomes 1-19 as in the VCF.
# No fallback release: a failed download must stop, never fall back to another assembly.
# The cache file is named after the assembly, so an old cached v4 download is never reused.
ENSEMBL_RELEASE="52"
GFF_URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${ENSEMBL_RELEASE}/gff3/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v3.${ENSEMBL_RELEASE}.chr.gff3.gz"
GFF_GZ="$CACHE_DIR/Populus_trichocarpa.Pop_tri_v3.gff3.gz"

if [[ ! -f "$GFF_GZ" ]]; then
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

# Count metadata rows
N_META=$(tail -n+2 "$META_OUT" | wc -l)
echo "    VCF samples: $N_SAMPLES"
echo "    Metadata rows: $N_META"

if [[ "$N_SAMPLES" -ne "$N_META" ]]; then
    echo "WARNING: Sample count mismatch! VCF=$N_SAMPLES, Metadata=$N_META"
fi

# Check for zero coordinates
N_ZERO=$(awk -F'\t' 'NR>1 && $3==0 && $4==0' "$META_OUT" | wc -l)
if [[ "$N_ZERO" -gt 0 ]]; then
    echo "WARNING: $N_ZERO samples have zero coordinates (may need manual fixing)"
fi

# GFF stats
N_GENES=$(grep -c 'mRNA' "$GFF_OUT" || echo "0")
echo "    GFF mRNA features: $N_GENES"

echo ""
echo "=== Populus balsamifera dataset preparation complete ==="
echo "Files:"
echo "  VCF:      $VCF_OUT"
echo "  Metadata: $META_OUT"
echo "  GFF:      $GFF_OUT"
echo ""
echo "Next: run pipeline with config_populus.yaml"
echo "  docker run ... snakemake --config mode=processing --configfile config_populus.yaml"
