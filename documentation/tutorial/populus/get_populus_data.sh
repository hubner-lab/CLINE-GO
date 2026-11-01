#!/usr/bin/env bash
# =============================================================================
# get_populus_data.sh <dest> — fetch the CLINE-GO tutorial inputs (Populus balsamifera)
# =============================================================================
# Writes into <dest> (the tutorial's configs read data/populus/, so run it from the repo
# root as `bash documentation/tutorial/populus/get_populus_data.sh data/populus`):
#
#   balsam_core336inds_42pops.vcf.gz       336 trees, 42 populations, 107,309 SNPs, 19 chr.
#                                          Fitzpatrick et al. 2021, Mol Ecol Resour,
#                                          doi:10.1111/1755-0998.13374; the authors' repo
#                                          github.com/stephenrkeller/Fitzpatrick_etal_MER_2021
#                                          at a pinned commit. The repo carries no license, so
#                                          CLINE-GO does not redistribute the VCF: you download
#                                          it from the authors. Kept gzipped: every pipeline
#                                          step that reads the input VCF (plink, bcftools,
#                                          vcftools) takes .vcf.gz.
#   balsam_core336inds_42pops_sorted.vcf.gz
#                                          the same VCF with chromosomes in numeric order:
#                                          the file the pipeline reads (see below for why).
#   Populus_trichocarpa.Pop_tri_v3.52.chr.gff3
#                                          P. trichocarpa gene models, assembly Pop_tri_v3
#                                          (JGI, GCA_000002775.3), Ensembl Plants release 52:
#                                          39,514 genes on chromosomes 1-19, named as in the
#                                          VCF. Unpacked, because the pipeline's GFF reader
#                                          (scripts/normalize_gff.py) opens it as plain text.
#                                          Not v4: the VCF's positions run past v4's
#                                          chromosome ends (a chromosome-1 SNP sits at
#                                          50,484,875; chromosome 1 is 49,788,581 bp in v4,
#                                          50,495,391 bp in v3).
#   populus_metadata.tsv                   site, sample, latitude, longitude: copied from this
#                                          directory (sources in SOURCES.md, sites.tsv).
#
# Both downloads are checked against a fixed sha256 and fetched again if a cached copy fails
# the check. Nothing else leaves this machine: no contact field, no identifier.
# =============================================================================
set -euo pipefail

if [[ $# -ne 1 || -z "${1:-}" ]]; then
    echo "usage: $0 <dest directory>   (e.g. data/populus)" >&2
    exit 2
fi
DEST="$1"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for tool in curl sha256sum gzip; do
    command -v "$tool" > /dev/null || { echo "ERROR: $tool is not installed" >&2; exit 1; }
done

VCF_COMMIT="7b39b0034e9224da5a62ee3f63f03e197c5735df"
VCF_NAME="balsam_core336inds_42pops.vcf.gz"
VCF_URL="https://raw.githubusercontent.com/stephenrkeller/Fitzpatrick_etal_MER_2021/${VCF_COMMIT}/${VCF_NAME}"
VCF_SHA256="053c486247aa9ea78c9fb9d447c9e727273f35e3c81647b4aa68d99a960109f6"

GFF_RELEASE="52"
GFF_GZ_NAME="Populus_trichocarpa.Pop_tri_v3.${GFF_RELEASE}.chr.gff3.gz"
GFF_URL="https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-${GFF_RELEASE}/gff3/populus_trichocarpa/${GFF_GZ_NAME}"
GFF_SHA256="3152ef21bb160998efc1e31367a867ed28df29be730b239f27642a5e9e2046b5"

META_NAME="populus_metadata.tsv"

mkdir -p "$DEST"

sha_ok() {  # sha_ok <file> <expected sha256>
    [[ -f "$1" ]] && [[ "$(sha256sum "$1" | cut -d' ' -f1)" == "$2" ]]
}

fetch() {  # fetch <url> <file> <expected sha256>
    local url="$1" out="$2" sha="$3"
    if sha_ok "$out" "$sha"; then
        echo "    $(basename "$out"): present, sha256 OK"
        return 0
    fi
    echo "    $(basename "$out"): downloading $url"
    curl -fL --progress-bar --retry 3 --retry-delay 5 -o "$out.part" "$url"
    if ! sha_ok "$out.part" "$sha"; then
        echo "ERROR: sha256 of $(basename "$out") is $(sha256sum "$out.part" | cut -d' ' -f1)," >&2
        echo "       expected $sha. The source changed or the download is corrupt." >&2
        rm -f "$out.part"
        exit 1
    fi
    mv "$out.part" "$out"
    echo "    $(basename "$out"): sha256 OK"
}

echo ">>> VCF (Fitzpatrick et al. 2021, commit ${VCF_COMMIT:0:7})"
fetch "$VCF_URL" "$DEST/$VCF_NAME" "$VCF_SHA256"

# The authors' VCF lists chromosomes 10-19 before 1-9. The pipeline's LD-decay step
# (scripts/ld_decay_pairs.R) pairs PLINK's variant table, which PLINK writes in numeric
# chromosome order, with the VCF's rows by position, and stops on any other order
# (filed 2026-10-31). So the pipeline reads a copy with the same header and the same records,
# sorted by chromosome then position; the downloaded file stays as it is.
SORTED_NAME="${VCF_NAME%.vcf.gz}_sorted.vcf.gz"
if [[ ! -s "$DEST/$SORTED_NAME" ]]; then
    echo "    writing $SORTED_NAME (chromosomes in numeric order, for the pipeline's LD-decay step)"
    { gzip -dc "$DEST/$VCF_NAME" | grep '^#'
      gzip -dc "$DEST/$VCF_NAME" | grep -v '^#' | LC_ALL=C sort -t "$(printf '\t')" -k1,1n -k2,2n -s
    } | gzip -n -c > "$DEST/$SORTED_NAME.part"
    n_in=$(gzip -dc "$DEST/$VCF_NAME" | grep -vc '^#')
    n_out=$(gzip -dc "$DEST/$SORTED_NAME.part" | grep -vc '^#')
    if [[ "$n_in" -ne "$n_out" ]]; then
        echo "ERROR: sorted VCF has $n_out records, the download $n_in" >&2
        rm -f "$DEST/$SORTED_NAME.part"; exit 1
    fi
    mv "$DEST/$SORTED_NAME.part" "$DEST/$SORTED_NAME"
else
    echo "    $SORTED_NAME: present"
fi

echo ">>> GFF (Ensembl Plants ${GFF_RELEASE}, Pop_tri_v3)"
fetch "$GFF_URL" "$DEST/$GFF_GZ_NAME" "$GFF_SHA256"
GFF_OUT="$DEST/${GFF_GZ_NAME%.gz}"
gzip -dc "$DEST/$GFF_GZ_NAME" > "$GFF_OUT.part"
mv "$GFF_OUT.part" "$GFF_OUT"
# The checked .gz stays beside it, so a rerun verifies it instead of downloading again.
echo "    unpacked to $(basename "$GFF_OUT") (the pipeline reads the GFF as plain text)"

echo ">>> Metadata"
cp "$HERE/$META_NAME" "$DEST/$META_NAME"

# Every metadata sample must be in the VCF header (the pipeline stops on one that is not).
missing=$(comm -23 <(tail -n +2 "$DEST/$META_NAME" | cut -f2 | sort) \
                   <(gzip -dc "$DEST/$VCF_NAME" | grep -m1 '^#CHROM' | cut -f10- | tr '\t' '\n' | sort))
if [[ -n "$missing" ]]; then
    echo "ERROR: metadata samples not in the VCF: $missing" >&2
    exit 1
fi

n_meta=$(($(wc -l < "$DEST/$META_NAME") - 1))
n_sites=$(tail -n +2 "$DEST/$META_NAME" | cut -f1 | sort -u | wc -l)
n_genes=$(awk -F'\t' '$3 == "gene"' "$GFF_OUT" | wc -l)
echo ""
echo "Done: $DEST"
echo "  VCF       $SORTED_NAME (the pipeline's input; $VCF_NAME is the authors' file)"
echo "  GFF       $(basename "$GFF_OUT") ($n_genes genes)"
echo "  metadata  $META_NAME ($n_meta trees, $n_sites sites)"
