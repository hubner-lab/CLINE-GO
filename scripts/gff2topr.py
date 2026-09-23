#!/usr/bin/env python
## modified script from https://github.com/totajuliusd/topr
#
# GFF3 -> topr gene-annotation TSV.
#
# Called only by the Shiny app (scripts/clinego.app/R/fct_regions.R, via
# processx::run) to build the regional-plot annotation on demand; no Snakemake
# rule invokes it.
#
# The conversion is TWO-PASS on purpose. The single-pass version attached an exon
# with `if parent_id in genes`, i.e. only when the parent's own line had already
# been read — GFF3 mandates no such ordering, and an exon listed before its parent
# was dropped with no message. Children are therefore buffered and resolved after
# the whole file is read.
import sys
import re

if len(sys.argv) != 6:
    print("Usage: {} <file> <feature_name> <field_name> <biotype_name> <output_filename>".format(sys.argv[0]))
    sys.exit(1)

file_name = sys.argv[1]
feature_name = sys.argv[2]
field_name = sys.argv[3]
biotype_name = sys.argv[4]
output_filename = sys.argv[5]

# stderr, not stdout: the caller captures both, but stdout is the place a future
# consumer would look for data.
print('\t'.join([file_name, feature_name, field_name, biotype_name, output_filename]),
      file=sys.stderr)


def attr(attrs, key):
    """Value of a GFF3 attribute, or None. Keys are matched exactly."""
    m = re.search(r'(?:^|;)\s*{}=([^;]+)'.format(re.escape(key)), attrs)
    return m.group(1) if m else None


genes = {}          # feature id -> row being built
parent_of = {}      # feature id -> its Parent id (any feature type)
children = {"exon": [], "CDS": []}   # (parent id, start, end), file order
no_name = 0

try:
    with open(file_name, 'r') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            # "\r\n", not "\n": a CRLF GFF must not leave "\r" on the last attribute.
            columns = line.rstrip("\r\n").split('\t')
            if len(columns) < 9:  # Skip malformed lines
                continue

            chrom, start, end, attrs = columns[0], columns[3], columns[4], columns[8]
            feat_id = attr(attrs, "ID")
            parent_id = attr(attrs, "Parent")
            if feat_id and parent_id:
                # First parent only; a multi-parent feature is not a gene model
                # this plot can draw.
                parent_of[feat_id] = parent_id.split(',')[0]

            if columns[2] in children:
                if parent_id:
                    children[columns[2]].append(
                        (parent_id.split(',')[0], start, end))

            if columns[2] == feature_name:
                if not feat_id:
                    continue
                name = attr(attrs, field_name)
                if name is None:
                    # Was: drop the whole gene. A missing display name is not a
                    # reason to omit a gene from the regional plot.
                    name = feat_id
                    no_name += 1
                # Honour the biotype attribute when the GFF carries it. It used
                # to be hardcoded because requiring it dropped every gene of a
                # GFF without one (H. spontaneum); defaulting does both jobs.
                gene_type = attr(attrs, biotype_name) or 'protein_coding'
                if feat_id not in genes:
                    genes[feat_id] = {
                        "gene_start": start,
                        "gene_end": end,
                        "chr": chrom,
                        "gene_name": name,
                        "biotype": gene_type,
                        "exon_chromstart": [],
                        "exon_chromend": []
                    }

except FileNotFoundError:
    print("File not found: {}".format(file_name))
    sys.exit(1)


def owner(child_parent_id):
    """Walk Parent links up to the selected feature (exon -> mRNA -> gene)."""
    seen = set()
    node = child_parent_id
    while node is not None and node not in genes:
        if node in seen:          # malformed cycle
            return None
        seen.add(node)
        node = parent_of.get(node)
    return node


# exon is the primary structure; CDS is the same fallback read_gff_exons()
# (scripts/R/lib/gff_parsing.R) applies, so a CDS-only GFF — SIMDATA's own, and
# every Ensembl-style "no exon lines" annotation — still draws gene structure.
blocks = children["exon"] if children["exon"] else children["CDS"]
if not children["exon"] and children["CDS"]:
    print("No exon features; using CDS", file=sys.stderr)

orphans = 0
for parent_id, start, end in blocks:
    gene_id = owner(parent_id)
    if gene_id is None:
        orphans += 1
        continue
    genes[gene_id]["exon_chromstart"].append(start)
    genes[gene_id]["exon_chromend"].append(end)

if orphans:
    print("WARNING: {} exon/CDS features have no '{}' ancestor and were dropped"
          .format(orphans, feature_name), file=sys.stderr)
if no_name:
    print("WARNING: {} '{}' features carry no '{}' attribute; using their ID as the symbol"
          .format(no_name, feature_name, field_name), file=sys.stderr)

# Print the output
with open(output_filename, 'w') as w:
    w.write("\t".join(["chrom", "gene_start", "gene_end", "gene_symbol", "biotype", "exon_chromstart", "exon_chromend"]) + '\n')
    for gene, data in genes.items():
        if not gene:
            continue
        exon_starts = ",".join(data.get("exon_chromstart", []))
        exon_ends = ",".join(data.get("exon_chromend", []))
        w.write("\t".join([
            data["chr"], 
            data["gene_start"], 
            data["gene_end"], 
            data["gene_name"], 
            data["biotype"], 
            exon_starts, 
            exon_ends
        ]) + '\n')
