# config.sh -- shared parameters for the Trifolium gene-identity re-run.
# Sourced by the numbered scripts; holds no secrets and prints nothing.

# Lab server (read-only use: every remote command only reads and prints to stdout).
LAB_HOST="lab"
LAB_RESULTS="/mnt/data/eugene/ADAPTOGENE/Trifolium86MainK5ERpc5_results"
LAB_GENOME_DIR="/mnt/data/eugene/Trifolium/data/Genome"
LAB_GTF="${LAB_GENOME_DIR}/Trifolium_genes.gtf"                       # BRAKER gene models (full exon structure)
LAB_FASTA="${LAB_GENOME_DIR}/Trifolium_Israeliticum_ref_11052025.fa"   # assembly, indexed (.fai)
LAB_GFFREAD="/usr/local/bin/gffread"
LAB_SAMTOOLS="/usr/local/bin/samtools"
LAB_CHR_PREFIX="Chr_"   # assembly/GTF use Chr_1..Chr_6, the pipeline tables use 1..6

# BLAST+ from a pinned nixpkgs revision (blast 2.17.0,
# /nix/store/3kb3321b1ncs3n4dnwwgzxz33wzqfksj-blast-2.17.0).
NIXPKGS_REF="github:NixOS/nixpkgs/901c54c08bdc20e00c2b3c73716b97497b1b3c74"

# UniProt REST (stream endpoint returns the full result set in one response).
UNIPROT_STREAM="https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_FIELDS="accession,id,reviewed,gene_names,gene_oln,gene_orf,protein_name,length,xref_panther,xref_ensemblplants,organism_id"

# Search parameters.
EVALUE="1e-5"            # BLASTP / TBLASTN E-value cutoff
MAX_TARGET_SEQS="500"    # kept by BLAST; the top 5 per query are cut afterwards (max_target_seqs acts inside the search)
TOP_N="5"                # hits kept per query and database in the long table
THREADS="8"

# Microsynteny (regions 1 and 12).
SYN_REGIONS="1_8649433-8911055 5_5557888-5883272"
SYN_BITSCORE_FRAC="0.90" # candidate Arabidopsis loci: bitscore >= 90% of the gene's best locus bitscore
SYN_AGI_WINDOW="100"     # same Arabidopsis chromosome, |AGI number difference| <= 100 (AGI codes step by 10 per gene)
SYN_MIN_SUPPORT="2"      # a candidate is syntenic if >= 2 other region genes have a nearby candidate (different locus)

# TBLASTN trigger and genomic windows.
TBLASTN_SCOV_MAX="0.60"  # run TBLASTN for a named gene whose best Arabidopsis hit covers < 60% of the Arabidopsis protein
WINDOW_FLANK="100000"    # genomic window = region +/- 100 kb

# NCBI Batch CD-Search.
CDSEARCH_URL="https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
CDSEARCH_POLL_S="10"
CDSEARCH_MAX_WAIT_S="3600"
