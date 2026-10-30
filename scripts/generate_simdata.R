#!/usr/bin/env Rscript
# Generate synthetic VCF and GFF data for pipeline testing
# Creates population-structured SNPs with climate-associated variants

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
# No default: the old relative "data/" regenerated the working SIMDATA fixtures in place when
# run from the repo root, dropping their injected relatedness and preGEA samples (b9be33).
if (length(args) != 1) stop("usage: Rscript generate_simdata.R <output_dir>")
OUTDIR <- args[1]

message("INFO: Generating SIMDATA synthetic dataset")
message(paste0("INFO: Output directory: ", OUTDIR))

#=============================================================================
# CONFIGURATION
#=============================================================================

# Samples per population
samples_per_pop <- 10
populations <- c("Negev", "TelAviv", "Galilee")
n_samples <- samples_per_pop * length(populations)

# Sample names
sample_names <- c(
  paste0("NEG", sprintf("%02d", 1:10)),
  paste0("TAV", sprintf("%02d", 1:10)),
  paste0("GAL", sprintf("%02d", 1:10))
)

# Population assignment (1=Negev/Desert, 2=TelAviv/Mediterranean, 3=Galilee/North)
pop_assign <- rep(1:3, each = samples_per_pop)

# Chromosomes
chromosomes <- paste0("chr", 1:5)
snps_per_chr <- 20
total_snps <- length(chromosomes) * snps_per_chr

# Sample with high missingness (to test filtering)
high_miss_sample <- "NEG10"
high_miss_rate <- 0.7  # 70% missing

#=============================================================================
# SNP CLUSTER DEFINITIONS (for sigSNPs)
#=============================================================================

# Cluster 1: chr1, 10-11Mb - Desert adaptation (GO:0009414 response to water deprivation)
# Cluster 2: chr2, 20-21Mb - Population structured but no enrichment
# Cluster 3: chr3, 15-16Mb - Cold adaptation (GO:0009409 response to cold)

clusters <- list(
  cluster1 = list(chr = "chr1", start = 10000000, end = 11000000, n_snps = 5,
                  pattern = "desert"),  # High alt freq in Negev
  cluster2 = list(chr = "chr2", start = 20000000, end = 21000000, n_snps = 5,
                  pattern = "mixed"),   # Population structured, mixed
  cluster3 = list(chr = "chr3", start = 15000000, end = 16000000, n_snps = 5,
                  pattern = "cold")     # High alt freq in Galilee
)

#=============================================================================
# GENERATE VCF
#=============================================================================

message("INFO: Generating VCF file")

# VCF header
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##fileDate=20260126",
  "##source=SIMDATA_generator",
  "##reference=synthetic",
  paste0("##contig=<ID=chr1,length=50000000>"),
  paste0("##contig=<ID=chr2,length=50000000>"),
  paste0("##contig=<ID=chr3,length=50000000>"),
  paste0("##contig=<ID=chr4,length=50000000>"),
  paste0("##contig=<ID=chr5,length=50000000>"),
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
         paste(sample_names, collapse = "\t"))
)

# Function to generate genotypes based on population and pattern
generate_genotypes <- function(pop_assign, pattern, n_samples) {
  gts <- character(n_samples)

  for (i in 1:n_samples) {
    pop <- pop_assign[i]

    # Set allele frequency based on pattern and population
    if (pattern == "desert") {
      # Desert adaptation: high alt freq in Negev (pop 1)
      p_alt <- switch(pop, 0.8, 0.3, 0.1)
    } else if (pattern == "cold") {
      # Cold adaptation: high alt freq in Galilee (pop 3)
      p_alt <- switch(pop, 0.1, 0.3, 0.8)
    } else if (pattern == "mixed") {
      # Population structured but not extreme
      p_alt <- switch(pop, 0.6, 0.4, 0.2)
    } else {
      # Random background
      p_alt <- runif(1, 0.2, 0.8)
    }

    # Generate diploid genotype
    allele1 <- rbinom(1, 1, p_alt)
    allele2 <- rbinom(1, 1, p_alt)
    gts[i] <- paste0(allele1, "/", allele2)
  }

  return(gts)
}

# Generate SNP records
vcf_records <- list()
snp_id <- 0

for (chr in chromosomes) {
  # Check if this chromosome has a cluster
  chr_cluster <- NULL
  for (cl_name in names(clusters)) {
    if (clusters[[cl_name]]$chr == chr) {
      chr_cluster <- clusters[[cl_name]]
      chr_cluster$name <- cl_name
    }
  }

  # Generate positions for this chromosome
  # Regular SNPs spread across chromosome
  regular_positions <- sort(sample(1000000:45000000, snps_per_chr -
                                     ifelse(is.null(chr_cluster), 0, chr_cluster$n_snps)))

  # Add cluster SNPs if present
  if (!is.null(chr_cluster)) {
    cluster_positions <- sort(sample(chr_cluster$start:chr_cluster$end, chr_cluster$n_snps))
    all_positions <- sort(c(regular_positions, cluster_positions))
  } else {
    all_positions <- regular_positions
    cluster_positions <- integer(0)
  }

  for (pos in all_positions) {
    snp_id <- snp_id + 1

    # Determine if this is a cluster SNP
    is_cluster_snp <- pos %in% cluster_positions
    pattern <- if (is_cluster_snp && !is.null(chr_cluster)) chr_cluster$pattern else "random"

    # Generate genotypes
    gts <- generate_genotypes(pop_assign, pattern, n_samples)

    # Add missing data
    # Random missingness (~5% for most samples)
    for (i in 1:n_samples) {
      if (sample_names[i] == high_miss_sample) {
        # High missingness sample
        if (runif(1) < high_miss_rate) gts[i] <- "./."
      } else {
        # Normal missingness
        if (runif(1) < 0.05) gts[i] <- "./."
      }
    }

    # Create VCF record
    ref <- sample(c("A", "C", "G", "T"), 1)
    alt <- sample(setdiff(c("A", "C", "G", "T"), ref), 1)

    record <- paste(
      chr, pos, paste0("snp", snp_id), ref, alt, ".", "PASS", "DP=100", "GT",
      paste(gts, collapse = "\t"),
      sep = "\t"
    )

    vcf_records[[length(vcf_records) + 1]] <- record
  }
}

#=============================================================================
# LD FILLER SNPS (non-significant, for LD decay testing)
#=============================================================================
# Creates blocks of correlated SNPs at increasing distances from anchor points.
# Different base flip rates per population create distinct LD decay curves:
#   Negev:   slow decay (~150-200kb half-decay, bottleneck simulation)
#   TelAviv: medium decay (~100-120kb)
#   Galilee: fast decay (~60-80kb)

message("INFO: Generating LD filler SNPs")

ld_distances <- c(2000, 5000, 10000, 25000, 50000, 100000, 200000, 350000, 500000)
ld_blocks_per_chr <- 5

# Excluded zones: cluster regions (avoid placing LD blocks in sig SNP clusters)
excluded_zones <- list(
  list(chr = "chr1", start = 9500000, end = 11500000),   # cluster1
  list(chr = "chr2", start = 19500000, end = 21500000),   # cluster2
  list(chr = "chr3", start = 14500000, end = 16500000)    # cluster3
)

# Population-specific flip rate: determines how fast LD decays
flip_rate <- function(distance_bp, population) {
  base_rate <- switch(as.character(population),
    "1" = 0.003,   # Negev: slow decay
    "2" = 0.005,   # TelAviv: medium
    "3" = 0.008    # Galilee: fast decay
  )
  pmin(0.5, base_rate * distance_bp / 1000)
}

# Flip genotype alleles with given probability
flip_genotypes <- function(anchor_gts, distance_bp, pop_assign, n_samples) {
  flipped <- anchor_gts
  for (i in 1:n_samples) {
    rate <- flip_rate(distance_bp, pop_assign[i])
    # Parse genotype
    alleles <- strsplit(anchor_gts[i], "/")[[1]]
    if (any(alleles == ".")) { flipped[i] <- "./."; next }
    a1 <- as.integer(alleles[1])
    a2 <- as.integer(alleles[2])
    # Flip each allele independently
    if (runif(1) < rate) a1 <- 1 - a1
    if (runif(1) < rate) a2 <- 1 - a2
    flipped[i] <- paste0(a1, "/", a2)
  }
  flipped
}

# Check if position is in an excluded zone
is_excluded <- function(chr, pos, zones) {
  for (z in zones) {
    if (z$chr == chr && pos >= z$start && pos <= z$end) return(TRUE)
  }
  FALSE
}

ld_snp_count <- 0
for (chr in chromosomes) {
  # Find safe anchor positions for this chromosome
  for (block in 1:ld_blocks_per_chr) {
    # Try to find a non-excluded anchor position
    max_anchor <- 45000000 - max(ld_distances) - 1000  # leave room for linked SNPs
    attempts <- 0
    repeat {
      anchor_pos <- sample(2000000:max_anchor, 1)
      if (!is_excluded(chr, anchor_pos, excluded_zones)) break
      attempts <- attempts + 1
      if (attempts > 100) { anchor_pos <- 2000000 + block * 5000000; break }
    }

    # Generate anchor genotype (random allele freq per population)
    anchor_gts <- generate_genotypes(pop_assign, "random", n_samples)
    # Add missingness to anchor (low rate, no high-miss sample treatment for LD SNPs)
    for (i in 1:n_samples) {
      if (sample_names[i] == high_miss_sample) {
        if (runif(1) < high_miss_rate) anchor_gts[i] <- "./."
      } else {
        if (runif(1) < 0.03) anchor_gts[i] <- "./."
      }
    }

    # Add anchor SNP
    snp_id <- snp_id + 1
    ld_snp_count <- ld_snp_count + 1
    ref <- sample(c("A", "C", "G", "T"), 1)
    alt <- sample(setdiff(c("A", "C", "G", "T"), ref), 1)
    record <- paste(chr, anchor_pos, paste0("snp", snp_id), ref, alt, ".", "PASS", "DP=100", "GT",
                    paste(anchor_gts, collapse = "\t"), sep = "\t")
    vcf_records[[length(vcf_records) + 1]] <- record

    # Add linked SNPs at each distance
    for (dist in ld_distances) {
      linked_pos <- anchor_pos + dist
      linked_gts <- flip_genotypes(anchor_gts, dist, pop_assign, n_samples)
      snp_id <- snp_id + 1
      ld_snp_count <- ld_snp_count + 1
      ref <- sample(c("A", "C", "G", "T"), 1)
      alt <- sample(setdiff(c("A", "C", "G", "T"), ref), 1)
      record <- paste(chr, linked_pos, paste0("snp", snp_id), ref, alt, ".", "PASS", "DP=100", "GT",
                      paste(linked_gts, collapse = "\t"), sep = "\t")
      vcf_records[[length(vcf_records) + 1]] <- record
    }
  }
}
message(paste0("INFO: Added ", ld_snp_count, " LD filler SNPs (",
               length(chromosomes) * ld_blocks_per_chr, " blocks x ",
               length(ld_distances) + 1, " SNPs)"))

#=============================================================================
# SORT VCF RECORDS AND WRITE
#=============================================================================
# VCF records need to be sorted by chromosome then position
# Parse chr and pos from each record to sort
parse_chr_pos <- function(record) {
  fields <- strsplit(record, "\t")[[1]]
  list(chr = fields[1], pos = as.integer(fields[2]))
}
record_info <- lapply(vcf_records, parse_chr_pos)
record_chrs <- sapply(record_info, `[[`, "chr")
record_pos <- sapply(record_info, `[[`, "pos")
sort_order <- order(match(record_chrs, chromosomes), record_pos)
vcf_records <- vcf_records[sort_order]

# Write VCF
vcf_file <- file.path(OUTDIR, "SIMDATA.vcf")
writeLines(c(vcf_header, unlist(vcf_records)), vcf_file)
message(paste0("INFO: Written VCF with ", length(vcf_records), " SNPs to ", vcf_file))

#=============================================================================
# GENERATE GFF
#=============================================================================

message("INFO: Generating GFF file")

# Gene definitions near SNP clusters
# Cluster 1 (chr1:10-11Mb): Genes with GO:0009414 (response to water deprivation)
# Cluster 2 (chr2:20-21Mb): Genes without relevant GO terms
# Cluster 3 (chr3:15-16Mb): Genes with GO:0009409 (response to cold)

genes <- list(
  # Cluster 1 genes - drought response
  list(chr = "chr1", start = 10050000, end = 10055000, strand = "+",
       id = "gene001", name = "DRY1", go = "GO:0009414"),
  list(chr = "chr1", start = 10150000, end = 10158000, strand = "-",
       id = "gene002", name = "DRY2", go = "GO:0009414"),
  list(chr = "chr1", start = 10350000, end = 10360000, strand = "+",
       id = "gene003", name = "DRY3", go = "GO:0009414"),
  list(chr = "chr1", start = 10550000, end = 10556000, strand = "+",
       id = "gene004", name = "DRY4", go = "GO:0009414"),
  list(chr = "chr1", start = 10750000, end = 10762000, strand = "-",
       id = "gene005", name = "DRY5", go = "GO:0009414"),

  # Cluster 2 genes - no enrichment (various GO terms)
  list(chr = "chr2", start = 20050000, end = 20058000, strand = "+",
       id = "gene006", name = "MIX1", go = "GO:0006412"),
  list(chr = "chr2", start = 20250000, end = 20260000, strand = "-",
       id = "gene007", name = "MIX2", go = "GO:0006396"),
  list(chr = "chr2", start = 20450000, end = 20455000, strand = "+",
       id = "gene008", name = "MIX3", go = "GO:0007049"),
  list(chr = "chr2", start = 20650000, end = 20662000, strand = "+",
       id = "gene009", name = "MIX4", go = "GO:0006811"),
  list(chr = "chr2", start = 20850000, end = 20858000, strand = "-",
       id = "gene010", name = "MIX5", go = "GO:0055085"),

  # Cluster 3 genes - cold response
  list(chr = "chr3", start = 15050000, end = 15060000, strand = "+",
       id = "gene011", name = "CLD1", go = "GO:0009409"),
  list(chr = "chr3", start = 15250000, end = 15258000, strand = "-",
       id = "gene012", name = "CLD2", go = "GO:0009409"),
  list(chr = "chr3", start = 15450000, end = 15462000, strand = "+",
       id = "gene013", name = "CLD3", go = "GO:0009409"),
  list(chr = "chr3", start = 15650000, end = 15655000, strand = "+",
       id = "gene014", name = "CLD4", go = "GO:0009409"),
  list(chr = "chr3", start = 15850000, end = 15860000, strand = "-",
       id = "gene015", name = "CLD5", go = "GO:0009409"),

  # Background genes across chromosomes (no cluster association)
  list(chr = "chr1", start = 5000000, end = 5010000, strand = "+",
       id = "gene016", name = "BKG1", go = "GO:0008150"),
  list(chr = "chr1", start = 25000000, end = 25008000, strand = "-",
       id = "gene017", name = "BKG2", go = "GO:0003674"),
  list(chr = "chr2", start = 5000000, end = 5012000, strand = "+",
       id = "gene018", name = "BKG3", go = "GO:0005575"),
  list(chr = "chr2", start = 35000000, end = 35007000, strand = "-",
       id = "gene019", name = "BKG4", go = "GO:0008150"),
  list(chr = "chr3", start = 5000000, end = 5009000, strand = "+",
       id = "gene020", name = "BKG5", go = "GO:0003674"),
  list(chr = "chr3", start = 35000000, end = 35011000, strand = "-",
       id = "gene021", name = "BKG6", go = "GO:0005575"),
  list(chr = "chr4", start = 10000000, end = 10008000, strand = "+",
       id = "gene022", name = "BKG7", go = "GO:0008150"),
  list(chr = "chr4", start = 25000000, end = 25010000, strand = "-",
       id = "gene023", name = "BKG8", go = "GO:0003674"),
  list(chr = "chr5", start = 10000000, end = 10012000, strand = "+",
       id = "gene024", name = "BKG9", go = "GO:0005575"),
  list(chr = "chr5", start = 30000000, end = 30008000, strand = "-",
       id = "gene025", name = "BKG10", go = "GO:0008150")
)

# GFF header
gff_header <- c(
  "##gff-version 3",
  "##species synthetic test organism",
  "##genome-build SIMDATA v1.0"
)

# Generate GFF records (gene and mRNA features)
gff_records <- list()

for (gene in genes) {
  # Gene feature
  gene_attrs <- paste0(
    "ID=", gene$id, ";",
    "Name=", gene$name, ";",
    "description=", gene$name, " protein;",
    "biotype=protein_coding;",
    "ontology=", gene$go
  )

  gene_record <- paste(
    gene$chr, "SIMDATA", "gene",
    gene$start, gene$end, ".", gene$strand, ".",
    gene_attrs,
    sep = "\t"
  )
  gff_records[[length(gff_records) + 1]] <- gene_record

  # mRNA feature
  mrna_id <- paste0(gene$id, ".1")
  mrna_attrs <- paste0(
    "ID=", mrna_id, ";",
    "Parent=", gene$id, ";",
    "Name=", gene$name, "-mRNA;",
    "description=", gene$name, " transcript;",
    "biotype=protein_coding;",
    "ontology=", gene$go
  )

  mrna_record <- paste(
    gene$chr, "SIMDATA", "mRNA",
    gene$start, gene$end, ".", gene$strand, ".",
    mrna_attrs,
    sep = "\t"
  )
  gff_records[[length(gff_records) + 1]] <- mrna_record

  # CDS feature (simplified: one CDS spanning most of the gene)
  cds_start <- gene$start + 100
  cds_end <- gene$end - 100
  cds_attrs <- paste0(
    "ID=", mrna_id, ".cds;",
    "Parent=", mrna_id
  )

  cds_record <- paste(
    gene$chr, "SIMDATA", "CDS",
    cds_start, cds_end, ".", gene$strand, "0",
    cds_attrs,
    sep = "\t"
  )
  gff_records[[length(gff_records) + 1]] <- cds_record
}

# Write GFF
gff_file <- file.path(OUTDIR, "SIMDATA.gff3")
writeLines(c(gff_header, unlist(gff_records)), gff_file)
message(paste0("INFO: Written GFF with ", length(genes), " genes to ", gff_file))

#=============================================================================
# SUMMARY
#=============================================================================

message("\n=== SIMDATA Summary ===")
message(paste0("Samples: ", n_samples, " (", samples_per_pop, " per population)"))
message(paste0("Populations: ", paste(populations, collapse = ", ")))
message(paste0("High-missingness sample: ", high_miss_sample, " (~", high_miss_rate*100, "% missing)"))
message(paste0("SNPs: ", length(vcf_records), " (", total_snps, " original + ", ld_snp_count, " LD filler)"))
message(paste0("Chromosomes: ", paste(chromosomes, collapse = ", ")))
message("\nSNP Clusters:")
message("  Cluster 1 (chr1:10-11Mb): Desert adaptation, GO:0009414")
message("  Cluster 2 (chr2:20-21Mb): Mixed pattern, no enrichment")
message("  Cluster 3 (chr3:15-16Mb): Cold adaptation, GO:0009409")
message(paste0("\nLD Blocks: ", length(chromosomes) * ld_blocks_per_chr, " blocks (",
               ld_blocks_per_chr, " per chromosome)"))
message("  Per-population LD decay: Negev=slow, TelAviv=medium, Galilee=fast")
message("\nGenes: ", length(genes))
message("  - 5 genes with GO:0009414 (response to water deprivation)")
message("  - 5 genes with GO:0009409 (response to cold)")
message("  - 5 genes with various GO terms (no enrichment)")
message("  - 10 background genes")
message("\nFiles created:")
message(paste0("  ", vcf_file))
message(paste0("  ", gff_file))
message("\nINFO: SIMDATA generation complete")
