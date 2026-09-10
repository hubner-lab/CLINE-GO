#!/usr/bin/env Rscript
# mvp_write_configs.R -- emit one config_MVP{seed}.yaml per row of benchmarks/mvp_seeds.tsv.
#
# Usage:  Rscript benchmarks/mvp_write_configs.R [--seeds=1232548,...] [--k-best=NULL]
#
# Every value that deviates from scripts/clinego.app/inst/config_default.yaml is deviating
# for a stated reason, and the reason is written into the emitted YAML as a comment -- these
# configs are benchmark artifacts that have to be defensible in a paper.
#
# The MVP genome is 20 linkage groups x 50,000 sites = 1 Mb total (Lotterhos 2023, Methods).
# Several pipeline defaults are sized for real genomes and would span an entire linkage group
# here; those are rescaled below.

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x
A <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
  a <- sub("^--", "", a)
  A[[gsub("-", "_", sub("=.*$", "", a))]] <- sub("^[^=]*=", "", a)
}

ROOT      <- Sys.getenv("PIPELINE_ROOT", "/pipeline")
SEEDS_TSV <- file.path(ROOT, "benchmarks/mvp_seeds.tsv")

man <- fread(SEEDS_TSV, colClasses = c(seed = "character"))
if (!is.null(A$seeds)) {
  want <- strsplit(A$seeds, ",")[[1]]
  man <- man[seed %in% want]
  if (!nrow(man)) stop("No manifest rows matched --seeds=", A$seeds)
}

template <- function(seed, arm, arch, fst, r2t, r2s, n_causal, k_best, k_authors) {
  p <- paste0("MVP", seed)
  sprintf('#=============================================================================
# CLINE-GO config -- MVP benchmark replicate, seed %s
#
# Source : Lotterhos 2023, PNAS 120(12):e2220313120, doi:10.1073/pnas.2220313120
#          Data: BCO-DMO doi:10.26008/1912/bco-dmo.889769.1
# Built by: benchmarks/convert_mvp.R (see that script for the scoring contract)
#
# Replicate profile (from the deposit\'s own summary table, see benchmarks/mvp_seeds.tsv):
#   arm          : %s
#   architecture : %s
#   landscape    : structure magnitude meanFst = %.4f
#   confounding  : R2(PC1~temperature) = %.4f   -> bio_1 is the PARTIALLY CONFOUNDED axis
#                  R2(PC1~Env2)        = %.4f   -> bio_2 is the ORTHOGONAL CONTROL axis
#   truth        : %s causal loci at MAF > 0.01
#
# Truth tables live beside the inputs: truth_temp.tsv (score with --traits=bio_1 only),
# truth_sal.tsv (--traits=bio_2 only), truth_any.tsv (--traits=all). Never mix the pairing.
#=============================================================================

Input:
  dir: "data/mvp/%s/"              # MUST end in "/" (common.smk concatenates, no os.path.join)
  vcf: "%s.vcf"
  metadata: "%s_metadata.tsv"
  gff: ""                                  # simulated genome -- annotation/enrichment inert

project_name: %s

#-----------------------------------------------------------------------------
# FILTER
# maf 0.01 (NOT the 0.05 real-data default): the deposit\'s genotype matrix and its
# mutations/truth table are BOTH already filtered at MAF > 0.01 in these same 1000
# individuals, so 0.01 is the no-op that keeps every causal locus in the truth set
# testable. Filtering harder here would recreate exactly the defect that disqualified the
# Laruson benchmark (truth defined at one MAF, VCF at another -> recall capped by
# construction). loci_freq.tsv carries a_freq_subset so any other denominator stays a
# column filter rather than a re-run.
#-----------------------------------------------------------------------------
Filter:
  maf: 0.01
  snp_miss: 0.1
  sample_miss: 0.5
  min_depth: null
  max_depth: null
  het_outlier_sd: null
  relatedness: 0.99
  relatedness_action: keep                 # simulated sibs are real signal, not artifacts

#-----------------------------------------------------------------------------
# LD -- rescaled to the simulated genome.
# Linkage groups are 50 kb, so the 100 kb default window would span an entire LG.
#-----------------------------------------------------------------------------
LD:
  window: 10
  step: 2
  r2: 0.2

#-----------------------------------------------------------------------------
# SNMF
# k_best comes from the deposit\'s own per-seed K column ("number of populations used in
# analyses"), here K = %s, floored at 2 because the pipeline\'s structure plots and
# pop_diff_test need at least two clusters.
#
# Why not the cross-entropy minimum: this landscape is 100 demes on a continuous
# stepping-stone grid, so there is no discrete K to find. Verified empirically on seed
# 1231288 -- sNMF cross-entropy decreases monotonically across the whole K=2..10 range
# (0.5416, 0.5204, 0.5061, 0.4987, 0.4943, 0.4889, 0.4848, 0.4807, 0.4770) with no interior
# minimum, so "pick the minimum" would just return k_end. Using the authors\' published K is
# the same principle the seed selection itself rests on: screen on the statistics the
# original study computed, not on our own re-derivation.
#
# k_best is in any case only the CENTRE of the PreGEA sweeps (LFMM K = k_best +/- k_offset,
# EMMAX/RDA PCs = 0..n_pcs_max); the benchmark exists to find the best setting by scoring
# against truth, not to trust this value.
#-----------------------------------------------------------------------------
sNMF:
  k_start: 2
  k_end: 10
  k_best: %s
  repeats: 10

Map:
  climate_extent: "auto"
  gap: 0.5
  resolution: 0.5                          # unused under Climate.source: custom
  zoom_extent: "NULL"

#-----------------------------------------------------------------------------
# CLIMATE -- custom is MANDATORY for a simulated landscape.
# The default worldclim path would extract real bio_1..bio_19 from real rasters keyed on
# lat/lon, silently attaching real-world climate unrelated to the environment the causal
# loci were simulated against -- a failure mode with no error message.
#   bio_1 = temp_opt (deme temperature optimum)  -- the confounded axis
#   bio_2 = sal_opt  (deme Env2/"salinity" optimum) -- the orthogonal control
# Climate.predictors and PreGEA.predictors MUST equal Climate.custom.columns: they are not
# cross-validated, and a mismatch fails mid-DAG with a column-not-found error.
# No future_table: MVP ships home-deme fitness only, never transplant fitness, so there is
# no honest offset target here. mode=maladaptation is out of scope for this dataset.
#-----------------------------------------------------------------------------
Climate:
  enabled: true
  source: custom
  predictors: "bio_1,bio_2"
  custom:
    present_table: "%s_env_present.tsv"
    columns: "bio_1,bio_2"
    key: site
    grid_resolution: 0.1                   # demes are 0.2 deg apart -- must stay below that
  Varpart:
    response: pcs
    response_var_cutoff: 0.8
    response_max_pcs: 20
    response_min_pcs: 2
    structure_table: qmatrix
    permutations: 999

Population:
  calc_stats: false

Piemap:
  alpha: 0.6
  show_labels: false
  label_size: 8
  pie_scale: 1.0
  use_points: false

LDdecay:
  group_by: "site"
  min_samples: 10
  max_distance: 50                         # LGs are 50 kb; 500 kb default is off-scale
  scope: "both"

#-----------------------------------------------------------------------------
# GEA -- the three methods the benchmark ladder scripts sweep
# (benchmarks/score_ladders.sh: LFMM-K, EMMAX-#PCs, RDA-Condition()-PCs).
# snp_clumping_distance rescaled: the 100 kb default would collapse an entire 50 kb linkage
# group into a single region.
#-----------------------------------------------------------------------------
GEA:
  configs:
    - method: "EMMAX"
      adjust: "bonf"
      threshold: \'0.05\'
    - method: "LFMM"
      adjust: "bonf"
      threshold: \'0.05\'
    - method: "RDA"
      adjust: "bonf"
      threshold: \'0.05\'
  snp_clumping_distance: 5000
  promoter_length: 1000

PreGEA:
  predictors: "bio_1,bio_2"
  k_offset: 2
  n_pcs_max: 10
  Advanced:
    collinearity_r: 0.7
    vif_max: 10.0
    axis_alpha: 0.05
    permutations: 999
  TransferGuard:
    enabled: no

#-----------------------------------------------------------------------------
# GWAS -- second validation arm, free: MVP ships the true simulated phenotypes
# (phen_temp, phen_sal) and home-deme fitness in metadata columns 5+, so phenotype-based
# detection is scorable against the SAME causal-locus truth tables as the GEA arm.
# No missing phenotype values exist, so missing_strategy never engages.
#-----------------------------------------------------------------------------
GWAS:
  configs:
    - method: "EMMAX"
      adjust: "bonf"
      threshold: \'0.05\'
    - method: "BLINK"
      adjust: "bonf"
      threshold: \'0.05\'
  missing_strategy: "MEAN"
  snp_clumping_distance: 5000
  promoter_length: 1000

GFF:
  feature: "mRNA"
  gene_name: "description"
  biotype: "biotype"
  go_field: "NULL"

Enrichment:
  top_terms: 20
  plot_width: 12
  plot_height: 10
  cnet_label: "gene_id"

GEAxGWAS:
  pairwise:
    window_size: 25000
    min_snps: 2

# Declares this project as a SIMULATED benchmark replicate, not a real dataset.
# The Shiny app hides simulated projects from the ordinary project selector and
# surfaces them only in benchmark mode. A declared flag rather than a name
# convention: "starts with MVP" is the implicit kind of rule that let two datasets
# share the name TEST in 2026-07 (CLAUDE.md).
Simulation:
  enabled: true
  source: "MVP / Lotterhos 2023"
  truth_table: "truth_any.tsv"
',
    seed, arm, arch, fst, r2t, r2s, n_causal,
    p, p, p, p,           # Input.dir, vcf, metadata, project_name
    k_authors, k_best,    # sNMF comment block + k_best
    p)                    # Climate.custom.present_table
}

for (i in seq_len(nrow(man))) {
  r <- man[i]
  out <- file.path(ROOT, paste0("config_MVP", r$seed, ".yaml"))
  writeLines(template(r$seed, r$arm, r$architecture, r$meanFst,
                      r$r2_pc1_temp, r$r2_pc1_sal, r$n_causal_maf01,
                      r$k_best, r$K_authors), out)
  message("INFO: wrote ", out, "  (k_best=", r$k_best, ")")
}
message("INFO: ", nrow(man), " config(s) written")
