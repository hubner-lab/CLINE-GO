library(dplyr)
library(data.table)
library(gradientForest)
library(stringr)
library(qs)
args = commandArgs(trailingOnly=TRUE)
##############
LFMM = args[1]
SIGSNPS = args[2]       # Selected_SNPs.tsv with SNPID column
VCFSNP = args[3]        # .vcfsnp file
REMOVED = args[4]       # .removed file
SAMPLES = args[5]
CLIM_PRESENT_SITE = args[6]
PREDICTORS_SELECTED = args[7] %>% str_split(',') %>% unlist
NTREE = args[8] %>% as.numeric
COR_THRESHOLD = args[9] %>% as.numeric
SPATIAL_CORRECTION = args[10]         # "with" or "without"
MODEL_TYPE = args[11]   # "adaptive" or "random"
OUTPUT = args[12]
# Appended after OUTPUT (positions 13-14), matching rda.R's "append after
# outputs, never insert mid-list" convention — inserting earlier would
# silently shift MODEL_TYPE/OUTPUT. "NULL" (literal string) when
# SPATIAL_CORRECTION == 'without' (workflow/rules/maladaptation.smk's
# nospatial branch never reads these).
DBMEM_VECTORS = args[13]   # PreGEA/tables/spatial/dbmem_vectors.tsv, or "NULL"
DBMEM_SELECTED = args[14]  # PreGEA/tables/varpart/dbmem_selected.tsv, or "NULL"
# Appended 2026-08-22, same "append after outputs" convention as 13-14 above.
RESPONSE_UNIT = if (length(args) >= 15) args[15] else 'individual'  # 'individual' | 'site_frequency'
FREQ_MIN_CALLS = if (length(args) >= 16) as.integer(args[16]) else 2L
##############

if (!RESPONSE_UNIT %in% c('individual', 'site_frequency')) {
  stop(sprintf("ERROR: RESPONSE_UNIT must be 'individual' or 'site_frequency'; got '%s'",
               RESPONSE_UNIT))
}

set.seed(42)
message(paste0('INFO: Building ', MODEL_TYPE, ' Gradient Forest model'))

# Load LFMM genotype matrix
lfmm_dt <- fread(LFMM)
message(paste0('INFO: LFMM matrix: ', nrow(lfmm_dt), ' samples x ', ncol(lfmm_dt), ' SNPs'))

# Load SNP IDs and removed SNPs
sigsnps <- fread(SIGSNPS, header = TRUE)$SNPID
vcfsnp <- fread(VCFSNP, header = FALSE) %>%
  dplyr::select(V1, V2) %>%
  setNames(c('chr', 'pos')) %>%
  dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
  .$chrpos

removed_dt <- fread(REMOVED, header = FALSE)
if (!is.null(removed_dt) && nrow(removed_dt) > 0) {
  removed <- removed_dt %>%
    dplyr::select(V1, V2) %>%
    setNames(c('chr', 'pos')) %>%
    dplyr::mutate(chrpos = paste0(chr, ':', pos)) %>%
    .$chrpos
  mask_adaptive <- (!vcfsnp %in% removed) & (vcfsnp %in% sigsnps)
} else {
  mask_adaptive <- vcfsnp %in% sigsnps
}

n_adaptive <- sum(mask_adaptive)
message(paste0('INFO: Adaptive SNPs: ', n_adaptive))

# Select SNPs based on model type
if (MODEL_TYPE == 'adaptive') {
  snp_subset <- lfmm_dt[, ..mask_adaptive]
  message(paste0('INFO: Using ', ncol(snp_subset), ' adaptive SNPs'))
} else if (MODEL_TYPE == 'random') {
  # [CHANGED 2026-08-10] Size-matched to the adaptive set. Was
  # min(max(n_adaptive, 300), ncol(lfmm_dt)) -- a 300-SNP floor, so a 185-SNP curated set was
  # "controlled" by a 300-SNP random model, 62% larger than the thing it controls for. The
  # neutral model only means anything as a like-for-like comparison, and nothing recorded the
  # size difference downstream. See docs/pipeline_improvement_requests.md request M.
  n_random <- min(n_adaptive, ncol(lfmm_dt))
  if (n_random < 50) {
    message(paste0('WARNING: random model has only ', n_random, ' SNPs (matched to the adaptive ',
                   'set). Gradient Forest turnover functions are unstable at this size; the ',
                   'previous 300-SNP floor existed to avoid this, and was removed because it ',
                   'broke the size match. Treat this neutral model as indicative only.'))
  }
  # [CHANGED 2026-08-10] Seed derived from the run label. Previously the global set.seed(42) at
  # the top of this script was the only seed, and the sampling pool is the full matrix regardless
  # of {run_label} -- so EVERY random_model.qs in a project was the identical draw. Verified on
  # Trifolium86FinalPC3: the three run_labels' response column sets were identical(), 300/300
  # overlap. The pipeline reported three neutral models and had one. The spatial tag is stripped
  # so that {label}_nospatial and {label}_spatial share a draw, keeping the spatial-correction
  # comparison like-for-like.
  run_label <- sub('_(no)?spatial$', '', basename(dirname(OUTPUT)))
  set.seed(42L + sum(utf8ToInt(run_label)))
  random_cols <- sample(names(lfmm_dt), n_random)
  snp_subset <- lfmm_dt[, ..random_cols]
  message(paste0('INFO: Using ', ncol(snp_subset), ' random SNPs (size-matched to the ',
                 n_adaptive, '-SNP adaptive set; run_label "', run_label, '", seed ',
                 42L + sum(utf8ToInt(run_label)), ')'))
} else {
  stop(paste0('ERROR: Unknown MODEL_TYPE: ', MODEL_TYPE))
}

# Guard: these columns become gradientForest() RESPONSE variables, and gradientForest() has no
# na.action. An unimputed LFMM file encodes a missing call as the literal 9, which would enter the
# regression forest as a genotype 3.5x further from the reference homozygote than the real
# alternate homozygote is -- so the forest learns which samples were CALLED, not which carry the
# allele. Checked on the subset actually passed to the model. Same guard as scripts/rda_offset.R:106.
# Keyed on the sentinel 9 rather than on `max > 2`, so a legitimately polyploid dataset
# (sNMF.ploidy > 2, dosages 3,4,...) is not falsely rejected. 9 is LFMM's reserved missing code,
# so it is never a valid dosage in this format regardless of ploidy.
subset_mat <- as.matrix(snp_subset)
if (RESPONSE_UNIT == 'individual' && any(subset_mat == 9, na.rm = TRUE)) {
  frac9 <- mean(subset_mat == 9, na.rm = TRUE)
  stop(sprintf(paste0('Gradient Forest: genotype matrix contains the LFMM missing code 9 in ',
                      '%.1f%% of cells - this is an UNIMPUTED LFMM file. gradientForest() has no ',
                      'na.action and would read 9 as a genotype, so the forest learns which ',
                      'samples were CALLED rather than which carry the allele. Requires an ',
                      'imputed matrix (*_imp*). Check the Snakemake input wiring.'), 100 * frac9))
}
rm(subset_mat)

# Compute maxLevel for gradient forest. Recomputed after site aggregation below when
# RESPONSE_UNIT == 'site_frequency' (n rows drops from samples to sites).
maxLevel <- log2(0.368 * nrow(lfmm_dt) / 2)

# Load samples (also the sample-ID key for the dbMEM join below)
samples <- fread(SAMPLES,
                 colClasses = c("site" = "character",
                                'sample' = 'character',
                                'latitude' = 'numeric',
                                'longitude' = 'numeric'))

# Load climate predictors
env.bio <- fread(CLIM_PRESENT_SITE) %>% dplyr::select(!!PREDICTORS_SELECTED)

# ── Site-allele-frequency aggregation (imputation-sensitivity control) ──────────────────
# Collapse samples to one row per site, with each SNP's value the observed allele frequency
# mean(x[x != 9]) / 2. Nothing is imputed: a site's frequency is computed from the calls it
# actually has, and a site with fewer than FREQ_MIN_CALLS observed calls at a SNP yields NA,
# which drops that SNP entirely (gradientForest() has no na.action). `site_rows` is reused
# below to collapse the dbMEM scores the same way.
site_rows <- NULL
if (RESPONSE_UNIT == 'site_frequency') {
  stopifnot(nrow(snp_subset) == nrow(samples))
  site_levels <- unique(samples$site)
  site_rows   <- lapply(site_levels, function(s) which(samples$site == s))
  names(site_rows) <- site_levels

  geno <- as.matrix(snp_subset)
  geno[geno == 9] <- NA_real_
  n_calls <- t(vapply(site_rows, function(ix) colSums(!is.na(geno[ix, , drop = FALSE])),
                      numeric(ncol(geno))))
  freq <- t(vapply(site_rows, function(ix) colMeans(geno[ix, , drop = FALSE], na.rm = TRUE) / 2,
                   numeric(ncol(geno))))
  freq[n_calls < FREQ_MIN_CALLS] <- NA_real_

  n_in       <- ncol(freq)
  keep_full  <- colSums(is.na(freq)) == 0
  # A SNP constant across every site carries no turnover signal and gradientForest() cannot
  # fit it. Note this is variance ACROSS sites, NOT within-site polymorphism: Fitzpatrick &
  # Keller's "polymorphic in >5 populations" filter assumes an outcrosser and removes 100% of
  # markers on a het-collapsed selfer, where every site is internally fixed but sites differ
  # from each other — which is exactly the informative case.
  keep_var   <- keep_full & (apply(freq, 2, function(v) var(v, na.rm = TRUE)) > 0)
  keep_var[is.na(keep_var)] <- FALSE
  freq <- freq[, keep_var, drop = FALSE]

  message(sprintf('INFO: site-frequency aggregation: %d samples -> %d sites', nrow(geno), nrow(freq)))
  message(sprintf('INFO: observed calls per SNP per site: min %.2f, mean %.2f, max %.2f',
                  min(rowMeans(n_calls)), mean(n_calls), max(rowMeans(n_calls))))
  message(sprintf('INFO: SNPs %d in -> %d kept (%d dropped: a site below %d observed calls; %d constant across sites)',
                  n_in, ncol(freq), sum(!keep_full), FREQ_MIN_CALLS, sum(keep_full & !keep_var)))
  if (ncol(freq) > 0) {
    message(sprintf('INFO: %.1f%% of site-frequency cells are exactly 0 or 1 (near-fixed)',
                    100 * mean(freq == 0 | freq == 1, na.rm = TRUE)))
  }
  if (ncol(freq) < 10) {
    # Do NOT stop(). This model exists only to feed the imputation-sensitivity check, and a
    # diagnostic must never take mode=maladaptation down with it. Write a sentinel the
    # checker recognises, and let every real offset product finish.
    msg <- sprintf(paste0('only %d SNPs survived observed-only aggregation (need >=10) — too ',
                          'few observed calls per site. Lower ',
                          'Maladaptation.methods.gradient_forest.freq_min_calls (currently %d), ',
                          'or accept that this panel cannot support the check.'),
                   ncol(freq), FREQ_MIN_CALLS)
    message('WARNING: Gradient Forest site-frequency model: ', msg)
    qsave(list(status = 'insufficient', reason = msg,
               n_snps = ncol(freq), n_sites = nrow(freq)), OUTPUT)
    message('INFO: wrote insufficient-data sentinel to: ', OUTPUT)
    quit(save = 'no', status = 0)
  }

  snp_subset <- as.data.frame(freq)
  # Predictors are already constant within a site, so take the first row of each.
  env.bio <- env.bio[vapply(site_rows, `[`, integer(1), 1), , drop = FALSE]
  maxLevel <- log2(0.368 * nrow(snp_subset) / 2)
  message(sprintf('INFO: maxLevel recomputed for %d site rows: %.2f', nrow(snp_subset), maxLevel))
}

# Build input matrix and predictor list. Spatial covariates are the climate
# module's forward-selected dbMEM vectors (adespatial::dbmem() + ordiR2step,
# climate.smk Block 2 — scripts/pregea_dbmem.R + scripts/pregea_varpart.R,
# moved out of PreGEA in the module split), NOT the old ad hoc
# vegan::pcnm(dist(raw lon/lat)) + "keep half the positive eigenvalues"
# heuristic this replaces. "Keep all positive MEMs" was considered and
# rejected: docs/rda_research.md's A20/C4 rows always pair dbMEM WITH forward
# selection, never recommend the raw set alone as a GF predictor block.
if (SPATIAL_CORRECTION == 'with') {
  dbmem_dt <- fread(DBMEM_VECTORS, colClasses = c(sample = "character"))
  sel_dt   <- fread(DBMEM_SELECTED)
  selected_mems <- sel_dt[selected == TRUE, mem]

  mem_cols <- grep("^MEM\\d+$", names(dbmem_dt), value = TRUE)
  if (length(mem_cols) == 0 || length(selected_mems) == 0) {
    stop("HARD STOP: dbMEM produced no usable spatial vectors for this project ",
         "(see PreGEA/tables/spatial/dbmem_diagnostics.tsv and ",
         "PreGEA/tables/varpart/dbmem_selected.tsv for why — e.g. too few ",
         "unique sample coordinates, or ordiR2step selected 0 MEMs). A ",
         "'spatial' Gradient Forest model cannot be built on this dataset. ",
         "Set Maladaptation.methods.gradient_forest.spatial_correction to ",
         "'without', or investigate why dbMEM/forward-selection degenerated.")
  }

  # dbmem_vectors.tsv is keyed by sample (preGEA's own row order); this
  # script's LFMM/env.bio matrices are positionally aligned to `samples`
  # (SAMPLES arg) — match() explicitly here rather than assuming row order,
  # since dbMEM comes from a different source file than the other two.
  dbmem_dt <- dbmem_dt[match(samples$sample, sample)]
  stopifnot(!anyNA(dbmem_dt$sample))  # every GF sample must resolve to a dbMEM row

  mem_scores <- as.data.frame(dbmem_dt[, ..selected_mems])
  colnames(mem_scores) <- selected_mems
  message(paste0('INFO: Using ', ncol(mem_scores), ' forward-selected dbMEM axes: ',
                 paste(selected_mems, collapse = ', ')))

  # Collapse dbMEM scores to site rows when the response was aggregated, so the predictor
  # block keeps the same row count as snp_subset. dbMEMs are built from coordinates, which
  # are identical within a site, so the site mean is the site's value.
  if (RESPONSE_UNIT == 'site_frequency') {
    mem_scores <- as.data.frame(t(vapply(
      site_rows,
      function(ix) colMeans(mem_scores[ix, , drop = FALSE]),
      numeric(ncol(mem_scores))
    )))
    colnames(mem_scores) <- selected_mems
  }

  input_matrix <- cbind(env.bio, mem_scores, snp_subset)
  predictor_vars <- c(colnames(env.bio), colnames(mem_scores))
} else {
  input_matrix <- cbind(env.bio, snp_subset)
  predictor_vars <- colnames(env.bio)
}

message(paste0('INFO: Input matrix: ', nrow(input_matrix), ' samples x ', ncol(input_matrix), ' variables'))
message(paste0('INFO: Predictor variables: ', paste(predictor_vars, collapse = ', ')))
message(paste0('INFO: Response variables (SNPs): ', ncol(snp_subset)))

# Run Gradient Forest
gf <- gradientForest(input_matrix,
                     predictor.vars = predictor_vars,
                     response.vars = colnames(snp_subset),
                     ntree = NTREE,
                     maxLevel = maxLevel,
                     trace = TRUE,
                     corr.threshold = COR_THRESHOLD)

# Save model
qsave(gf, OUTPUT)
message(paste0('INFO: Saved ', MODEL_TYPE, ' GF model to: ', OUTPUT))
