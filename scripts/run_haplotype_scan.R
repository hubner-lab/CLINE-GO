#!/usr/bin/env Rscript
# run_haplotype_scan.R
# Run crosshap haplotyping across multiple epsilon values for selected regions
# Produces HapObject .qs files and clustree visualization plots
#
# Args:
#   1  selected_regions  - TSV with region_id, chr, start, end, snp_count
#   2  filtered_vcf      - Full filtered VCF (real genotypes)
#   3  ld_vcf            - Imputed VCF for LD computation (or same as filtered)
#   4  metadata_file     - Aligned metadata TSV (site, sample, lat, lon, traits...)
#   5  metadata_type     - 'site' or 'cluster_K{N}'
#   6  clusters_file     - Cluster assignments TSV (or 'NULL')
#   7  mgmin             - DBSCAN MinPts
#   8  minhap            - Minimum haplotype group size
#   9  epsilon_range     - Comma-separated epsilon values
#  10  min_snps          - Minimum SNPs per region (skip if fewer)
#  11  intermediate_dir  - Output dir for .qs HapObject files
#  12  plots_dir         - Output dir for clustree plots

args <- commandArgs(trailingOnly = TRUE)
selected_regions <- args[1]
filtered_vcf     <- args[2]
ld_vcf           <- args[3]
metadata_file    <- args[4]
metadata_type    <- args[5]
clusters_file    <- args[6]
mgmin            <- as.integer(args[7])
minhap           <- as.integer(args[8])
epsilon_range    <- as.numeric(strsplit(args[9], ",")[[1]])
min_snps         <- as.integer(args[10])
intermediate_dir <- args[11]
plots_dir        <- args[12]

library(data.table)
library(crosshap)
library(qs)
library(ggplot2)

message("=== Haplotype Scan ===")
message("Regions: ", selected_regions)
message("Filtered VCF: ", filtered_vcf)
message("LD VCF: ", ld_vcf)
message("Metadata type: ", metadata_type)
message("MGmin: ", mgmin, "  minHap: ", minhap)
message("Epsilon range: ", paste(epsilon_range, collapse = ", "))
message("Min SNPs threshold: ", min_snps)

# Read regions
regions <- fread(selected_regions)
message("Processing ", nrow(regions), " regions")

# Read metadata for phenotype and grouping
meta <- fread(metadata_file, colClasses = c("site" = "character", "sample" = "character"))
sample_col <- names(meta)[2]  # sample column
site_col   <- names(meta)[1]  # site column

# Phenotype: use first numeric trait (column 5+)
if (ncol(meta) >= 5) {
  trait_cols <- names(meta)[5:ncol(meta)]
  first_trait <- trait_cols[1]
} else {
  trait_cols <- character(0)
  first_trait <- NULL
}
message("Using trait '", ifelse(is.null(first_trait), "NONE (dummy)", first_trait), "' for crosshap phenotype input")

# Prepare metadata grouping. crosshap's metadata argument is a CATEGORICAL label per
# individual whose per-haplotype frequencies are contrasted, so a grouping must be a hard
# assignment, never a continuous value.
if (metadata_type == "site") {
  meta_groups <- meta[, .(ind = get(sample_col), group = get(site_col))]
} else if (grepl("^cluster", metadata_type)) {
  if (clusters_file == "NULL") stop("metadata_type '", metadata_type, "' needs a clusters table (arg 6)")
  # clusters_K{k}.tsv is the sNMF Q-matrix: sample, site, C1..CK (extract_clusters.R). It
  # carries ancestry PROPORTIONS only — there is no assignment column — so the assignment is
  # derived here as the argmax over the Q columns, the same rule the app's
  # cluster_pop_summary() uses. (Before 2026-09-13 a grep for "assigned|cluster" matched
  # nothing and the fallback took the LAST column, C{K}: every sample became its own group,
  # labelled by an ancestry fraction — audit B22.)
  clusters <- fread(clusters_file, colClasses = c(sample = "character", site = "character"))
  q_cols <- grep("^C[0-9]+$", names(clusters), value = TRUE)
  if (length(q_cols) == 0) stop("clusters table ", clusters_file, " has no C1..CK ancestry columns")
  if (!"sample" %in% names(clusters)) stop("clusters table ", clusters_file, " has no 'sample' column")
  q_mat <- as.matrix(clusters[, ..q_cols])
  meta_groups <- data.table(
    ind   = clusters$sample,
    group = q_cols[max.col(q_mat, ties.method = "first")]
  )
  message("Cluster grouping from ", length(q_cols), " Q columns: ",
          paste(sprintf("%s=%d", names(table(meta_groups$group)), as.integer(table(meta_groups$group))), collapse = ", "))
} else {
  stop("Unknown metadata_type '", metadata_type, "': expected 'site' or 'cluster_K<k>'")
}

# Initialize region status tracking
region_status <- data.table(
  region_id = character(),
  chr = character(),
  start = integer(),
  end = integer(),
  snp_count = integer(),
  status = character(),
  failure_reason = character(),
  epsilon_succeeded = character()
)

# Helper function to add status row
add_status <- function(rid, chr, start_pos, end_pos, n_snps, status_val, reason = "", eps_succeeded = "") {
  region_status <<- rbind(region_status, data.table(
    region_id = rid,
    chr = as.character(chr),
    start = as.integer(start_pos),
    end = as.integer(end_pos),
    snp_count = as.integer(n_snps),
    status = status_val,
    failure_reason = reason,
    epsilon_succeeded = eps_succeeded
  ))
}

# Helper function to extract region VCF
extract_region_vcf <- function(input_vcf, output_vcf, chr, start, end) {
  cmd <- paste0(
    "awk -v chr='", chr, "' -v start=", start, " -v end=", end,
    " 'BEGIN{OFS=\"\\t\"} /^#/{print; next} $1==chr && $2>=start && $2<=end{print}' ",
    input_vcf, " > ", output_vcf
  )
  system(cmd)
}

# Helper function to count SNPs in VCF
count_vcf_snps <- function(vcf_file) {
  if (!file.exists(vcf_file)) return(0)
  vcf_lines <- readLines(vcf_file, warn = FALSE)
  data_lines <- vcf_lines[!grepl("^#", vcf_lines)]
  length(data_lines)
}

# Process each region
for (i in seq_len(nrow(regions))) {
  reg <- regions[i]
  rid <- reg$region_id
  chr <- reg$chr
  start_pos <- reg$start
  end_pos <- reg$end

  message("\n--- Region ", rid, ": ", chr, ":", start_pos, "-", end_pos, " ---")

  # Create temp directory for this region
  tmp_dir <- file.path(intermediate_dir, paste0("tmp_region_", rid))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  region_vcf <- file.path(tmp_dir, "region.vcf")
  region_ld_vcf <- file.path(tmp_dir, "region_ld.vcf")

  # Extract region from VCF
  message("Extracting region VCF: ", chr, ":", start_pos, "-", end_pos)
  tryCatch({
    extract_region_vcf(filtered_vcf, region_vcf, chr, start_pos, end_pos)
  }, error = function(e) {
    message("ERROR: VCF extraction failed for region ", rid, ": ", e$message)
    add_status(rid, chr, start_pos, end_pos, 0, 'failed_vcf_extraction', e$message, '')
    unlink(tmp_dir, recursive = TRUE)
    return(NULL)
  })

  # Count SNPs in extracted VCF
  n_snps <- count_vcf_snps(region_vcf)
  message("Extracted ", n_snps, " SNPs")

  # Check minimum SNP threshold
  if (n_snps < min_snps) {
    message("WARNING: Region ", rid, " has only ", n_snps, " SNPs (minimum ", min_snps, " required), skipping")
    add_status(rid, chr, start_pos, end_pos, n_snps, 'skipped_low_snps',
               paste0('Insufficient SNPs: ', n_snps, ' < ', min_snps, ' minimum'), '')
    unlink(tmp_dir, recursive = TRUE)
    next
  }

  # Extract region from LD VCF (imputed, for LD computation)
  extract_region_vcf(ld_vcf, region_ld_vcf, chr, start_pos, end_pos)

  # Compute pairwise LD with plink
  ld_prefix <- file.path(tmp_dir, "region_ld")
  cmd_plink <- paste0(
    "plink --vcf ", region_ld_vcf,
    " --const-fid 0",
    " --r2 square --out ", ld_prefix,
    " --allow-extra-chr 2>&1 | grep -v 'Warning\\|written'"
  )
  message("Computing LD matrix...")
  system(cmd_plink, ignore.stdout = TRUE, ignore.stderr = TRUE)

  ld_file <- paste0(ld_prefix, ".ld")
  if (!file.exists(ld_file)) {
    message("WARNING: LD computation failed for region ", rid, ", skipping")
    add_status(rid, chr, start_pos, end_pos, n_snps, 'failed_ld_computation',
               'Plink LD computation failed', '')
    unlink(tmp_dir, recursive = TRUE)
    next
  }

  # Prepare phenotype file (two columns, no header: ind | value)
  pheno_file <- file.path(tmp_dir, "pheno.txt")
  if (!is.null(first_trait)) {
    pheno_data <- meta[, .(ind = get(sample_col), value = get(first_trait))]
    pheno_data <- pheno_data[!is.na(value)]
  } else {
    # No phenotype columns — use dummy values for crosshap
    pheno_data <- meta[, .(ind = get(sample_col), value = 0)]
  }
  fwrite(pheno_data, pheno_file, sep = "\t", col.names = FALSE, quote = FALSE)

  # Prepare metadata file (two columns, no header: ind | group)
  meta_grp_file <- file.path(tmp_dir, "metadata.txt")
  fwrite(meta_groups, meta_grp_file, sep = "\t", col.names = FALSE, quote = FALSE)

  # Run crosshap with comprehensive error handling
  hap_result <- tryCatch({
    message("Reading VCF for crosshap...")
    vcf_obj <- read_vcf(region_vcf)

    # Ensure unique SNP IDs (crosshap requires this)
    if (any(duplicated(vcf_obj$ID)) || all(vcf_obj$ID == ".")) {
      vcf_obj <- dplyr::mutate(vcf_obj, ID = paste0("SNP_", POS))
      message("Renamed SNP IDs to SNP_<POS> for uniqueness")
    }

    message("Reading LD matrix...")
    ld_obj <- read_LD(ld_file, vcf_obj)

    message("Reading phenotype...")
    pheno_obj <- read_pheno(pheno_file)

    message("Reading metadata...")
    meta_obj <- read_metadata(meta_grp_file)

    # Adjust MGmin if needed (don't exceed n_snps / 3)
    local_mgmin <- min(mgmin, max(3L, as.integer(n_snps / 3)))
    if (local_mgmin != mgmin) {
      message("Adjusted MGmin from ", mgmin, " to ", local_mgmin, " (region has ", n_snps, " SNPs)")
    }

    message("Running haplotyping with ", length(epsilon_range), " epsilon values...")
    hap_obj <- run_haplotyping(
      vcf = vcf_obj,
      LD = ld_obj,
      pheno = pheno_obj,
      metadata = meta_obj,
      epsilon = epsilon_range,
      MGmin = local_mgmin,
      minHap = minhap
    )

    hap_obj  # Return result

  }, error = function(e) {
    full_error <- conditionMessage(e)
    if (grepl("must be an array|subscript out of bounds|Varfile", full_error)) {
      message("WARNING: No haplotype structure found for region ", rid,
              " (zero marker groups — SNPs may be in linkage equilibrium)")
      add_status(rid, chr, start_pos, end_pos, n_snps, 'failed_haplotyping',
                 'No haplotype structure (zero marker groups); try higher epsilon or lower MGmin', '')
    } else {
      message("ERROR: Haplotyping failed for region ", rid, ": ", full_error)
      add_status(rid, chr, start_pos, end_pos, n_snps, 'failed_haplotyping',
                 substr(full_error, 1, 200), '')
    }
    return(NULL)
  })

  # Clean up temp files
  unlink(tmp_dir, recursive = TRUE)

  # Check if haplotyping succeeded
  if (is.null(hap_result)) {
    next
  }

  # crosshap HapObject is a named list of per-epsilon results
  # Check which epsilons produced valid results
  if (length(hap_result) == 0) {
    message("WARNING: No marker groups found for region ", rid, " (SNPs may be too sparse or not in LD)")
    add_status(rid, chr, start_pos, end_pos, n_snps, 'failed_haplotyping',
               'No Marker Groups identified at any epsilon value', '')
    next
  }

  # Extract epsilon names that succeeded
  eps_names <- names(hap_result)
  eps_succeeded <- eps_names[grepl('Haplotypes_', eps_names)]

  if (length(eps_succeeded) == 0) {
    message("WARNING: No epsilon values produced results for region ", rid)
    add_status(rid, chr, start_pos, end_pos, n_snps, 'failed_haplotyping',
               'No Marker Groups identified at any epsilon value', '')
    next
  }

  message("SUCCESS: Region ", rid, " haplotyped successfully (", length(eps_succeeded), "/", length(epsilon_range), " epsilons)")

  # Record success
  add_status(rid, chr, start_pos, end_pos, n_snps, 'success', '',
             paste(eps_succeeded, collapse = ','))

  # Save HapObject
  hap_file <- file.path(intermediate_dir, paste0("Region_", rid, "_HapObject.qs"))
  qsave(hap_result, hap_file)
  message("Saved HapObject: ", hap_file)

  # Generate clustree visualization (requires >= 2 epsilon results)
  if (length(eps_succeeded) >= 2) {
    tryCatch({
      message("Generating clustree visualizations...")

      # Marker Group clustree (default: first trait Pheno coloring)
      p_mg <- clustree_viz(hap_result, type = "MG") +
                ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 30))
      clustree_mg_file <- file.path(plots_dir, paste0("Region_", rid, "_clustree_MG"))
      qsave(p_mg, paste0(clustree_mg_file, ".qs"))
      ggsave(paste0(clustree_mg_file, ".png"), p_mg, width = 10, height = 8, dpi = 150)
      ggsave(paste0(clustree_mg_file, ".svg"), p_mg, width = 10, height = 8)
      message("Saved MG clustree: ", clustree_mg_file)

      # Haplotype clustree (default)
      p_hap <- clustree_viz(hap_result, type = "hap") +
                 ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 30))
      clustree_hap_file <- file.path(plots_dir, paste0("Region_", rid, "_clustree_hap"))
      qsave(p_hap, paste0(clustree_hap_file, ".qs"))
      ggsave(paste0(clustree_hap_file, ".png"), p_hap, width = 10, height = 8, dpi = 150)
      ggsave(paste0(clustree_hap_file, ".svg"), p_hap, width = 10, height = 8)
      message("Saved hap clustree: ", clustree_hap_file)

      # Per-trait clustree variants (patch Pheno for each trait)
      if (length(trait_cols) > 1) {
        for (trait in trait_cols) {
          tryCatch({
            message("  Per-trait clustree for '", trait, "'...")
            hap_trait <- hap_result
            for (eps_name in names(hap_trait)) {
              if (!is.null(hap_trait[[eps_name]]$Indfile)) {
                ind_names  <- gsub("\\.", "-", hap_trait[[eps_name]]$Indfile$Ind)
                hap_trait[[eps_name]]$Indfile$Pheno <- meta[[trait]][match(ind_names, meta[[sample_col]])]
              }
            }
            pt_mg <- clustree_viz(hap_trait, type = "MG") +
                       ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 30))
            pt_hap <- clustree_viz(hap_trait, type = "hap") +
                        ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 30))
            ggsave(file.path(plots_dir, paste0("Region_", rid, "_clustree_MG_",  trait, ".png")),
                   pt_mg,  width = 10, height = 8, dpi = 150)
            ggsave(file.path(plots_dir, paste0("Region_", rid, "_clustree_hap_", trait, ".png")),
                   pt_hap, width = 10, height = 8, dpi = 150)
            message("  Saved per-trait clustree for '", trait, "'")
          }, error = function(e) {
            message("  WARNING: Per-trait clustree failed for '", trait, "': ", e$message)
          })
        }
      }
    }, error = function(e) {
      message("WARNING: Clustree visualization failed for region ", rid, ": ", e$message)
    })
  } else {
    message("NOTE: Only ", length(eps_succeeded), " epsilon(s) produced results; ",
            "clustree requires >= 2. Skipping clustree visualization.")
  }
}

# Write status summary file
message("\n=== Writing scan status summary ===")
status_file <- file.path(intermediate_dir, 'scan_status.tsv')

# Compute summary statistics
n_total <- nrow(regions)
n_success <- sum(region_status$status == 'success')
n_skipped <- sum(region_status$status == 'skipped_low_snps')
n_failed_ld <- sum(region_status$status == 'failed_ld_computation')
n_failed_hap <- sum(region_status$status == 'failed_haplotyping')
n_failed_vcf <- sum(region_status$status == 'failed_vcf_extraction')
n_failed_other <- sum(grepl('failed', region_status$status) &
                      !region_status$status %in% c('failed_ld_computation', 'failed_haplotyping', 'failed_vcf_extraction'))

# Add summary row (using special format similar to Overlap_summary.tsv)
summary_row <- data.table(
  region_id = 'SUMMARY',
  chr = 'total',
  start = n_total,
  end = n_success,
  snp_count = n_skipped,
  status = as.character(n_failed_hap + n_failed_ld + n_failed_vcf + n_failed_other),
  failure_reason = '',
  epsilon_succeeded = ''
)

# Write full status table (detail + summary)
full_status <- rbind(region_status, summary_row)
fwrite(full_status, status_file, sep = '\t', quote = FALSE)
message("Wrote status summary: ", status_file)

# Print summary to log
message("\n=== Haplotype scan summary ===")
message("Total regions: ", n_total)
message("Succeeded: ", n_success)
message("Skipped (low SNPs): ", n_skipped)
message("Failed (LD computation): ", n_failed_ld)
message("Failed (haplotyping): ", n_failed_hap)
message("Failed (VCF extraction): ", n_failed_vcf)
if (n_failed_other > 0) message("Failed (other): ", n_failed_other)

# Note: scan_done.flag is managed by Snakemake's touch() — always created when rule completes.
# We write scan_failed.flag as an additional diagnostic when all regions fail.
# The Shiny app checks for actual HapObject files (not just the flag) to determine success.
fail_flag_file <- file.path(intermediate_dir, 'scan_failed.flag')
if (file.exists(fail_flag_file)) unlink(fail_flag_file)

if (n_success > 0) {
  message("\n=== SUCCESS: ", n_success, "/", n_total, " regions completed successfully ===")
} else {
  # Write failure flag with helpful message for human inspection
  failure_msg <- paste0(
    "All ", n_total, " regions failed. Check scan_status.tsv for details.\n\n",
    "Common issues:\n",
    "- Regions have too few SNPs (current threshold: ", min_snps, " SNPs)\n",
    "  -> Increase region size, reduce HAPLOTYPE_SCAN_MIN_SNPS, or use different source\n",
    "- MGmin too high for SNP density (current: ", mgmin, ")\n",
    "  -> Reduce HAPLOTYPE_SCAN_MGMIN\n",
    "- Epsilon range doesn't match LD structure\n",
    "  -> Try wider range or different values\n\n",
    "Failure breakdown:\n",
    "  Skipped (low SNPs): ", n_skipped, "\n",
    "  Failed (LD computation): ", n_failed_ld, "\n",
    "  Failed (haplotyping): ", n_failed_hap, "\n",
    "  Failed (VCF extraction): ", n_failed_vcf, "\n"
  )
  writeLines(failure_msg, fail_flag_file)
  message("\n=== FAILURE: 0/", n_total, " regions succeeded ===")
  message("Wrote failure details to: ", fail_flag_file)

  # Print most common failure reasons
  if (nrow(region_status) > 0) {
    message("\nMost common failure reasons:")
    failure_summary <- region_status[status != 'success', .N, by = failure_reason][order(-N)]
    for (i in seq_len(min(3, nrow(failure_summary)))) {
      reason <- failure_summary$failure_reason[i]
      count <- failure_summary$N[i]
      message("  (", count, "x) ", substr(reason, 1, 80))
    }
  }
}

message("\n=== Haplotype scan complete ===")
