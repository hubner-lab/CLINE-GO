#!/usr/bin/env Rscript
# convert_mvp.R -- convert ONE MVP simulation replicate (Lotterhos 2023, PNAS
# 120(12):e2220313120; data BCO-DMO doi:10.26008/1912/bco-dmo.889769.1) into CLINE-GO's
# standard input set, using Climate.source: custom.
#
# Usage:
#   Rscript convert_mvp.R --seed=1232548 --raw-dir=data/mvp/raw --out-dir=data/mvp/MVP1232548 \
#       [--lon0=34.0] [--lat0=31.0] [--grid-step=0.2]
#
# Inputs (extracted by benchmarks/fetch_mvp.sh; RAW_DIR is never modified):
#   {raw}/genotypes/{seed}_Rout_Gmat_sample.txt.gz   SNPs x 1000 individuals, entries 0/1/2
#                                                    (derived-allele counts), MAF>0.01 in sample
#   {raw}/mutations/{seed}_Rout_muts_full.txt.gz     one row per SNP, ROWS ALIGNED TO Gmat ROWS
#   {raw}/individuals/{seed}_Rout_ind_subset.txt.gz  one row per sampled individual
#
# Outputs (all written to OUT_DIR):
#   MVP{seed}.vcf                 synthesized VCF: CHROM=LG (1-20), POS=per-LG position,
#                                 ID=mutname, REF=A/ALT=T, GT only, no missing data
#   MVP{seed}_metadata.tsv        site, sample, latitude, longitude, phen_temp, phen_sal, fitness
#   MVP{seed}_env_present.tsv     site, bio_1 (=temp_opt), bio_2 (=sal_opt) -- one row per deme
#   truth_temp.tsv                chr, pos, category  <- score with --traits=bio_1 ONLY
#   truth_sal.tsv                 chr, pos, category  <- score with --traits=bio_2 ONLY
#   truth_any.tsv                 chr, pos, category  <- score with --traits=all
#   loci_freq.tsv                 per-locus allele frequency + causal flags, so any MAF
#                                 denominator is a column filter rather than a re-run
#   conversion_provenance.tsv     one row: counts, checksums, and every branch taken
#
# ---------------------------------------------------------------------------------------
# SCORING CONTRACT -- read before using the truth tables.
#
# `causal_temp` and `causal_sal` are INDEPENDENT per-axis classifications: a locus can be
# causal for temperature and neutral-linked for Env2. benchmarks/lib_detection.R's
# load_truth() takes a single `category` column, and call_by_threshold() unions calls across
# traits -- so scoring bio_2's hits against temperature-only truth silently converts real
# detections into false positives. The pairing is fixed and must never be mixed. It also
# differs by method family, because the trait column names do:
#
#   LFMM / EMMAX (one column per predictor, "bio_1" / "bio_2"):
#     --truth=truth_temp.tsv  --traits=bio_1     confounded axis    (R^2(PC1~temp) 0.31-0.57)
#     --truth=truth_sal.tsv   --traits=bio_2     orthogonal control (R^2(PC1~Env2) ~ 0.00)
#     --truth=truth_any.tsv   --traits=all       union
#
#   RDA (multivariate; rda.R writes ONE "climate_multivariate" column, so there is no
#   per-axis p-value and --traits=bio_1 would match nothing):
#     --truth=truth_any.tsv   --traits=all       the only correct pairing
#
#   GWAS (trait columns are the metadata phenotype names):
#     --truth=truth_temp.tsv  --traits=phen_temp
#     --truth=truth_sal.tsv   --traits=phen_sal
#     --truth=truth_any.tsv   --traits=fitness   home-deme fitness is polygenic across BOTH
#                                                axes, so single-axis truth is a category error
#
# Category vocabulary matches lib_detection.R:load_truth(): causal / linked_neutral /
# background_neutral. Only background_neutral loci count as false positives -- which is also
# the convention of the source paper: "the performance statistics were calculated by only
# including truly neutral loci unaffected by selection on linkage groups 11 to 20 and the
# QTNs" (Lotterhos 2023, Methods). Headline precision is therefore precision_strict.
#
# NO future environment table is written. MVP ships home-deme fitness only, not transplant
# fitness, so a synthetic "future" shift would be a fabricated target, not data. The genomic
# offset arm belongs to Capblancq et al. 2026 -- see docs/gea_simulation_benchmarks.md Sec 4.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------ args ----
parse_kv_args <- function(argv) {
  kv <- list()
  for (a in argv) {
    if (!grepl("^--", a)) stop("Unrecognized argument (expected --key=value): ", a)
    a <- sub("^--", "", a)
    k <- sub("=.*$", "", a)
    v <- sub("^[^=]*=", "", a)
    kv[[gsub("-", "_", k)]] <- v
  }
  kv
}
A <- parse_kv_args(commandArgs(trailingOnly = TRUE))

`%||%` <- function(x, y) if (is.null(x) || is.na(x) || identical(x, "")) y else x
SEED      <- A$seed
RAW_DIR   <- A$raw_dir   %||% "data/mvp/raw"
OUT_DIR   <- A$out_dir   %||% file.path("data/mvp", paste0("MVP", SEED))
# Origin is offset by half a raster cell (0.05) so demes land at CELL CENTRES rather than on
# cell boundaries: stage_custom_climate.R:160 hard-errors if two demes resolve to the same
# raster cell, and boundary-landing points are where that collision comes from.
LON0      <- as.numeric(A$lon0      %||% "34.05")
LAT0      <- as.numeric(A$lat0      %||% "31.05")
GRID_STEP <- as.numeric(A$grid_step %||% "0.2")

if (is.null(SEED) || SEED == "") stop("Usage: convert_mvp.R --seed=N [--raw-dir=] [--out-dir=]")

PROJECT   <- paste0("MVP", SEED)
SITES_N   <- 100L   # 10 x 10 landscape grid (Lotterhos 2023 Methods)
INDS_N    <- 1000L  # 10 individuals per deme
N_LG      <- 20L    # "The genome consisted of 20 linkage groups each with 50,000 sites."
LG_LEN    <- 50000L
SEL_LG    <- 10L    # "QTNs could evolve on the first 10 linkage groups"

gmat_path <- file.path(RAW_DIR, "genotypes",   paste0(SEED, "_Rout_Gmat_sample.txt.gz"))
muts_path <- file.path(RAW_DIR, "mutations",   paste0(SEED, "_Rout_muts_full.txt.gz"))
ind_path  <- file.path(RAW_DIR, "individuals", paste0(SEED, "_Rout_ind_subset.txt.gz"))
for (f in c(gmat_path, muts_path, ind_path)) {
  if (!file.exists(f)) stop("ERROR: expected input not found: ", f)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("INFO: converting seed ", SEED, " -> ", OUT_DIR)

prov <- list(seed = SEED, project = PROJECT)

# ------------------------------------------------------------------ read ----
# All three MVP tables were written by R's write.table() with row names: the header line
# carries one FEWER field than every data row, and strings are double-quoted, space-separated.
# fread's own header handling misaligns every column by one on such a file, so detect the
# row-name column explicitly and name the columns ourselves.
read_mvp_gz <- function(path, what) {
  con <- gzfile(path, open = "rt")
  head2 <- readLines(con, n = 2L)
  close(con)
  if (length(head2) < 2L) stop("ERROR: ", what, " has no data rows: ", path)
  hdr  <- scan(text = head2[1], what = "", quiet = TRUE)
  row1 <- scan(text = head2[2], what = "", quiet = TRUE)
  has_rowname <- length(row1) == length(hdr) + 1L
  if (!has_rowname && length(row1) != length(hdr)) {
    stop("ERROR: ", what, " header has ", length(hdr), " fields but row 1 has ", length(row1),
         " -- cannot determine the layout of ", path)
  }
  cn <- if (has_rowname) c(".rowname", hdr) else hdr
  dt <- fread(cmd = paste("zcat", shQuote(path)), header = FALSE, skip = 1L,
              sep = " ", col.names = cn, showProgress = FALSE)
  # .rowname is kept, not dropped: for Gmat it carries the mutname, which is what lets the
  # genotype matrix be joined to the mutations table by KEY rather than by row order.
  attr(dt, "had_rowname") <- has_rowname
  dt
}

message("INFO: reading mutations table")
muts <- read_mvp_gz(muts_path, "muts_full")
message("INFO:   ", nrow(muts), " loci x ", ncol(muts), " cols (row-name column: ",
        attr(muts, "had_rowname"), ")")

message("INFO: reading individuals table")
ind <- read_mvp_gz(ind_path, "ind_subset")
message("INFO:   ", nrow(ind), " individuals x ", ncol(ind), " cols (row-name column: ",
        attr(ind, "had_rowname"), ")")

message("INFO: reading genotype matrix")
gmat <- read_mvp_gz(gmat_path, "Gmat_sample")
message("INFO:   ", nrow(gmat), " rows x ", ncol(gmat), " cols (row-name column: ",
        attr(gmat, "had_rowname"), ")")

need <- function(dt, cols, what) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) {
    stop("ERROR: ", what, " missing required column(s): ", paste(miss, collapse = ", "),
         "\n  present: ", paste(names(dt), collapse = ", "))
  }
}
need(muts, c("LG", "pos_pyslim", "mutname", "a_freq_subset", "causal",
             "causal_temp", "causal_sal"), "muts_full")
need(ind,  c("subpopID", "indID", "x", "y", "temp_opt", "sal_opt",
             "phen_temp", "phen_sal", "fitness"), "ind_subset")

# --------------------------------------------- contract 2: row alignment ----
# The deposit states the mutations rows "correspond to rows in the Genotypes file", but row
# order alone is unverifiable and a silent misalignment would key every truth table to the
# wrong loci. The real file gives us something much stronger: Gmat's row names ARE the
# mutnames and its column names ARE the indIDs, so both joins can be done by key.
if (!attr(gmat, "had_rowname")) {
  stop("ERROR: Gmat has no row-name column -- expected it to carry the mutname")
}
gmat_mutname <- as.character(gmat$.rowname)
gmat[, ".rowname" := NULL]
gmat_indid <- names(gmat)

if (nrow(gmat) != nrow(muts)) {
  stop("ERROR: Gmat rows (", nrow(gmat), ") != muts rows (", nrow(muts), ")")
}
if (ncol(gmat) != nrow(ind)) {
  stop("ERROR: Gmat columns (", ncol(gmat), ") != individuals (", nrow(ind), ")")
}
if (nrow(ind) != INDS_N) warning("Expected ", INDS_N, " sampled individuals, got ", nrow(ind))
prov$n_loci_raw <- nrow(muts)
prov$n_individuals <- nrow(ind)

# --- loci: join Gmat rows to muts rows by mutname -------------------------------------
muts[, mutname := as.character(mutname)]
if (setequal(gmat_mutname, muts$mutname) && !anyDuplicated(gmat_mutname)) {
  locus_join <- "by_mutname"
  muts <- muts[match(gmat_mutname, muts$mutname)]     # muts row i now == Gmat row i
} else if (identical(length(gmat_mutname), nrow(muts))) {
  locus_join <- "by_row_order"
  message("NOTE: Gmat row names do not form the same set as muts$mutname -- falling back to ",
          "the deposit's stated row-order correspondence")
} else {
  stop("ERROR: cannot align Gmat rows to muts rows")
}
message("INFO: locus join mode: ", locus_join)
prov$locus_join_mode <- locus_join

# --- drop loci that fall outside the documented genome --------------------------------
# The simulated genome is N_LG * LG_LEN = 20 * 50,000 = 1,000,000 sites, and pos_pyslim is
# 1-based (contract 1 below proves that against both the LG column and mutname). So the
# valid range is 1 .. 1,000,000, and a locus outside it has no linkage group under the
# (k-1)*50000+1 .. k*50000 layout. Contract 1 refuses to guess a coordinate system, so
# without this guard the whole replicate fails to convert.
#
# Measured across all 92 benchmark seeds, exactly TWO such loci exist, one per end, each in
# a different replicate and each NON-CAUSAL:
#   seed 1232620  "1-0"          pos 0          -> LG 0, categorised neutral-linked
#   seed 1233058  "21-1000001"   pos 1,000,001  -> LG 21, categorised neutral
# Both are single off-by-one boundary artefacts, 1 locus out of 25,576 and 11,247.
#
# They are DROPPED rather than remapped. Remapping either one (to LG 1 pos 1, or LG 20 pos
# 50,000) would change the coordinate contract that every truth join rests on, for the sake
# of one non-causal locus in one replicate -- and in the pos-1000001 case it would also
# collide with a real locus if one already sits at that position.
#
# The drop happens WHILE muts and gmat are still row-aligned, so the identical mask goes to
# both -- the same discipline contract 3 uses for duplicate collapse. The counts are written
# to the provenance for every seed (0 almost everywhere), so any count table can see that
# these two replicates carry one locus fewer than the deposit reports.
GENOME_LEN  <- N_LG * LG_LEN
out_of_range <- muts$pos_pyslim < 1L | muts$pos_pyslim > GENOME_LEN
prov$n_loci_below_start <- sum(muts$pos_pyslim < 1L)
prov$n_loci_past_end    <- sum(muts$pos_pyslim > GENOME_LEN)
if (any(out_of_range)) {
  print(muts[out_of_range, .(mutname, LG, pos_pyslim, causal)])
  message("NOTE: dropping ", sum(out_of_range), " locus/loci outside the genome ",
          "(1 .. ", GENOME_LEN, ") -- no linkage group exists for them under the ",
          "1-based layout. This seed's SNP count will be ", sum(out_of_range),
          " below the deposit's published figure.")
  if (any(as.logical(muts$causal[out_of_range]))) {
    stop("ERROR: a CAUSAL locus sits outside the genome -- refusing to drop it silently.")
  }
  muts <- muts[!out_of_range]
  gmat <- gmat[!out_of_range]
  if (nrow(gmat) != nrow(muts)) stop("ERROR: muts/Gmat misaligned after the out-of-range drop")
}

# --- individuals: join Gmat columns to ind rows by indID ------------------------------
ind_ids <- as.character(ind$indID)
if (setequal(gmat_indid, ind_ids) && !anyDuplicated(gmat_indid)) {
  individual_join <- "by_indID"
  ind <- ind[match(gmat_indid, ind_ids)]              # ind row j now == Gmat column j
} else {
  individual_join <- "by_column_order"
  message("NOTE: Gmat column names are not indIDs -- falling back to positional order")
}
message("INFO: individual join mode: ", individual_join)
prov$individual_join_mode <- individual_join

# ------------------------------- contract 1: coordinate system of the loci --
# `mutname` is documented as "a unique ID for the mutation based on linkage group and
# position" and is written as "<LG>-<pos>". It is therefore the primary evidence for the
# coordinate convention, rather than inferring one from the range of pos_pyslim. Cross-check
# it against pos_pyslim, which muts_metadata.md defines as "position as output from pyslim,
# which is 1 less than SliM output".
if (!is.numeric(muts$LG)) muts[, LG := as.integer(LG)]
lg_rng <- range(muts$LG)
if (lg_rng[1] < 1L || lg_rng[2] > N_LG) {
  stop("ERROR: LG outside the documented 1..", N_LG, " range: [", lg_rng[1], ", ", lg_rng[2], "]")
}

# mutname is "<LG>-<position>", optionally with a "_<n>" suffix marking a recurrent mutation
# that arose at the same map position at a different time -- the exact case muts_metadata.md
# warns about, and the one contract 3 below has to collapse.
MUTNAME_RE <- "^([0-9]+)-([0-9]+)(_[0-9]+)?$"
if (!all(grepl(MUTNAME_RE, muts$mutname))) {
  stop("ERROR: mutname does not parse as '<LG>-<pos>[_<n>]': ",
       paste(head(muts$mutname[!grepl(MUTNAME_RE, muts$mutname)], 3), collapse = ", "))
}
lg_mn  <- as.integer(sub(MUTNAME_RE, "\\1", muts$mutname))
pos_mn <- as.numeric(sub(MUTNAME_RE, "\\2", muts$mutname))
n_recurrent <- sum(nzchar(sub(MUTNAME_RE, "\\3", muts$mutname)))
prov$n_recurrent_mutname <- n_recurrent
message("INFO: ", n_recurrent, " locus name(s) carry a recurrent-mutation suffix")

if (!identical(lg_mn, muts$LG)) {
  stop("ERROR: the LG encoded in mutname disagrees with the LG column -- ",
       sum(lg_mn != muts$LG), " row(s) differ. Refusing to guess the coordinate system.")
}
off <- unique(pos_mn - muts$pos_pyslim)
prov$mutname_minus_pyslim <- paste(sort(off), collapse = ",")
if (length(off) != 1L) {
  stop("ERROR: mutname position minus pos_pyslim is not constant (values: ",
       prov$mutname_minus_pyslim, ") -- the two columns are not the same coordinate.")
}
message("INFO: mutname position = pos_pyslim + ", off, " (same coordinate, cross-validated)")

print(muts[, .(n = .N, pos_min = min(pos_pyslim), pos_max = max(pos_pyslim)), by = LG][order(LG)])

# Determine the coordinate scale, and cross-check it against the LG column rather than
# trusting a range heuristic alone.
#
# pos_pyslim is 1-BASED within the 1 Mb genome, so linkage group k spans
# (k-1)*50000+1 .. k*50000 and the boundary locus belongs to the LOWER group. Confirmed by
# the deposit itself: seed 1232908 carries mutname "5-250000" with LG 5, and seed 1232548
# carries "13-650000" with LG 13 -- both exactly on a boundary, both assigned to the lower
# group. A 0-based reading (pos %/% LG_LEN + 1) puts those loci one group too high, and
# adding 1 to build the VCF position pushes them to 50001, outside the contig.
if (max(pos_mn) > LG_LEN) {
  pos_scale <- "genome_global"
  lg_from_pos <- as.integer((muts$pos_pyslim - 1L) %/% LG_LEN) + 1L
  if (!identical(lg_from_pos, muts$LG)) {
    print(muts[lg_from_pos != LG][1:min(5, .N), .(mutname, LG, pos_pyslim)])
    stop("ERROR: LG derived from the global position disagrees with the LG column for ",
         sum(lg_from_pos != muts$LG), " row(s) -- the ", LG_LEN,
         "-site linkage-group layout does not hold for this seed.")
  }
  muts[, pos := as.integer(pos_pyslim - (LG - 1L) * LG_LEN)]
} else {
  pos_scale <- "per_lg"
  muts[, pos := as.integer(pos_pyslim)]
}
message("INFO: source positions are ", pos_scale,
        " -> VCF POS emitted per linkage group, 1..", LG_LEN)
prov$pos_scale <- pos_scale
if (!all(muts$pos >= 1L & muts$pos <= LG_LEN)) {
  bad <- muts[pos < 1L | pos > LG_LEN][1:min(5, .N)]
  print(bad[, .(LG, mutname, pos_pyslim, pos)])
  stop("ERROR: derived POS outside 1..", LG_LEN, " under the '", pos_scale,
       "' interpretation -- see rows above.")
}
muts[, chr := LG]

# ------------------------------- the decisive alignment proof ---------------
# Counting rows only proves the tables are the same length. This proves row i of the genotype
# matrix really is locus i of the mutations table: a_freq_subset is documented as the derived
# allele frequency in exactly these 1000 individuals, and the matrix holds derived-allele
# counts, so the two must agree to rounding. If they don't, every chr:pos in every truth
# table points at the wrong locus -- silently, with nothing downstream to catch it.
G <- as.matrix(gmat)
storage.mode(G) <- "integer"
if (anyNA(G)) stop("ERROR: genotype matrix contains NA -- the source is documented as complete")
if (!all(G %in% 0:2)) {
  stop("ERROR: genotype matrix holds values outside {0,1,2}: ",
       paste(head(sort(unique(as.vector(G))), 10), collapse = ", "))
}
rm(gmat)
af_from_G <- rowSums(G) / (2 * ncol(G))
af_dev <- max(abs(af_from_G - muts$a_freq_subset))
prov$af_max_deviation <- signif(af_dev, 4)
message(sprintf("INFO: allele-frequency alignment check -- max |AF(Gmat) - a_freq_subset| = %.3g",
                af_dev))
AF_TOL <- 1e-3
if (af_dev > AF_TOL) {
  worst <- order(abs(af_from_G - muts$a_freq_subset), decreasing = TRUE)[1:5]
  print(data.table(mutname = muts$mutname[worst], from_gmat = af_from_G[worst],
                   a_freq_subset = muts$a_freq_subset[worst]))
  stop("ERROR: genotype matrix and mutations table are misaligned (max deviation ", af_dev,
       " > ", AF_TOL, "). Refusing to emit a truth table keyed to the wrong loci.")
}

# ------------------------------------ contract 3: duplicate (chr,pos) -------
# muts_metadata.md: "a few rare mutations may arise on the same 'position' on the SLiM
# genetic map but at different times." The scoring harness keys on chr:pos
# (lib_detection.R:28), and this pipeline is biallelic-only, so duplicates must collapse to
# one row. Keep the causal row wherever a duplicate group contains one, and apply the
# identical row mask to Gmat so alignment survives.
muts[, .row := .I]
muts[, .is_causal := as.logical(causal)]
setorder(muts, chr, pos, -.is_causal)
keep_rows <- muts[, .SD[1L], by = .(chr, pos)]$.row
n_dropped <- nrow(muts) - length(keep_rows)
if (n_dropped > 0L) {
  message("INFO: collapsing ", n_dropped, " duplicate (chr,pos) row(s) -- keeping the causal ",
          "row within each group")
}
prov$n_duplicate_collapsed <- n_dropped

keep_rows <- sort(keep_rows)
muts <- muts[order(.row)][.row %in% keep_rows]
G <- G[keep_rows, , drop = FALSE]
stopifnot(nrow(G) == nrow(muts))

# Final ordering for the VCF: sorted by (CHROM, POS), applied to both tables together.
ord  <- order(muts$chr, muts$pos)
muts <- muts[ord]
G <- G[ord, , drop = FALSE]
if (anyDuplicated(muts[, .(chr, pos)])) stop("ERROR: duplicate (chr,pos) survived the collapse")
prov$n_loci <- nrow(muts)
message("INFO: ", nrow(muts), " loci retained after de-duplication")

# ------------------------------------ contract 4/5: metadata + coordinates --
# subpopID / indID are numeric, which is exactly the fread colClasses footgun documented in
# CLAUDE.md. Prefixing makes every downstream ID unambiguously character.
ind[, site   := sprintf("deme_%03d", as.integer(subpopID))]
ind[, sample := sprintf("ind_%05d", as.integer(indID))]
if (anyDuplicated(ind$sample)) stop("ERROR: individual IDs are not unique after prefixing")

x_rng <- range(ind$x); y_rng <- range(ind$y)
message(sprintf("INFO: deme grid x [%g, %g], y [%g, %g]", x_rng[1], x_rng[2], y_rng[1], y_rng[2]))
ind[, longitude := LON0 + GRID_STEP * (x - x_rng[1])]
ind[, latitude  := LAT0 + GRID_STEP * (y - y_rng[1])]
if (max(ind$longitude) > 180 || min(ind$longitude) < -180 ||
    max(ind$latitude)  >  90 || min(ind$latitude)  <  -90) {
  stop("ERROR: derived coordinates leave the valid WGS84 range -- lower --grid-step or --lon0/--lat0")
}
prov$coord_box <- sprintf("lon [%.3f, %.3f] lat [%.3f, %.3f]",
                         min(ind$longitude), max(ind$longitude),
                         min(ind$latitude),  max(ind$latitude))
prov$grid_step_deg <- GRID_STEP

# Env values are the deme optima that drove selection. They must come from this table and
# never from a raster lookup -- see docs/gea_simulation_benchmarks.md Sec 5.
site_env <- unique(ind[, .(site, longitude, latitude, bio_1 = temp_opt, bio_2 = sal_opt)])
if (anyDuplicated(site_env$site)) {
  stop("ERROR: temp_opt/sal_opt/coordinates are not constant within a deme -- cannot derive a ",
       "single per-site environment row. Inspect ", ind_path, ".")
}
setorder(site_env, site)
message("INFO: ", nrow(site_env), " distinct demes")
if (nrow(site_env) != SITES_N) warning("Expected ", SITES_N, " demes, got ", nrow(site_env))
prov$n_sites <- nrow(site_env)

# stage_custom_climate.R:160 errors when two demes land in one raster cell, so deme spacing
# along each axis must exceed Climate.custom.grid_resolution.
axis_sep <- function(v) { u <- sort(unique(v)); if (length(u) < 2L) Inf else min(diff(u)) }
min_sep <- min(axis_sep(site_env$longitude), axis_sep(site_env$latitude))
message(sprintf("INFO: minimum deme separation %.3f deg per axis -- must exceed %s",
                min_sep, "Climate.custom.grid_resolution"))
prov$min_coord_sep_deg <- round(min_sep, 4)

# Traits land in metadata columns 5+ (prepare_phenotypes.R:26 reads them positionally).
metadata <- ind[, .(site, sample, latitude, longitude, phen_temp, phen_sal, fitness)]
setorder(metadata, sample)
meta_path <- file.path(OUT_DIR, paste0(PROJECT, "_metadata.tsv"))
fwrite(metadata, meta_path, sep = "\t")
message("INFO: wrote ", meta_path)

env_path <- file.path(OUT_DIR, paste0(PROJECT, "_env_present.tsv"))
fwrite(site_env[, .(site, bio_1, bio_2)], env_path, sep = "\t")
message("INFO: wrote ", env_path)

# ------------------------------------------------------------------ VCF -----
# Synthesize a minimal, valid, biallelic diploid VCF from the 0/1/2 derived-allele counts
# (G was built, validated and de-duplicated above). There is no missing data in the source
# matrix; vcf2lfmm.R's zero-missingness guard handles the downstream consequence
# (LEA::impute errors on a matrix with nothing to impute).
GT_LOOKUP <- c("0/0", "0/1", "1/1")
gt_chr <- matrix(GT_LOOKUP[G + 1L], nrow = nrow(G), ncol = ncol(G))

# Gmat columns are in ind row order: either natively (positional join) or because the
# by-ID branch above reordered them to match.
sample_order <- ind$sample
vcf_path <- file.path(OUT_DIR, paste0(PROJECT, ".vcf"))
hdr <- c(
  "##fileformat=VCFv4.2",
  sprintf("##source=convert_mvp.R seed=%s (Lotterhos 2023, BCO-DMO doi:10.26008/1912/bco-dmo.889769.1)", SEED),
  sprintf("##contig=<ID=%d,length=%d>", seq_len(N_LG), LG_LEN),
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
          sample_order), collapse = "\t")
)
fixed <- paste(muts$chr, muts$pos, muts$mutname, "A", "T", ".", "PASS", ".", "GT", sep = "\t")
geno  <- do.call(paste, c(as.data.frame(gt_chr, stringsAsFactors = FALSE), sep = "\t"))
con <- file(vcf_path, open = "wt")
writeLines(hdr, con)
writeLines(paste(fixed, geno, sep = "\t"), con)
close(con)
message("INFO: wrote ", vcf_path, " (", nrow(muts), " variants x ", length(sample_order), " samples)")
prov$vcf <- basename(vcf_path)

# ---------------------------------------------------------------- truth -----
# Source levels -> harness vocabulary. "neutral-linked" means a neutral locus sitting on a
# linkage group where QTNs could evolve (LG 1-10), NOT a distance window.
LEVEL_MAP <- c("causal" = "causal",
               "neutral-linked" = "linked_neutral",
               "neutral" = "background_neutral")

map_truth <- function(col, label) {
  v <- as.character(muts[[col]])
  # 1-trait seeds carry no causal loci on the unused axis. Whether the deposit encodes that
  # as "neutral" or as NA varies, and an NA here must not become an unclassified locus: a
  # locus with no effect on this trait IS neutral for it, and its linkage-group membership
  # is what decides linked vs background -- the same rule truth_any uses.
  if (anyNA(v)) {
    n_na <- sum(is.na(v))
    v[is.na(v)] <- ifelse(muts$chr[is.na(v)] <= SEL_LG, "neutral-linked", "neutral")
    message("NOTE: ", col, " had ", n_na, " NA value(s) -- classified as neutral for this ",
            "axis, linked vs background by linkage group (LG 1-", SEL_LG, " carry QTNs)")
  }
  unknown <- setdiff(unique(v), names(LEVEL_MAP))
  if (length(unknown)) {
    stop("ERROR: unexpected level(s) in ", col, ": ", paste(unknown, collapse = ", "),
         " -- expected only ", paste(names(LEVEL_MAP), collapse = ", "))
  }
  out <- data.table(chr = muts$chr, pos = muts$pos, category = unname(LEVEL_MAP[v]))
  p <- file.path(OUT_DIR, paste0("truth_", label, ".tsv"))
  fwrite(out, p, sep = "\t")
  tab <- table(factor(out$category, levels = unname(LEVEL_MAP)))
  message(sprintf("INFO: wrote %s -- %d causal, %d linked_neutral, %d background_neutral",
                  p, tab[["causal"]], tab[["linked_neutral"]], tab[["background_neutral"]]))
  as.list(tab)
}
t_temp <- map_truth("causal_temp", "temp")
t_sal  <- map_truth("causal_sal",  "sal")

# truth_any: causal on EITHER axis; otherwise linked (selected LGs 1-10) vs background.
any_cat <- ifelse(as.logical(muts$causal), "causal",
                  ifelse(muts$chr <= SEL_LG, "linked_neutral", "background_neutral"))
truth_any <- data.table(chr = muts$chr, pos = muts$pos, category = any_cat)
fwrite(truth_any, file.path(OUT_DIR, "truth_any.tsv"), sep = "\t")
t_any <- as.list(table(factor(any_cat, levels = unname(LEVEL_MAP))))
message(sprintf("INFO: wrote truth_any.tsv -- %d causal, %d linked_neutral, %d background_neutral",
                t_any[["causal"]], t_any[["linked_neutral"]], t_any[["background_neutral"]]))

if (t_temp[["causal"]] == 0L) warning("truth_temp.tsv has zero causal loci for seed ", SEED)
if (t_sal[["causal"]]  == 0L) {
  message("NOTE: truth_sal.tsv has zero causal loci -- expected for a 1-trait seed, where Env2 ",
          "is a pure negative control")
}

prov$causal_temp <- t_temp[["causal"]]; prov$linked_temp <- t_temp[["linked_neutral"]]
prov$background_temp <- t_temp[["background_neutral"]]
prov$causal_sal  <- t_sal[["causal"]];  prov$linked_sal  <- t_sal[["linked_neutral"]]
prov$background_sal <- t_sal[["background_neutral"]]
prov$causal_any  <- t_any[["causal"]];  prov$linked_any  <- t_any[["linked_neutral"]]
prov$background_any <- t_any[["background_neutral"]]

# ------------------------------------------------------------ loci_freq -----
# a_freq_subset is the derived-allele frequency in the same 1000 individuals the Gmat
# describes. Carrying it (and MAF) through means any alternative MAF denominator -- e.g. the
# 0.05 used on real data -- is a column filter on this table, never a pipeline re-run.
loci <- muts[, .(chr, pos, mutname, a_freq_subset,
                 maf = pmin(a_freq_subset, 1 - a_freq_subset),
                 causal = as.logical(causal), causal_temp, causal_sal)]
fwrite(loci, file.path(OUT_DIR, "loci_freq.tsv"), sep = "\t")
prov$maf_min <- round(min(loci$maf), 5)
prov$n_causal_maf05 <- loci[causal & maf >= 0.05, .N]
message("INFO: wrote loci_freq.tsv -- min MAF ", prov$maf_min,
        "; causal loci surviving MAF>=0.05: ", prov$n_causal_maf05,
        " of ", loci[causal == TRUE, .N])

# ----------------------------------------------------------- provenance -----
prov$gmat_md5 <- tools::md5sum(gmat_path)[[1]]
prov$muts_md5 <- tools::md5sum(muts_path)[[1]]
prov$ind_md5  <- tools::md5sum(ind_path)[[1]]
prov$converted_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
fwrite(as.data.table(prov), file.path(OUT_DIR, "conversion_provenance.tsv"), sep = "\t")
message("INFO: wrote conversion_provenance.tsv")
message("INFO: convert_mvp.R complete for seed ", SEED)
