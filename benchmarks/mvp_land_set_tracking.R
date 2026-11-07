#!/usr/bin/env Rscript
# =============================================================================
# mvp_land_set_tracking.R -- what each method's SNP set follows, on all 1,800 replicates of the three MVP
# landscapes (plan ~/.claude/plans/let-s-run-the-lanpscape-dynamic-pony.md, Phase 1). Base: the per-replicate loop
# of mvp_land_set_profile.R (deme allele frequencies from the VCF, neutral PC1, per-SNP |r|), unchanged.
#
# Per replicate x set (LFMM2, EMMAX, RDA, the 1/3, 2/3 and 3/3 sets, causal loci) x part:
#   all                  the whole set
#   class parts          the deposit's truth (loci_freq.tsv) and the paper's linked rule (option B of
#                        mvp_panel_tables.R, as detection_per_seed.tsv counts TP / FP): causal temperature, causal
#                        salinity, causal both (pleiotropy replicates), causal unassigned (causal in truth_any but in
#                        neither trait list, e.g. MVP1231110 3:23169 at MAF 0.0105), linked (non-causal within the
#                        replicate's LD distance d_bp of a causal locus on its chromosome), neutral (the rest = FP)
#   linked by trait      linked temperature / linked salinity: every causal locus within d_bp is causal for that trait
#                        only; linked mixed: loci of both traits (or pleiotropic, or unassigned) within d_bp
#   neutral by place     neutral, causal chromosomes (linkage groups 1-10, beyond d_bp) / neutral, neutral chromosomes
#                        (11-20)
#   test parts           LFMM2 and EMMAX only: picked by the temperature test only, the salinity test only, or both.
#                        The per-predictor top lists are rebuilt from the p tables by the toppct rule of
#                        compute_pval_threshold() (scripts/R/utils/pval_threshold.R): N = ceiling(share * rows / 100
#                        - 1e-9), NA rows counted, inclusive cutoff, share 0.25 per predictor (as mvp_land_lfmm_check.py).
# n and the median over the part's SNPs of |r| across demes between allele frequency and temperature (bio_1),
# salinity (bio_2) and neutral PC1 (first PC of the deme frequencies of SNPs neutral for both traits); and the PCA of
# the part (pca1_*): |r| of the PC1 of the part's deme frequencies with temperature and salinity, and its variance share.
# Per replicate: r2 of PC1 with each predictor, selection regime, demographic scenario, K; and per replicate x set
# the outcome under the median of the three offset methods (rank_long.tsv): tie-split weight in ranks 4-6
# (mvp_rank_ties.R), loss against the replicate's best set.
#
# Asserts: set-level values = set_profile.tsv (the earlier run of the same code); causal + linked = TP and neutral =
# FP of detection_per_seed.tsv; union of the test lists = the set file; % in ranks 4-6 per landscape = t_ranks.tsv;
# 600 replicates per landscape.
#
# Usage: Rscript benchmarks/mvp_land_set_tracking.R <out_dir>
#        TRACK_N=<k>: pilot on the first k replicates per landscape (asserts on the full counts skipped)
#        TRACK_CORES (default 24)
# =============================================================================
suppressMessages({library(data.table); library(parallel)})
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
OUT <- args[1]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
setDTthreads(1)
NPILOT <- as.integer(Sys.getenv("TRACK_N", "0"))
CORES  <- as.integer(Sys.getenv("TRACK_CORES", "24"))
EVAL <- "/pipeline/benchmarks/mvp_eval"
LMS  <- file.path(EVAL, "landscape_ms")
source("/pipeline/benchmarks/mvp_rank_ties.R")

LANDS <- data.table(
    land = c("SS-Clines", "SS-Mtn", "Est-Clines"),
    gal  = c("figures_ld_gallery14", "figures_land_ssmtn_gallery", "figures_land_estclines_gallery"),
    rank = c("figures_ld_offset14/rank", "figures_land_ssmtn_offset15_ssmtn/rank",
             "figures_land_estclines_offset15_estclines/rank"))
PS <- fread(file.path(EVAL, "landscape_compare/per_seed.tsv"), colClasses = c(seed = "character"))
SEEDS <- unique(PS[family == "TP", .(landscape, seed, arch)])
stopifnot(all(SEEDS[, .N, by = landscape]$N == 600L), nrow(SEEDS) == 1800L)
if (NPILOT > 0L) SEEDS <- SEEDS[, .SD[seq_len(NPILOT)], by = landscape]
MAN <- fread("/pipeline/benchmarks/mvp_seeds.tsv", colClasses = c(seed = "character"))[seed %in% SEEDS$seed]
MAN[, regime := fifelse(grepl("unequal-S$", architecture), "unequal-S", fifelse(grepl("equal-S$", architecture), "equal-S", NA_character_))]
stopifnot(nrow(MAN) == nrow(SEEDS), !anyNA(MAN$regime), !anyDuplicated(MAN$seed))
SEEDS <- merge(SEEDS, MAN[, .(seed, regime, demog = demog_level_sub, k_best)], by = "seed", sort = FALSE)

LDBP <- rbindlist(lapply(c("work/ld_window_lab/p4/seed_bp_out600.tsv", "work/mvp_land/results/seed_ld_bp_ssmtn.tsv",
                           "work/mvp_land/results/seed_ld_bp_estclines.tsv"),
                         function(f) fread(file.path("/pipeline", f), colClasses = c(seed = "character"))))
stopifnot(!anyDuplicated(LDBP$seed), all(SEEDS$seed %in% LDBP$seed), !anyNA(LDBP$d_bp))
SEEDS <- merge(SEEDS, LDBP, by = "seed", sort = FALSE)

SETS  <- c(LFMM = "gea.solo_lfmm", EMMAX = "gea.solo_emmax", RDA = "gea.solo_rda", truth = "truth",
           `1/3 methods` = "gea.union", `2/3 methods` = "gea.best", `3/3 methods` = "gea.intersect3")
PANELS <- c("LFMM", "EMMAX", "RDA", "1/3 methods", "2/3 methods", "3/3 methods")   # the six sets of the rank tables
TESTS <- c("LFMM", "EMMAX")
CLASSES <- c("causal temperature", "causal salinity", "causal both", "causal unassigned", "linked", "neutral")
LSPLIT  <- c("linked temperature", "linked salinity", "linked mixed")
NSPLIT  <- c("neutral, causal chromosomes", "neutral, neutral chromosomes")
TPARTS  <- c("temperature test only", "salinity test only", "both tests")
tree <- function(L, s) if (L == "SS-Clines")
    sprintf("/pipeline/_archive/mvp_ssclines/results/MVP%shd_results", s) else sprintf("/pipeline/MVP%shd_results", s)
toppct_keep <- function(p, share) {
    n <- ceiling(share * length(p) / 100 - 1e-9)
    ok <- sort(p[!is.na(p)])
    stopifnot(length(ok) >= n)
    !is.na(p) & p <= ok[n]
}
med <- function(x) if (length(x)) median(x, na.rm = TRUE) else NA_real_

one <- function(i) {
    L <- SEEDS$landscape[i]; s <- SEEDS$seed[i]; D <- sprintf("/pipeline/data/mvp/MVP%s", s)
    # ---- unchanged from mvp_land_set_profile.R
    meta <- fread(file.path(D, sprintf("MVP%s_metadata.tsv", s)), colClasses = c(site = "character", sample = "character"))
    env  <- fread(file.path(D, sprintf("MVP%s_env_present.tsv", s)), colClasses = c(site = "character"))
    vcf  <- fread(cmd = sprintf("grep -v '^##' %s", file.path(D, sprintf("MVP%s.vcf", s))), colClasses = "character")
    ids  <- paste0(vcf[["#CHROM"]], ":", vcf$POS)
    G <- as.matrix(vcf[, -(1:9)]); G <- (substr(G, 1, 1) == "1") + (substr(G, 3, 3) == "1")
    site <- meta$site[match(colnames(vcf)[-(1:9)], meta$sample)]
    stopifnot(!anyNA(site))
    demes <- sort(unique(site))
    F <- sapply(demes, function(d) rowMeans(G[, site == d, drop = FALSE]) / 2)   # SNP x deme
    e <- env[match(demes, site)]
    tr <- fread(file.path(D, "loci_freq.tsv"))
    stopifnot(identical(paste0(tr$chr, ":", tr$pos), ids))
    bg <- which(tr$causal_temp == "neutral" & tr$causal_sal == "neutral")
    bgF <- F[bg, ]; bgF <- bgF[apply(bgF, 1, sd) > 0, ]
    pc1 <- prcomp(t(bgF))$x[, 1]
    ok <- apply(F, 1, sd) > 0
    rt <- rs <- rp <- rep(NA_real_, nrow(F))
    rt[ok] <- abs(cor(t(F[ok, ]), e$bio_1)); rs[ok] <- abs(cor(t(F[ok, ]), e$bio_2)); rp[ok] <- abs(cor(t(F[ok, ]), pc1))
    # ---- new: class of every SNP: causal by trait from the deposit's truth, linked by option B (within d_bp of a
    # causal locus on the same chromosome, inclusive, as near_causal_bp() of mvp_panel_tables.R)
    ct <- tr$causal_temp == "causal"; cs <- tr$causal_sal == "causal"; cz <- as.logical(tr$causal)
    ctr <- fcase(ct & cs, "causal both", ct, "causal temperature", cs, "causal salinity",
                 cz, "causal unassigned")                                      # causal in truth_any, in neither trait list
    stopifnot(identical(!is.na(ctr), cz))
    d_bp <- SEEDS$d_bp[i]
    chr <- as.character(tr$chr)
    near_t <- near_s <- near_o <- rep(FALSE, nrow(tr))
    for (ch in unique(chr[cz])) {
        ic <- which(cz & chr == ch); ia <- which(chr == ch)
        for (q in ic) {
            w <- ia[abs(tr$pos[ia] - tr$pos[q]) <= d_bp]
            if (ctr[q] == "causal temperature") near_t[w] <- TRUE
            else if (ctr[q] == "causal salinity") near_s[w] <- TRUE
            else near_o[w] <- TRUE
        }
    }
    lk <- !cz & (near_t | near_s | near_o)
    cls <- fifelse(cz, ctr, fifelse(lk, "linked", "neutral"))
    lsp <- fifelse(lk, fifelse(near_t & !near_s & !near_o, "linked temperature",
                               fifelse(near_s & !near_t & !near_o, "linked salinity", "linked mixed")), NA_character_)
    nsp <- fifelse(cls == "neutral", fifelse(tr$causal_temp == "neutral" & tr$causal_sal == "neutral",
                                             "neutral, neutral chromosomes", "neutral, causal chromosomes"), NA_character_)
    stopifnot(!anyNA(cls))
    # ---- new: the per-predictor top lists of the two per-predictor methods
    T <- tree(L, s)
    test_of <- lapply(setNames(TESTS, TESTS), function(m) {
        fs <- Sys.glob(file.path(T, "GEA/tables/methods", m, sprintf("%s_pvalues_K*.tsv", m)))
        fs <- fs[!grepl("_sig_", fs)]
        stopifnot(length(fs) == 1L)
        pv <- fread(fs, colClasses = c(SNPID = "character", chr = "character"))
        data.table(SNPID = pv$SNPID, t = toppct_keep(pv$bio_1, 0.25), s = toppct_keep(pv$bio_2, 0.25))[t | s]
    })
    # PCA of a part's SNPs (user, 2026-11-07): PC1 of the deme x SNP allele frequencies (centred, not scaled, as the
    # neutral PC1 above), its |r| across demes with each predictor, and its share of the part's frequency variance
    pca1_of <- function(jj) {
        jj <- jj[ok[jj]]
        if (!length(jj)) return(c(NA_real_, NA_real_, NA_real_))
        X <- t(F[jj, , drop = FALSE])
        if (length(jj) == 1L) return(c(abs(cor(X[, 1], e$bio_1)), abs(cor(X[, 1], e$bio_2)), 1))
        pc <- prcomp(X, center = TRUE, scale. = FALSE)
        c(abs(cor(pc$x[, 1], e$bio_1)), abs(cor(pc$x[, 1], e$bio_2)), pc$sdev[1]^2 / sum(pc$sdev^2))
    }
    part_row <- function(nm, p, jj) {
        q <- pca1_of(jj)
        data.table(set = nm, part = p, n = length(jj), r_temp = med(rt[jj]), r_sal = med(rs[jj]), r_pc1 = med(rp[jj]),
                   pca1_r_temp = q[1], pca1_r_sal = q[2], pca1_pve = q[3])
    }
    out <- rbindlist(lapply(names(SETS), function(nm) {
        x <- fread(file.path(T, "_intermediate/snp_sets", SETS[[nm]], "selected_snps.tsv"), colClasses = c(SNPID = "character"))
        j <- match(x$SNPID, ids)
        stopifnot(!anyNA(j), !anyDuplicated(j))
        R <- list(part_row(nm, "all", j))
        for (k in CLASSES) R[[length(R) + 1L]] <- part_row(nm, k, j[cls[j] == k])
        for (k in LSPLIT)  R[[length(R) + 1L]] <- part_row(nm, k, j[!is.na(lsp[j]) & lsp[j] == k])
        for (k in NSPLIT)  R[[length(R) + 1L]] <- part_row(nm, k, j[!is.na(nsp[j]) & nsp[j] == k])
        if (nm %in% TESTS) {
            tt <- test_of[[nm]]
            if (!setequal(tt$SNPID, x$SNPID)) stop(sprintf("%s %s: test lists %d, set %d, common %d", s, nm,
                                                         nrow(tt), nrow(x), length(intersect(tt$SNPID, x$SNPID))))
            w <- tt[match(x$SNPID, tt$SNPID)]
            R[[length(R) + 1L]] <- part_row(nm, TPARTS[1], j[w$t & !w$s])
            R[[length(R) + 1L]] <- part_row(nm, TPARTS[2], j[w$s & !w$t])
            R[[length(R) + 1L]] <- part_row(nm, TPARTS[3], j[w$t & w$s])
        }
        rbindlist(R)
    }))
    out[, `:=`(landscape = L, seed = s, arch = SEEDS$arch[i], regime = SEEDS$regime[i], demog = SEEDS$demog[i],
               k_best = SEEDS$k_best[i], r2_pc1_temp = cor(pc1, e$bio_1)^2, r2_pc1_sal = cor(pc1, e$bio_2)^2)]
    out[]
}
t0 <- Sys.time()
RES <- mclapply(seq_len(nrow(SEEDS)), function(i) tryCatch(one(i), error = function(e) {
    message("FAIL ", SEEDS$seed[i], ": ", conditionMessage(e)); NULL }), mc.cores = CORES)
R <- rbindlist(RES)
message(sprintf("replicates done: %d of %d in %.1f min", uniqueN(R$seed), nrow(SEEDS),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (uniqueN(R$seed) != nrow(SEEDS)) stop("only ", uniqueN(R$seed), " replicates done")
setcolorder(R, c("landscape", "seed", "arch", "regime", "demog", "k_best", "r2_pc1_temp", "r2_pc1_sal", "set", "part"))

# ------------------------------------------------------------------ asserts against existing tables
# 1. set level = the set profile (same code, earlier run)
SP <- fread(file.path(LMS, "set_profile.tsv"), colClasses = c(seed = "character"))
A <- R[part == "all"]
CNT <- dcast(R[part %in% CLASSES], landscape + seed + set ~ part, value.var = "n")
A <- merge(A, CNT, by = c("landscape", "seed", "set"))
x <- merge(A[set %in% unique(SP$set)], SP, by = c("landscape", "seed", "set"), suffixes = c("", ".sp"))
stopifnot(nrow(x) == nrow(A[set %in% unique(SP$set)]))
dif <- x[, .(n = max(abs(n - n.sp)), r_temp = max(abs(r_temp - r_temp.sp)), r_sal = max(abs(r_sal - r_sal.sp)),
             r_pc1 = max(abs(r_pc1 - r_pc1.sp)), r2t = max(abs(r2_pc1_temp - r2_pc1_temp.sp)),
             ctemp = max(abs((`causal temperature` + `causal both`) / n - sh_ctemp)),
             csal = max(abs((`causal salinity` + `causal both`) / n - sh_csal)))]
print(dif)
stopifnot(all(unlist(dif) < 1e-12))
message("assert 1: set-level n, |r| and causal composition equal set_profile.tsv in ", nrow(x), " rows")
# 2. classes = detection_per_seed.tsv (TP = causal + linked, FP = neutral)
DET <- rbindlist(lapply(seq_len(nrow(LANDS)), function(k)
    fread(file.path(EVAL, LANDS$gal[k], "detection_per_seed.tsv"), colClasses = c(seed = "character"))[
        label %in% PANELS, .(landscape = LANDS$land[k], seed, set = label, n_det = n, n_causal, n_linked, FP)]))
y <- merge(A[set != "truth"], DET, by = c("landscape", "seed", "set"))
stopifnot(nrow(y) == nrow(A[set != "truth"]))
stopifnot(y[, all(n == n_det)], y[, all(`causal temperature` + `causal salinity` + `causal both` + `causal unassigned` == n_causal)],
          y[, all(linked == n_linked)], y[, all(neutral == FP)])
message("assert 2: causal, linked and neutral counts equal detection_per_seed.tsv in ", nrow(y), " rows")
# 3. the test lists: union = the set file is checked per replicate inside one(); here the parts add up
z <- merge(A[set %in% TESTS, .(landscape, seed, set, n)],
           dcast(R[part %in% TPARTS], landscape + seed + set ~ part, value.var = "n"), by = c("landscape", "seed", "set"))
stopifnot(z[, all(`temperature test only` + `salinity test only` + `both tests` == n)])
message("assert 3: the three test parts add up to the set in ", nrow(z), " rows")
u <- merge(A[, .(landscape, seed, set, linked, neutral)],
           dcast(R[part %in% c(LSPLIT, NSPLIT)], landscape + seed + set ~ part, value.var = "n"), by = c("landscape", "seed", "set"))
stopifnot(u[, all(`linked temperature` + `linked salinity` + `linked mixed` == linked)],
          u[, all(`neutral, causal chromosomes` + `neutral, neutral chromosomes` == neutral)])
message("assert 3b: linked by trait and neutral by place add up in ", nrow(u), " rows")

# ------------------------------------------------------------------ outcome, median of the three offset methods
RL <- rbindlist(lapply(seq_len(nrow(LANDS)), function(k)
    fread(file.path(EVAL, LANDS$rank[k], "rank_long.tsv"), colClasses = c(seed = "character"))[
        engine == "median of 3" & seed %in% SEEDS$seed][, landscape := LANDS$land[k]]))
RX <- rank_split(RL, c("seed", "engine"))
W456 <- RX[, .(w456 = sum(w[r >= 4L])), by = .(seed, panel)]
ACC <- RL[empty == FALSE & !is.na(accuracy), .(seed, panel, accuracy)][, loss := max(accuracy) - accuracy, by = seed]
OC <- merge(RL[, .(landscape, seed, panel, rank, empty)], W456, by = c("seed", "panel"))
OC <- merge(OC, ACC[, .(seed, panel, accuracy, loss)], by = c("seed", "panel"), all.x = TRUE)
OC <- OC[panel %in% PANELS, .(landscape, seed, set = panel, rank, w456, accuracy, loss, empty)]
if (NPILOT == 0L) {
    KEY <- c("LFMM2" = "LFMM", "EMMAX" = "EMMAX", "RDA" = "RDA", "1/3 set" = "1/3 methods", "2/3 set" = "2/3 methods",
             "3/3 set" = "3/3 methods")
    TRK <- fread(file.path(LMS, "tables/t_ranks.tsv"))[offset_method == "Median of the three" & statistic == "ranks 4–6" &
                                                       set %in% names(KEY) & landscape != "all"]
    TRK[, set := unname(KEY[set])]
    v <- merge(OC[, .(pct = 100 * mean(w456)), by = .(landscape, set)], TRK[, .(landscape, set, value)], by = c("landscape", "set"))
    stopifnot(nrow(v) == 18L, max(abs(v$pct - v$value)) < 1e-9)
    message("assert 4: % in ranks 4-6 per landscape and set equals t_ranks.tsv (18 rows)")
    stopifnot(all(R[part == "all", uniqueN(seed), by = landscape]$V1 == 600L))
    message("assert 5: 600 replicates per landscape")
}
OC <- merge(OC, unique(R[, .(landscape, seed, arch, regime, demog, k_best, r2_pc1_temp, r2_pc1_sal)]), by = c("landscape", "seed"))

# ------------------------------------------------------------------ write
fwrite(R, file.path(OUT, "set_tracking.tsv"), sep = "\t")
fwrite(OC, file.path(OUT, "set_outcome.tsv"), sep = "\t")
FR <- merge(R, A[, .(landscape, seed, set, n_all = n)], by = c("landscape", "seed", "set"))
SUM <- FR[, .(replicates = .N, replicates_nonempty = sum(n > 0), median_n = median(n),
              median_pct_of_set = median(100 * n / n_all),
              r_temp = median(r_temp, na.rm = TRUE), r_temp_q25 = quantile(r_temp, .25, na.rm = TRUE, names = FALSE),
              r_temp_q75 = quantile(r_temp, .75, na.rm = TRUE, names = FALSE),
              r_sal = median(r_sal, na.rm = TRUE), r_sal_q25 = quantile(r_sal, .25, na.rm = TRUE, names = FALSE),
              r_sal_q75 = quantile(r_sal, .75, na.rm = TRUE, names = FALSE),
              r_pc1 = median(r_pc1, na.rm = TRUE),
              pca1_r_temp = median(pca1_r_temp, na.rm = TRUE), pca1_r_sal = median(pca1_r_sal, na.rm = TRUE),
              pca1_pve = median(pca1_pve, na.rm = TRUE)), by = .(landscape, set, part)]
fwrite(SUM, file.path(OUT, "set_tracking_summary.tsv"), sep = "\t")
print(SUM[part %in% c("all", TPARTS, "linked", LSPLIT, "neutral") & set != "truth",
          .(landscape, set, part, median_pct_of_set, r_temp, r_sal, r_pc1)], digits = 3)
message("written: set_tracking.tsv (", nrow(R), " rows), set_outcome.tsv (", nrow(OC), "), set_tracking_summary.tsv (", nrow(SUM), ")")
