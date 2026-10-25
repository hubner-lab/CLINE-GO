# =============================================================================
# fig3_style.R -- colour roles of the Trifolium figures (Fig 3, S9-S13, offset pack), all from ltc "minou".
# Plan ~/.claude/plans/scalable-doodling-allen.md (2026-10-16; user: the whole paper in minou, gradient maps unchanged).
# source() after benchmarks/ms_style.R is found: this file sources it, then adds only Trifolium roles.
# Every colour here is a MINOU value, white / MS_INK, or a mix of two of them (asserted at the bottom).
# =============================================================================
local({
    here <- if (dir.exists("/pipeline")) "/pipeline" else getwd()
    source(file.path(here, "benchmarks/ms_style.R"), local = FALSE)
})

mix <- function(a, b, t) grDevices::colorRampPalette(c(a, b))(101)[round(t * 100) + 1]

# sNMF K = 5 ancestry clusters: the five minou hues; grey stays for references.
FIG3_CLUSTERS <- unname(MINOU[c("teal", "red", "amber", "sage", "navy")])

# EMMAX predictors, all 16 (user 2026-10-17: draw bio_11 and bio_19 as their own series, not "bio_8 = bio_11"; they
# are identical to bio_8 / bio_16 at every plant, which the captions state). Was 14 colours before 2026-10-17.
# Temperature: red -> orange -> amber; precipitation: teal -> navy. Class is the meaning, shade only separates.
FIG3_TEMP <- paste0("bio_", 1:11)
FIG3_PREC <- paste0("bio_", c(12, 13, 15, 16, 19))
FIG3_PRED_COL <- c(
    setNames(grDevices::colorRampPalette(c(MINOU[["red"]], MS_ORANGE, MINOU[["amber"]]))(length(FIG3_TEMP)), FIG3_TEMP),
    setNames(grDevices::colorRampPalette(c(MINOU[["teal"]], MINOU[["navy"]]))(length(FIG3_PREC)), FIG3_PREC))
FIG3_CLASS_COL <- c(temperature = MINOU[["red"]], precipitation = MINOU[["teal"]])

# Manhattan: chromosomes alternate two light tints; thresholds per method; agreement-set status.
FIG3_CHR      <- c(mix("#FFFFFF", MINOU[["grey"]], 0.55), mix("#FFFFFF", MINOU[["navy"]], 0.45))
# RDA is sage everywhere (points, threshold, "RDA only"): teal and navy already carry the precipitation predictors.
FIG3_THR      <- c(EMMAX = MINOU[["red"]], RDA = MINOU[["sage"]])
FIG3_STATUS   <- c("both methods" = MINOU[["navy"]], "EMMAX only" = MINOU[["red"]], "RDA only" = MINOU[["sage"]])
# Miami v2 classes (fancy-scribbling-garden.md Phase 2): the agreement set split into the identical SNP and the linked
# partner, plus significant SNPs with no partner in the other method (smaller, orange: not amber, the region shading).
FIG3_AGREE    <- c("same SNP, both methods" = MINOU[["navy"]], "linked, within 130,811 bp" = MINOU[["teal"]],
                   "one method, no partner" = MS_ORANGE)
# Miami v3 (idempotent-conjuring-noodle.md Phase 2; user review v2: "two groups, one-method sig and agreement SNPs,
# two colours, two sizes, significant in red"): agreement = selected_snps.tsv (identical + linked), red and larger.
FIG3_AGREE2   <- c("both methods (within 130,811 bp)" = MINOU[["red"]], "one method" = MS_ORANGE)
FIG3_STRIPE   <-c("in a repeated vector" = MINOU[["red"]], "locus-specific" = MINOU[["teal"]])
FIG3_K        <- c("K = 3" = MINOU[["amber"]], "K = 5" = MINOU[["navy"]])
FIG3_REGION   <- mix("#FFFFFF", MINOU[["amber"]], 0.35)   # region shading behind the points
FIG3_GENE     <- MINOU[["navy"]]                          # one gene colour (deep-dive design)

# Offset consensus categories (Fig 1 and S1-S4; gradient maps keep their own scale).
FIG3_RISK <- c("higher risk" = MINOU[["red"]], "methods disagree" = MINOU[["amber"]], "lower risk" = MINOU[["teal"]])

# ---- assertion: every role is minou, white, ink, or a two-colour mix of them --------------------------------------
local({
    base <- c(unname(MINOU), MS_ORANGE, "#FFFFFF", MS_INK)
    ok <- unique(toupper(unlist(lapply(base, function(a) lapply(base, function(b)
        grDevices::colorRampPalette(c(a, b))(101))))))
    ok <- c(ok, toupper(grDevices::colorRampPalette(c(MINOU[["red"]], MS_ORANGE, MINOU[["amber"]]))(length(FIG3_TEMP))),
            toupper(grDevices::colorRampPalette(c(MINOU[["teal"]], MINOU[["navy"]]))(length(FIG3_PREC))))
    roles <- c(FIG3_CLUSTERS, FIG3_PRED_COL, FIG3_CLASS_COL, FIG3_CHR, FIG3_THR, FIG3_STATUS, FIG3_STRIPE, FIG3_K,
               FIG3_REGION, FIG3_GENE, FIG3_AGREE, FIG3_AGREE2,
               FIG3_RISK, MS_SEQUENTIAL, MS_DIVERGING)
    bad <- roles[!toupper(roles) %in% ok]
    if (length(bad)) stop("fig3_style.R: colours outside minou: ", paste(bad, collapse = ", "))
})
