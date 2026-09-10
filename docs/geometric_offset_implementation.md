# Geometric Genetic Offset — Implementation Guide

> **Status:** research + design spec. **No code written yet.** This document is the brief for a next-session
> code agent implementing `genetic_offset` as a second maladaptation method beside `gradient_forest`.
> Confidence tags: `[verified]` primary source · `[verify-container]` confirm against pinned **LEA 3.22**
> (`?genetic.gap` inside Docker) before coding · `[codebase]` grounded in this repo · `[literature]` author
> recommendation · `[design]` design judgment.

## 1. TL;DR

Geometric genetic offset (the "genetic gap") is a single LEA function call. `genetic.gap()` fits an `lfmm2`
model internally, derives the covariance matrix of SNP effect sizes, and returns a **squared (quadratic)
distance in environmental space** between present and future climate — one value per supplied location.
`[verified]`

**The whole data contract already exists in the pipeline and is row-aligned.** `[codebase]` The only genuinely
new code is one R script (`scripts/geometric_offset.R`). The integration cost is in plumbing: the maladaptation
"method registry" advertises a factory pattern but is **single-method-hardcoded** — rules, summary, and Shiny
all assume `gradient_forest`. Adding a method touches ~6 areas (Section 7). `[codebase]`

**Key discovery that simplifies the design:** `genetic.gap()` returns offset *per row of `new.env`/`pred.env`*,
not tied to genotype rows. Feed per-sample climate → per-sample offset (piemap). Feed per-raster-cell climate →
per-cell offset (landscape raster). Both input tables already exist (`_site` and `_all` climate tables that GF
already consumes). So Geometric GO maps **almost 1:1 onto GF's output set** — same `offset_raster.tif`,
`genetic_offset_map.tsv`, `genetic_offset_site.tsv`, and piemap. `[verified]` + `[codebase]`

## 2. Method summary `[literature]`

Gain et al. 2023 (*MBE* 40(6):msad140, DOI 10.1093/molbev/msad140) give a quantitative theory unifying genomic
offset statistics. The **geometric genomic offset** is implemented as `genetic.gap()` in LEA. Theory: under
**local Gaussian stabilizing selection** with traits at **local equilibrium**, the offset is a quadratic
distance in env space weighted by the covariance matrix of LFMM effect sizes — it has a dual interpretation as
a quadratic distance in environmental *and* genetic space, and predicts **relative fitness loss** after
environmental change. The function also returns a `distance` field = a confounder-corrected, multi-predictor
**modified RONA** statistic.

## 3. `genetic.gap()` API reference `[verified]` (LEA man page `genetic_gap.Rd`, repo `bcm-uga/LEA` branch `devel`)

> Tag every signature/field below `[verify-container]` before coding: run `?genetic.gap` and
> `str(genetic.gap(...))` on the offset_example data inside the `cline-go:latest` container (LEA 3.22).

### Signature
```r
genetic.gap(input, env, new.env, pred.env, K, scale, candidate.loci)
```

### Arguments
| arg | meaning |
|---|---|
| `input` | Genotype matrix or path, **lfmm format, NO missing values (no 9/NA)**. Use `impute()` (NMF) to complete. |
| `env` | Env covariates at the **genotype sampling locations** (the fitting env). `env` format, numeric, no missing. |
| `new.env` | Present env at the **prediction locations** (may differ from `env`). **Defaults to `env`.** |
| `pred.env` | Predicted (future) env at the **same locations as `new.env`**, same dimensions. Typically WorldClim-derived. |
| `K` | Number of latent factors (integer or sequence; a sequence returns the averaged prediction). From the PCA-screeplot elbow of the genotype matrix. |
| `scale` | Logical. If `TRUE`, **all env matrices are centered/scaled using the columnwise mean & sd of the `env` matrix.** Use it to read variable importance off `eigenvalues`/`vectors`. **Does NOT change the offset values.** |
| `candidate.loci` | Vector of column labels/indices of loci to include. Default = all loci. |

### Return object (a list)
| field | meaning |
|---|---|
| `offset` | Vector of genomic offset (genetic gap) values, **one per location in `new.env`/`pred.env`**. The maladaptation score. |
| `distance` | Vector of modified-RONA env distances (confounder-corrected, multi-predictor), one per location. |
| `eigenvalues` | Eigenvalues of the effect-size covariance matrix — relative importance of env-variable combinations (use with `scale=TRUE`). |
| `vectors` | Eigenvectors (env-variable combinations) sorted by importance. The **importance-plot analog** to GF. |

### Canonical example (from the man page) `[verified]`
```r
data("offset_example")
Y <- offset_example$geno; X <- offset_example$env; X.pred <- offset_example$env.pred
plot(prcomp(Y))                       # PCA screeplot suggests K
g.gap <- genetic.gap(input = Y, env = X, pred.env = X.pred, K = 2)   # new.env defaults to env
round(g.gap$offset, 3)                # per-sample offset
g.gap.scaled <- genetic.gap(input = Y, env = X, pred.env = X.pred, scale = TRUE, K = 2)
# "Scaling does not change genetic gaps":  plot(g.gap$offset, g.gap.scaled$offset)  -> identity line
barplot(g.gap.scaled$eigenvalues)     # variable-importance axes
g.gap.scaled$vectors[,1:2]            # variable loadings
```

## 4. The two design issues — RESOLVED

### 4.1 Per-sample vs landscape raster — BOTH are native `[verified]`
`offset` is computed *per row of `new.env`/`pred.env`*, decoupled from genotype rows. So:
- **Per-sample (piemap):** `new.env = climate_site` (present per-sample), `pred.env = climate_future_site`
  (future per-sample) → one offset per sample → piemap. `[codebase]` both tables exist.
- **Per-cell (landscape raster):** `new.env = climate_all` (present per raster cell), `pred.env =
  climate_future_all` (future per raster cell) → one offset per cell → reshape to a `terra` raster. `[codebase]`
  both `_all` tables exist and are already consumed by `gradient_forest_offset.R`.

> Earlier assumption (Explore agents + first advisor pass) was that the raster needs custom projection and
> should be deferred. **Primary source overrides this:** the raster is a native `genetic.gap()` call on the
> per-cell climate tables. The François *Spatial_Prediction_Genomic_Offset* tutorial demonstrates the identical
> computation done by hand via `lfmm2(..., effect.sizes=TRUE)` → `B <- mod@B` →
> `mean(((env.new - env.pred)[cell,] %*% t(B[candidates,]))^2)` per cell → `terra::rast`. `[verified]` Either
> route works; **prefer the single `genetic.gap()` call** feeding the `_all` tables (less custom code,
> identical math). `[design]`

**Consequence:** Geometric GO can fill GF's full output set — `offset_raster.tif`, `genetic_offset_map.tsv`,
`genetic_offset_site.tsv`, and the piemap — so it reuses the same `mala_*` path helpers and the same
`plot_piemap.R` and raster-write logic. `[codebase]`

### 4.2 Scaling — pass RAW, don't pre-scale, don't double-scale `[verified]`
`scale=TRUE` standardizes all env matrices by the **`env` (training) mean/sd**, and the man-page example proves
**scaling does not change the offset** — it only affects `eigenvalues`/`vectors`. The correctness requirement is
that `new.env` and `pred.env` are on the **same scale as each other** and consistent with `env`.

**Recommendation `[design]`:** pass the **RAW** present and future climate tables (`O['climate_site']` /
`O['climate_future_site']`, and `O['climate_all']` / `O['climate_future_all']`) — they are already on the same
raw scale (future is RAW, same units as present). Set `scale=TRUE` to get interpretable `eigenvalues`/`vectors`.
**Do NOT** pass the pre-scaled `O['climate_site_scaled']` (it scales present only → mismatched with raw future),
and **do NOT** hand-scale future separately. The tutorial's manual route scales `env.new`/`env.pred` with the
*training env's* mean/sd precisely to preserve this consistency — `genetic.gap(scale=TRUE)` does it internally.
`[verify-container]` confirm via `?genetic.gap` that `scale` uses `env`'s parameters for all three matrices.

## 5. Author recommendations & interpretation caveats `[literature]` (Gain et al. 2023)

- **K:** estimate from the PCA-screeplot elbow of the genotype matrix. Authors note latent-factor correction
  **increases predictive performance**. `[design]` Use the pipeline's `K_BEST` (same K as the LFMM2 GEA) for
  consistency and comparability; optionally pass a short sequence (genetic.gap averages). `[codebase]`
- **Candidate vs whole-genome:** authors **recommend candidate loci** from the GEA — *"Loci identified in the
  GEA study increased the performance of GO statistics compared with using whole genomic data, and the
  improvements were substantial."* They also note whole-genome predictions stay **close** to candidate-based
  ones. → Default to the curated SNP set (`candidate.loci` = GEA-significant indices); whole-genome ("all") is a
  valid fallback. This dovetails with the pipeline's existing `snp_sets` fan-out. `[codebase]` (Note tension
  with Lind & Lotterhos 2025 cross-method benchmark in `maladaptation_methods.md`, which found <3% median
  advantage for adaptive sets — record both; the authors' own guidance is candidate-first.)
- **Unitless squared distance:** the offset is a quadratic distance with **no absolute interpretation**; not
  magnitude-comparable across datasets or methods. **Compare by rank** (Spearman) against GF, never by raw
  value. `[literature]` + `[design]`
- **Equilibrium / Gaussian stabilizing selection:** assumes allele frequencies at local equilibrium under
  selection for an intermediate optimum. Traits out of equilibrium add an intercept term. `[literature]`
- **Novel / non-analog climates:** linear models make the offset **invariant under translation in the niche** —
  predictions are most relevant near the center of the species range, **less reliable at range margins / novel
  climates**. `[literature]`
- **Imputation:** `input` must be complete (no 9/NA). `[verified]`

## 6. Data contract — all inputs exist, row-aligned `[codebase]`

| Geometric-GO input | Pipeline path | Notes |
|---|---|---|
| `input` (genotypes, imputed) | `W['lfmm_imp_full']` = `{WORK_FILT}{VCF_BASE}_K{K_BEST}imp.lfmm` | space-delim, no header, rows=samples, cols=SNPs, 0/1/2, **no missing** ✓ |
| sample order | `W['samples_order']` = `{INTER}samples/samples_order.list` | VCF column order |
| `env` (fitting) = present per-sample | `O['climate_site']` | cols `sample`+`bio_1..bio_19`, **RAW** |
| `new.env`/`pred.env` per-sample | `O['climate_site']` / `O['climate_future_site']` | RAW present / RAW future, same row order |
| `new.env`/`pred.env` per-cell (raster) | `O['climate_all']` / `O['climate_future_all']` | RAW present / RAW future per raster cell (GF already uses these) |
| present raster (geometry to write back) | GF's `present_raster` input | for `terra::rast` cell layout |
| `K` | `K_BEST` (config `sNMF.k_best`) | same K as LFMM2 GEA |
| predictors | `PREDICTORS_SELECTED` | use the **same set as GF** for comparability |
| `candidate.loci` | `snp_set_file(set_name)` → `selected_snps.tsv` (cols `SNPID,chr,pos,min_pvalue`) | map SNPIDs → genotype column indices |

**Row-alignment invariant (no joins):** LFMM row i ↔ metadata row i ↔ `climate_site` row i ↔
`climate_future_site` row i — all forced to VCF column order. The `_all` present/future tables are cell-aligned
to each other in the same order GF consumes. `[codebase]`

> **`candidate.loci` — CRITICAL, silent-wrong-output trap** `[verified]` + `[design]`: `genetic.gap` fits
> `lfmm2` on the **whole** genotype matrix regardless of `candidate.loci`; the candidate set restricts **only
> the offset aggregation** (the man page: loci "included in the computation of the genetic offset"; the tutorial
> confirms `B <- mod@B` is computed for all loci, then the offset uses `B[candidates,]`). Therefore **always
> pass the FULL imputed matrix as `input`** and restrict via the `candidate.loci` index vector. **Do NOT subset
> the genotype matrix before fitting.**
>
> This is the opposite of GF: `gradient_forest_model.R` *subsets the matrix* (`snp_subset <- lfmm_dt[,
> ..mask_adaptive]`, lines ~31–55) because GF refits per SNP set. For `genetic.gap`, reuse GF's SNP-ID→column
> logic **only to build the integer index vector** (`which(vcfsnp %in% selected_snps$SNPID)`), then hand that
> vector to `candidate.loci` — never to subset `input`. Passing a subsetted matrix would refit the latent
> factors on candidates only → wrong `B` → **plausible-but-wrong offsets that fail silently**.
>
> `genetic.gap` wants **column positions in the genotype matrix** in the same SNP order as `vcfsnp` (GF arg 3).
> The `"all"` SNP-set sentinel → omit `candidate.loci` (default = all loci). `[codebase]`

## 7. Integration map — 6 touch-points `[codebase]`

The registry's `engine`/`*_script` fields are **vestigial** (only `.keys()` and `supports_random_model` are
read). Rules hardcode script paths. Mirror GF's hardcode pattern; do **not** build a real factory (out of
scope). `[design]`

1. **Registry** — `workflow/methods/maladaptation.py`, add a `geometric_offset` entry. Set
   `supports_random_model=False` (lfmm2 is deterministic → no neutral/random model). Add a flag such as
   `supports_cumulative_importance=False` (no turnover concept). `supports_spatial=False` (no PCNM — env-only).
2. **`workflow/rules/maladaptation.smk`** — `wildcard_constraints: method` is pinned to `r"gradient_forest"`
   (~line 12). Widen to `r"gradient_forest|geometric_offset"`. Add **one** rule that runs
   `scripts/geometric_offset.R` producing offset_raster + map + site TSV in a single `genetic.gap()`-based call
   (no separate model/offset split — GF needs that only because random forest is fit once then projected).
3. **`workflow/rules/common.smk` `get_targets()` maladaptation branch** (~1309–1350) — loops
   `set_name × spatial_tag × method`. Currently requests `mala_cumimp` + `mala_importance` **unconditionally**.
   Gate per-method: skip `mala_cumimp` for `geometric_offset`; for `mala_importance` emit an
   **eigenvalues/vectors** plot instead (variable importance). Force `spatial_tag = nospatial` only for
   `geometric_offset` (no PCNM). `[verify-container]` confirm genetic.gap has no spatial-correction arg (it does
   not — env-only).
4. **`mala_*` path helpers** (`common.smk` ~744–824) — already method-parameterized (`method` is first arg).
   **No change.** Reuse `mala_offset_raster/_map_values/_site_values(method,…)` and
   `mala_offset_piemap(method,…,size_trait)`.
5. **`workflow/rules/summary.smk` (~105–146) + `scripts/write_summary.R` (~285–324)** — both **hardcode
   `'gradient_forest'`** and do not loop methods. `write_summary.R` keys rows by
   `tag <- basename(dirname(p))` = `{run_label}_{spatial_tag}`, which **drops `{method}`** → two methods would
   **collide** on metric names (e.g. `offset_mean_{run_label}_{spatial_tag}`). Fix: loop `ACTIVE_MALA_METHODS`
   and include method in the row key (e.g. `basename(dirname(dirname(p)))`).
6. **Shiny** — `R/fct_paths.R` is already `method=`-parameterized (no change). `R/mod_maladaptation.R` has **no
   method selector** and hardcodes `method="gradient_forest"`; add a selector and thread `method` through
   `find_gf_suffixes()` (and gate the cumulative-importance card off for `geometric_offset`).

## 8. New R script spec — `scripts/geometric_offset.R` `[design]`

One script, one `genetic.gap()` call set. Mirror GF script conventions: positional `commandArgs(trailingOnly=
TRUE)`, `message()` logging, `data.table::fread(..., colClasses=c(sample="character"))` for ID columns
(per CLAUDE.md fread rule), libraries loaded at top, `tryCatch` around optional plots.

Proposed positional args (finalize against the rule):
```
1  LFMM_IMP_FULL      W['lfmm_imp_full']           # imputed genotypes
2  VCFSNP             SNP order file (GF arg 3)     # to map candidate SNPIDs -> column indices
3  CANDIDATE_SNPS     snp_set_file(set_name)        # selected_snps.tsv; "all" sentinel => all loci
4  ENV_SITE_PRESENT   O['climate_site']  (RAW)      # env (fitting) + new.env (per-sample present)
5  ENV_SITE_FUTURE    O['climate_future_site'] RAW  # pred.env (per-sample future)
6  ENV_ALL_PRESENT    O['climate_all']   (RAW)      # new.env (per-cell present) for raster
7  ENV_ALL_FUTURE     O['climate_future_all'] RAW   # pred.env (per-cell future) for raster
8  PRESENT_RASTER     GF present_raster             # cell geometry for terra write-back
9  PREDICTORS         PREDICTORS_SELECTED           # subset/order bio_* columns
10 K                  K_BEST
11 OUT_SITE_TSV       mala_offset_site_values(...)  # per-sample offset (+ RONA distance)
12 OUT_MAP_TSV        mala_offset_map_values(...)   # per-cell offset
13 OUT_RASTER_TIF     mala_offset_raster(...)       # terra raster of per-cell offset
14 OUT_IMPORTANCE     mala_importance(...)          # eigenvalues/vectors barplot
```

Behavior — **ONE `genetic.gap()` call, not three** `[design]` + `[verified]`:

`offset` is per-row of `new.env`/`pred.env`, each row independent, and **scale-invariant** (man-page example:
`plot(g.gap$offset, g.gap.scaled$offset)` is the identity). Stack the per-sample rows on top of the per-cell
rows and run a single fit with `scale=TRUE` — yields valid offsets for both AND `eigenvalues`/`vectors` for
free. One `lfmm2` refit on the full barley matrix instead of three.

```r
# cand = integer column indices of CANDIDATE_SNPS against VCFSNP order; omit (NULL) if "all"
# present_site/future_site: per-sample RAW env;  present_all/future_all: per-cell RAW env
# all four already subset to PREDICTORS in one fixed column order
new.env  <- rbind(present_site, present_all)     # n_site rows, then n_cell rows
pred.env <- rbind(future_site,  future_all)
g <- genetic.gap(input = LFMM, env = present_site,            # env = per-sample FITTING env only
                 new.env = new.env, pred.env = pred.env,
                 K = K_BEST, scale = TRUE, candidate.loci = cand)
n_site      <- nrow(present_site)
site_offset <- g$offset[seq_len(n_site)]                      # -> OUT_SITE_TSV (sample, offset, rona=g$distance[seq_len(n_site)])
cell_offset <- g$offset[(n_site + 1):length(g$offset)]        # -> OUT_MAP_TSV + reshape onto PRESENT_RASTER -> OUT_RASTER_TIF
# g$eigenvalues / g$vectors -> barplot OUT_IMPORTANCE (no separate scale=TRUE call needed)
```

- Build `cand` = column indices of `CANDIDATE_SNPS` against `VCFSNP` order (or omit for "all"). **Pass the FULL
  `input` matrix — never subset it** (see §6 trap).
- `env` stays the per-sample fitting env; only `new.env`/`pred.env` are stacked. Slice `g$offset` by row range
  to recover site vs cell — no second fit, no scale inconsistency. `scale=TRUE` everywhere (does not move the
  offset; only makes `eigenvalues`/`vectors` interpretable).
- Piemap produced by the existing `plot_gf_offset_piemap` rule reusing `plot_piemap.R` on `OUT_SITE_TSV`
  (no new piemap script). `[codebase]`

`[verify-container]` Before wiring the args, confirm in Docker: (a) `genetic.gap` accepts stacked
`new.env`/`pred.env` row counts ≠ genotype row count (per-cell grids — the tutorial's 36k-cell example shows it
does); (b) it tolerates the large stacked `new.env` (n_cells × n_predictors) without excessive memory — if not,
fall back to the tutorial's manual `lfmm2(effect.sizes=TRUE)` + `B` per-cell route, which streams cell-by-cell;
(c) whether `eigenvalues`/`vectors` reflect the **candidate** subset or all loci when `candidate.loci` is set —
this determines what the importance plot means.

## 9. Config block `[design]`

```yaml
Maladaptation:
  snp_sets: "all"            # existing SNP-set fan-out also feeds candidate.loci
  methods:
    gradient_forest: { ... }            # unchanged
    geometric_offset:
      candidate_only: TRUE              # TRUE: candidate.loci = GEA SNP set; FALSE: all loci
      k: ""                             # blank => use K_BEST
      # NO ntree / cor_threshold / spatial_correction / random_model — not applicable
```

## 10. Verification (next session, not now) `[design]`

1. Container sanity first: `Rscript -e 'library(LEA); data(offset_example); str(genetic.gap(offset_example$geno, offset_example$env, pred.env=offset_example$env.pred, K=2))'` — confirm fields and per-row `offset` on LEA 3.22.
2. Run `mode=maladaptation` on **TEST** (`config.yaml`) with `geometric_offset` configured. Confirm
   `Maladaptation/{plots,tables}/geometric_offset/{run_label}_nospatial/` gets `genetic_offset_site.tsv`,
   `genetic_offset_map.tsv`, `offset_raster.tif`, the offset piemap, and the eigenvalues importance plot.
3. Confirm `Pipeline_summary.tsv` reports **both** methods with no metric-name collision.
4. Spearman rank-correlate per-site geometric offset vs GF offset (expect positive, not identical).
5. Shiny: method selector switches GF ↔ geometric_offset piemaps; cumulative-importance card hidden for
   geometric_offset.

## 11. RDA Offset — implemented (2026-07-29 correction) `[design]` + `[literature]`

**Correction, dated 2026-07-29:** this section originally said RDA offset was contingent on RDA-GEA
landing first and would be Shiny-driven, not a batch method. Both premises changed: RDA-GEA landed
(`scripts/rda.R`), and RDA offset shipped as a **standard batch maladaptation method**
(`workflow/methods/maladaptation.py` → `rda_offset`, `scripts/rda_offset.R`) — the human-judgment
steps (which axes, structure correction) are resolved once via config (`Maladaptation.methods.rda_offset.axes`
/ `.condition_pcs`) plus the existing curated-SNP-set mechanism shared with `geometric_offset`, not through
a bespoke interactive UI. RDA-based genomic offset (Capblancq & Forester 2021, *Methods Ecol. Evol.*
12(12):2298–2309, DOI 10.1111/2041-210X.13722) projects present and future climate into the RDA-constrained
adaptive space (via `method="loadings"`, never `predict()` — broken in the published source) and measures
eigenvalue-weighted Euclidean distance between them per site/cell. Full decision table, gotchas, and
citations: `docs/rda_research.md` Part B; condensed summary: `docs/maladaptation_methods.md`'s RDA Genetic
Offset section. Do not re-research — Part B is the primary spec.

## 12. Sources
- LEA man page `genetic_gap.Rd` — `github.com/bcm-uga/LEA` (branch `devel`, `man/genetic_gap.Rd`); pinned **LEA
  3.22** in `Dockerfile`. `[verified]`
- LEA man page `lfmm2.Rd` (same repo) — `lfmm2(input, env, K, lambda=1e-5, effect.sizes=FALSE)` → `@U,@V`, and
  `@B` with `effect.sizes=TRUE`. `[verified]`
- François tutorial *Visualizing geographic predictions of genetic offset measures*
  (`membres-timc.imag.fr/Olivier.Francois/LEA/files/Spatial_Prediction_Genomic_Offset.html`) — manual raster
  route via `lfmm2(effect.sizes=TRUE)` + `B`. `[verified]`
- Gain, Decap, Lindo, François et al. (2023). *A quantitative theory for genomic offset statistics.* **MBE**
  40(6):msad140, DOI 10.1093/molbev/msad140 (PMC10306404). `[literature]`
- Gain & François (2021). *LEA 3.* Mol. Ecol. Resour. 21(8):2738. `[literature]`
- Capblancq & Forester (2021). RDA genomic offset. DOI 10.1111/2041-210X.13722. `[literature]`
