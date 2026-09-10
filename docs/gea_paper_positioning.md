# Positioning CLINE-GO against Lotterhos 2023 — what to claim, and what not to

Derived from the journal-07 benchmark (`work/journal/07_mvp_gea_crossseed_portfolio.Rmd`),
14 MVP replicates × 2 MAF arms × 5 structure-correction rungs × 11 methods.
Source study: **Lotterhos 2023**, PNAS 120(12):e2220313120. Data: BCO-DMO
doi:10.26008/1912/bco-dmo.889769.1. Scoring contract: `benchmarks/MVP_README.md`.

Every number below is reproducible from `benchmarks/mvp_eval/` — the per-axis tables live in
`axis07/` (MAF 0.01) and `axis07_m05/` (MAF 0.05), the cross-seed roll-up in `report07/`.

---

## 1. The one-paragraph argument

Lotterhos 2023 evaluated each GEA method **alone, at one threshold**. Every method they tested
lands at one of two unusable extremes: conservative methods recover ~8 % of causal loci, and
liberal ones reach useful power only at 89–99 % false discovery. Nothing sits in between.
CLINE-GO's contribution is the step they did not test — **combining methods with per-method
thresholds** — which recovers 2.5× the causal loci at the same (zero, median) false-discovery
rate, using a configuration fixed on an independent replicate and applied without any tuning.

---

## 2. Claim ladder — strongest first

| # | Claim | Evidence | Verdict |
|---|---|---|---|
| 1 | Our implementation reproduces the source study | LFMM @ q<0.05, per replicate: median \|Δ TPR\| = **0.000** on both axes, 14/20 cells exact to 3 d.p. | **strong** |
| 2 | No published single method has a usable operating point | their own published TPR/FDR, medians over our 12 replicates (§4) | **strong** — uses *their* numbers |
| 3 | Naive all-method union — the common default — inherits the worst of it | TPR 0.181 at FDR **0.969** | **strong** |
| 4 | A structured combination dominates, with no truth table | TPR 2.5× / 1.5× at median FDR 0.000, higher on **19/20** cells | **strong, headline** |
| 5 | The winning structure transfers between datasets | selection bias median **0.011–0.017** F1 | **strong, and novel** |
| 6 | Optimal structure depends on trait architecture | 2-method union for moderately polygenic, `at_least_2`@5 kb for highly polygenic; 4/4 per stratum | moderate (4 seeds/stratum) |
| 7 | MAF 0.05 improves ranking | 70 % of method×seed pairs, median ×1.11 | moderate — direction replicates, magnitude does not |

**Claim 5 is the methodological novelty.** Anyone can report "our combination scored best on our
data". Showing that a configuration fixed on one dataset still works on datasets that took no
part in choosing it is what makes it a *recommendation* rather than a fitted result.

---

## 3. Headline table — the achievable configuration

**LFMM + EMMAX, plain union, `Filter.maf: 0.05`, per-method thresholds imported verbatim from an
independent replicate. No truth table used at any point.**

| Axis | Our TPR | Published TPR | Ratio | Our FDR (median) | Published FDR | TPR ≥ pub. | Strict dominance |
|---|---:|---:|---:|---:|---:|---:|---:|
| Temp (confounded, R²(PC1~env) 0.31–0.57) | **0.211** | 0.084 | **2.5×** | 0.000 | 0.000 | **12/12** | 6/12 |
| Env2 (orthogonal control) | **0.522** | 0.350 | **1.5×** | 0.000 | 0.000 | **8/8** | 5/8 |

**Write it as:** *uniformly higher power at a median false-discovery rate of zero, matching the
published method's; strict per-replicate dominance in roughly half of cases.*

Do **not** write "strictly dominant" unqualified. TPR is higher on 19/20 replicate×axis cells and
never lower, but our FDR is non-zero on 6 cells (worst 0.667, on a 324-causal replicate).

---

## 4. Their paper has no usable operating point

Published values, medians over the same 12 replicates (`benchmarks/mvp_published_baseline.tsv`):

| Method (published) | TPR | FDR | Precision |
|---|---:|---:|---:|
| LFMM, temp | 0.084 | 0.000 | 1.000 |
| LFMM, Env2 | 0.350 | 0.000 | 1.000 |
| Kendall AF~env | 0.373 / 0.300 | — | — |
| RDA, uncorrected | 0.197 | 0.892 | 0.108 |
| RDA, structure-corrected | 0.140 | 0.922 | 0.078 |
| Phenotype GWAS, temp | 0.678 | **0.989** | 0.011 |
| Phenotype GWAS, Env2 | 0.963 | **0.993** | 0.007 |

⚠ **`pub_lfmm_aucpr_temp_neut` is not citable.** It is not reproducible from the deposit's own
temperature p-values (we re-score 0.203 against a published 0.0783, and the published figure
matches the *salinity* column within the offset every other cell shows — consistent with a column
swap). Flagged as `lfmm_temp_aucpr_suspect` in the baseline table. Use the re-scored value from
`benchmarks/mvp_score_published.R`, or omit.

---

## 5. Table hierarchy — what each row is for

| Configuration | Truth table? | Temp TPR | Temp FDR | Role in paper |
|---|---|---:|---:|---|
| Published LFMM @ q<0.05 | — | 0.084 | 0.000 | baseline |
| Ours, identical config | no | 0.060 | 0.000 | **validation** (median \|Δ\| = 0.000) |
| Naive all-method union @ bonf 0.05 | no | 0.181 | **0.969** | **strawman** = pipeline's current default |
| **LFMM+EMMAX union, MAF 0.05, imported thresholds** | **no** | **0.211** | **0.000** | **headline** |
| LFMM+EMMAX+GLM `at_least_2`@5 kb, oracle-calibrated | **yes** | 0.174 | 0.000 | **headroom only — label as such** |

**The last row must never be presented as achievable.** Its per-method thresholds are chosen by
best-F1 against the truth table. A real user has none. The gap between it and the row above is the
value of a truth-free calibrator (`mode=pregea`), which journal 07 did not measure.

---

## 6. Where the gain comes from — a two-factor result

| Step | Temp TPR | Temp FDR |
|---|---:|---:|
| LFMM alone, imported threshold, MAF 0.01 | 0.060 | 0.000 |
| + combine with EMMAX (MAF 0.01) | 0.167 | 0.500 |
| + `Filter.maf: 0.05` | **0.211** | **0.000** |

Combination raises power but degrades precision at MAF 0.01; the MAF 0.05 filter is what makes the
*imported* calibration land correctly. Neither ingredient alone reaches the headline. This also
explains why journal 06 saw the MAF effect so strongly on its single replicate.

---

## 7. Secondary result — transferability (the novel part)

Fixing the subset/rule/window on one replicate and re-using it elsewhere costs a median of
**0.011–0.017 F1** versus searching each dataset from scratch. The gap between the best scheme and
the pipeline default is **0.22** — an order of magnitude larger.

Practical claim: **users should not search method subsets on their own data.** Take the published
structure, tune only the thresholds. Tail risk exists (max observed bias 0.24), so present the
median as a recommendation, not a guarantee.

---

## 8. Architecture-dependent recommendation

| Trait architecture | Recommended structure | Median gain over default | Wins |
|---|---|---:|---:|
| Moderately polygenic (28–70 causal) | `LFMM+RDA` or `LFMM+EMMAX`, plain union, no window | +0.272 | 4/4 |
| Highly polygenic (268–395 causal) | `LFMM+RDA+EMMAX`, `at_least_2`, 5 kb window | +0.182 | 4/4 |

LFMM appears in every recommendation; none contains more than three methods. **Both all-11 schemes
rank last on both arms** — more methods is actively worse, which contradicts the pipeline's own
current behaviour.

---

## 9. Actionable pipeline finding to state in the discussion

`SIGSNPS_METHOD = 'All'` is hardcoded at `workflow/rules/common.smk:192`, so batch `mode=gea` can
only union. That union is **the worst of the seven schemes tested**: median 1 242 SNPs called, and
it loses to a tuned single method on 8/8 replicates. Replacing it with `LFMM+RDA+EMMAX` at
`at_least_2`/5 kb gives **+0.13 F1 at 4.5× fewer calls**. `combine_sigsnps.R` already implements
the needed rules; only the config key is missing.

---

## 10. Reviewer objections — pre-written answers

| Objection | Answer |
|---|---|
| "You tuned on the test set" | Transfer/re-search split is explicit and the bias is quantified (0.011–0.017 F1). The headline row uses thresholds imported from a replicate not in the evaluation. |
| "F1 favours your operating point" | Headline is reported as TPR/FDR at a fixed threshold, using **their** `*_neutSNPs` convention verbatim (only `background_neutral` counts as a false positive). |
| "Cherry-picked replicates" | All 14 from the deposit's own selection, including 2 controls predicted *a priori* to fail. |
| "RDA underperforms in your hands" | Stated as a data limit, not a method limit: MVP has 2 predictors, `rda.R`'s `K_sel` floors at 2, so RDA cannot show the advantage it is recommended for. Our RDA reproduces theirs (0.1353 vs published 0.1373). |
| "Why is your LFMM temp TPR lower than published?" | 0.060 vs 0.084 — median \|Δ\| is 0.000; the difference is two small-denominator replicates (2 and 6 causal loci) where one locus flips TPR wholesale. |

---

## 11. Honest limitations — put these in the paper, not in a reviewer's mouth

1. **Two predictors.** MVP ships `bio_1`/`bio_2` only. Real projects in this repo use **19**.
   This is the largest external-validity gap, and it specifically handicaps RDA.
2. **One simulation model.** All 14 replicates share one SLiM model, Fst 0.171–0.217, a 10×10 deme
   grid. Transferability is demonstrated *between replicates*, not between models, and not on real
   data.
3. **No genetic-offset arm.** MVP has home-deme fitness only, never transplant fitness, so there is
   no honest offset target. That belongs to Capblancq et al. 2026 (Dryad doi:10.5061/dryad.jdfn2z3m6).
4. **Oligogenic architectures untested.** All four oligogenic replicates have <15 causal loci,
   where one locus moves F1 more than the differences being compared. Reported, excluded from votes.
5. **Recall denominators.** The cross-seed F1 tables use `truth_any`; the head-to-head tables use
   per-axis truth, because RDA has no per-axis p-value (`rda.R` writes one
   `climate_multivariate` column). Never mix the two.
6. **Per-seed `TOP_MAX`.** Rank-based rules are capped per replicate (2× testable causal, floor
   100), so absolute F1 is **not** comparable between replicates — only within-replicate deltas and
   ranks are.
7. **WZA excluded.** ~10k SNPs over 20 linkage groups gives order 10¹–10² windows, an order of
   magnitude short of what window statistics need.

---

## 12. Open work that would strengthen the paper

| Item | Why it matters | Cost |
|---|---|---|
| Score combinations at **`mode=pregea`-nominated** thresholds | Closes the gap between the imported-threshold row and the oracle row; would let us claim the pipeline calibrates itself without truth | ~1 h, same harness |
| Per-axis scoring at rungs other than c3 | Current head-to-head uses the default rung only | ~2 h |
| Cross-validated anchor (each replicate as anchor in turn) | Removes the single-anchor dependency in claim 5 | ~14× current scoring cost |
| A second simulation model | Addresses limitation 2 directly | large |

---

## 13. Suggested figures

1. **Precision–recall plane**, all published methods (their values) + our configurations, one panel
   per axis. Makes "no usable operating point" visible in one glance.
2. **Transfer gain box plot** — `report07/transfer_gain.png`, already generated.
3. **Selection-bias box plot** — `report07/selection_bias.png`, already generated.
4. **Two-factor decomposition** (§6) as a three-point cascade.
