# GEA (Genotype-Environment Association) method registry.
#
# Each entry defines one association method. Adding a new method requires:
#   1. A dict entry here.
#   2. An R script at the declared 'script' path that reads the standard
#      input formats for its engine and writes pvalue/qvalue TSVs.
#
# engine:              dispatch key used by rule factories ("emmax"|"lfmm"|"gapit"|"rda")
# script:              Rscript path relative to /pipeline/
# gapit_model:         model name passed to gapit.R (gapit engine only)
# supports_phenotypes: method can be used in phenotype_association (GWAS) mode
# supports_drop:       method supports DROP mode (per-trait sample subsetting)
# multivariate:        method emits ONE p-value column for the whole predictor set
#                      (a pseudo-trait), not one column per predictor
# pseudo_trait:        the single column name when multivariate=True, else None
# supports_wza:        method's p-values can be fed through compute_wza.R
# adjust_default:      default significance rule name for this method's
#                      GEA.configs row ("bonf"|"qval"|"top"|"toppct"|"custom"). Seeds
#                      the Shiny config sidebar's method editor and the GEA
#                      tab's per-method threshold table.
# threshold_default:   default value paired with adjust_default, as a STRING —
#                      the config YAML stores thresholds as strings and the
#                      value is a filename/wildcard component (see
#                      _assoc_downstream.smk's `adjust` wildcard constraint).
# significance_family: what the method's p-value column MEANS.
#                      "univariate_pvalue" = one test per (SNP, trait).
#                      "multivariate_pvalue" = one joint test per SNP over the
#                      whole predictor set (RDA's rdadapt). A future
#                      Bayes-factor method (e.g. BayPass) would declare
#                      "bayes_factor" here so the threshold layer can dispatch
#                      on family instead of special-casing the method name.
# params:              declarative per-method hyperparameters, surfaced in the
#                      Shiny GEA sidebar's method editor and resolved at
#                      config-parse time by resolve_method_params() in common.smk.
#                      Each spec is built by P(type, default, **kw) below.
#                      default=K_BEST_SENTINEL resolves to sNMF.k_best — this is
#                      the sole mechanism that keeps existing configs' behavior
#                      byte-identical (LFMM K / EMMAX #PCs / GAPIT PCA.total all
#                      already equal k_best today; the sentinel preserves that
#                      without hard-coding a literal that would silently diverge
#                      from a project's own k_best on upgrade).

K_BEST_SENTINEL = "@k_best"


def P(type, default, **kw):
    """One hyperparameter spec. kw: min, max, choices, help."""
    return {"type": type, "default": default, **kw}


def _gapit(model):
    return {
        "engine": "gapit",
        "script": "scripts/gapit.R",
        "gapit_model": model,
        "supports_phenotypes": True,
        "supports_drop": True,
        "multivariate": False,
        "pseudo_trait": None,
        "supports_wza": True,
        "adjust_default": "bonf",
        "threshold_default": "0.05",
        "significance_family": "univariate_pvalue",
        "params": {
            "n_pcs": P("int", K_BEST_SENTINEL, min=0, max=20,
                       help="PCA.total passed to GAPIT. Default: sNMF k_best."),
        },
    }


GEA_METHODS = {
    "EMMAX": {
        "engine": "emmax",
        "script": "scripts/emmax.R",
        "supports_phenotypes": True,
        "supports_drop": True,
        "multivariate": False,
        "pseudo_trait": None,
        "supports_wza": True,
        "adjust_default": "bonf",
        "threshold_default": "0.05",
        "significance_family": "univariate_pvalue",
        "params": {
            "n_pcs": P("int", K_BEST_SENTINEL, min=0, max=20,
                       help="Number of PCA covariates used as fixed effects. "
                            "Default: sNMF k_best."),
            "kinship": P("enum", "BN", choices=["BN", "IBS"],
                         help="emmax-kin estimator: BN (Balding-Nichols, default) "
                              "or IBS (emmax-kin -s flag)."),
        },
    },
    "LFMM": {
        "engine": "lfmm",
        "script": "scripts/lfmm.R",
        "supports_phenotypes": False,
        "supports_drop": False,
        "multivariate": False,
        "pseudo_trait": None,
        "supports_wza": True,
        "adjust_default": "bonf",
        "threshold_default": "0.05",
        "significance_family": "univariate_pvalue",
        "params": {
            "K": P("int", K_BEST_SENTINEL, min=1, max=20,
                   help="Number of latent factors in lfmm2(). Default: sNMF "
                        "k_best. Read the preGEA K-ladder before overriding."),
            "genomic_control": P("bool", True,
                   help="lfmm2.test(genomic.control=...). Default TRUE — "
                        "matches production behavior. PreGEA's own K-ladder "
                        "always fits with this OFF so lambda_GC stays "
                        "informative (docs/rda_research.md C.0); this "
                        "switch only affects the production GEA run."),
        },
    },
    "RDA": {
        "engine": "rda",
        "script": "scripts/rda.R",
        # False keeps RDA out of GWAS_METHODS (gwas.py filters GEA_METHODS by
        # this flag) — RDA is GEA-only in this pipeline.
        "supports_phenotypes": False,
        "supports_drop": False,
        "multivariate": True,
        "pseudo_trait": "climate_multivariate",
        "supports_wza": True,
        "adjust_default": "bonf",
        # Capblancq & Forester 2021 (MEE 12:2298-2309) apply Bonferroni 0.01/m
        # on the full marker count — docs/rda_research.md A.0/A.6. NOT 0.05:
        # rdadapt emits ONE joint test per SNP, so the univariate 0.05
        # convention has no standing here. The 2018 companion paper (Capblancq
        # et al.) uses q<0.1 instead; the disagreement between the two papers
        # is unresolved in the literature (rda_research.md A.4 item 1) — rda.R
        # keeps a 4-rule sensitivity ladder in its diagnostics table so both
        # published rules stay visible alongside whichever one is applied.
        "threshold_default": "0.01",
        "significance_family": "multivariate_pvalue",
        "params": {
            "axes": P("str", "auto",
                      help="K = retained constrained axes. 'auto' = axes with "
                           "anova.cca(by='axis') p < axis_alpha, capped at "
                           "RDA$CCA$rank, floored at 2. Or a fixed integer >= 2."),
            "axis_alpha": P("float", 0.05, min=0.001, max=0.5,
                            help="Per-axis significance cutoff when axes='auto'."),
            "condition_pcs": P("int", K_BEST_SENTINEL, min=0, max=20,
                               help="Number of LD-pruned PCs in Condition(). "
                                    "0 disables structure correction (uncorrected "
                                    "scan — not the literature default)."),
            "predictor_set": P("str", "auto",
                               help="'auto' = Climate.predictors. Or a comma-"
                                    "separated subset (e.g. from a preGEA "
                                    "ordiR2step ladder). Never auto-applied."),
            "permutations": P("int", 999, min=99, max=9999,
                              help="anova.cca permutations (literature: 999)."),
            "fit_mode": P("enum", "auto", choices=["auto", "full", "pruned"],
                          help="'full' = literature default (Capblancq & Forester "
                               "2021 explicitly do NOT pre-prune). 'pruned' is an "
                               "ENGINEERING FALLBACK, not a scientific alternative. "
                               "'auto' uses full unless the size gate trips."),
            "full_fit_max_snps": P("int", 2000000, min=10000,
                                   help="SNP-count gate for fit_mode='auto'."),
            "full_fit_max_gb": P("float", 64.0, min=1.0,
                                 help="Estimated dense-matrix budget (GB) for "
                                      "fit_mode='auto'. Estimate = 8*n*m*4 bytes."),
            "seed": P("int", 42, help="RNG seed for anova.cca permutations."),
        },
    },
    "BLINK":   _gapit("BLINK"),
    "FarmCPU": _gapit("FarmCPU"),
    "MLM":     _gapit("MLM"),
    "MLMM":    _gapit("MLMM"),
    "GLM":     _gapit("GLM"),
    "CMLM":    _gapit("CMLM"),
    "ECMLM":   _gapit("ECMLM"),
    "SUPER":   _gapit("SUPER"),
}
