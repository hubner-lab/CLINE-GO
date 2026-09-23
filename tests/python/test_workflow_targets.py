"""workflow/rules/ — the config parse and the target graph.

The largest untested surface in the repository. ``common.smk`` (2 354 lines)
turns a config dict into the list of files each mode is asked to build; the
``.smk`` files turn rules into the files that can actually be built. Nothing
compared the two except ``snakemake -n``, which needs Docker, a VCF, a metadata
file and a real ``_results/`` tree — so it has never been part of the merge gate.

Five assertion families, each with a failure mode the R suites cannot reach:

  * TestTargetsAreProducible — every path ``get_targets(mode)`` returns matches
    some rule's ``output:``. An orphaned target is a ``MissingRule`` error today,
    raised by Snakemake only when someone runs that mode.
  * TestTargetSetShape — no duplicates, nothing outside ``{PROJECT}_results/``.
  * TestRegimeGating / TestNumericSpellingIsStable — two invariants that already
    hold and are worth pinning: ``gwas_only`` must not reach climate, and a
    config value spelled ``7.0`` must give the same paths as ``7`` (the bug
    common.smk:80-89 documents at length, which produced a second ``_work/``
    tree and ``cross_entropy_K2-7.0.png`` beside ``...K2-6.png``).
  * TestConfigValidation — configs the parser must refuse (inverted K sweep,
    K<1, k_best outside the sweep, an unknown ``adjust``, a non-numeric
    ``threshold``) together with the ones it must keep accepting. Every reject
    case here was a quarantined skip until common.smk grew
    ``_validate_k_sweep`` / ``_validate_adjust_threshold``.
  * TestConfigValidationGaps — what is LEFT of that quarantine. Correct-behaviour
    assertions for configs the parser still accepts and the run then cannot
    honour, each behind a ``skipTest`` naming where it is filed. Same convention
    as tests/testthat/test-known-bugs.R: fixing one means deleting a skip, and
    these must never be weakened to match current behaviour.

The harness that makes all of this possible is documented in _smk_harness.py.
"""
import unittest

import _smk_harness as H

FILED = ("filed in docs/pipeline_improvement_requests.md — config accepted at "
         "parse time, rejected (or silently ignored) at run time")


class _Loaded(unittest.TestCase):
    """Default config loaded once; every mode's target list resolvable."""

    @classmethod
    def setUpClass(cls):
        cls.ns = H.load()
        H.save_snp_set(cls.ns["PROJECT"], "curated")   # mode=maladaptation reads the FS
        cls.patterns = H.producible_patterns(cls.ns)
        cls.targets = {m: H.targets(cls.ns, m) for m in H.MODES}


class TestRuleGraphIsFullyModelled(_Loaded):
    """Guards the harness itself: an output it cannot evaluate is a blind spot."""

    def test_every_output_expression_evaluates(self):
        failures = [(f, r, e, err) for f, r, e, _v, err in H.rule_outputs(self.ns) if err]
        self.assertEqual(
            failures, [],
            "rule output expressions that could not be evaluated — every one of "
            "them is a rule whose products are invisible to TestTargetsAreProducible")

    def test_rule_graph_is_not_empty(self):
        self.assertGreater(len(self.patterns), 200)


class TestTargetsAreProducible(_Loaded):
    def test_every_target_matches_a_rule_output(self):
        for mode in H.MODES:
            with self.subTest(mode=mode):
                orphans = sorted({t for t in self.targets[mode]
                                  if not any(rx.match(t) for rx, _f, _r in self.patterns)})
                self.assertEqual(orphans, [],
                                 f"mode={mode} asks for files no rule can build")

    def test_unknown_mode_is_rejected(self):
        with self.assertRaises(ValueError):
            H.targets(self.ns, "not_a_mode")


class TestTargetSetShape(_Loaded):
    def test_no_duplicate_targets(self):
        for mode in H.MODES:
            with self.subTest(mode=mode):
                seen, dupes = set(), []
                for t in self.targets[mode]:
                    (dupes.append(t) if t in seen else seen.add(t))
                self.assertEqual(sorted(set(dupes)), [])

    def test_targets_stay_inside_the_results_tree(self):
        outdir = self.ns["OUTDIR"]
        for mode in H.MODES:
            with self.subTest(mode=mode):
                escaped = sorted({t for t in self.targets[mode] if not t.startswith(outdir)})
                self.assertEqual(escaped, [], f"mode={mode} writes outside {outdir}")

    def test_every_mode_requests_the_summary_flag(self):
        # pipeline_summary.tsv is what makes a mode a mode (summary.smk:66).
        for mode in H.MODES:
            with self.subTest(mode=mode):
                self.assertIn(self.ns["W"]["summary_done"], self.targets[mode])


class TestRegimeGating(unittest.TestCase):
    """gwas_only declares a project with no coordinates: climate must be unreachable."""

    @classmethod
    def setUpClass(cls):
        config = H.make_config(Regime__mode="gwas_only")
        config["GEA"]["configs"] = []      # RDA rejects Climate.enabled: false at parse
        cls.ns = H.load(config)

    def test_climate_is_forced_off(self):
        self.assertFalse(self.ns["CLIMATE_ENABLED"])

    def test_climate_dependent_modes_refuse_to_run(self):
        for mode in ("climate", "pregea", "gea", "maladaptation"):
            with self.subTest(mode=mode):
                with self.assertRaises(ValueError):
                    H.targets(self.ns, mode)

    def test_runnable_modes_emit_no_climate_targets(self):
        for mode in ("processing", "prestructure", "structure", "traits",
                     "gwas", "gea_x_gwas"):
            with self.subTest(mode=mode):
                leaked = [t for t in H.targets(self.ns, mode)
                          if "/climate/" in t or "/GEA/" in t or "/Maladaptation/" in t]
                self.assertEqual(leaked, [])

    def test_unknown_regime_is_rejected_at_parse(self):
        with self.assertRaises(ValueError):
            H.load(H.make_config(Regime__mode="gwas-only"))


class TestNumericSpellingIsStable(unittest.TestCase):
    """`k_end: 7.0` and `k_end: 7` must address the same files.

    The Shiny sidebar round-trips numeric fields through R's as.numeric(), so
    yaml::write_yaml() stores counts as doubles. Several reach output PATHS.
    """

    def test_float_and_int_spellings_give_identical_targets(self):
        as_int = H.make_config()
        as_int["sNMF"].update(k_end=7, k_best=3)
        as_int["Maladaptation"]["methods"]["gradient_forest"]["ntree"] = 1000
        as_int["Maladaptation"]["methods"]["geometric_offset"]["k"] = 3
        as_int["GWAS"]["promoter_length"] = 10000
        float_ns, int_ns = H.load(), H.load(as_int)
        for mode in [m for m in H.MODES if m != "maladaptation"]:
            with self.subTest(mode=mode):
                self.assertEqual(sorted(H.targets(float_ns, mode)),
                                 sorted(H.targets(int_ns, mode)))


class TestConfigValidation(unittest.TestCase):
    """Configs the parser must refuse, and the ones it must keep accepting.

    Every case here was a QUARANTINE skip until the parse-time checks landed
    (common.smk `_validate_k_sweep` / `_validate_adjust_threshold`): each one
    was accepted at parse time and then either ignored or rejected mid-run,
    after the expensive rules had already burned compute. The accept-side tests
    are not padding — a validator that is too strict fails the same way a
    missing one does, only louder.
    """

    def test_k_start_above_k_end_is_rejected(self):
        # k_range(9, 3) == [] while snmf.R's 9:3 is valid in R: sNMF runs, the
        # cross-entropy plot is drawn, no per-K output is extracted, exit 0.
        config = H.make_config()
        config["sNMF"].update(k_start=9, k_end=3)
        with self.assertRaises(ValueError):
            H.load(config)

    def test_k_below_one_is_rejected(self):
        # K=0 is not an ancestry model; the targets are generated anyway.
        config = H.make_config()
        config["sNMF"].update(k_start=0)
        with self.assertRaises(ValueError):
            H.load(config)

    def test_k_best_outside_the_swept_range_is_rejected(self):
        # K_BEST=99 while the sweep is 2..7: every downstream path keys on
        # K_BEST, so mode=gea asks for EMMAX_pvalues_K99.tsv from a Q-matrix
        # prestructure was never asked to produce.
        config = H.make_config()
        config["sNMF"].update(k_start=2, k_end=7, k_best=99)
        with self.assertRaises(ValueError):
            H.load(config)

    def test_unknown_adjust_is_rejected_at_parse(self):
        # compute_pval_threshold() stop()s on an unknown adjustment
        # (scripts/R/utils/pval_threshold.R:147) — but only after the method's
        # genome scan has run.
        config = H.make_config()
        config["GEA"]["configs"][0] = {"method": "EMMAX", "adjust": "bonferoni",
                                       "threshold": "0.05"}
        with self.assertRaises(ValueError):
            H.load(config)

    def test_non_numeric_threshold_is_rejected_at_parse(self):
        config = H.make_config()
        config["GEA"]["configs"][0] = {"method": "EMMAX", "adjust": "bonf",
                                       "threshold": "abc"}
        with self.assertRaises(ValueError):
            H.load(config)

    def test_top_n_below_one_is_rejected(self):
        # adjust='top' takes a SNP count; the app's
        # threshold_value_valid_for_type() requires >= 1.
        config = H.make_config()
        config["GEA"]["configs"][0] = {"method": "EMMAX", "adjust": "top",
                                       "threshold": "0.5"}
        with self.assertRaises(ValueError):
            H.load(config)

    def test_probability_thresholds_above_one_are_rejected(self):
        for adjust in ("bonf", "qval"):
            config = H.make_config()
            config["GEA"]["configs"][0] = {"method": "EMMAX", "adjust": adjust,
                                           "threshold": "5"}
            with self.subTest(adjust=adjust):
                with self.assertRaises(ValueError):
                    H.load(config)

    def test_gwas_configs_are_validated_too(self):
        config = H.make_config()
        config["GWAS"]["configs"][0] = {"method": "EMMAX", "adjust": "bonferoni",
                                        "threshold": "0.05"}
        with self.assertRaisesRegex(ValueError, "GWAS.configs"):
            H.load(config)

    def test_every_supported_rule_still_parses(self):
        # The four rules compute_pval_threshold() implements must all survive
        # the new check, and must still reach the target list verbatim.
        for adjust, threshold in (("bonf", "0.05"), ("qval", "0.1"),
                                  ("top", "50"), ("custom", "1e-6")):
            config = H.make_config()
            config["GEA"]["configs"] = [{"method": "EMMAX", "adjust": adjust,
                                         "threshold": threshold}]
            with self.subTest(rule=f"{adjust}_{threshold}"):
                ns = H.load(config)
                self.assertEqual(ns["GEA_CONFIGS"], {"EMMAX": f"{adjust}_{threshold}"})
                self.assertTrue(any(f"sig_snps_{adjust}_{threshold}.tsv" in t
                                    for t in H.targets(ns, "gea")))

    def test_every_app_valid_threshold_still_parses(self):
        # Boundary values the app's threshold_value_valid_for_type()
        # (fct_combine.R) accepts. Refusing them here would make the Shiny
        # sidebar write configs the pipeline cannot parse.
        for adjust, threshold in (("bonf", "1"), ("qval", "1"), ("top", "1"),
                                  ("custom", "5")):
            config = H.make_config()
            config["GEA"]["configs"] = [{"method": "EMMAX", "adjust": adjust,
                                         "threshold": threshold}]
            with self.subTest(rule=f"{adjust}_{threshold}"):
                ns = H.load(config)
                self.assertEqual(ns["GEA_CONFIGS"], {"EMMAX": f"{adjust}_{threshold}"})

    def test_k_start_one_is_still_allowed(self):
        # K=1 is LEA's cross-entropy baseline, and the app's config schema
        # (fct_config_schema.R) sets min = 1 for sNMF.k_start.
        config = H.make_config()
        config["sNMF"].update(k_start=1)
        ns = H.load(config)
        self.assertEqual(ns["K_START"], 1)

    def test_k_best_null_is_still_allowed(self):
        # k_best is unset until the cross-entropy plot has been read; every
        # pre-structure mode has to parse without it.
        config = H.make_config()
        config["sNMF"]["k_best"] = None
        ns = H.load(config)
        self.assertIsNone(ns["K_BEST"])

    def test_the_shipped_k_sweep_still_parses(self):
        ns = H.load()
        self.assertEqual(ns["K_START"], 2)
        self.assertEqual(ns["K_END"], 7)
        self.assertEqual(ns["K_BEST"], 3)


class TestConfigValidationGaps(unittest.TestCase):
    """QUARANTINE — correct behaviour for configs the parser wrongly accepts.

    Every skip names a filing. Deleting a skip is how one gets fixed. None of
    these may be rewritten to assert what the code does today.
    """

    def test_out_of_range_frequencies_are_rejected(self):
        # check_float() only tests castability, so maf=-1 reaches plink AND the
        # _work/ directory tag: _work/maf-1_miss0.2_smiss0.5_rel0.99/.
        self.skipTest("Filter.maf / Filter.snp_miss / LD.r2 have no range check — " + FILED)
        for group, key, value in (("Filter", "maf", -1.0), ("Filter", "maf", 0.9),
                                  ("Filter", "snp_miss", 5.0), ("LD", "r2", 7.0)):
            config = H.make_config()
            config[group][key] = value
            with self.subTest(setting=f"{group}.{key}={value}"):
                with self.assertRaises(ValueError):
                    H.load(config)


if __name__ == "__main__":
    unittest.main()
