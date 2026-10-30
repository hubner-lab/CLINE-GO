"""Release Phase 4b: the "missing tables" are targeted, each by its own new rule.

Seven numbers lived only in plots or logs (00d7a8 lambda_GC, 6bf499 Tracy-Widom, bffcfc the
correlograms' matrices, f59f46 cross-entropy, bb27e1 the per-trait thresholds, 922fa9 the SNP-set
summary; 1d2490 is a column, tested in R). Each table is written by a rule whose output is new,
so an existing results tree gains the file without rerunning anything: the plot and scan rules
those numbers come from keep their outputs, inputs and shell text (pca_plot's projections feed
EMMAX; the cross-entropy and trait correlogram PNGs are inputs of write_summary; a new output on
assoc_find_sig_snps would rewrite every sig table and rerun everything downstream of it).

These tests pin, for every mode, which file is asked for and which rule makes it, plus the two
helpers whose output is text: method_pvalue_correction() and snp_sets_summary_tsv().
"""
import json
import os
import re
import shutil
import unittest

import _smk_harness as H


def producer(patterns, path):
    """Names of the rules whose output pattern matches `path`."""
    return sorted({rule for regex, _smk, rule in patterns if regex.match(path)})


class _Loaded(unittest.TestCase):
    config = None

    @classmethod
    def setUpClass(cls):
        cls.ns = H.load(cls.config)
        H.save_snp_set(cls.ns["PROJECT"], "curated")
        cls.patterns = H.producible_patterns(cls.ns)
        cls.targets = {m: H.targets(cls.ns, m) for m in H.MODES}
        cls.k = cls.ns["K_BEST"]

    def assertMadeBy(self, mode, suffix_regex, rule):
        hits = [t for t in self.targets[mode] if re.search(suffix_regex, t)]
        self.assertEqual(len(hits), 1, f"mode={mode}: expected one target matching "
                                       f"{suffix_regex}, got {hits}")
        self.assertEqual(producer(self.patterns, hits[0]), [rule],
                         f"{hits[0]} must come from {rule} alone")


class TestNewTablesAreTargeted(_Loaded):

    def test_prestructure_tables(self):
        self.assertMadeBy("prestructure", r"/PreStructure/tables/tracy_widom\.tsv$",
                          "tracy_widom_table")
        self.assertMadeBy("prestructure", r"/PreStructure/tables/cross_entropy_K2-7\.tsv$",
                          "cross_entropy_table")

    def test_correlation_tables(self):
        self.assertMadeBy("climate", r"/climate/tables/present/climate_correlation\.tsv$",
                          "correlation_table")
        self.assertMadeBy("traits", r"/Traits/tables/trait_correlation\.tsv$",
                          "trait_correlation_table")
        self.assertMadeBy("traits", r"/Traits/tables/trait_climate_correlation\.tsv$",
                          "trait_climate_correlation_table")

    def test_structure_writes_the_invariant_predictor_table(self):
        # Release A2 decision (2026-10-30): mode=climate refuses at parse time on an invariant
        # predictor, so the table the app strikes chips from is written one mode earlier.
        self.assertMadeBy("structure", r"/climate/tables/present/climate_invariant_predictors\.tsv$",
                          "check_climate_variance")

    def test_lambda_and_threshold_tables_beside_every_pvalue_table(self):
        for mode, source, methods in (("gea", "GEA", ("EMMAX", "LFMM", "RDA")),
                                      ("gwas", "GWAS", ("EMMAX", "BLINK"))):
            for m in methods:
                base = rf"/{source}/tables/methods/{m}/{m}"
                self.assertMadeBy(mode, rf"{base}_lambda_K{self.k}\.tsv$", "assoc_genomic_inflation")
                self.assertMadeBy(mode, rf"{base}_thresholds_K{self.k}\.tsv$", "assoc_thresholds")

    def test_snp_sets_summary(self):
        self.assertMadeBy("maladaptation", r"/Maladaptation/tables/snp_sets_summary\.tsv$",
                          "snp_sets_summary")

    def test_new_names_escape_the_discovery_globs(self):
        # The app and check_invariants.R find p-value and sig tables by these globs; a new table
        # caught by one would be read as a p-value or sig table.
        new = [t for m in ("gea", "gwas") for t in self.targets[m]
               if re.search(r"_(lambda|thresholds)_K\d+\.tsv$", t)]
        self.assertEqual(len(new), 10)
        for t in new:
            name = os.path.basename(t)
            self.assertNotRegex(name, r"_pvalues_K")
            self.assertNotRegex(name, r"_sig_snps_")
            self.assertNotRegex(name, r"_wza_K")


class TestGwasOnlyAsksForNoClimateTable(unittest.TestCase):
    """gwas_only has no climate: the trait table stays, the climate-joined ones go.
    (test_workflow_targets.py's regime test already keeps /climate/ out of these modes.)"""

    @classmethod
    def setUpClass(cls):
        config = H.make_config(Regime__mode="gwas_only")
        config["GEA"]["configs"] = []      # RDA rejects Climate.enabled: false at parse
        cls.ns = H.load(config)

    def test_traits_keeps_its_own_table_only(self):
        traits = H.targets(self.ns, "traits")
        self.assertTrue(any(t.endswith("/Traits/tables/trait_correlation.tsv") for t in traits))
        self.assertFalse(any(t.endswith("/Traits/tables/trait_climate_correlation.tsv")
                             for t in traits))

    def test_structure_asks_for_no_invariant_predictor_table(self):
        self.assertFalse(any(t.endswith("climate_invariant_predictors.tsv")
                             for t in H.targets(self.ns, "structure")))


class TestPvalueCorrection(_Loaded):

    def test_lfmm_rda_and_the_rest(self):
        f = self.ns["method_pvalue_correction"]
        self.assertEqual(f("GEA", "LFMM"), "genomic_control")
        self.assertEqual(f("GEA", "RDA"), "rdadapt_gif")
        self.assertEqual(f("GEA", "EMMAX"), "none")
        self.assertEqual(f("GWAS", "BLINK"), "none")

    def test_lfmm_without_genomic_control(self):
        cfg = H.make_config(GEA__configs=[
            {"method": "LFMM", "adjust": "bonf", "threshold": "0.05",
             "params": {"genomic_control": False}}])
        self.assertEqual(H.load(cfg)["method_pvalue_correction"]("GEA", "LFMM"), "none")


class TestSnpSetsSummary(_Loaded):

    HEADER = ("set\torigin\tsource_module\tstrategy\tclumping_distance\tn_snps\t"
              "usable_for_offset\tused\treason")

    def test_all_reports_used_and_skipped_sets(self):
        records = {"gea.best": ({"origin": "config", "source_module": "GEA",
                                 "strategy": "Cross-method >=2", "clumping_distance": 3055}, 187),
                   "tiny": ({"name": "tiny", "source_module": "GEA", "strategy": "Union",
                             "clumping_distance": 1000000}, 2),
                   "bare": ({}, 5)}
        text = self.ns["snp_sets_summary_tsv"](["bare", "gea.best"], records)
        self.assertEqual(text.splitlines(), [
            self.HEADER,
            "bare\tunknown\t\t\t\t5\tTRUE\tTRUE\tused",
            "gea.best\tconfig\tGEA\tCross-method >=2\t3055\t187\tTRUE\tTRUE\tused",
            "tiny\tshiny\tGEA\tUnion\t1000000\t2\tFALSE\tFALSE\tbelow_min_offset_snps",
        ])

    def test_an_explicit_list_reports_the_rest_as_not_listed(self):
        ns = H.load(H.make_config(Maladaptation__snp_sets=["a"]))
        text = ns["snp_sets_summary_tsv"](["a"], {"a": ({}, 4), "b": ({}, 9)})
        self.assertEqual(text.splitlines()[1:], ["a\tunknown\t\t\t\t4\tTRUE\tTRUE\tused",
                                                 "b\tunknown\t\t\t\t9\tTRUE\tFALSE\tnot_listed"])

    def test_store_count_and_manifest_record(self):
        # The count is the resolver's own: data lines of selected_snps.tsv (save_snp_set: 3).
        with H.patched():
            self.assertEqual(self.ns["snp_set_snp_count"]("curated"), 3)
            self.assertIn("curated", self.ns["stored_snp_sets"]())
            self.assertEqual(self.ns["read_snp_set_record"]("curated"), {})
        older = H.save_snp_set(self.ns["PROJECT"], "older")
        store = os.path.dirname(older)
        with open(os.path.join(store, "manifest.json"), "w") as fh:
            json.dump([{"name": "older", "n_snps": 3, "source_module": "GEA",
                        "strategy": "Union", "clumping_distance": 1000000}], fh)
        try:
            with H.patched():
                rec = self.ns["read_snp_set_record"]("older")
            self.assertEqual(rec["clumping_distance"], 1000000)
            line = self.ns["snp_sets_summary_tsv"](["older"], {"older": (rec, 3)}).splitlines()[1]
            self.assertEqual(line, "older\tshiny\tGEA\tUnion\t1000000\t3\tTRUE\tTRUE\tused")
        finally:
            # The sandbox store is shared by every test module in this process.
            os.remove(os.path.join(store, "manifest.json"))
            shutil.rmtree(older)


if __name__ == "__main__":
    unittest.main()
