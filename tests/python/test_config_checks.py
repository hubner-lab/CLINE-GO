"""Config checks: a config no run can honour is refused at parse time, by name.

Parse-time checks run on every load of common.smk. The data-dependent ones run in
get_targets() for the mode that would fail, read a table an earlier mode wrote, and
skip while it is not on disk (the consuming script then still stops at run time).

Pinned, one class per finding:
  * 99f9ed  counts are integers: a fraction or a boolean is refused, never truncated;
            an integral float or numeric string ("7.0", '10') passes, and the value
            that reaches paths and shell stays as written.
  * 74f9f7  Filter.maf in [0, 0.5], snp_miss / sample_miss / relatedness and LD.r2 in
            [0, 1] (the app sidebar's bounds); the work-dir tag is unchanged.
  * 173795  a boolean is a YAML bool, 0/1 or an explicit word: 'off' is false, and any
            other word is refused (it used to read as true).
  * b5fcd3  Map.resolution is a WorldClim resolution; Map.climate_extent is 'auto',
            'world' or one "min_lon,max_lon,min_lat,max_lat" string (a YAML list
            expanded to four argv words). WorldClim only.
  * 5dad05  custom climate columns are named bio_<N>.
  * de3737  an unreadable, empty or non-tab metadata file is an error, not "no traits".
  * 2eb456  geometric_offset.k: null means K_BEST, as '' does.
  * 0c4bcc  a gwas_only project with GEA.configs (config_default.yaml ships them)
            builds mode=gea_x_gwas from GWAS alone.
  * a3a9ed  custom future_tables with the single-future plot or sensitivity check.
  * 130e0e + 07991d  a fitted predictor that does not vary across the samples (zero
            variance, one value, all NA) is refused by name in climate, pregea, gea
            and maladaptation; an invariant predictor nobody fits is not.
  * 8dbf84  an auto_* distance whose LD-decay rule the table cannot give is refused in
            mode=gea / mode=gwas; a declared distance is not checked. The decision
            agrees with region_distance.R's on the same tables (R is run).
  * 895bac  a spatial Gradient Forest with Climate.predictors >= sites - 1.
  * 932a8d + e83e54 wiring: the invariance check reads the predictor universe by
            name; the climate summary counts the selected predictors.
Phase 2b (input formats):
  * 1d4be2  metadata columns 1-4 are site, sample, latitude, longitude by name.
  * 634e7c  a GWAS trait the {pheno_trait} wildcard cannot match is refused by name.
  * 899581  sNMF.ploidy other than 2 is refused (sNMF runs diploid only).
"""
import os
import shutil
import subprocess
import tempfile
import unittest

import _smk_harness as H
from test_emmax_imputed_input import rule_graph

GEA_NO_RDA = [{"method": "EMMAX", "adjust": "bonf", "threshold": "0.05"},
              {"method": "LFMM", "adjust": "bonf", "threshold": "0.05"}]


def load(**overrides):
    return H.load(H.make_config(**overrides))


def write(path, text, binary=False):
    """Write a /pipeline/... file where the harness reads it (its sandbox)."""
    real = os.path.join(H._SANDBOX, path[len("/pipeline/"):])
    os.makedirs(os.path.dirname(real), exist_ok=True)
    with open(real, "wb" if binary else "w") as handle:
        handle.write(text)
    return real


def tsv(*rows):
    return "".join("\t".join(str(v) for v in row) + "\n" for row in rows)


class TestIntegers(unittest.TestCase):          # 99f9ed
    def test_fractions_and_booleans_are_refused(self):
        for overrides in ({"cpu": 0.5}, {"cpu": True}, {"sNMF__repeats": 2.9},
                          {"sNMF__repeats": False}, {"sNMF__k_best": 3.5},
                          {"sNMF__k_end": "7.5"}):
            with self.subTest(overrides=overrides):
                with self.assertRaisesRegex(ValueError, "must be an integer"):
                    load(**overrides)

    def test_counts_below_one_are_refused(self):
        with self.assertRaisesRegex(ValueError, "cpu must be >= 1"):
            load(cpu=0)
        with self.assertRaisesRegex(ValueError, "repeats must be an integer >= 1"):
            load(sNMF__repeats=0)

    def test_integral_spellings_pass_and_stay_as_written(self):
        ns = load(cpu=4.0, sNMF__repeats="10", sNMF__k_best="3.0")
        self.assertEqual(ns["CPU"], 4.0)
        self.assertEqual(ns["REPEAT"], 10)
        self.assertEqual(ns["K_BEST"], 3)
        self.assertEqual(ns["K_END"], 7)            # the fixture's k_end: 7.0


class TestRanges(unittest.TestCase):            # 74f9f7
    def test_out_of_range_values_are_refused(self):
        for key, value in (("Filter__maf", -1), ("Filter__maf", 0.9),
                           ("Filter__snp_miss", 5), ("Filter__sample_miss", 1.5),
                           ("Filter__relatedness", 1.2), ("LD__r2", 7), ("LD__r2", -0.1)):
            with self.subTest(key=key, value=value):
                name = key.replace("__", ".").lower()
                with self.assertRaisesRegex(ValueError, rf"{name} must be a number in"):
                    load(**{key: value})

    def test_a_boolean_is_not_a_number(self):
        with self.assertRaisesRegex(ValueError, "filter.maf must be float"):
            load(Filter__maf=True)

    def test_bounds_are_inclusive(self):
        for overrides in ({"Filter__maf": 0}, {"Filter__snp_miss": 1}, {"LD__r2": 1},
                          {"LD__r2": 0}, {"Filter__relatedness": 0},
                          {"Filter__maf": 0.5, "sNMF__imputed_maf": 0.1}):
            with self.subTest(overrides=overrides):
                load(**overrides)

    def test_path_tags_keep_the_written_values(self):
        ns = load()
        self.assertEqual(ns["FILT_TAG"], "maf0.05_miss0.2_smiss0.5_rel0.99")
        self.assertEqual(ns["LD_TAG"], "ld0.2_win100_step20")


class TestBooleans(unittest.TestCase):          # 173795
    def test_words(self):
        for value, want in (("off", "F"), ("on", "T"), ("Yes", "T"), ("no", "F"),
                            ("n", "F"), (1, "T"), (0, "F"), (True, "T"), ("", "F")):
            with self.subTest(value=value):
                self.assertEqual(load(Piemap__show_labels=value)["PIEMAP_SHOW_LABELS"], want)

    def test_climate_enabled_off_is_off(self):
        ns = load(Climate__enabled="off", GEA__configs=GEA_NO_RDA)
        self.assertFalse(ns["CLIMATE_ENABLED"])

    def test_null_reads_as_false(self):
        self.assertFalse(load(Population__calc_stats=None)["CALC_POP_STATS"])

    def test_other_words_are_refused(self):
        for key, value in (("Population__calc_stats", "none"),
                           ("Climate__enabled", "disabled"), ("Piemap__show_labels", 2)):
            with self.subTest(key=key, value=value):
                with self.assertRaisesRegex(ValueError, "must be true or false"):
                    load(**{key: value})

    def test_method_bool_params_share_the_rule(self):
        def lfmm(gc):
            return [{"method": "LFMM", "adjust": "bonf", "threshold": "0.05",
                     "params": {"genomic_control": gc}}]
        self.assertIs(load(GEA__configs=lfmm("off"))["GEA_PARAMS"]["LFMM"]["genomic_control"], False)
        with self.assertRaisesRegex(ValueError, "genomic_control"):
            load(GEA__configs=lfmm("none"))


class TestMap(unittest.TestCase):               # b5fcd3
    def test_resolution_outside_worldclim_is_refused(self):
        for value in (1.0, 7, "30s", True):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "Map.resolution"):
                    load(Map__resolution=value)

    def test_worldclim_resolutions_pass_as_written(self):
        self.assertEqual(load(Map__resolution="2.5")["RESOLUTION"], "2.5")
        load(Map__resolution=10)

    def test_extent_must_be_one_string_of_four_numbers(self):
        for value in ([34, 36, 29, 34], "36,34,29,34", "34,36", "a,b,c,d", "34,36,34,29"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "Map.climate_extent"):
                    load(Map__climate_extent=value)
        load(Map__climate_extent="34,36,29,34")
        load(Map__climate_extent="world")

    def test_not_checked_without_worldclim(self):
        load(Regime__mode="gwas_only", GEA__configs=GEA_NO_RDA, Map__resolution=1.0)
        load(Climate__source="custom", Map__resolution=1.0,
             Climate__custom={"present_table": "env.tsv", "future_table": "fut.tsv",
                              "columns": "bio_1,bio_2,bio_3"})


class TestCustomColumns(unittest.TestCase):     # 5dad05
    def custom(self, columns):
        return load(Climate__source="custom",
                    Climate__custom={"present_table": "env.tsv", "future_table": "fut.tsv",
                                     "columns": columns})

    def test_non_bio_names_are_refused(self):
        with self.assertRaisesRegex(ValueError, r"bio_<N>.*tmin_jan"):
            self.custom("tmin_jan,bio_2")

    def test_bio_names_pass(self):
        self.assertEqual(self.custom("bio_1,bio_2,bio_3")["ALL_BIO"], ["bio_1", "bio_2", "bio_3"])


class TestMetadata(unittest.TestCase):          # de3737
    def config(self, name, content, binary=False):
        folder = f"cc_meta_{name}/"
        write(f"/pipeline/{folder}x.vcf", "##fileformat=VCFv4.2\n")
        write(f"/pipeline/{folder}meta.tsv", content, binary=binary)
        return H.make_config(Input={"dir": folder, "vcf": "x.vcf", "metadata": "meta.tsv",
                                    "gff": "testdata.gff3"})

    def test_comma_separated_file_is_named(self):
        cfg = self.config("comma", "site,sample,latitude,longitude,height\nA,s1,30,34,1.2\n")
        with self.assertRaisesRegex(ValueError, "tab-separated.*comma-separated"):
            H.load(cfg)

    def test_empty_file(self):
        with self.assertRaisesRegex(ValueError, "is empty"):
            H.load(self.config("empty", ""))

    def test_unreadable_file(self):
        cfg = self.config("bytes", b"site\tsample\tlatitude\tlongitude\n\xff\xfe\t\x80\n",
                          binary=True)
        with self.assertRaisesRegex(ValueError, "cannot be read"):
            H.load(cfg)

    def test_too_few_columns(self):
        with self.assertRaisesRegex(ValueError, "header has 3"):
            H.load(self.config("three", "site\tsample\theight\nA\ts1\t1.2\n"))

    def test_four_columns_mean_no_traits(self):
        ns = H.load(self.config("four", "site\tsample\tlatitude\tlongitude\nA\ts1\t30\t34\n"))
        self.assertFalse(ns["META_HAS_PHENO"])


class TestGeometricK(unittest.TestCase):        # 2eb456
    KEY = "Maladaptation__methods__geometric_offset__k"

    def test_null_is_k_best(self):
        self.assertEqual(load(**{self.KEY: None})["GO_K"], "")

    def test_an_integer_stays_as_written(self):
        self.assertEqual(load(**{self.KEY: 3.0})["GO_K"], "3.0")

    def test_a_fraction_or_word_is_refused(self):
        for value in (2.5, "abc"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "geometric_offset.k must be an integer"):
                    load(**{self.KEY: value})


class TestGwasOnlyOverlap(unittest.TestCase):   # 0c4bcc
    def rules(self, ns):
        rules = rule_graph(ns)
        pair = next(r for r in rules if r["name"] == "compute_pairwise_overlaps")
        summary = next(r for r in rules if r["name"] == "write_summary"
                       and "overlap_miami" in r["raw"].get("input", ""))
        return pair, summary

    def test_gwas_only_reads_no_gea_table(self):
        ns = load(Regime__mode="gwas_only", GEA__configs=GEA_NO_RDA)
        self.assertFalse(ns["GEA_RUNS"])
        pair, summary = self.rules(ns)
        self.assertEqual(pair["input"]["selected"], [ns["O"]["pheno_selected_snps"]])
        self.assertEqual(pair["params"]["gea_snps"], "NULL")
        self.assertEqual(summary["input"]["_0"], [])
        self.assertNotIn(ns["O"]["overlap_miami"], H.targets(ns, "gea_x_gwas"))

    def test_standard_reads_both(self):
        ns = load()
        self.assertTrue(ns["GEA_RUNS"])
        pair, summary = self.rules(ns)
        self.assertEqual(pair["input"]["selected"],
                         [ns["O"]["selected_snps"], ns["O"]["pheno_selected_snps"]])
        self.assertEqual(summary["input"]["_0"], [ns["O"]["overlap_miami"]])


class TestCustomFutureScenarios(unittest.TestCase):   # a3a9ed
    def targets(self, project, future, **overrides):
        custom = {"present_table": "env.tsv", "columns": "bio_1,bio_2,bio_3", **future}
        ns = load(project_name=project, Climate__source="custom", Climate__custom=custom,
                  **overrides)
        H.save_snp_set(project, "setA")
        return H.targets(ns, "maladaptation")

    def test_future_tables_with_the_single_future_steps_are_refused(self):
        with self.assertRaisesRegex(ValueError, r"future_tables.*fut_a, fut_b.*emit_plots"):
            self.targets("cc_fut_refused", {"future_tables": ["fut_a.tsv", "fut_b.tsv"]})

    def test_future_tables_pass_without_them(self):
        self.targets("cc_fut_off", {"future_tables": ["fut_a.tsv", "fut_b.tsv"]},
                     Maladaptation__emit_plots=False,
                     Maladaptation__methods__gradient_forest__sensitivity_check=False)

    def test_one_future_table_passes(self):
        self.targets("cc_fut_one", {"future_table": "fut.tsv"})


# A climate_present_site.tsv: bio_1 and bio_4 vary; bio_2 and bio_14 are constant,
# bio_3 has one value, bio_5 none.
SITE = tsv(("sample", "bio_1", "bio_2", "bio_3", "bio_4", "bio_5", "bio_14"),
           ("s1", 1.5, 5, 7, 4, "", 0),
           ("s2", 2.5, 5, "", 5, "", 0),
           ("s3", 3.5, 5, "NA", 6, "", 0))


class TestInvariantPredictors(unittest.TestCase):    # 130e0e, 07991d
    def targets(self, project, mode, site=SITE, saved_set=False, **overrides):
        ns = load(project_name=project, **overrides)
        if site is not None:
            write(ns["O"]["climate_site"], site)
        if saved_set:
            H.save_snp_set(project, "setA")
        return H.targets(ns, mode)

    def test_climate_refuses_a_zero_variance_predictor(self):
        with self.assertRaisesRegex(ValueError, r"Climate.predictors: bio_2 \(zero_variance\)"):
            self.targets("cc_inv_climate", "climate", Climate__predictors="bio_1,bio_2")

    def test_gea_refuses_a_one_value_predictor(self):
        with self.assertRaisesRegex(ValueError, r"bio_3 \(single_value\)"):
            self.targets("cc_inv_gea", "gea", Climate__predictors="bio_1,bio_3")

    def test_gea_checks_an_explicit_rda_predictor_set(self):
        rda = {"method": "RDA", "adjust": "bonf", "threshold": "0.01",
               "params": {"predictor_set": "bio_1,bio_2"}}
        with self.assertRaisesRegex(ValueError, r"GEA.configs RDA params.predictor_set: bio_2"):
            self.targets("cc_inv_rda", "gea", Climate__predictors="bio_1,bio_4",
                         GEA__configs=GEA_NO_RDA + [rda])

    def test_pregea_refuses_an_all_na_predictor(self):
        with self.assertRaisesRegex(ValueError, r"PreGEA.predictors: bio_5 \(all_na\)"):
            self.targets("cc_inv_pregea", "pregea", PreGEA__predictors="bio_1,bio_5")

    def test_maladaptation_names_its_own_list(self):
        with self.assertRaisesRegex(ValueError, r"Maladaptation.predictors: bio_2"):
            self.targets("cc_inv_mala", "maladaptation", saved_set=True,
                         Climate__predictors="bio_1,bio_2,bio_4",
                         Maladaptation__predictors="bio_1,bio_2")

    def test_an_invariant_predictor_nobody_fits_passes(self):
        # SIMDATA: bio_14 is invariant and not in Climate.predictors.
        self.targets("cc_inv_unused", "climate", Climate__predictors="bio_1,bio_4")
        self.targets("cc_inv_unused", "gea", Climate__predictors="bio_1,bio_4")

    def test_skipped_before_mode_structure(self):
        self.targets("cc_inv_fresh", "gea", site=None, Climate__predictors="bio_1,bio_2")


# LD-decay tables in the columns ld_decay_analyze.R writes. Reliable rows carry the
# reason field as fwrite quotes an empty string.
LD_HEADER = ("group", "scope", "short_r2", "background_r2", "decay_critical_bp",
             "decay_half_bp", "decay_r2thr_bp", "r2_threshold", "decay_method",
             "ld_status", "ld_status_reason")
WHY = "short-range r2 0.3 is 1.4x the critical r2 of unlinked markers (0.2), below 2"
LD_TABLES = {
    # group All: reliable genome-wide, chr 2 unreliable, chr 3 misses r2_threshold.
    # C2: unreliable everywhere. C3: genome-wide misses critical, chromosomes fine.
    "mixed": tsv(LD_HEADER,
                 ("All", "genome_wide", 0.44, 0.09, 33036, 13141, 14974, 0.2, "loess_log", "reliable", '""'),
                 ("All", "1", 0.45, 0.09, 32978, 9320, 10202, 0.2, "loess_log", "reliable", '""'),
                 ("All", "2", 0.17, 0.09, 32521, 34209, "", 0.2, "loess_log", "unreliable", WHY),
                 ("All", "3", 0.43, 0.09, 31953, 13403, "", 0.2, "loess_log", "reliable", '""'),
                 ("C2", "genome_wide", 0.32, 0.23, 2407, 11571, 3382, 0.2, "loess_log", "unreliable", WHY),
                 ("C2", "1", 0.24, 0.23, 2159, 7363, 2858, 0.2, "loess_log", "unreliable", WHY),
                 ("C3", "genome_wide", 0.58, 0.23, "", 15874, 27940, 0.2, "loess_log", "reliable", '""'),
                 ("C3", "1", 0.56, 0.23, 22852, 12620, 25508, 0.2, "loess_log", "reliable", '""')),
    # written before the reliability check: no ld_status column, so each row's status is
    # derived from short_r2 / background_r2 (96e4fd): All 4.9x reliable, C2 1.39x not.
    "no_status": tsv(LD_HEADER[:9],
                     ("All", "genome_wide", 0.44, 0.09, 33036, 13141, 14974, 0.2, "loess_log"),
                     ("All", "1", 0.45, 0.09, "", 9320, 10202, 0.2, "loess_log"),
                     ("C2", "genome_wide", 0.32, 0.23, 2407, 11571, 3382, 0.2, "loess_log"),
                     ("C2", "1", 0.24, 0.23, 2159, 7363, 2858, 0.2, "loess_log")),
    "empty": tsv(LD_HEADER),
    # the PopLDdecay-era table: no rule columns
    "old": tsv(("group", "scope", "half_decay_bp"), ("All", "genome_wide", 1500)),
}


class TestLdRuleDecisionMatchesR(unittest.TestCase):    # 8dbf84
    """The parse-time mirror and region_distance.R decide alike on every case."""

    def test_python_and_r_agree(self):
        rscript = shutil.which("Rscript")
        self.assertIsNotNone(rscript, "Rscript is in the image; the agreement needs it")
        ns = load()
        tmp = tempfile.mkdtemp(prefix="clinego_ldrule_")
        try:
            cases = []
            for name, text in LD_TABLES.items():
                path = os.path.join(tmp, f"{name}.tsv")
                with open(path, "w") as handle:
                    handle.write(text)
                header, rows = ns["_read_tsv"](path)
                for spec in ("auto_genome_wide", "auto_per_chromosome"):
                    for ld_rule in ("critical", "half_decay", "r2_threshold"):
                        for group in ("All", "C2", "C3", "X"):
                            stops = ns["_ld_auto_refusal"](header, rows, spec, ld_rule, group)
                            cases.append((path, spec, ld_rule, group, stops is not None))
            manifest = os.path.join(tmp, "cases.tsv")
            with open(manifest, "w") as handle:
                handle.write(tsv(*[c[:4] for c in cases]))
            script = (
                "suppressPackageStartupMessages(library(data.table)); "
                f"source('{H.REPO_ROOT}/scripts/R/lib/region_distance.R'); "
                "cases <- read.delim(commandArgs(TRUE)[1], header = FALSE, colClasses = 'character'); "
                "for (i in seq_len(nrow(cases))) cat(tryCatch({ suppressMessages("
                "resolve_clumping_distance(cases[i, 2], cases[i, 1], cases[i, 3], cases[i, 4])); "
                "'OK' }, error = function(e) 'STOP'), '\\n', sep = '')")
            out = subprocess.run([rscript, "-e", script, manifest], capture_output=True,
                                 text=True, timeout=120)
            self.assertEqual(out.returncode, 0, out.stderr)
            r_stops = [line == "STOP" for line in out.stdout.split()]
            self.assertEqual(len(r_stops), len(cases), out.stdout)
            disagree = [(os.path.basename(c[0]), c[1], c[2], c[3], "python stops" if c[4] else "R stops")
                        for c, r in zip(cases, r_stops) if c[4] != r]
            self.assertEqual(disagree, [])
            self.assertTrue(any(r_stops) and not all(r_stops), "the cases must cover both outcomes")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestLdRuleRefusal(unittest.TestCase):     # 8dbf84
    def targets(self, project, mode, table="mixed", **overrides):
        ns = load(project_name=project, **overrides)
        write(ns["O"]["ld_decay_table"], LD_TABLES[table])
        return H.targets(ns, mode)

    def test_an_unreliable_genome_wide_curve_is_refused(self):
        with self.assertRaisesRegex(ValueError, r"GEA.snp_clumping_distance = auto_genome_wide.*unreliable"):
            self.targets("cc_ld_unrel", "gea", GEA__snp_clumping_distance="auto_genome_wide",
                         GEA__ld_decay_group="C2")

    def test_a_rule_not_reached_is_refused(self):
        with self.assertRaisesRegex(ValueError, r"rule 'critical' is not reached"):
            self.targets("cc_ld_reach", "gea", GEA__snp_clumping_distance="auto_genome_wide",
                         GEA__ld_decay_group="C3")

    def test_per_chromosome_falls_back_to_a_reliable_genome_wide_value(self):
        self.targets("cc_ld_perchr", "gea", GEA__snp_clumping_distance="auto_per_chromosome",
                     GEA__clumping_ld_rule="r2_threshold")

    def test_a_reliable_curve_passes_in_both_modes(self):
        self.targets("cc_ld_ok", "gea", GEA__snp_clumping_distance="auto_genome_wide")
        self.targets("cc_ld_ok", "gwas", GEA__snp_clumping_distance="auto_genome_wide")

    def test_a_declared_distance_is_not_checked(self):
        # d150: a declared 150000 bp beside an explicit, unused rule on an unreliable curve.
        self.targets("cc_ld_fixed", "gea", GEA__snp_clumping_distance=150000,
                     GEA__clumping_ld_rule="critical", GEA__ld_decay_group="C2")

    def test_gwas_inherits_and_is_checked(self):
        with self.assertRaisesRegex(ValueError, r"GWAS.snp_clumping_distance"):
            self.targets("cc_ld_gwas", "gwas", GEA__snp_clumping_distance="auto_genome_wide",
                         GEA__ld_decay_group="C2")

    def test_a_declared_snp_set_with_an_auto_gap_is_checked(self):
        sets = [{"name": "ovl", "strategy": "Cross-method", "clumping_distance": "auto_genome_wide"}]
        with self.assertRaisesRegex(ValueError, r"GEA.snp_sets gea.ovl clumping_distance"):
            self.targets("cc_ld_set", "gea", GEA__snp_clumping_distance=1000000,
                         GEA__ld_decay_group="C2", GEA__snp_sets=sets)

    def test_skipped_before_mode_structure(self):
        ns = load(project_name="cc_ld_fresh", GEA__snp_clumping_distance="auto_genome_wide",
                  GEA__ld_decay_group="C2")
        H.targets(ns, "gea")

    def test_a_table_without_ld_status_is_judged_on_its_own_curve(self):    # 96e4fd
        # It used to read as reliable; C2's 0.32 / 0.23 is refused like a written status.
        with self.assertRaisesRegex(ValueError, r"unreliable: short-range r2 0.32 is 1.39x"):
            self.targets("cc_ld_legacy", "gea", table="no_status",
                         GEA__snp_clumping_distance="auto_genome_wide", GEA__ld_decay_group="C2")
        self.targets("cc_ld_legacy_ok", "gea", table="no_status",
                     GEA__snp_clumping_distance="auto_genome_wide")


class TestSaturatedSpatialGF(unittest.TestCase):    # 895bac
    def targets(self, project, n_sites, **overrides):
        ns = load(project_name=project, **overrides)
        if n_sites is not None:
            rows = [("site", "sample", "latitude", "longitude")]
            rows += [(f"S{i % n_sites}", f"s{i}", 30 + i, 34) for i in range(2 * n_sites)]
            write(ns["W"]["metadata_climate_valid"], tsv(*rows))
        H.save_snp_set(project, "setA")
        return H.targets(ns, "maladaptation")

    def test_predictors_at_sites_minus_one_are_refused(self):
        for spatial in ("both", "with"):
            with self.subTest(spatial=spatial):
                with self.assertRaisesRegex(ValueError, r"3 Climate.predictors on 4 sites"):
                    self.targets(f"cc_gf_{spatial}", 4,
                                 Maladaptation__methods__gradient_forest__spatial_correction=spatial)

    def test_without_spatial_correction_passes(self):
        self.targets("cc_gf_without", 4,
                     Maladaptation__methods__gradient_forest__spatial_correction="without")

    def test_enough_sites_pass(self):
        self.targets("cc_gf_five", 5)

    def test_metadata_counts_sites_before_mode_structure(self):
        # testdata's metadata has 9 sites: 8 predictors leave no residual df.
        many = ",".join(f"bio_{i}" for i in range(1, 9))
        with self.assertRaisesRegex(ValueError, r"8 Climate.predictors on 9 sites"):
            self.targets("cc_gf_meta", None, Climate__predictors=many)
        self.targets("cc_gf_meta_ok", None)


class TestInvariantTableWiring(unittest.TestCase):  # 932a8d, e83e54
    def rules(self, ns):
        return {(r["name"], r["smk"]): r for r in rule_graph(ns)}

    def test_the_variance_check_names_the_predictor_universe(self):
        for overrides, want in (({}, ",".join(f"bio_{i}" for i in range(1, 20))),
                                ({"Climate__source": "custom",
                                  "Climate__custom": {"present_table": "env.tsv",
                                                      "future_table": "fut.tsv",
                                                      "columns": "bio_2,bio_7"}}, "bio_2,bio_7")):
            with self.subTest(overrides=overrides):
                rule = self.rules(load(**overrides))[("check_climate_variance", "climate.smk")]
                self.assertEqual(rule["params"]["columns"], want)

    def test_the_climate_summary_counts_the_selected_predictors(self):
        ns = load(Climate__predictors="bio_1,bio_12")
        summary = next(r for r in rule_graph(ns) if r["name"] == "write_summary"
                       and "climate {params.summary_tsv}" in r["shell"])
        self.assertEqual(summary["params"]["predictors"], "bio_1,bio_12")
        self.assertIn("{params.predictors}", summary["shell"])


# ------------------------------------------------------------------ Phase 2b: input formats
def meta_config(name, content, **overrides):
    """A config whose Input.metadata is `content` (TestMetadata.config, plus overrides)."""
    folder = f"cc_meta2b_{name}/"
    write(f"/pipeline/{folder}x.vcf", "##fileformat=VCFv4.2\n")
    write(f"/pipeline/{folder}meta.tsv", content)
    return H.make_config(Input={"dir": folder, "vcf": "x.vcf", "metadata": "meta.tsv",
                                "gff": "testdata.gff3"}, **overrides)


class TestMetadataColumnNames(unittest.TestCase):   # 1d4be2
    def test_columns_one_to_four_are_checked_by_name(self):
        # gwasonly is the finding's case: GWAS ran on t3 alone, t1 and t2 sliced off as
        # "latitude"/"longitude".
        for name, header in (("gwasonly", "site\tsample\tt1\tt2\tt3"),
                             ("latlon", "site\tsample\tlat\tlon\theight"),
                             ("order", "sample\tsite\tlatitude\tlongitude\theight")):
            with self.subTest(header=header):
                with self.assertRaisesRegex(
                        ValueError, "columns 1-4 must be named site, sample, latitude, longitude"):
                    H.load(meta_config(name, header + "\nA\ts1\t1\t2\t3\n"))

    def test_the_canonical_header_passes_and_traits_follow_it(self):
        ns = H.load(meta_config("canon", "site\tsample\tlatitude\tlongitude\theight\n"
                                         "A\ts1\t30\t34\t1.2\n"))
        self.assertEqual(ns["PHENO_TRAITS"], ["height"])

    def test_empty_coordinates_pass(self):
        # A gwas_only project keeps the two columns, empty (filter_coord_samples.R reads them by name).
        ns = H.load(meta_config("nocoord", "site\tsample\tlatitude\tlongitude\theight\n"
                                           "A\ts1\t\t\t1.2\n"))
        self.assertEqual(ns["PHENO_TRAITS"], ["height"])


class TestGwasTraitNames(unittest.TestCase):        # 634e7c
    HEADER = "site\tsample\tlatitude\tlongitude\t100_seed_weight\tplant.height\tyield-2020\theight\n"
    ROW = "A\ts1\t30\t34\t1\t2\t3\t4\n"

    def test_names_the_wildcard_cannot_take_are_refused_at_parse_time(self):
        cfg = meta_config("badtraits", self.HEADER + self.ROW, GWAS__traits="")
        with self.assertRaisesRegex(ValueError,
                                    "'100_seed_weight', 'plant.height', 'yield-2020'"):
            H.load(cfg)

    def test_traits_left_out_by_gwas_traits_are_not_checked(self):
        ns = H.load(meta_config("filtered", self.HEADER + self.ROW, GWAS__traits="height"))
        self.assertEqual(ns["PHENO_TRAITS"], ["height"])

    def test_not_checked_without_gwas(self):
        ns = H.load(meta_config("nogwas", self.HEADER + self.ROW, GWAS__configs=[]))
        self.assertEqual(ns["PHENO_TRAITS"], [])

    def test_the_check_and_every_wildcard_share_one_pattern(self):
        # A constraint spelled out again in gwas.smk could drift from the parse-time check.
        with open(os.path.join(H.REPO_ROOT, "workflow", "rules", "gwas.smk")) as handle:
            text = handle.read()
        self.assertNotIn('pheno_trait = r"', text)
        self.assertGreater(text.count("pheno_trait = PHENO_TRAIT_RE"), 0)
        ns = load()
        self.assertEqual(ns["ASSOC_SOURCES"]["GWAS"]["trait_regex"], ns["PHENO_TRAIT_RE"])


class TestPloidy(unittest.TestCase):                 # 899581
    def test_another_ploidy_is_refused(self):
        with self.assertRaisesRegex(ValueError, "diploid only"):
            load(sNMF__ploidy=4)
        with self.assertRaisesRegex(ValueError, "must be an integer"):
            load(sNMF__ploidy=2.5)

    def test_absent_or_two_passes(self):
        for overrides in ({}, {"sNMF__ploidy": 2}, {"sNMF__ploidy": "2.0"}):
            with self.subTest(overrides=overrides):
                self.assertEqual(load(**overrides)["PLOIDY"], 2)


if __name__ == "__main__":
    unittest.main()
