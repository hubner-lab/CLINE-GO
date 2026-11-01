"""project_root: where a run reads its inputs and writes its results (release v1.0).

common.smk builds INDIR, OUTDIR and LOGDIR from the optional config key
`project_root`. Pinned here:
  * absent, the three strings are byte-identical to the pre-v1.0 form
    ('/pipeline/' + Input.dir, '/pipeline/<p>_results/', '/pipeline/<p>_logs/') for
    Input.dir 'data/', 'data' and 'test_data/', so no existing tree reruns;
  * an absolute existing directory puts all three under it;
  * a relative or a missing directory is refused at parse time, naming the key;
  * the baked-code guard: with BUILD_INFO_PATH on disk and no project_root the load
    stops; without that file it is quiet.
"""
import os
import shutil
import tempfile
import unittest

import _smk_harness as H

TEST_DATA = os.path.join(H.REPO_ROOT, "test_data")
INPUTS = ("testdata.vcf", "testdata_metadata.tsv", "testdata.gff3")


def _sandbox_inputs(test, input_dir):
    """Put the test inputs where the harness reads /pipeline/<input_dir>/, removed
    again after the test (the sandbox is shared by every test in the run)."""
    real = os.path.join(H._SANDBOX, input_dir)
    if not os.path.exists(real):
        test.addCleanup(shutil.rmtree, real, True)
    os.makedirs(real, exist_ok=True)
    for name in INPUTS:
        shutil.copy(os.path.join(TEST_DATA, name), os.path.join(real, name))


class DefaultRootIsTodaysForm(unittest.TestCase):

    def test_strings_unchanged_for_every_input_dir_spelling(self):
        for input_dir in ("data/", "data", "test_data/"):
            with self.subTest(input_dir=input_dir):
                if input_dir != "test_data/":
                    _sandbox_inputs(self, input_dir.rstrip("/"))
                ns = H.load(H.make_config(Input__dir=input_dir))
                self.assertEqual(ns["INDIR"], "/pipeline/" + input_dir)
                self.assertEqual(ns["OUTDIR"], "/pipeline/testdata_results/")
                self.assertEqual(ns["LOGDIR"], "/pipeline/testdata_logs/")

    def test_empty_string_reads_as_absent(self):
        ns = H.load(H.make_config(project_root=""))
        self.assertEqual(ns["OUTDIR"], "/pipeline/testdata_results/")


class AbsoluteRoot(unittest.TestCase):

    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="clinego_root_")
        shutil.copytree(TEST_DATA, os.path.join(self.root, "test_data"))

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def test_all_three_under_the_root(self):
        for spelled in (self.root, self.root + "/"):
            with self.subTest(project_root=spelled):
                ns = H.load(H.make_config(project_root=spelled))
                self.assertEqual(ns["INDIR"], self.root + "/test_data/")
                self.assertEqual(ns["OUTDIR"], self.root + "/testdata_results/")
                self.assertEqual(ns["LOGDIR"], self.root + "/testdata_logs/")
                self.assertTrue(os.path.isdir(ns["OUTDIR"]))

    def test_inputs_are_read_from_the_root_not_from_pipeline(self):
        os.remove(os.path.join(self.root, "test_data", "testdata.vcf"))
        with self.assertRaisesRegex(ValueError, "input.vcf"):
            H.load(H.make_config(project_root=self.root))


class RefusedRoots(unittest.TestCase):

    def test_relative_root(self):
        for rel in ("data", "./proj", "proj/"):
            with self.subTest(project_root=rel), \
                 self.assertRaisesRegex(ValueError, "project_root.*absolute"):
                H.load(H.make_config(project_root=rel))

    def test_missing_root(self):
        missing = os.path.join(tempfile.gettempdir(), "clinego_no_such_root_7f3a")
        self.assertFalse(os.path.exists(missing))
        with self.assertRaisesRegex(ValueError, "project_root.*not an existing directory"):
            H.load(H.make_config(project_root=missing))


class BakedCodeGuard(unittest.TestCase):
    """BUILD_INFO_PATH exists only in an image with the code baked in."""

    def setUp(self):
        self.ns = H.load()
        self.tmp = tempfile.mkdtemp(prefix="clinego_buildinfo_")
        self.info = os.path.join(self.tmp, "BUILD_INFO.json")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_constant_is_the_image_path(self):
        self.assertEqual(self.ns["BUILD_INFO_PATH"], "/pipeline/BUILD_INFO.json")

    def test_fires_with_build_info_and_no_root(self):
        with open(self.info, "w") as handle:
            handle.write('{"version": "1.0"}')
        with self.assertRaisesRegex(ValueError, "the code runs from the image: give "
                                    "--config project_root=<dir>, or use clinego run"):
            self.ns["resolve_project_root"]({}, build_info_path=self.info)

    def test_quiet_without_build_info(self):
        self.assertEqual(self.ns["resolve_project_root"]({}, build_info_path=self.info),
                         "/pipeline/")

    def test_root_given_passes_the_guard(self):
        with open(self.info, "w") as handle:
            handle.write("{}")
        self.assertEqual(
            self.ns["resolve_project_root"]({"project_root": self.tmp},
                                            build_info_path=self.info),
            self.tmp + "/")

    def test_a_load_reads_the_module_constant(self):
        """Through the harness, /pipeline/BUILD_INFO.json resolves into its sandbox."""
        sandboxed = os.path.join(H._SANDBOX, "BUILD_INFO.json")
        self.assertFalse(os.path.exists(os.path.join(H.REPO_ROOT, "BUILD_INFO.json")))
        with open(sandboxed, "w") as handle:
            handle.write("{}")
        try:
            with self.assertRaisesRegex(ValueError, "the code runs from the image"):
                H.load()
        finally:
            os.remove(sandboxed)
        H.load()


if __name__ == "__main__":
    unittest.main()
