"""scripts/design_adequacy.py — the two pure functions, the two I/O ones, and main().

Why this file is worth more than its size: `jacobi_eigenvalues()` is a hand-rolled
cyclic-Jacobi eigensolver with no library behind it, and the numbers it produces
(numerical_rank, effective_dimensionality, condition_number) are exactly the ones a
user quotes when deciding whether a design can support a GEA at all.
"""
import contextlib
import io
import os
import tempfile
import unittest

import _support

da = _support.load_script("design_adequacy")


def approx(got, want, tol=1e-9):
    return all(abs(g - w) <= tol for g, w in zip(got, want)) and len(got) == len(want)


class TestJacobiEigenvalues(unittest.TestCase):
    """Known spectra. Each matrix here has eigenvalues obtainable by hand."""

    def test_identity(self):
        self.assertTrue(approx(da.jacobi_eigenvalues([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                               [1.0, 1.0, 1.0]))

    def test_diagonal_returns_descending(self):
        # Also pins the documented "Descending" contract: input order is 3,1,2.
        self.assertTrue(approx(da.jacobi_eigenvalues([[3, 0, 0], [0, 1, 0], [0, 0, 2]]),
                               [3.0, 2.0, 1.0]))

    def test_two_by_two(self):
        self.assertTrue(approx(da.jacobi_eigenvalues([[2, 1], [1, 2]]), [3.0, 1.0]))

    def test_repeated_eigenvalue(self):
        # [[4,1,1],[1,4,1],[1,1,4]] = 3I + J -> eigenvalues 6, 3, 3.
        self.assertTrue(approx(da.jacobi_eigenvalues([[4, 1, 1], [1, 4, 1], [1, 1, 4]]),
                               [6.0, 3.0, 3.0]))

    def test_one_by_one(self):
        # n == 1: the `for p in range(n - 1)` rotation loop never executes.
        self.assertTrue(approx(da.jacobi_eigenvalues([[5.0]]), [5.0]))

    def test_trace_is_preserved(self):
        m = [[4.0, 1.0, 0.5], [1.0, 3.0, -0.2], [0.5, -0.2, 2.0]]
        ev = da.jacobi_eigenvalues(m)
        self.assertAlmostEqual(sum(ev), sum(m[i][i] for i in range(3)), places=9)

    def test_does_not_mutate_its_argument(self):
        m = [[2.0, 1.0], [1.0, 2.0]]
        before = [row[:] for row in m]
        da.jacobi_eigenvalues(m)
        self.assertEqual(m, before)

    def test_singular_matrix_smallest_eigenvalue_may_be_negative(self):
        """Do NOT assert strict PSD here.

        design_adequacy.py:336-339 documents that on a singular correlation matrix
        lambda_min lands at the noise floor and CAN come back slightly negative —
        which is the whole reason the caller guards `cond` with `lam_min > 1e-12`
        instead of comparing raw. Assert the noise floor, not non-negativity.
        """
        # col3 == col1 + col2 exactly -> the correlation matrix is singular.
        c1 = [1.0, 2.0, 3.0, 4.0, 5.0, 7.0]
        c2 = [2.0, 1.0, 4.0, 3.0, 8.0, 5.0]
        c3 = [a + b for a, b in zip(c1, c2)]
        R, _ = da.corr_matrix([c1, c2, c3])
        ev = da.jacobi_eigenvalues(R)
        self.assertGreaterEqual(ev[-1], -1e-9)
        self.assertLess(ev[-1], 1e-9)
        # ... and that this is exactly the case the caller routes to inf.
        cond = (ev[0] / ev[-1]) if ev[-1] > 1e-12 else float("inf")
        self.assertEqual(cond, float("inf"))

    def test_non_convergence_is_silent(self):
        """PINNED DEFECT, not desired behaviour.

        With max_sweeps exhausted the function returns whatever it has, with no
        warning, no exception and no flag — a caller cannot tell a converged
        result from an abandoned one. Filed in docs/pipeline_improvement_requests.md.
        """
        m = [[2.0, 1.0], [1.0, 2.0]]
        got = da.jacobi_eigenvalues(m, max_sweeps=0)
        self.assertTrue(approx(got, [2.0, 2.0]))   # the untouched diagonal


class TestCorrMatrix(unittest.TestCase):

    def test_returns_matrix_and_sd(self):
        out = da.corr_matrix([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0]])
        self.assertIsInstance(out, tuple)
        self.assertEqual(len(out), 2)

    def test_perfect_and_anti_correlation(self):
        R, _ = da.corr_matrix([[1.0, 2.0, 3.0, 4.0],
                               [2.0, 4.0, 6.0, 8.0],
                               [4.0, 3.0, 2.0, 1.0]])
        self.assertAlmostEqual(R[0][1], 1.0, places=12)
        self.assertAlmostEqual(R[0][2], -1.0, places=12)

    def test_diagonal_is_one_and_matrix_is_symmetric(self):
        R, _ = da.corr_matrix([[1.0, 5.0, 2.0, 9.0],
                               [3.0, 1.0, 4.0, 1.0],
                               [2.0, 7.0, 1.0, 8.0]])
        for i in range(3):
            self.assertAlmostEqual(R[i][i], 1.0, places=12)
            for j in range(3):
                self.assertAlmostEqual(R[i][j], R[j][i], places=12)

    def test_sd_uses_n_minus_one_denominator(self):
        # var = ((1-3)^2 + (3-3)^2 + (5-3)^2) / 2 = 4 -> sd 2
        _, sd = da.corr_matrix([[1.0, 3.0, 5.0]])
        self.assertAlmostEqual(sd[0], 2.0, places=12)

    def test_zero_variance_column_gets_a_synthetic_identity_row(self):
        """Deliberate branch (design_adequacy.py:101-103), and main() depends on it.

        An invariant predictor gets 1.0 on its own diagonal and 0.0 everywhere
        else — which is why main() screens invariant predictors out BEFORE the
        spectrum (:264-274): left in, each would add an eigenvalue of exactly 1
        to numerical_rank and effective_dimensionality.
        """
        R, sd = da.corr_matrix([[1.0, 2.0, 3.0], [7.0, 7.0, 7.0]])
        self.assertEqual(sd[1], 0.0)
        self.assertEqual(R[1][1], 1.0)
        self.assertEqual(R[0][1], 0.0)
        self.assertEqual(R[1][0], 0.0)

    def test_single_row_takes_the_sd_zero_path(self):
        _, sd = da.corr_matrix([[4.0], [9.0]])
        self.assertEqual(sd, [0.0, 0.0])


class TestReadTsv(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, text):
        p = os.path.join(self.dir.name, "t.tsv")
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def test_splits_header_and_rows(self):
        hdr, rows = da.read_tsv(self.write("a\tb\n1\t2\n3\t4\n"))
        self.assertEqual(hdr, ["a", "b"])
        self.assertEqual(rows, [["1", "2"], ["3", "4"]])

    def test_blank_lines_are_skipped(self):
        hdr, rows = da.read_tsv(self.write("a\tb\n1\t2\n\n\n3\t4\n\n"))
        self.assertEqual(len(rows), 2)

    def test_short_row_is_not_padded(self):
        """Pinned: a malformed TSV yields a short list, not a header-width one."""
        _, rows = da.read_tsv(self.write("a\tb\tc\n1\t2\n"))
        self.assertEqual(rows[0], ["1", "2"])

    def test_header_only_file_yields_no_rows(self):
        hdr, rows = da.read_tsv(self.write("a\tb\n"))
        self.assertEqual(hdr, ["a", "b"])
        self.assertEqual(rows, [])


class TestMain(unittest.TestCase):
    """main() writes its flags and its usage text to stderr; swallow it so a green
    run stays readable. Failures still surface through the assertions."""

    @staticmethod
    def run_main(argv):
        with contextlib.redirect_stderr(io.StringIO()):
            return da.main(argv)


    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.d = self.dir.name
        self.out = os.path.join(self.d, "design_adequacy.tsv")

    def test_wrong_argument_count_returns_2(self):
        self.assertEqual(self.run_main(["design_adequacy.py", "a", "b"]), 2)
        self.assertEqual(self.run_main(["design_adequacy.py", "a", "b", "c", "d", "e"]), 2)

    def test_metadata_without_site_column_returns_1(self):
        meta = os.path.join(self.d, "m.tsv")
        clim = os.path.join(self.d, "c.tsv")
        with open(meta, "w") as fh:
            fh.write("sample\tlatitude\nID1\t30\n")
        with open(clim, "w") as fh:
            fh.write("sample\tbio_1\nID1\t1.0\n")
        self.assertEqual(self.run_main(["x", meta, clim, "all", self.out]), 1)

    def test_unknown_predictor_returns_1(self):
        meta, clim = _support.make_fixture(self.d, n_sites=5, n_pred=3, per_site=3,
                                           n_latent=2, noise=0.5, seed=11)
        self.assertEqual(self.run_main(["x", meta, clim, "bio_1,bio_99", self.out]), 1)

    def test_adequate_design_is_full_rank_and_unflagged(self):
        meta, clim = _support.make_fixture(self.d, n_sites=60, n_pred=4, per_site=5,
                                           n_latent=4, noise=1.6, seed=2)
        self.assertEqual(self.run_main(["x", meta, clim, "all", self.out]), 0)
        t = _support.read_table(self.out)
        self.assertEqual(t["n_sites"][0], "60")
        self.assertEqual(t["n_samples"][0], "300")
        self.assertEqual(t["n_predictors"][0], "4")
        self.assertEqual(t["numerical_rank"][0], "4")
        self.assertEqual(t["numerical_rank_ceiling"][0], "4")
        self.assertEqual(t["n_invariant_predictors"][0], "0")
        self.assertEqual(t["n_flags_fail"][0], "0")

    def test_saturated_design_fails_on_residual_df(self):
        """18 predictors over 11 sites: rank-deficient BY CONSTRUCTION."""
        meta, clim = _support.make_fixture(self.d, n_sites=11, n_pred=18, per_site=6,
                                           n_latent=3, noise=0.05, seed=1,
                                           singleton_sites=2)
        self.assertEqual(self.run_main(["x", meta, clim, "all", self.out]), 0)
        t = _support.read_table(self.out)
        self.assertEqual(t["n_sites"][0], "11")
        self.assertEqual(t["site_level_residual_df"][0], str(11 - 18 - 1))
        self.assertEqual(t["site_level_residual_df"][1], "FAIL")
        # rank is capped by the design (n_sites - 1 = 10), not by the data.
        self.assertEqual(t["numerical_rank_ceiling"][0], "10")
        self.assertEqual(t["numerical_rank"][1], "WARN")

    def test_messy_fixture_accounts_for_every_dropped_row(self):
        """The counting contract, on the one fixture that exercises all four drops.

        6 sites x 4 samples in the metadata = 24. The climate table omits 2 rows
        (S05), carries an ORPHAN with no metadata entry, blanks bio_2 on all 4
        S03 rows and on 1 S00 row, and has an invariant bio_12.
        """
        meta, clim = _support.make_messy_fixture(self.d)
        self.assertEqual(self.run_main(["x", meta, clim, "all", self.out]), 0)
        t = _support.read_table(self.out)
        self.assertEqual(t["n_samples_in_metadata"][0], "24")
        self.assertEqual(t["n_sites_in_metadata"][0], "6")
        # S03 loses every row -> it is not a site any more.
        self.assertEqual(t["n_sites"][0], "5")
        # 23 climate rows - 1 orphan - 4 (S03) - 1 (S00) = 17 that carried values.
        self.assertEqual(t["n_samples"][0], "17")
        self.assertEqual(t["samples_missing_from_climate"][0], "7")
        self.assertEqual(t["n_invariant_predictors"][0], "1")
        self.assertIn("bio_12", t["n_invariant_predictors"][2])
        self.assertEqual(t["n_predictors_in_spectrum"][0], "3")

    def test_output_is_a_four_column_tsv(self):
        meta, clim = _support.make_fixture(self.d, n_sites=8, n_pred=3, per_site=4,
                                           n_latent=3, noise=0.8, seed=5)
        self.assertEqual(self.run_main(["x", meta, clim, "all", self.out]), 0)
        with open(self.out) as fh:
            lines = fh.read().rstrip("\n").split("\n")
        self.assertEqual(lines[0], "metric\tvalue\tflag\tnote")
        for ln in lines[1:]:
            self.assertEqual(len(ln.split("\t")), 4)


if __name__ == "__main__":
    unittest.main()
