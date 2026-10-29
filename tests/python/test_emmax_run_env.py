"""scripts/emmax_run.sh runs EMMAX in a fixed environment (Phase 4a, 2026-10-29).

The bundled emmax-intel64 links a static Intel MKL whose code path follows memory alignment,
and the alignment followed the size of the process environment: one more variable in the
container moved SIMDATA's GWAS EMMAX p-values by up to 7e-5 relative, with byte-identical
inputs. The output also changed with the thread count, and MKL_CBWR alone did not remove the
environment dependence (measured). The wrapper therefore execs EMMAX under `env -i` with
exactly MKL_CBWR=COMPATIBLE, OMP_NUM_THREADS=1 and MKL_NUM_THREADS=1.

The behavioural test hands the wrapper /usr/bin/printenv as its "binary" and asks for named
variables only, so it prints their values and never the whole environment. The wrapper's
preflight refuses anything but a native x86-64 host, so that test is skipped elsewhere.
"""
import os
import platform
import re
import subprocess
import unittest

from _support import REPO_ROOT

WRAPPER = os.path.join(REPO_ROOT, "scripts", "emmax_run.sh")
FIXED = {"MKL_CBWR": "COMPATIBLE", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}


def _native_x86_64():
    if platform.machine() != "x86_64":
        return False
    try:
        with open("/proc/cpuinfo") as handle:
            return any(re.match(r"(flags|vendor_id)\s*:", line) for line in handle)
    except OSError:
        return False


class FixedEnvironmentTest(unittest.TestCase):
    def test_the_wrapper_execs_emmax_under_env_i_with_the_three_settings(self):
        with open(WRAPPER) as handle:
            text = handle.read()
        self.assertRegex(text, r'(?m)^exec env -i \$FIXED_ENV "\$BIN" "\$@"$')
        fixed = re.search(r'(?m)^FIXED_ENV="([^"]*)"$', text).group(1).split()
        for name, value in FIXED.items():
            self.assertIn(f"{name}={value}", fixed)
        # $BIN is checked first, so the exec can never run a bare `env`.
        self.assertLess(text.index('[ -x "$BIN" ]'), text.index("exec env -i"))

    @unittest.skipUnless(_native_x86_64(), "emmax_run.sh refuses a non-native x86-64 host")
    def test_the_callers_environment_never_reaches_the_binary(self):
        env = dict(os.environ, CLINEGO_TEST_SENTINEL="leak", MKL_CBWR="AVX2",
                   MKL_NUM_THREADS="8", OMP_NUM_THREADS="8")
        res = subprocess.run([WRAPPER, "/usr/bin/printenv", "MKL_CBWR", "OMP_NUM_THREADS",
                              "MKL_NUM_THREADS", "CLINEGO_TEST_SENTINEL"],
                             env=env, capture_output=True, text=True)
        # printenv exits 1 when a named variable is unset: the sentinel must be.
        self.assertEqual(res.returncode, 1, res.stderr)
        self.assertEqual(res.stdout.splitlines(), ["COMPATIBLE", "1", "1"])


if __name__ == "__main__":
    unittest.main()
