"""Registry hyperparameters must reach the rules that run the method.

WHY THIS IS A SEPARATE TEST
---------------------------
`workflow/methods/gea.py` declares each method's hyperparameters (`params`), and
`common.smk`'s resolve_method_params() defaults, coerces, range-checks and REJECTS
unknown names for BOTH GEA.configs and GWAS.configs. Everything up to the point of
use is therefore validated — and nothing checks the point of use. A param can be
accepted, bounds-checked, written into the project config by the Shiny sidebar,
shown with a PreGEA "recommended" badge, and then never reach the script: the run
silently uses the registry default instead of the value the user chose.

That is what happened to every GWAS hyperparameter. GWAS_PARAMS was parsed at
common.smk and exported to the app's config bundle, but no rule in gwas.smk read
it: EMMAX ran with n_pcs = sNMF k_best whatever the config said, and both GAPIT
call sites stopped at argument 14, so gapit.R's `N_PCS <- if (length(args) >= 15)
... else K_BEST` fallback took over. mode=gea honoured the same param on the same
registry entry, so a GEA scan and a GWAS scan of one project were corrected for
structure differently.

A positional-arity check CANNOT see this: arg 15 is optional, so a 14-token
call site is a legal invocation of gapit.R. The contract broken here is not
positional arity but "a declared param has a consumer", which is why it lives in
its own file.

CONSUMED means a `.get("<param>")` lookup rooted at the module's own params dict
(GEA_PARAMS / GWAS_PARAMS) exists somewhere under workflow/ — either in a rule
file directly, or in the body of a helper that a rule file calls WITH that dict
(common.smk's emmax_kinship_climate_path() and gapit_shared_npcs() are both that
shape). It is not "the param name appears in the file": `kinship` occurs a dozen
times in gwas.smk as a rule name and an input key while nothing reads the param.
REFUSED is the escape hatch, and it demands that the param fail LOUDLY at parse
time rather than being quietly dropped.
"""
import os
import re
import unittest

import _support
from _support import REPO_ROOT

gea = _support.import_methods("gea")
gwas = _support.import_methods("gwas")

RULES = os.path.join(REPO_ROOT, "workflow", "rules")

# module label -> (registry, rule file that runs the method)
MODULES = {
    "GEA":  (gea.GEA_METHODS,  "gea.smk"),
    "GWAS": (gwas.GWAS_METHODS, "gwas.smk"),
}

# (module, method, param) -> why the pipeline cannot honour it. Every entry MUST
# be rejected at config-parse time; test_refused_params_are_rejected_loudly
# asserts common.smk raises with the param named in the message.
REFUSED = {
    ("GWAS", "EMMAX", "kinship"):
        "Every GWAS kinship rule declares a .aBN.kinf output and calls "
        "emmax-kin without -s, so the IBS estimator has no rule on the GWAS "
        "side. GEA has one (kinship_gea_climate_ibs). Accepting 'IBS' and "
        "running BN would be silent wrong science, so common.smk raises.",
}


def read(path):
    with open(os.path.join(RULES, path), errors="replace") as fh:
        return fh.read()


SYMBOL = {"GEA": "GEA_PARAMS", "GWAS": "GWAS_PARAMS"}

# GEA_PARAMS.get("EMMAX", {}).get("n_pcs", K_BEST)  /  ...get(_method, {}).get("n_pcs")
_DIRECT = re.compile(r"\b(GEA_PARAMS|GWAS_PARAMS)\b\s*\.get\([^)]*\)\s*"
                     r"\.get\(\s*['\"](\w+)['\"]")
# gapit_shared_npcs(GEA_GAPIT_CONFIGS, GEA_PARAMS, 'GEA')  — helper called with the dict
_INDIRECT = re.compile(r"\b(\w+)\s*\(([^()]*\b(?:GEA_PARAMS|GWAS_PARAMS)\b[^()]*)\)")
_KEYS = re.compile(r"\.get\(\s*['\"](\w+)['\"]")
# RDA reads its params by subscript off a local alias: _rda_p = GEA_PARAMS["RDA"]
# then _rda_p['axes'] — no .get() anywhere.
_SUBSCRIPTS = re.compile(r"\[\s*['\"](\w+)['\"]\s*\]")
_ALIAS = re.compile(r"\b(\w+)\s*=\s*(GEA_PARAMS|GWAS_PARAMS)\s*\[")


def _smk_sources():
    for fn in sorted(os.listdir(RULES)):
        if fn.endswith(".smk"):
            yield fn, read(fn)


def _helper_body(name, text):
    """Source of `def name(...)`: its header plus the indented block under it.

    Bounded by the first subsequent line that starts in column 0 — NOT by the
    next `def`. common.smk puts module-level guard code between functions, and a
    def-to-def span swallowed 400 lines of unrelated config parsing, so every
    `.get("...")` key in the file looked like a consumed hyperparameter.
    """
    m = re.search(r"^def\s+" + re.escape(name) + r"\s*\(", text, re.M)
    if not m:
        return ""
    lines = text[m.start():].splitlines()
    body = [lines[0]]
    for line in lines[1:]:
        if line.strip() and not line[0].isspace():
            break
        body.append(line)
    return "\n".join(body)


def _strip_refusal_guards(text):
    """Blank out `if <PARAMS>...:` blocks whose body raises.

    Reading a param in order to REFUSE it is not consuming it — without this the
    parse-time guard for GWAS EMMAX kinship='IBS' would itself make the param
    look wired, and the refusal would silently become the evidence that no
    refusal is needed.
    """
    lines = text.splitlines()
    out = list(lines)
    for i, line in enumerate(lines):
        m = re.match(r"^(\s*)if\b.*\b(?:GEA_PARAMS|GWAS_PARAMS)\b.*:\s*$", line)
        if not m:
            continue
        indent = len(m.group(1))
        block = [i]
        for j in range(i + 1, len(lines)):
            if lines[j].strip() and (len(lines[j]) - len(lines[j].lstrip())) <= indent:
                break
            block.append(j)
        if any("raise" in lines[j] for j in block):
            for j in block:
                out[j] = ""
    return "\n".join(out)


def consumed_params():
    """{module: {param names a rule actually reads out of that module's dict}}.

    Restricted to names the registry declares: a helper body may read unrelated
    keys, and only declared names can answer the question this file asks.
    """
    declared = {p for _m, _meth, p in declared_params()}
    out = {mod: set() for mod in MODULES}
    inv = {v: k for k, v in SYMBOL.items()}
    sources = dict(_smk_sources())
    common = sources.get("common.smk", "")
    for _fn, raw in sources.items():
        text = _strip_refusal_guards(raw)
        for sym, param in _DIRECT.findall(text):
            out[inv[sym]].add(param)
        # `_rda_p = GEA_PARAMS["RDA"]` … `_rda_p['axes']`
        for alias, sym in _ALIAS.findall(text):
            alias_lines = "\n".join(l for l in text.splitlines()
                                    if re.search(r"\b" + re.escape(alias) + r"\s*\[", l))
            out[inv[sym]].update(_SUBSCRIPTS.findall(alias_lines))
        for helper, argtext in _INDIRECT.findall(raw):
            for sym in inv:
                if sym not in argtext:
                    continue
                body = _helper_body(helper, common) or _helper_body(helper, text)
                out[inv[sym]].update(_KEYS.findall(body) + _SUBSCRIPTS.findall(body))
    for mod in out:
        out[mod] &= declared
    return out


def declared_params():
    """[(module, method, param), ...] over both registries."""
    out = []
    for module, (registry, _rule_file) in MODULES.items():
        for method, cfg in registry.items():
            for param in (cfg.get("params") or {}):
                out.append((module, method, param))
    return out


def _raise_messages(text):
    """Every `raise ...(...)` argument body in a .smk/py source, paren-matched."""
    msgs = []
    for m in re.finditer(r"\braise\s+\w+\(", text):
        depth, i = 0, m.end() - 1
        while i < len(text):
            if text[i] == "(":
                depth += 1
            elif text[i] == ")":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        msgs.append(text[m.end():i])
    return msgs


class TestEveryDeclaredParamHasAConsumer(unittest.TestCase):
    def test_declared_params_are_consumed_or_refused(self):
        consumed = consumed_params()
        orphans = []
        for module, method, param in declared_params():
            if (module, method, param) in REFUSED:
                continue
            if param not in consumed[module]:
                orphans.append(f"{module}.configs {method} params.{param}")
        self.assertEqual(
            sorted(orphans), [],
            "declared hyperparameter(s) that no rule reads — the run silently "
            "uses the registry default instead of the configured value. Wire "
            "them into the rule, or add a REFUSED row plus a parse-time raise: "
            f"{sorted(orphans)}")

    def test_no_stale_refusal(self):
        """A param that got wired must lose its REFUSED row."""
        consumed = consumed_params()
        stale = [k for k in REFUSED if k[2] in consumed[k[0]]]
        self.assertEqual(stale, [], f"now referenced in the rule file: {stale}")

    def test_refused_params_are_rejected_loudly(self):
        common = read("common.smk")
        msgs = _raise_messages(common)
        for module, method, param in REFUSED:
            hit = [m for m in msgs
                   if param in m and f"{module}.configs" in m]
            self.assertTrue(
                hit,
                f"{module}.configs {method} params.{param} is listed as REFUSED "
                f"but common.smk never raises with '{module}.configs' and "
                f"'{param}' in the message — it is being dropped silently, "
                "which is the defect this file exists to prevent.")

    def test_the_consumption_scan_resolves_both_shapes(self):
        """Direct (GEA_PARAMS.get(...).get("n_pcs")) and indirect (a helper
        called with the dict) must both resolve, or every param looks orphaned
        and the test degrades to noise."""
        consumed = consumed_params()
        self.assertIn("n_pcs", consumed["GEA"])    # direct, gea.smk EMMAX rule
        self.assertIn("kinship", consumed["GEA"])  # indirect, emmax_kinship_climate_path()
        self.assertIn("n_pcs", consumed["GWAS"])
        self.assertNotIn("kinship", consumed["GWAS"])

    def test_the_scan_actually_sees_the_registry(self):
        """Guard against the registries importing empty."""
        found = declared_params()
        self.assertIn(("GEA", "EMMAX", "n_pcs"), found)
        self.assertIn(("GWAS", "BLINK", "n_pcs"), found)
        self.assertIn(("GEA", "LFMM", "K"), found)
        self.assertGreater(len(found), 15)


class TestGapitCallSitesAgreeOnArity(unittest.TestCase):
    """gapit.R's n_pcs is positional argument 15 and OPTIONAL, so a call site
    that omits it exits 0 and writes every declared output while silently using
    K_BEST as the PCA covariate count. Both GEA and GWAS run the same script
    through the same registry entry, so the two must pass the same arity.

    Self-contained on purpose: every gapit.R call site is a plain
    `Rscript /pipeline/scripts/gapit.R a b c ... > {log} 2>&1` block, so joining
    backslash continuations and splitting on whitespace up to the redirect is an
    exact tokenisation for this one script."""

    _RULE = re.compile(r"^\s*rule\s+(\w+)\s*:", re.M)
    _CALL = re.compile(r"/pipeline/scripts/gapit\.R\s+([^>]*)>")

    def _gapit_calls(self):
        out = {}
        for fn in sorted(os.listdir(RULES)):
            if not fn.endswith(".smk"):
                continue
            text = re.sub(r"\\\n", " ", read(os.path.join(RULES, fn)))
            rules = [(m.start(), m.group(1)) for m in self._RULE.finditer(text)]
            for m in self._CALL.finditer(text):
                name = [n for pos, n in rules if pos < m.start()][-1]
                out[f"{fn}:{name}"] = m.group(1).split()
        return out

    def test_every_gapit_call_site_passes_n_pcs(self):
        calls = self._gapit_calls()
        self.assertGreaterEqual(len(calls), 3, f"call sites not found: {sorted(calls)}")
        short = {k: len(v) for k, v in calls.items() if len(v) < 15}
        self.assertEqual(
            short, {},
            "gapit.R call site(s) stop before argument 15 (n_pcs), so GAPIT runs "
            f"with N_PCS = K_BEST no matter what the config says: {short}")
        wrong = {k: v[14] for k, v in calls.items() if v[14] != "{params.n_pcs}"}
        self.assertEqual(wrong, {}, f"argument 15 is not {{params.n_pcs}}: {wrong}")


if __name__ == "__main__":
    unittest.main()
