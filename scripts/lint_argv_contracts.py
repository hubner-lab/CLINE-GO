#!/usr/bin/env python3
"""lint_argv_contracts.py — cross-check Snakemake call sites against the positional
argv contracts of the Rscript wrappers they invoke.

WHY THIS EXISTS
---------------
Every `scripts/*.R` wrapper binds its arguments BY POSITION
(`args <- commandArgs(trailingOnly = TRUE)`; `IN <- args[1]`; `OUT <- args[9]`),
and every call site is a free-form `shell:` string in `workflow/rules/*.smk`.
Nothing connects the two. The contract is 76 scripts wide, up to 27 positions
deep (`scripts/rda_offset.R`), and it is checked by nobody:

  * Add an argument to a script and forget one of its call sites -> every
    argument after the insertion point silently shifts. The script still exits 0
    and still writes its declared outputs, so Snakemake is satisfied and the
    tests pass; only the CONTENT is wrong. An output path arriving in an input
    slot means a script overwrites one of its own inputs.
  * A `params.*_str` that space-joins a file list expands to N shell tokens.
    Quoted, it is one argument; unquoted, it shifts everything after it by N-1.
    Quoting is the only thing standing between the pipeline and that shift, and
    it is a per-call-site convention with no enforcement.
  * A bare `{output}` in a rule with k outputs expands to k tokens, not one.

This linter is a STATIC check: it never runs a script. It reads the argv indices
each script actually dereferences, reads the token sequence each call site
passes, and reports where the two cannot agree.

CHECKS
------
  arity_mismatch         call site passes a token count the script cannot consume
                         (fewer than its required positions, or more than its
                         maximum, for a script with a fixed contract)
  unquoted_multi_token   a space-joined `params` value, or a bare `{output}` /
                         `{input}` in a multi-output/multi-input rule, is passed
                         unquoted into a fixed positional slot
  slot_role_mismatch     the call site passes {output} into a slot the script
                         binds to an input variable, or {input.x} into a slot it
                         binds to an OUT_* variable — the signature of a shift
                         that arity alone cannot see, because a trailing optional
                         argument absorbs the token count
  sparse_argv            a script reads args[n] but never args[n-1] — either a
                         dead position or an off-by-one
  no_arity_guard         a script with >= MIN_GUARDED_ARGS positions never
                         inspects `length(args)`, so a shifted call site cannot
                         be detected at run time

Exit status is 1 if any error-severity violation is found, 0 otherwise —
matching `scripts/check_invariants.R`'s validator convention.

USAGE
    python3 scripts/lint_argv_contracts.py [REPO_ROOT] [--tsv OUT.tsv] [--quiet]
"""

from __future__ import annotations

import os
import re
import sys

# A script with at least this many positional arguments is expected to assert
# something about length(args). Below it, a shift is usually caught by an
# immediate file-not-found.
MIN_GUARDED_ARGS = 6

# ---------------------------------------------------------------------------
# R side: what positions does a script dereference, and is any of it optional?
# ---------------------------------------------------------------------------

_ARG_IDX = re.compile(r"\bargs\s*\[\s*(\d+)\s*\]")
_ARGS_LEN = re.compile(r"length\s*\(\s*args\s*\)")
# `length(args) >= 7`, `length(args) > 6`, `n >= 7` after `n <- length(args)`
_OPT_GUARD = re.compile(r"length\s*\(\s*args\s*\)\s*(>=|>|==)\s*(\d+)")
# Variadic consumption: args[seq_len(n-1)], args[-1], args[n], args[2:n], ...
_VARIADIC = re.compile(
    r"args\s*\[\s*(?:seq_len|seq_along|-|[0-9]+\s*:\s*n\b|n\b|c\()", re.VERBOSE
)


_ASSIGN_ARG1 = re.compile(r"([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)\s*args\s*\[\s*1\s*\]")
# `OUT_FILE = args[8]`, `OUTPUT <- args[5]`,
# `OUT_DIAGNOSTICS <- if (length(args) >= 19 && nzchar(args[19])) args[19] else NA`
_ASSIGN_LHS = re.compile(r"^\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*(?:<-|=)(?!=)")
_OUT_NAME = re.compile(r"(^|_)(OUT|OUTPUT|OUTFILE)($|_)|^OUT", re.IGNORECASE)


def _is_mode_dispatch(text: str) -> bool:
    """True if args[1] selects between >= 3 alternative argv contracts."""
    for m in _ASSIGN_ARG1.finditer(text):
        var = re.escape(m.group(1))
        eqs = re.findall(rf"\b{var}\s*==\s*['\"]", text)
        if len(eqs) >= 3:
            return True
    return False


class ScriptContract:
    """The positional contract of one `scripts/*.R` wrapper."""

    def __init__(self, path: str, text: str):
        self.path = path
        idx = sorted({int(m.group(1)) for m in _ARG_IDX.finditer(text)})
        self.indices = idx
        self.max_index = idx[-1] if idx else 0
        self.missing = [i for i in range(1, self.max_index + 1) if i not in idx]
        self.has_len_check = bool(_ARGS_LEN.search(text))
        self.variadic = bool(_VARIADIC.search(text))
        # A position is OPTIONAL only when EVERY place the script dereferences it
        # also tests length(args) on that line — the `if (length(args) >= 14)
        # args[14] else "NULL"` idiom this repo uses throughout. Anything read
        # unconditionally is required, however high its index, so appending a new
        # mandatory argument above an existing optional one is still caught.
        guarded = {}
        for line in text.splitlines():
            hits = [int(m.group(1)) for m in _ARG_IDX.finditer(line)]
            if not hits:
                continue
            is_guard_line = bool(_ARGS_LEN.search(line))
            for k in hits:
                guarded[k] = guarded.get(k, True) and is_guard_line
        optional = {k for k, g in guarded.items() if g}
        req = [k for k in idx if k not in optional]
        self.optional = sorted(optional)
        self.required = max(req) if req else 0

        # MODE-DISPATCHING scripts (write_summary.R) hold one contract PER mode
        # selected by args[1], so no single arity applies and the guards of one
        # mode say nothing about another's. Arity is unenforceable for these;
        # only the upper bound (nothing passed is ignored) still holds.
        self.multi_contract = _is_mode_dispatch(text)
        if self.multi_contract:
            self.required = 1

        # Slot ROLE, read off the variable the script binds each position to.
        # Arity alone cannot see a one-token deletion that a trailing optional
        # argument absorbs (emmax.R's args[11] does exactly that), but a shift
        # always lands an output path in a slot the script treats as an input, or
        # the reverse — and that IS visible.
        # A slot counts as an output slot only when EVERY binding of it names an
        # output. A script with two mode branches can bind one position to
        # OUT_RASTER in one branch and PRESENT_ALL in the other
        # (stage_custom_climate.R:59,65) — ambiguous, so unenforceable, so skipped
        # rather than guessed at.
        roles: dict[int, set[bool]] = {}
        for line in text.splitlines():
            hits = [int(m.group(1)) for m in _ARG_IDX.finditer(line)]
            if not hits:
                continue
            m = _ASSIGN_LHS.match(line)
            if not m:
                continue
            is_out = bool(_OUT_NAME.search(m.group(1)))
            for k in hits:
                roles.setdefault(k, set()).add(is_out)
        self.out_slots = {k for k, v in roles.items() if v == {True}}
        self.in_slots = {k for k, v in roles.items() if v == {False}}

    @property
    def name(self) -> str:
        return os.path.basename(self.path)


def load_contracts(repo: str) -> dict[str, ScriptContract]:
    out: dict[str, ScriptContract] = {}
    sdir = os.path.join(repo, "scripts")
    for fn in sorted(os.listdir(sdir)):
        if not fn.endswith(".R"):
            continue
        p = os.path.join(sdir, fn)
        with open(p, errors="replace") as fh:
            text = fh.read()
        if "commandArgs" not in text:
            continue
        out[fn] = ScriptContract(os.path.join("scripts", fn), text)
    return out


# ---------------------------------------------------------------------------
# Snakemake side: rules, their params/input/output shapes, and their shell text
# ---------------------------------------------------------------------------

# The name is optional: rules generated in a loop (one per GAPIT model) are
# anonymous `rule:` blocks named by a `name:` directive, and requiring a name
# here silently skipped all four of them.
_RULE = re.compile(r"^(\s*)rule(?:\s+([A-Za-z_][A-Za-z0-9_]*))?\s*:")
_SECTION = re.compile(r"^(\s*)(name|input|output|params|log|shell|run|benchmark|resources"
                      r"|threads|wildcard_constraints|priority|message)\s*:(.*)$")
_KEY = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=")


class Rule:
    def __init__(self, name: str, file: str, line: int):
        self.name = name
        self.file = file
        self.line = line
        self.sections: dict[str, list[str]] = {}

    # -- shapes ------------------------------------------------------------
    def _keys(self, section: str) -> list[str]:
        return [m.group(1) for l in self.sections.get(section, [])
                for m in [_KEY.match(l)] if m]

    def n_outputs(self) -> int:
        body = self.sections.get("output", [])
        keys = self._keys("output")
        if keys:
            return len(keys)
        # keyless: count comma-separated top-level entries
        return max(1, _count_toplevel_entries(body))

    def n_inputs(self) -> int:
        keys = self._keys("input")
        if keys:
            return len(keys)
        return max(1, _count_toplevel_entries(self.sections.get("input", [])))

    def multi_token_params(self) -> set[str]:
        """params keys whose value is a SPACE-joined string (multi shell token)."""
        out = set()
        text = "\n".join(self.sections.get("params", []))
        # ' '.join(...) / " ".join(...) — comma/semicolon joins stay one token
        for m in re.finditer(
            r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(?:lambda[^:]*:\s*)?[\"']\s[\"']\s*\.join",
            text,
        ):
            out.add(m.group(1))
        return out

    def list_inputs(self) -> set[str]:
        """input keys whose value is plainly a list (expands to many tokens)."""
        out = set()
        cur = None
        buf: list[str] = []
        for l in self.sections.get("input", []):
            m = _KEY.match(l)
            if m:
                if cur:
                    out.update(_maybe_list(cur, " ".join(buf)))
                cur, buf = m.group(1), [l[m.end():]]
            elif cur:
                buf.append(l)
        if cur:
            out.update(_maybe_list(cur, " ".join(buf)))
        return out


def _maybe_list(key: str, value: str) -> set[str]:
    v = value.strip()
    if v.startswith("[") or re.search(r"lambda[^:]*:\s*\[", v) or ".join" in v:
        return {key}
    if re.search(r"\bfor\b.+\bin\b", v) and "[" in v:
        return {key}
    return set()


def _count_toplevel_entries(lines: list[str]) -> int:
    depth = 0
    n = 1
    for ch in "".join(lines):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            n += 1
    return n


def parse_rules(path: str, rel: str) -> list[Rule]:
    with open(path, errors="replace") as fh:
        lines = fh.read().splitlines()
    rules: list[Rule] = []
    cur: Rule | None = None
    cur_indent = 0
    section: str | None = None
    for i, raw in enumerate(lines, start=1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        m = _RULE.match(raw)
        if m:
            cur = Rule(m.group(2) or "<anonymous>", rel, i)
            cur_indent = len(m.group(1))
            rules.append(cur)
            section = None
            continue
        if cur is None:
            continue
        indent = len(raw) - len(raw.lstrip())
        if indent <= cur_indent:
            cur = None
            section = None
            continue
        s = _SECTION.match(raw)
        # Any depth below the rule line, not a fixed +4: gea_x_gwas.smk indents
        # its sections by 2, and an exact +4 silently dropped both its rules.
        if s and len(s.group(1)) > cur_indent:
            section = s.group(2)
            cur.sections.setdefault(section, [])
            rest = s.group(3)
            if section == "name" and cur.name == "<anonymous>":
                cur.name = rest.strip()  # the f-string source, e.g. f"gwas_a_{...}"
            if rest.strip():
                cur.sections[section].append(rest)
            continue
        if section:
            cur.sections[section].append(raw)
    return rules


# ---------------------------------------------------------------------------
# Extract Rscript invocations out of a shell: block
# ---------------------------------------------------------------------------

_RSCRIPT = re.compile(r"Rscript\s+(?:\S*/)?scripts/([A-Za-z0-9_]+\.R)")


_STR_PREFIX = re.compile(r"^[A-Za-z]*")


def _literal_body(lit: str) -> str:
    """Strip a Python string literal's prefix and delimiters, keep its content."""
    lit = _STR_PREFIX.sub("", lit, count=1)
    for q in ('"' * 3, "'" * 3, '"', "'"):
        if lit.startswith(q) and lit.endswith(q) and len(lit) >= 2 * len(q):
            return lit[len(q):-len(q)]
    return lit


def shell_text(rule: Rule) -> str:
    """Concatenate a shell: block's PYTHON string literals into one command.

    Extracting the literals with `tokenize` rather than deleting quote characters
    is load-bearing. The pipeline's only protection against a space-joined file
    list shifting every later argument IS a pair of shell-level double quotes
    (`"{params.sigsnps_str}"`), so a naive quote strip erases exactly the thing
    this linter has to see. Python delimiters go; shell quotes stay.
    """
    body = "\n".join(rule.sections.get("shell", []))
    parts: list[str] = []
    import io
    import tokenize as tk

    try:
        # Wrap in parentheses so an indented, implicitly-concatenated multi-line
        # literal is a syntactically complete expression.
        for tok in tk.generate_tokens(io.StringIO("(" + body + "\n)").readline):
            if tok.type == tk.STRING:
                parts.append(_literal_body(tok.string))
    except (tk.TokenError, IndentationError, SyntaxError):
        parts = [body]
    return " ".join(parts).replace("\\\n", " ")


def tokenize_invocation(cmd: str) -> list[str]:
    """Tokens passed to the script, up to the first redirect. Quotes preserved."""
    toks: list[str] = []
    i = 0
    n = len(cmd)
    while i < n:
        c = cmd[i]
        if c.isspace():
            i += 1
            continue
        if c in ">|;&":
            break
        if c == "2" and cmd[i:i + 4] == "2>&1":
            break
        start = i
        quoted = False
        while i < n and not cmd[i].isspace():
            if cmd[i] == '"':
                quoted = True
                i += 1
                while i < n and cmd[i] != '"':
                    i += 1
            i += 1
        tok = cmd[start:i]
        if tok.startswith(">") or tok.startswith("2>"):
            break
        toks.append(tok)
        del quoted
    return toks


class Violation:
    def __init__(self, check, severity, where, detail):
        self.check = check
        self.severity = severity
        self.where = where
        self.detail = detail

    def __str__(self):
        return f"[{self.severity}] {self.check}  {self.where}\n        {self.detail}"


def token_width(tok: str, rule: Rule, mt_params: set[str], list_inputs: set[str]):
    """(min_tokens, max_tokens_or_None, reason) this token expands to at run time."""
    bare = tok.strip()
    is_quoted = bare.startswith('"') and bare.endswith('"') and len(bare) >= 2
    inner = bare[1:-1] if is_quoted else bare
    if is_quoted:
        return (1, 1, None)
    m = re.fullmatch(r"\{params\.([A-Za-z0-9_]+)\}", inner)
    if m and m.group(1) in mt_params:
        return (1, None, f"params.{m.group(1)} is a space-joined list")
    if re.fullmatch(r"\{output\}", inner):
        k = rule.n_outputs()
        return (k, k, f"bare {{output}} in a {k}-output rule" if k > 1 else None)
    if re.fullmatch(r"\{input\}", inner):
        k = rule.n_inputs()
        return (k, k, f"bare {{input}} in a {k}-input rule" if k > 1 else None)
    m = re.fullmatch(r"\{input\.([A-Za-z0-9_]+)\}", inner)
    if m and m.group(1) in list_inputs:
        return (1, None, f"input.{m.group(1)} is a list")
    return (1, 1, None)


def lint(repo: str) -> tuple[list[Violation], dict]:
    contracts = load_contracts(repo)
    violations: list[Violation] = []
    stats = {"scripts": len(contracts), "rules": 0, "call_sites": 0}

    # ---- script-intrinsic checks ----
    for fn, c in sorted(contracts.items()):
        if c.missing:
            violations.append(Violation(
                "sparse_argv", "error", f"{c.path}",
                f"reads args[{c.max_index}] but never args["
                + ",".join(str(i) for i in c.missing) + "]"))
        if c.max_index >= MIN_GUARDED_ARGS and not c.has_len_check:
            violations.append(Violation(
                "no_arity_guard", "warning", f"{c.path}",
                f"{c.max_index} positional arguments, no length(args) check — a "
                f"shifted call site cannot be detected at run time"))

    # ---- call-site checks ----
    rdir = os.path.join(repo, "workflow", "rules")
    for fn in sorted(os.listdir(rdir)):
        if not fn.endswith(".smk"):
            continue
        rules = parse_rules(os.path.join(rdir, fn), f"workflow/rules/{fn}")
        stats["rules"] += len(rules)
        for rule in rules:
            cmd = shell_text(rule)
            mt = rule.multi_token_params()
            li = rule.list_inputs()
            for m in _RSCRIPT.finditer(cmd):
                script = m.group(1)
                stats["call_sites"] += 1
                c = contracts.get(script)
                if c is None:
                    violations.append(Violation(
                        "unknown_script", "error",
                        f"{rule.file}:{rule.line} rule {rule.name}",
                        f"invokes scripts/{script}, which does not read commandArgs()"))
                    continue
                toks = tokenize_invocation(cmd[m.end():])
                lo = hi = 0
                unbounded = []
                for t in toks:
                    a, b, why = token_width(t, rule, mt, li)
                    lo += a
                    if b is None:
                        unbounded.append((t, why))
                        hi = None
                    elif hi is not None:
                        hi += b
                where = f"{rule.file}:{rule.line} rule {rule.name} -> {script}"
                if c.variadic:
                    continue  # contract is "N inputs then an output"; arity is free
                for t, why in unbounded:
                    violations.append(Violation(
                        "unquoted_multi_token", "error", where,
                        f"{t} is passed unquoted into a fixed positional slot "
                        f"({why}); it expands to >1 shell token and shifts every "
                        f"later argument"))
                if not unbounded:
                    for slot, t in enumerate(toks, start=1):
                        inner = t.strip().strip('"')
                        is_out = bool(re.fullmatch(r"\{output(\.[A-Za-z0-9_]+)?\}", inner))
                        is_in = bool(re.fullmatch(r"\{input(\.[A-Za-z0-9_]+)?\}", inner))
                        if slot in c.out_slots and is_in:
                            violations.append(Violation(
                                "slot_role_mismatch", "error", where,
                                f"position {slot} is an OUTPUT slot in {script} but "
                                f"the call site passes {t} — the script would "
                                f"overwrite one of its own inputs"))
                        elif is_out and slot in c.in_slots:
                            violations.append(Violation(
                                "slot_role_mismatch", "error", where,
                                f"position {slot} is an INPUT slot in {script} "
                                f"(its output slots are {sorted(c.out_slots)}) but "
                                f"the call site passes {t} — arguments are shifted"))
                    if lo < c.required:
                        violations.append(Violation(
                            "arity_mismatch", "error", where,
                            f"call site passes {lo} tokens; script requires "
                            f"{c.required} (max {c.max_index})"))
                    elif lo > c.max_index:
                        violations.append(Violation(
                            "arity_mismatch", "error", where,
                            f"call site passes {lo} tokens; script reads at most "
                            f"args[{c.max_index}] — {lo - c.max_index} ignored"))
    return violations, stats


def main(argv: list[str]) -> int:
    repo = "."
    tsv = None
    quiet = False
    rest = []
    i = 0
    while i < len(argv):
        if argv[i] == "--tsv":
            tsv = argv[i + 1]
            i += 2
        elif argv[i] == "--quiet":
            quiet = True
            i += 1
        else:
            rest.append(argv[i])
            i += 1
    if rest:
        repo = rest[0]

    violations, stats = lint(repo)
    errors = [v for v in violations if v.severity == "error"]

    if tsv:
        with open(tsv, "w") as fh:
            fh.write("check\tseverity\twhere\tdetail\n")
            for v in violations:
                fh.write(f"{v.check}\t{v.severity}\t{v.where}\t"
                         f"{v.detail.replace(chr(10), ' ')}\n")
    if not quiet:
        for v in violations:
            print(v)
        print(f"\n{stats['scripts']} argv-reading scripts, {stats['rules']} rules, "
              f"{stats['call_sites']} Rscript call sites")
        print(f"{len(errors)} error(s), {len(violations) - len(errors)} warning(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
