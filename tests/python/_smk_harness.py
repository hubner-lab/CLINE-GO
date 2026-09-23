"""Load ``workflow/rules/common.smk`` as ordinary Python, and model the rule graph.

``common.smk`` is 2 354 lines — the largest file in the repository — and it is the
only place that decides, for every mode, WHICH FILES the pipeline is asked to
build. Until now nothing tested it: the R suites cannot see it, and the only
check that ever exercised it is ``snakemake -n``, which needs Docker, a VCF and a
metadata file, so it is not something the merge gate can run.

It is testable, because it is plain Python:

  * ``workdir: OUTDIR`` (:1996) is a bare variable *annotation*, not a Snakemake
    directive as far as the interpreter is concerned, and there is no ``rule``
    block in the file — so ``exec()`` runs it to completion.
  * The only two globals it needs are ``config`` (a dict) and ``workflow`` (used
    once, for ``.basedir``).
  * Three filesystem calls have to be redirected, all at load time:
    ``os.path.exists`` and ``open`` (the input-file checks and the metadata
    header read, which both address ``/pipeline/...``) and ``os.makedirs``
    (:1994 creates the output tree). ``patched()`` does that; ``makedirs`` is
    sent to a tempdir so importing the module never writes into the repository.

``rule_outputs()`` then reads the OTHER ``.smk`` files as text, pulls each rule's
``output:`` block, and evaluates it in that same namespace — giving the set of
paths the workflow can actually produce, which is what a target list has to be
checked against.  Three shapes have to be handled or a producible target reads as
unreachable:

  * rules nested in ``for <var> in ...:`` factory loops (the GAPIT / GWAS model
    blocks) — ``_loop_env()`` replays the binding;
  * ``_out_path``-style locals assigned in the loop body — ``_local_assigns()``;
  * anonymous ``rule:`` blocks with a ``name:`` field, which have no rule name.

Configs are built as dicts by ``make_config()`` rather than read from YAML: PyYAML
is present in the image but the fixture wants to be perturbable field by field,
and a dict literal makes the perturbation the readable part of each test.  The
values mirror ``test_data/config_testdata.yaml`` (including its deliberately
float-spelled ``k_end: 7.0`` / ``ntree: 1000.0``, which exist to exercise the
integer casting at :89) with ONE deviation: ``GEA.snp_clumping_distance``
replaces the shipped file's ``region_distance``/``combine_gap`` pair, which
common.smk:276-287 accepts but answers with a DeprecationWarning per load. The
quick gate's warning baseline is zero, and the fixture should not be the thing
that breaks it.
"""
import ast
import contextlib
import copy
import glob
import os
import re
import tempfile
from unittest import mock

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COMMON_SMK = os.path.join(REPO_ROOT, "workflow", "rules", "common.smk")

# Everything common.smk writes goes here instead of into the repository.
_SANDBOX = tempfile.mkdtemp(prefix="clinego_smk_")

MODES = ["processing", "prestructure", "structure", "climate", "traits",
         "pregea", "gea", "gwas", "gea_x_gwas", "maladaptation"]


class _Workflow:
    """The single attribute common.smk reads off Snakemake's `workflow` object."""
    basedir = REPO_ROOT


_real_exists = os.path.exists
_real_makedirs = os.makedirs
_real_open = open
_real_glob = glob.glob


def _under_pipeline(path):
    path = str(path)
    return path[len("/pipeline/"):] if path.startswith("/pipeline/") else None


def _exists(path):
    tail = _under_pipeline(path)
    if tail is None:
        return _real_exists(path)
    return _real_exists(os.path.join(REPO_ROOT, tail)) or _real_exists(os.path.join(_SANDBOX, tail))


def _makedirs(path, **kw):
    tail = _under_pipeline(path)
    return _real_makedirs(os.path.join(_SANDBOX, tail) if tail is not None else path, **kw)


def _glob(pattern, **kw):
    """resolve_active_snp_sets() globs /pipeline/... for saved SNP sets."""
    tail = _under_pipeline(pattern)
    if tail is None:
        return _real_glob(pattern, **kw)
    return (_real_glob(os.path.join(REPO_ROOT, tail), **kw)
            + _real_glob(os.path.join(_SANDBOX, tail), **kw))


def _open(file, *a, **kw):
    tail = _under_pipeline(file)
    if tail is None:
        return _real_open(file, *a, **kw)
    in_repo = os.path.join(REPO_ROOT, tail)
    return _real_open(in_repo if _real_exists(in_repo) else os.path.join(_SANDBOX, tail), *a, **kw)


@contextlib.contextmanager
def patched():
    """Redirect the three filesystem calls common.smk makes against /pipeline.

    Needed around get_targets() too, not only around load(): several modes
    re-check input files (check_file_exists for the GFF) when their target list
    is built.
    """
    with mock.patch("os.path.exists", _exists), \
         mock.patch("os.makedirs", _makedirs), \
         mock.patch("glob.glob", _glob), \
         mock.patch("builtins.open", _open):
        yield


# --------------------------------------------------------------- config fixture
# Mirrors test_data/config_testdata.yaml. Float spellings are deliberate.
BASE_CONFIG = {
    "project_name": "testdata",
    "Input": {"dir": "test_data/", "vcf": "testdata.vcf",
              "metadata": "testdata_metadata.tsv", "gff": "testdata.gff3"},
    "Filter": {"maf": 0.05, "snp_miss": 0.2, "sample_miss": 0.5,
               "relatedness": 0.99, "relatedness_action": "remove"},
    "LD": {"window": 100, "step": 20, "r2": 0.2},
    "sNMF": {"k_start": 2, "k_end": 7.0, "k_best": 3.0, "repeats": 10},
    "Map": {"climate_extent": "auto", "gap": 0.5, "resolution": 2.5,
            "zoom_extent": "34.5,35.5,31,32.5"},
    "Climate": {"enabled": True, "predictors": "bio_1,bio_2,bio_3",
                "Varpart": {"response": "pcs", "response_var_cutoff": 0.8,
                            "response_max_pcs": 20, "response_min_pcs": 2,
                            "structure_table": "qmatrix", "permutations": 999}},
    "Population": {"calc_stats": False, "window_size": 3000000,
                   "custom_trait_file": "NULL"},
    "Piemap": {"alpha": 0.6, "show_labels": False, "label_size": 8,
               "pie_scale": 1.0, "use_points": False},
    "LDdecay": {"group_by": "cluster", "min_samples": 10,
                "max_distance": 500, "scope": "both"},
    "GEA": {"configs": [
                {"method": "EMMAX", "adjust": "bonf", "threshold": "0.05"},
                {"method": "LFMM", "adjust": "bonf", "threshold": "0.05"},
                {"method": "RDA", "adjust": "bonf", "threshold": "0.01",
                 "params": {"condition_pcs": 3, "axes": "auto"}}],
            "snp_clumping_distance": 1000000, "promoter_length": 10000},
    "PreGEA": {"predictors": "bio_1,bio_2,bio_3", "k_offset": 2, "n_pcs_max": 10,
               "Advanced": {"collinearity_r": 0.7, "vif_max": 10.0,
                            "axis_alpha": 0.05, "permutations": 199},
               "TransferGuard": {"enabled": False, "lfmm_k": "auto",
                                 "emmax_n_pcs": "auto"}},
    "GFF": {"feature": "mRNA", "gene_name": "description",
            "biotype": "biotype", "go_field": "ontology"},
    "Enrichment": {"top_terms": 20, "plot_width": 12, "plot_height": 10,
                   "cnet_label": "gene_id"},
    "Future": {"ssp": "585", "year": "2061-2080",
               "models": "MPI-ESM1-2-HR,UKESM1-0-LL,IPSL-CM6A-LR"},
    "Maladaptation": {"methods": {
        "gradient_forest": {"ntree": 1000.0, "cor_threshold": 0.5,
                            "spatial_correction": "both", "random_model": True,
                            "extrap": True},
        "geometric_offset": {"scale": True, "k": 3.0},
        "rda_offset": {"axes": "auto", "axis_alpha": 0.05, "permutations": 199,
                       "condition_pcs": 0, "seed": 42}},
        "snp_sets": "all"},
    "GEAxGWAS": {"pairwise": {"window_size": 500000, "min_snps": 2}},
    "GWAS": {"missing_strategy": "DROP",
             "traits": "height,flowering_time,disease_score",
             "configs": [{"method": "EMMAX", "adjust": "bonf", "threshold": "0.05"},
                         {"method": "BLINK", "adjust": "bonf", "threshold": "0.05"}],
             "promoter_length": 10000.0},
}


def make_config(**overrides):
    """Deep copy of BASE_CONFIG with `Group__key=value` overrides applied."""
    cfg = copy.deepcopy(BASE_CONFIG)
    for dotted, value in overrides.items():
        parts = dotted.split("__")
        node = cfg
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        node[parts[-1]] = value
    return cfg


def load(config=None):
    """exec common.smk with `config` and return its module namespace."""
    namespace = {
        "config": copy.deepcopy(BASE_CONFIG if config is None else config),
        "workflow": _Workflow(),
        "__name__": "clinego_common_smk",
        "__file__": COMMON_SMK,
    }
    with _real_open(COMMON_SMK) as handle:
        source = handle.read()
    with patched():
        exec(compile(source, COMMON_SMK, "exec"), namespace)
    return namespace


def targets(namespace, mode):
    with patched():
        return namespace["get_targets"](mode)


# ----------------------------------------------------------- rule-graph model
_RULE = re.compile(r"^\s*(?:rule|checkpoint)\s*([A-Za-z_][A-Za-z0-9_]*)?\s*:")
_SECTION = re.compile(r"^(\s*)(input|output|params|log|threads|shell|run|resources|"
                      r"wildcard_constraints|benchmark|priority|message|conda|group|"
                      r"retries|container|envmodules|localrule|default_target|"
                      r"notebook|script|name)\s*:")
_FOR = re.compile(r"^(\s*)for\s+([A-Za-z_][A-Za-z0-9_]*)\s+in\s+(.+?):\s*$")
_ASSIGN = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$")


def _indent(line):
    return len(line) - len(line.lstrip())


def _install_stubs(namespace):
    """Snakemake's output-flag helpers, which are not defined in common.smk."""
    def identity(value=None, *_a, **_kw):
        return value
    for flag in ("directory", "temp", "protected", "touch", "ancient",
                 "unpack", "report", "pipe"):
        namespace.setdefault(flag, identity)
    if "expand" not in namespace:
        import itertools

        def expand(pattern, **kw):
            if not kw:
                return [pattern]
            keys = list(kw)
            vals = [kw[k] if isinstance(kw[k], (list, tuple)) else [kw[k]] for k in keys]
            return [pattern.format(**dict(zip(keys, combo)))
                    for combo in itertools.product(*vals)]
        namespace["expand"] = expand
    return namespace


def _loop_env(lines, idx, namespace):
    """Bindings from every enclosing `for X in EXPR:` above the rule at `idx`.

    Returns [] when a loop iterates an empty sequence: the rule is not generated
    at all under this config, so it contributes no producible path.
    """
    envs = [{}]
    inner = _indent(lines[idx])
    for j in range(idx - 1, -1, -1):
        line = lines[j]
        if not line.strip():
            continue
        match = _FOR.match(line)
        if match and _indent(line) < inner:
            try:
                values = list(eval(match.group(3), namespace, {}))
            except Exception:
                values = []
            envs = [dict(env, **{match.group(2): v}) for env in envs for v in values]
            inner = _indent(line)
            if not envs:
                return []
    return envs


def _local_assigns(lines, idx, namespace, env):
    """Replay `name = expr` statements in the rule's own block (e.g. `_out_path`)."""
    env = dict(env)
    rule_indent = _indent(lines[idx])
    for j in range(max(0, idx - 40), idx):
        line = lines[j]
        if not line.strip() or _indent(line) > rule_indent:
            continue
        match = _ASSIGN.match(line)
        if not match or match.group(3).rstrip().endswith((",", "(", "[", "{", "\\")):
            continue
        try:
            env[match.group(2)] = eval(match.group(3), namespace, env)
        except Exception:
            pass                      # multi-line or wildcard-dependent: skipped
    return env


def _split_output_block(text):
    """Expressions in an `output:` block, with any `name =` labels dropped."""
    try:
        node = ast.parse("_f(" + text + "\n)", mode="eval")
    except SyntaxError:
        return []
    call = node.body
    return ([ast.unparse(a) for a in call.args]
            + [ast.unparse(k.value) for k in call.keywords])


def rule_outputs(namespace):
    """Every rule output in workflow/rules/*.smk, evaluated.

    -> list of (smk_basename, rule_name, expression, value, error)
    """
    _install_stubs(namespace)
    results = []
    for path in sorted(glob.glob(os.path.join(REPO_ROOT, "workflow", "rules", "*.smk"))):
        if os.path.basename(path) == "common.smk":
            continue
        with _real_open(path) as handle:
            lines = handle.read().splitlines()
        i, rule, rule_indent, envs = 0, None, 0, [{}]
        while i < len(lines):
            rule_match = _RULE.match(lines[i])
            if rule_match:
                rule = rule_match.group(1) or "<anonymous>"
                rule_indent = _indent(lines[i])
                envs = [_local_assigns(lines, i, namespace, e)
                        for e in _loop_env(lines, i, namespace)]
                i += 1
                continue
            section = _SECTION.match(lines[i])
            if section and section.group(2) == "output" and rule:
                indent = len(section.group(1))
                body = [lines[i][section.end():]]
                j = i + 1
                while j < len(lines):
                    line = lines[j]
                    if not line.strip():
                        body.append("")
                        j += 1
                        continue
                    here = _indent(line)
                    if here <= indent and (_SECTION.match(line) or _RULE.match(line)):
                        break
                    if here <= indent and line.lstrip().startswith('"""'):
                        break
                    if here <= rule_indent:
                        break
                    body.append(line)
                    j += 1
                block = "\n".join(body).strip().rstrip(",")
                for expression in _split_output_block(block):
                    for env in envs:
                        try:
                            results.append((os.path.basename(path), rule, expression,
                                            eval(expression, namespace, dict(env)), None))
                        except Exception as exc:
                            results.append((os.path.basename(path), rule, expression,
                                            None, f"{type(exc).__name__}: {exc}"))
                i = j
                continue
            i += 1
    return results


def producible_patterns(namespace):
    """Compiled regexes for every producible path; `{wildcard}` matches one segment."""
    patterns = []
    for smk, rule, _expr, value, error in rule_outputs(namespace):
        if error is not None:
            continue
        for item in (value if isinstance(value, (list, tuple)) else [value]):
            if isinstance(item, str):
                parts = re.split(r"(\{[^{}]*\})", item)
                regex = "".join(r"[^/]+" if p.startswith("{") else re.escape(p)
                                for p in parts)
                patterns.append((re.compile("^" + regex + "$"), smk, rule))
    return patterns


def save_snp_set(project, name):
    """Create a curated SNP set in the sandbox so mode=maladaptation resolves.

    resolve_active_snp_sets() (common.smk:1467) reads the filesystem, so the
    maladaptation target list is a function of run state, not of config alone.
    """
    directory = os.path.join(_SANDBOX, f"{project}_results", "_intermediate",
                             "snp_sets", name)
    _real_makedirs(directory, exist_ok=True)
    with _real_open(os.path.join(directory, "selected_snps.tsv"), "w") as handle:
        handle.write("chr\tpos\n1\t100\n")
    return directory
