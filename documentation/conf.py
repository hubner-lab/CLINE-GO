"""Sphinx configuration for the CLINE-GO documentation site.

Sphinx root: documentation/. Build locally with
    nix shell nixpkgs#python312 -c documentation/build_docs.sh <outdir> <venv dir>
Read the Docs reads .readthedocs.yaml at the repository root (warnings fail the build there too).

How pages use what this file defines
------------------------------------
Substitutions (all in ``myst_substitutions``):
  run_<config>_<mode>  the whole command for one pipeline mode, built by run_cmd();
                       <config> is a key of CONFIGS, <mode> one of MODES
  app                  the command that starts the Shiny app
  qc_* / tut_*         numbers from verified/numbers.tsv and tutorial/verified/numbers.tsv

  In text, write MyST ``{{key}}``. In a code block, write ``|key|``; it is replaced only inside a
  ``{code-block}`` directive (sphinx-substitution-extensions), never in a plain ``` fence:

      ```{code-block} text
      :class: clinego-run
      |run_populus_processing|
      ```

Runnable marker: every block a reader is meant to run carries ``:class: clinego-run`` (the HTML
gets ``<div class="clinego-run highlight-text ...">``), in page order, so a script can pull the
commands from the built pages. Output listings and file excerpts never carry it. Runnable blocks
use the ``text`` lexer, not ``bash``: bash highlighting wraps every word in its own <span>, so
the command would no longer appear verbatim in the HTML (``grep 'docker run'`` finds nothing).

In development: ``:::{in-development}`` ... ``:::`` (``_ext/clinego_docs.py``) marks what is not
built yet, with one title and one look everywhere. ``:::{slot} <title>`` is a dashed placeholder
for a figure, number or run time still to come; ``grep -rn '{slot}'`` lists what is left.

Install-page commands (git clone, docker build) carry no clinego-run class on purpose: a check
script that runs the marked commands must clone and build in its own way (never retag the image).

Tutorial chapters follow _templates/chapter.md (excluded from the build).

PACKAGED switches every command at once, from today's form (the repo mounted over /pipeline)
to the packaged image (code inside the image, one ``clinego`` command).
"""

import csv
import re
import sys
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(DOCS_DIR / "_ext"))

# -- Project -------------------------------------------------------------------------------------

project = "CLINE-GO"
author = "The CLINE-GO authors"
copyright = "2026, The CLINE-GO authors"

REPO_URL = "https://github.com/hubner-lab/CLINE-GO"

# -- Commands ------------------------------------------------------------------------------------

# False: the image holds no pipeline code, so every command mounts the repo over /pipeline and is
# run from the repo root. True: the packaged image (a later plan) with its one `clinego` command.
PACKAGED = False

# Packaged image tag; a placeholder until the packaging plan publishes one.
RELEASE = "<release>"

# Config name in a substitution -> its path, relative to the repo root.
CONFIGS = {
    "populus": "documentation/tutorial/populus/populus.yaml",
    "gardens": "documentation/tutorial/populus/populus_gardens.yaml",
    "testdata": "test_data/config_testdata.yaml",
}

MODES = (
    "processing",
    "prestructure",
    "structure",
    "climate",
    "traits",
    "pregea",
    "gea",
    "gwas",
    "gea_x_gwas",
    "maladaptation",
)

MEMORY = "20g"
USER_FLAG = "--user $(id -u):$(id -g)"
CONTINUATION = " \\\n  "  # one command, wrapped for reading; still pastes as one line of shell


def run_cmd(mode: str, config: str) -> str:
    """The whole command that runs one pipeline mode on one config, in the PACKAGED form."""
    if mode not in MODES:
        raise ValueError(f"run_cmd: unknown mode {mode!r}; expected one of {MODES}")
    if config not in CONFIGS:
        raise ValueError(f"run_cmd: unknown config {config!r}; expected one of {sorted(CONFIGS)}")
    config_path = CONFIGS[config]
    if PACKAGED:
        parts = [
            f"docker run --rm --name clinego-{mode} --memory={MEMORY} {USER_FLAG}",
            f"-v $PWD:/data ghcr.io/hubner-lab/cline-go:{RELEASE}",
            f"clinego run {mode} --config {config_path}",
        ]
    else:
        parts = [
            f"docker run {USER_FLAG} --rm --name clinego-{mode}",
            f"--memory={MEMORY} -e USER=pipeline -v $PWD:/pipeline cline-go:latest",
            # --scheduler sits before --config so the --configfile line stays short (flag order
            # is free); the interface's one-line form ends with it.
            f"snakemake -c4 -s Snakefile --scheduler greedy --config mode={mode}",
            f"--configfile {config_path}",
        ]
    return CONTINUATION.join(parts)


def app_cmd() -> str:
    """The command that starts the Shiny app on port 3838, in the PACKAGED form."""
    if PACKAGED:
        parts = [
            f"docker run --rm --name clinego-app --memory={MEMORY} {USER_FLAG} -p 3838:3838",
            f"-v $PWD:/data ghcr.io/hubner-lab/cline-go:{RELEASE}",
            "clinego app",
        ]
    else:
        parts = [
            f"docker run {USER_FLAG} --rm --name clinego-app",
            f"--memory={MEMORY} -e USER=pipeline -p 3838:3838 -v $PWD:/pipeline cline-go:latest",
            "Rscript /pipeline/scripts/clinego.app/dev.R",
        ]
    return CONTINUATION.join(parts)


# -- Verified numbers ----------------------------------------------------------------------------

# Each file: a header row, then tab-separated `key` and `value` (further columns are ignored).
# A missing file is no error; a malformed one, or a key defined twice, fails the build.
NUMBERS_FILES = (
    ("verified/numbers.tsv", "qc_"),           # the quick check (lane site)
    ("tutorial/verified/numbers.tsv", "tut_"),  # the Populus tutorial (lane populus)
)
NUMBER_KEY = re.compile(r"^[a-z][a-z0-9_]*$")


def load_numbers(files=NUMBERS_FILES, root=DOCS_DIR) -> dict:
    numbers: dict[str, str] = {}
    origin: dict[str, str] = {}
    for rel, prefix in files:
        path = root / rel
        if not path.is_file():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None or header[:2] != ["key", "value"]:
                raise ValueError(f"{rel}: the header must start with 'key<TAB>value', got {header}")
            count = 0
            for line_no, row in enumerate(reader, start=2):
                if not row or all(cell.strip() == "" for cell in row):
                    continue
                if len(row) < 2:
                    raise ValueError(f"{rel}:{line_no}: expected key<TAB>value, got {row}")
                key, value = row[0], row[1]
                if not NUMBER_KEY.match(key):
                    raise ValueError(f"{rel}:{line_no}: key {key!r} must match {NUMBER_KEY.pattern}")
                if not key.startswith(prefix):
                    raise ValueError(f"{rel}:{line_no}: key {key!r} must start with {prefix!r}")
                if value.strip() == "":
                    raise ValueError(f"{rel}:{line_no}: key {key!r} has an empty value")
                if key in numbers:
                    raise ValueError(
                        f"{rel}:{line_no}: key {key!r} is already defined in {origin[key]}"
                    )
                numbers[key] = value
                origin[key] = rel
                count += 1
        print(f"numbers: {count} keys from {rel}")
    return numbers


def build_substitutions() -> dict:
    subs = {f"run_{config}_{mode}": run_cmd(mode, config) for config in CONFIGS for mode in MODES}
    subs["app"] = app_cmd()
    numbers = load_numbers()
    clash = sorted(set(numbers) & set(subs))
    if clash:
        raise ValueError(f"numbers keys clash with command substitutions: {clash}")
    subs.update(numbers)
    return subs


# -- Extensions ----------------------------------------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_substitution_extensions",
    "clinego_docs",
]

myst_enable_extensions = [
    "attrs_block",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "substitution",
]
myst_heading_anchors = 3
myst_substitutions = build_substitutions()

# Replace |key| in every code-block (and literalinclude) without a per-block option.
substitutions_default_enabled = True

source_suffix = {".md": "markdown"}
root_doc = "index"

# Lane populus's tutorial inputs and both verified/ tables are data, never pages.
exclude_patterns = [
    "_build",
    "_templates/**",  # chapter.md is a template for writers, not a page
    "tutorial/populus/**",
    "tutorial/verified/**",
    "verified/**",
]

# -- HTML output ---------------------------------------------------------------------------------

html_theme = "sphinx_book_theme"
html_title = "CLINE-GO"
html_logo = "../assets/cline-go-logo.png"
html_static_path = ["_static"]
html_css_files = ["clinego.css"]

# Light by default (pipeline figures are drawn on white); the theme switcher stays in the header.
html_context = {"default_mode": "light"}

html_theme_options = {
    "repository_url": REPO_URL,
    "repository_branch": "main",
    "path_to_docs": "documentation",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": False,
    "use_download_button": False,
    "home_page_in_toc": True,
    "show_toc_level": 2,
    # pydata-sphinx-theme 0.20 otherwise adds a second search field in a header strip above
    # the page; the sidebar already has one.
    "navbar_persistent": [],
}
