"""Sphinx configuration for the CLINE-GO documentation site.

Sphinx root: documentation/. Build locally with
    nix shell nixpkgs#python312 -c documentation/build_docs.sh <outdir> <venv dir>
GitHub Pages: .github/workflows/docs.yml runs build_docs.sh on every push to main (warnings fail it there too).

How pages use what this file defines
------------------------------------
Substitutions (all in ``myst_substitutions``):
  run_<config>_<mode>  the whole command for one pipeline mode, built by run_cmd();
                       <config> is a key of CONFIGS, <mode> one of MODES
  app                  the command that starts the Shiny app
  build / example /    the launcher commands: build the image, copy the bundled test dataset,
  init_own             write config_myproject.yaml (PACKAGED only)
  config_<config>      the config file a reader passes to --config (config_testdata.yaml, ...)
  qc_* / tut_*         numbers from verified/numbers.tsv and tutorial/verified/numbers.tsv
                       (tut_* only while TUTORIAL is True)

  In text, write MyST ``{{key}}``. In a code block, write ``|key|``; it is replaced only inside a
  ``{code-block}`` directive (sphinx-substitution-extensions), never in a plain ``` fence:

      ```{code-block} text
      :class: clinego-run
      |run_testdata_processing|
      ```

Runnable marker: every block a reader is meant to run carries ``:class: clinego-run`` (the HTML
gets ``<div class="clinego-run highlight-text ...">``), in page order, so a script can pull the
commands from the built pages. Output listings and file excerpts never carry it. Runnable blocks
use the ``text`` lexer, not ``bash``: bash highlighting wraps every word in its own <span>, so
the command would no longer appear verbatim in the HTML (``grep 'docker run'`` finds nothing).

In development: ``:::{in-development}`` ... ``:::`` (``_ext/clinego_docs.py``) marks what is not
built yet, with one title and one look everywhere. ``:::{slot} <title>`` is a dashed placeholder
for a figure, number or run time still to come; ``grep -rn '{slot}'`` lists what is left.

Only the example, run and app blocks carry clinego-run. Install-page commands (git clone,
./clinego build, the PATH lines) and the shell steps around them (mkdir, cd) carry no class on
purpose: a check script that runs the marked commands must clone and build in its own way (never
retag the image) and choose its own folder. FAQ snippets (--show, unlock, docker stop) are never
marked either. The class also marks only commands that run as written: the quick check's
(testdata) blocks. The own-data and manual pages' commands (`clinego init myproject`,
run_own_<mode>) point at a project the reader builds first, so they stay unmarked.

Tutorial chapters follow _templates/chapter.md (excluded from the build).

PACKAGED switches every command at once, between the developer form (the repo mounted over
/pipeline, run from the repo root) and the packaged image (code inside the image, one ``clinego``
command run from the user's project folder). The released site is PACKAGED = True, and the page
prose is written for it: it names `clinego` literally, and the launcher substitutions (clone,
build, example, ...) exist only when PACKAGED is True, so the developer form is a fallback for the
run/app commands alone, not a second site.

TUTORIAL = False hides the parked Populus tutorial: tutorial/** joins exclude_patterns, the
populus and gardens configs leave CONFIGS, and tutorial/verified/numbers.tsv is not loaded. The
tutorial's files stay in the repository. To bring it back, set TUTORIAL = True and restore the
links that were removed when it was hidden (each as it stood at commit a3c1cb9, so
``git show a3c1cb9:documentation/<file>`` shows the old text):
  - index.md: the "Tutorial: balsam poplar" card (link tutorial/index), the "Does the offset mean
    anything? Chapter 7 ..." paragraph below the cards, and the toctree entry tutorial/index
  - quickstart.md: the last line's "Go on to the tutorial" link
  - manual/index.md: the link to the tutorial for the main settings in context
  - citing.md: the "Tutorial data" section (Fitzpatrick et al. 2021 and its data repository)
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

# True: the packaged image, code inside it, one `clinego` command run from the project folder.
# False: the developer form, the repo mounted over /pipeline and every command run from its root.
PACKAGED = True

# False hides the parked Populus tutorial (see the header for how to bring it back).
TUTORIAL = False

RELEASE_TAG = "v1.0"
# The branch or tag the clone command checks out. None clones main: the RELEASE_TAG tag does not
# exist on GitHub yet, and `git clone --branch v1.0` fails without it. Set it to RELEASE_TAG once
# the tag is pushed.
CLONE_BRANCH = None
CLONE_URL = "https://github.com/hubner-lab/CLINE-GO.git"

# Config name in a substitution -> (packaged path, developer path). The packaged path is relative
# to the project folder (`clinego example` writes config_testdata.yaml there, `clinego init
# myproject` config_myproject.yaml); the developer path is relative to the repo root.
CONFIGS = {
    "testdata": ("config_testdata.yaml", "test_data/config_testdata.yaml"),
    "own": ("config_myproject.yaml", "config_myproject.yaml"),
}
if TUTORIAL:
    CONFIGS["populus"] = (None, "documentation/tutorial/populus/populus.yaml")
    CONFIGS["gardens"] = (None, "documentation/tutorial/populus/populus_gardens.yaml")

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
    if PACKAGED:
        config_path = config_file(config)
        return f"clinego run {mode} --config {config_path}"
    else:
        config_path = CONFIGS[config][1]
        parts = [
            f"docker run {USER_FLAG} --rm --name clinego-{mode}",
            f"--memory={MEMORY} -e USER=pipeline -v $PWD:/pipeline cline-go:latest",
            # --scheduler sits before --config so the --configfile line stays short (flag order
            # is free); the interface's one-line form ends with it.
            f"snakemake -c4 -s Snakefile --scheduler greedy --config mode={mode}",
            f"--configfile {config_path}",
        ]
    return CONTINUATION.join(parts)


def config_file(config: str) -> str:
    """The file a reader passes to --config, as the current form names it."""
    packaged, developer = CONFIGS[config]
    path = packaged if PACKAGED else developer
    if path is None:
        raise ValueError(f"config {config!r} has no packaged form")
    return path


def launcher_cmds() -> dict:
    """The launcher's own commands (PACKAGED only): build, example, init, clone, PATH lines."""
    if not PACKAGED:
        return {}
    clone = f"git clone {CLONE_URL}"
    if CLONE_BRANCH:
        clone = f"git clone --branch {CLONE_BRANCH} {CLONE_URL}"
    return {
        "clone": clone,
        "build": "./clinego build",
        # The repository folder on PATH: `path` for this terminal; `path_keep`, run in that folder,
        # writes the same line with the folder's absolute path to ~/.bashrc for every new terminal.
        "path": 'export PATH="$PWD:$PATH"',
        "path_keep": r'echo "export PATH=\"$PWD:\$PATH\"" >> ~/.bashrc',
        "version": "clinego version",
        "example": "clinego example",
        "init_own": "clinego init myproject",
        "release_tag": RELEASE_TAG,
    }


def app_cmd() -> str:
    """The command that starts the Shiny app on port 3838, in the PACKAGED form."""
    if PACKAGED:
        return "clinego app"
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
NUMBERS_FILES = (("verified/numbers.tsv", "qc_"),)  # the quick check
if TUTORIAL:
    NUMBERS_FILES += (("tutorial/verified/numbers.tsv", "tut_"),)  # the Populus tutorial
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
    configs = [c for c in CONFIGS if not PACKAGED or CONFIGS[c][0] is not None]
    subs = {f"run_{config}_{mode}": run_cmd(mode, config) for config in configs for mode in MODES}
    subs.update({f"config_{config}": config_file(config) for config in configs})
    subs["app"] = app_cmd()
    subs.update(launcher_cmds())
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
if not TUTORIAL:
    exclude_patterns.append("tutorial/**")  # parked: files stay, pages are not built

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
