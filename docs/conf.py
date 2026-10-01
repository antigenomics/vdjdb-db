"""Sphinx configuration for the VDJdb specification site.

Served as a project page under the org custom domain, so it lives at
``docs.isalgo.dev/vdjdb-db/``. GitHub Pages is already provisioned for this repository with
``build_type: workflow``; nothing here has to create it.

**The schema tables are directives, not generated files.** ``docs/_ext/vdjdb_schema.py`` imports
:mod:`vdjdb.schema` at doc-build time, so there is no generated ``.rst`` in the tree that could go
stale between a schema change and someone remembering to re-run a script. That is the same reason
``docs/tuning/report.py`` generates the clustering tables rather than anyone transcribing them.
"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE / "_ext"))
sys.path.insert(0, str(HERE.parent / "src"))

project = "VDJdb"
author = "Mikhail Shugay"
copyright = "2026, Mikhail Shugay"
# Parsed from the package source rather than written out: arda's docs said 2.10.0 through twelve
# releases because the literal was never touched.
release = re.search(r'__version__ = "([^"]+)"',
                    (HERE.parent / "src/vdjdb/__init__.py").read_text()).group(1)
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "myst_parser",
    "vdjdb_schema",
]

# `docs/denoising.md` and `docs/clustering.md` are normative and carry `$...$` maths. They are
# included as Markdown rather than rewritten, so the specification has one copy, not two.
myst_enable_extensions = ["dollarmath", "colon_fence", "deflist"]
myst_heading_anchors = 3

napoleon_google_docstring = True
napoleon_numpy_docstring = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
# The docs job installs sphinx and polars, not the build's whole dependency tree.
autodoc_mock_imports = ["typer", "arda", "vdjtools", "mir", "sklearn", "igraph", "pandas",
                        "seqtree", "numpy", "clustereval"]

# Links into the Python docs are a convenience, and the build runs with `-W`, so an unreachable
# docs.python.org (a 503 on 2026-10-01 failed two docs runs and skipped the Pages deploy) must not
# fail it. Probe once with a short timeout; where it answers the mapping is used as before, and where
# it does not the mapping is dropped and the build says so, which loses only the links to Python's
# own documentation. Nothing in the specification is read from the inventory.
_PYTHON_INV = "https://docs.python.org/3/objects.inv"


def _reachable(url: str, timeout: float = 10.0) -> bool:
    import urllib.request

    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            return response.status == 200
    except OSError:
        return False


if _reachable(_PYTHON_INV):
    intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}
else:
    intersphinx_mapping = {}
    print(f"docs/conf.py: {_PYTHON_INV} is unreachable, building without links to the Python "
          "documentation")

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "tuning"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

html_baseurl = "https://docs.isalgo.dev/vdjdb-db/"
html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_title = f"VDJdb {release}"
html_theme_options = {
    "github_url": "https://github.com/antigenomics/vdjdb-db",
    "show_prev_next": False,
    "show_nav_level": 2,
    "navigation_depth": 3,
    "collapse_navigation": False,
    "header_links_before_dropdown": 4,
    "show_toc_level": 2,
}
# The stock `sidebar-nav-bs` renders only the children of the current top-level page; every page
# here is a sibling, so it would emit an empty heading. `site-nav.html` renders the whole tree.
html_sidebars = {"**": ["site-nav"], "index": []}
