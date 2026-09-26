"""The documentation site's directives, and whether the pages still name things that exist.

A full ``sphinx-build -W`` is the real gate and runs in CI. What is asserted here is the failure
that a docs build would only catch *after* someone renamed a table: every ``vdjdb-schema`` and
``vdjdb-vocabulary`` directive in the tree naming a table or vocabulary the registry still has.
The directives exist so the tables cannot go stale; this makes their *arguments* unable to go stale
too.
"""
from __future__ import annotations

import importlib.util
import re
from importlib import import_module
from pathlib import Path

import pytest

DOCS = Path("docs")
EXT = DOCS / "_ext/vdjdb_schema.py"


def _ext():
    spec = importlib.util.spec_from_file_location("vdjdb_schema_ext", EXT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _directive_args(name: str) -> list[tuple[Path, str]]:
    """Every ``(page, argument)`` for a MyST fence invoking ``name``."""
    pat = re.compile(rf"^```\{{{name}\}}\s+(\S+)", re.M)
    return [(p, m) for p in sorted(DOCS.rglob("*.md"))
            for m in pat.findall(p.read_text())]


def test_extension_imports_and_registers_three_directives():
    mod = _ext()
    assert {"VdjdbSchema", "VdjdbVocabulary", "VdjdbScoreRules"} <= set(dir(mod))
    assert mod.DEFAULT_COLUMNS == ("name", "title", "comment")


def test_every_schema_directive_names_a_declared_table():
    registry = import_module("vdjdb.schema.fields")
    used = _directive_args("vdjdb-schema")
    assert used, "no `vdjdb-schema` directives found -- the column reference lost its tables"
    for page, table in used:
        assert table in registry.TABLES, f"{page} names unknown table {table!r}"


def test_every_vocabulary_directive_names_a_known_vocabulary():
    known = {"species", "tables", "airr"}
    used = _directive_args("vdjdb-vocabulary")
    assert used
    for page, name in used:
        assert name in known, f"{page} names unknown vocabulary {name!r}"


def test_every_declared_table_is_documented_somewhere():
    """A new table that no page renders is a table nobody can look up."""
    registry = import_module("vdjdb.schema.fields")
    documented = {t for _, t in _directive_args("vdjdb-schema")}
    missing = set(registry.TABLES) - documented
    assert not missing, f"declared but not documented: {sorted(missing)}"


@pytest.mark.parametrize("page", ["index.rst", "getting-started.md", "builds.md",
                                  "dashboard.md", "submission.md", "outputs.md",
                                  "denoising.md", "clustering.md"])
def test_toctree_pages_exist(page):
    assert (DOCS / page).exists()
