"""Sphinx directives that read the field registry at doc-build time.

    .. vdjdb-schema:: vdjdb
       :columns: name, title, comment

    .. vdjdb-vocabulary:: species

    .. vdjdb-score-rules::

**There is no generated ``.rst`` in the tree, and that is the whole point.** A column table written
out by a script is correct until someone changes a column and forgets to re-run it; a table built
from :mod:`vdjdb.schema` while the page renders cannot be stale, because there is nothing to
refresh. The README and the code have already disagreed about the schema once -- the README's
complex-columns table omits ``d.beta`` and lists ``submitter`` -- which is the drift this prevents.

The directives raise rather than emitting an empty table when a name is unknown, so a typo in a
table or vocabulary name fails ``sphinx-build -W`` instead of publishing a blank.
"""
from __future__ import annotations

from importlib import import_module
from typing import ClassVar

from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.util.docutils import SphinxDirective

#: The :class:`~vdjdb.schema.Field` attributes a page may ask for, and their column headings.
ATTRS = {
    "name": "Column",
    "ships_as": "Column",
    "title": "Title",
    "comment": "Description",
    "type": "Type",
    "data_type": "Data type",
    "airr": "AIRR field",
    "visible": "Visible",
    "searchable": "Searchable",
    "autocomplete": "Autocomplete",
}
DEFAULT_COLUMNS = ("name", "title", "comment")


def _registry():
    """The field-registry *module*.

    ``vdjdb.schema`` re-exports a function named ``fields`` alongside the submodule of the same
    name, so ``from vdjdb.schema import fields`` binds the function. Import the submodule by path.
    """
    return import_module("vdjdb.schema.fields")


def _table(headings: list[str], rows: list[list[str]], widths: list[int]) -> nodes.table:
    table = nodes.table(classes=["vdjdb-schema"])
    group = nodes.tgroup(cols=len(headings))
    table += group
    for w in widths:
        group += nodes.colspec(colwidth=w)
    head = nodes.thead()
    group += head
    head += _row(headings)
    body = nodes.tbody()
    group += body
    for r in rows:
        body += _row(r)
    return table


def _inline(text: str) -> nodes.paragraph:
    """A paragraph with rst ``literal`` spans turned into real literal nodes.

    Built directly rather than handed to a parser: these strings are rst docstrings, but the pages
    that use this directive are MyST, so ``parse_text_to_nodes`` would run the Markdown parser over
    rst markup and publish the backticks verbatim -- measured, it does.
    """
    para = nodes.paragraph()
    for n, part in enumerate(text.split("``")):
        if not part:
            continue
        para += nodes.literal("", part) if n % 2 else nodes.Text(part)
    return para


def _row(cells: list[str]) -> nodes.row:
    row = nodes.row()
    for text in cells:
        entry = nodes.entry()
        # A paragraph rather than raw text: docutils writers expect block-level content in a cell,
        # and `literal` keeps column names in a monospace face without a role per cell.
        entry += nodes.paragraph("", "", nodes.Text(text))
        row += entry
    return row


class VdjdbSchema(SphinxDirective):
    """One table of the field registry, in its positional order."""

    required_arguments = 1
    optional_arguments = 0
    option_spec: ClassVar[dict] = {"columns": directives.unchanged}

    def run(self) -> list[nodes.Node]:
        registry = _registry()

        table = self.arguments[0].strip()
        if table not in registry.TABLES:
            raise self.error(f"unknown table {table!r}; known: "
                             f"{', '.join(sorted(registry.TABLES))}")
        cols = [c.strip() for c in self.options.get("columns", "").split(",") if c.strip()]
        cols = cols or list(DEFAULT_COLUMNS)
        if bad := [c for c in cols if c not in ATTRS]:
            raise self.error(f"unknown column(s) {bad}; known: {', '.join(ATTRS)}")
        rows = [[str(getattr(f, c) or "") for c in cols] for f in registry.fields(table)]
        widths = [40 if c == "comment" else 15 for c in cols]
        return [_table([ATTRS[c] for c in cols], rows, widths)]


class VdjdbVocabulary(SphinxDirective):
    """A controlled vocabulary the build enforces, as a bullet list."""

    required_arguments = 1

    #: Name -> a callable returning the sorted members. Only vocabularies the *code* enforces; a
    #: list maintained only for the docs would be exactly the drift this module exists to stop.
    def _sources(self) -> dict[str, list[str]]:
        registry = _registry()
        return {
            "species": sorted(registry.SPECIES),
            "tables": sorted(registry.TABLES),
            "airr": sorted(f"{k} -> {v}" for k, v in registry.AIRR_MAP.items()),
        }

    def run(self) -> list[nodes.Node]:
        name = self.arguments[0].strip()
        sources = self._sources()
        if name not in sources:
            raise self.error(f"unknown vocabulary {name!r}; known: {', '.join(sorted(sources))}")
        lst = nodes.bullet_list()
        for item in sources[name]:
            para = nodes.paragraph()
            para += nodes.literal("", item)
            entry = nodes.list_item()
            entry += para
            lst += entry
        return [lst]


class VdjdbScoreRules(SphinxDirective):
    """The confidence score, from the scoring module's own docstrings.

    Rendered from source rather than restated, for the same reason as the schema tables: the score
    is defined in :mod:`vdjdb.score.confidence` and a prose copy of it drifts.
    """

    required_arguments = 0

    def run(self) -> list[nodes.Node]:
        from vdjdb.score import confidence

        out: list[nodes.Node] = []
        dl = nodes.definition_list()
        for name in ("frequency", "cell_count", "sequencing_score", "row_score"):
            fn = getattr(confidence, name)
            item = nodes.definition_list_item()
            term = nodes.term()
            term += nodes.literal("", f"{name}()")
            item += term
            defn = nodes.definition()
            summary = " ".join((fn.__doc__ or "").strip().split("\n\n")[0].split())
            defn += _inline(summary)
            item += defn
            dl += item
        out.append(dl)
        sig = nodes.paragraph()
        sig += nodes.Text("The score is a maximum over the sample signature ")
        sig += nodes.literal("", ", ".join(confidence.SCORE_SIGNATURE))
        sig += nodes.Text(", so the same clonotype assayed twice takes the better of the two.")
        out.append(sig)
        return out


def setup(app):
    app.add_directive("vdjdb-schema", VdjdbSchema)
    app.add_directive("vdjdb-vocabulary", VdjdbVocabulary)
    app.add_directive("vdjdb-score-rules", VdjdbScoreRules)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
