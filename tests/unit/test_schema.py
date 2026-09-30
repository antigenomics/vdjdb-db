"""The field registry against its three sources of truth.

``attic/BuildDatabase.groovy`` is parsed **at test time**, never copied: it is the one correct
specification for the metadata and for the header-from-metadata derivation, and a copy would drift
exactly the way the nine duplicated column lists did.
"""
from __future__ import annotations

import re

import pytest

from vdjdb.config import repo_root
from vdjdb.schema import (
    ALL_COLUMNS,
    CHUNK_DEDUP_KEY,
    CLUSTER_MEMBERS_COLUMNS,
    EVIDENCE_COLUMNS,
    FIELDS,
    FULL_COLUMNS,
    MOTIF_PWMS_COLUMNS,
    SLIM_COLUMNS,
    TABLES,
    VDJDB_COLUMNS,
    VDJDB_WEB_COLUMNS,
    fields,
    header,
    render_meta,
    render_slim_meta,
)

GROOVY = repo_root() / "attic" / "BuildDatabase.groovy"


def _groovy_list(name: str) -> list[str]:
    """The string literals of a Groovy list assignment, in order."""
    src = GROOVY.read_text()
    start = src.index(f"{name} = [")
    depth, i = 0, src.index("[", start)
    for end in range(i, len(src)):
        depth += (src[end] == "[") - (src[end] == "]")
        if depth == 0:
            break
    body = src[i:end + 1]
    return [re.sub(r'\\t', "\t", m) for m in re.findall(r'"((?:[^"\\]|\\.)*)"', body)]


# --------------------------------------------------------------------------------------------
# Against the Groovy
# --------------------------------------------------------------------------------------------

def test_groovy_parser_finds_the_metadata() -> None:
    assert len(_groovy_list("METADATA_LINES")) == 22          # header + 21 rows
    assert len(_groovy_list("SLIM_METADATA_LINES")) == 17     # header + 16 rows


def test_chunk_columns_match_the_groovy_exactly() -> None:
    """``ALL_COLS`` is the chunk contract. It must not move without a migration."""
    groovy = _groovy_list("COMPLEX_COLUMNS") + _groovy_list("METHOD_COLUMNS") + _groovy_list("META_COLUMNS")
    assert list(ALL_COLUMNS) == groovy


def test_dedup_key_matches_the_groovy_signature_cols() -> None:
    assert list(CHUNK_DEDUP_KEY) == _groovy_list("SIGNATURE_COLS")


def test_metadata_differs_from_the_groovy_only_in_the_declared_ways() -> None:
    """Every difference from the shipped metadata is one of the three declared fixes.

    This is the test that makes ``rules/expected_diffs.toml``'s metadata rules honest: if the
    registry starts differing in some *other* way, this fails rather than the comparison silently
    absorbing it.
    """
    groovy = _groovy_list("METADATA_LINES")
    shipped = [line.split("\t")[0] for line in groovy[1:]]
    generated = list(VDJDB_COLUMNS)

    # 1. TCR_hash was missing entirely.
    assert set(generated) - set(shipped) == {"TCR_hash"}
    assert not set(shipped) - set(generated)
    # 2. vdjdb.score sat after cdr3fix; the data has it before method.
    assert shipped.index("vdjdb.score") > shipped.index("cdr3fix")
    assert generated.index("vdjdb.score") < generated.index("method")
    # 3. everything else keeps its order.
    assert [c for c in generated if c not in {"vdjdb.score", "TCR_hash"}] == \
           [c for c in shipped if c != "vdjdb.score"]


def test_the_four_web_rows_are_the_only_attribute_fix() -> None:
    """Only the ``web.*`` rows change attributes; every other row is reproduced value-for-value."""
    groovy = {line.split("\t")[0]: line for line in _groovy_list("METADATA_LINES")[1:]}
    changed = {f.name for f in fields("vdjdb")
               if f.name in groovy and f.meta_row() != groovy[f.name]}
    # vdjdb.score's title was updated to "Info" in production; the web.* rows are the field fix.
    assert changed == {"web.method", "web.method.seq", "web.cdr3fix.nc", "web.cdr3fix.unmp",
                       "vdjdb.score"}


def test_web_rows_get_a_real_data_type() -> None:
    """The shipped rows put ``0`` in ``data.type`` and ``factor`` in ``title``."""
    for name in ("web.method", "web.method.seq", "web.cdr3fix.nc", "web.cdr3fix.unmp"):
        f = FIELDS[name]
        assert f.data_type == "factor"
        assert f.visible == 0, "web.* fields are internal and must never be displayed"


def test_slim_metadata_gains_tcr_hash_and_fixes_the_geometry_order() -> None:
    shipped = [line.split("\t")[0] for line in _groovy_list("SLIM_METADATA_LINES")[1:]]
    assert set(SLIM_COLUMNS) - set(shipped) == {"TCR_hash"}
    # shipped has v.end before j.start, two positions too early; the data has j.start then v.end last
    assert shipped.index("v.end") < shipped.index("j.start")
    assert SLIM_COLUMNS.index("j.start") < SLIM_COLUMNS.index("v.end")
    assert SLIM_COLUMNS[-2:] == ("j.start", "v.end")


# --------------------------------------------------------------------------------------------
# Internal consistency -- the invariants vdjdb-web depends on
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("table", ["vdjdb", "vdjdb-web", "slim", "full",
                                   "cluster_members", "motif_pwms"])
def test_every_column_is_declared(table: str) -> None:
    for name in TABLES[table]:
        assert name in FIELDS, f"{table}: {name} has no Field declaration"


@pytest.mark.parametrize("table", sorted(TABLES))
def test_no_table_repeats_a_column(table: str) -> None:
    names = TABLES[table]
    assert len(set(names)) == len(names)


def test_header_is_the_name_column_of_the_metadata() -> None:
    """``BuildDatabase.groovy:411``'s invariant, restored. Its loss is why the files drifted."""
    for table in ("vdjdb", "vdjdb-web", "slim"):
        meta = render_meta(table).splitlines()
        assert header(table).split("\t") == [row.split("\t")[0] for row in meta[1:]]


def test_metadata_rows_are_rectangular() -> None:
    """``vdjdb-web`` splits on tab and indexes by position; a short row mistypes the column."""
    lines = render_meta("vdjdb-web").splitlines()
    assert all(len(line.split("\t")) == 8 for line in lines), "every meta row has 8 fields"
    assert len(lines) == 1 + len(VDJDB_WEB_COLUMNS)


def test_slim_metadata_is_two_columns() -> None:
    lines = render_slim_meta().splitlines()
    assert all(len(line.split("\t")) == 2 for line in lines)
    assert len(lines) == 1 + len(SLIM_COLUMNS)


def test_rendered_metadata_ends_with_a_newline() -> None:
    assert render_meta().endswith("\n")
    assert render_slim_meta().endswith("\n")


def test_positional_contract_widths() -> None:
    """The widths ``Motifs.scala`` hands Tablesaw as a fixed ``Array[ColumnType]``."""
    assert len(VDJDB_COLUMNS) == 22
    assert len(VDJDB_WEB_COLUMNS) == 27
    assert len(SLIM_COLUMNS) == 17
    assert len(FULL_COLUMNS) == 35
    assert len(CLUSTER_MEMBERS_COLUMNS) == 19
    assert len(MOTIF_PWMS_COLUMNS) == 27


def test_full_is_the_chunk_columns_plus_four_derived() -> None:
    assert FULL_COLUMNS[:len(ALL_COLUMNS)] == ALL_COLUMNS
    assert FULL_COLUMNS[len(ALL_COLUMNS):] == ("cdr3fix.alpha", "cdr3fix.beta",
                                               "vdjdb.score", "TCR_hash")


def test_web_table_is_the_release_table_plus_evidence() -> None:
    assert VDJDB_WEB_COLUMNS[:len(VDJDB_COLUMNS)] == VDJDB_COLUMNS
    assert VDJDB_WEB_COLUMNS[len(VDJDB_COLUMNS):] == EVIDENCE_COLUMNS


def test_no_metadata_value_contains_a_tab_or_newline() -> None:
    """A tab inside a comment would shift every later field of that row."""
    for f in FIELDS.values():
        for value in (f.name, f.type, f.data_type, f.title, f.comment):
            assert "\t" not in value and "\n" not in value, f"{f.name}: {value!r}"


def test_unknown_table_names_its_alternatives() -> None:
    with pytest.raises(KeyError, match="unknown table"):
        fields("nope")


def test_the_submission_template_is_a_projection_of_the_registry():
    """`template.tsv` replaced a 2016 `.xls` that nobody could regenerate.

    That spreadsheet was one of the nine places this module's docstring lists as having already
    drifted, and it had: its example rows wrote `antigen.species = HIV` where the vocabulary and
    that paper's own chunk now say `HIV-1`. Generated from `TABLES["chunk"]`, it cannot describe a
    column set the build does not read.
    """
    from vdjdb.schema import TABLES, header

    lines = (repo_root() / "template.tsv").read_text().rstrip("\n").split("\n")
    assert lines[0] == header("chunk"), (
        "template.tsv is stale; regenerate with `vdjdb schema --table chunk --format header`")
    assert len(lines) > 1, "a template with no example row teaches nothing"
    for i, line in enumerate(lines[1:], start=1):
        assert len(line.split("\t")) == len(TABLES["chunk"]), f"example row {i} is ragged"


def test_the_submission_template_passes_the_checks_a_submission_must_pass():
    """A template that would fail QC is worse than none: a submitter copies its shape."""
    from vdjdb.qc.lint import lint_file

    assert lint_file(repo_root() / "template.tsv") == []


def test_chunk_columns_are_an_order_and_not_a_gate():
    """`qc/lint.py` checks column *membership*, and must keep doing so.

    A submitted chunk may order its columns however it likes and may omit the optional ones, so
    declaring an order for the template must not turn that order into a requirement.
    """
    from vdjdb.schema import ALL_COLUMNS, CHUNK_COLUMNS, KEPT_CURATION_COLUMNS

    assert set(CHUNK_COLUMNS) <= set(ALL_COLUMNS) | set(KEPT_CURATION_COLUMNS)
    # the three curation columns a chunk may carry and the template does not name
    assert set(KEPT_CURATION_COLUMNS) - set(CHUNK_COLUMNS) == {
        "submitter", "comment", "method.pairing"}
