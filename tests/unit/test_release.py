"""Bundle assembly: what ships, what it is called, and whether it can be built twice.

The legacy release script gathered a bundle with ``cp *.txt``, which would have swept in seven side
tables the shipped zip does not contain -- so the script that produced the release was not the one
that defined it. These tests assert the declaration instead: every member named, the basenames
``vdjmatch`` looks for present, and the same inputs producing the same bytes.
"""
from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pytest

from vdjdb.release import bundle, manifest

#: `vdjmatch` looks inside a bundle for these by BASENAME, not by path.
VDJMATCH_NEEDS = {"vdjdb.txt", "vdjdb.slim.txt", "vdjdb_full.txt"}


@pytest.fixture
def build_dir(tmp_path):
    """A build directory with every required member of every bundle, and nothing else."""
    root = tmp_path / "out"
    for b in manifest.BUNDLES:
        for m in b.members:
            if m.optional:
                continue
            p = root / m.source
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f"contents of {m.source}\n")
    return root


def test_tag_scheme_is_enforced():
    assert bundle.version_of("v2026.09.1") == "2026.09.1"
    for bad in ("2026-06-03", "v2026.9.1", "v2026.09", "2026.09.1", "v2026.09.1-ZENODO"):
        with pytest.raises(ValueError, match="not v<YYYY>"):
            bundle.version_of(bad)


def test_latest_version_is_idempotent_and_keeps_history(tmp_path):
    p = tmp_path / "latest-version.txt"
    p.write_text("https://example.invalid/old.zip\n")
    first = bundle.prepare_latest("v2026.09.1", p)
    again = bundle.prepare_latest("v2026.09.1", p)
    assert first == again, "re-running a failed release must not push the URL twice"
    lines = p.read_text().splitlines()
    assert lines[0] == bundle.legacy_url("v2026.09.1")
    assert "https://example.invalid/old.zip" in lines, "history must be kept, not replaced"
    assert lines.count(lines[0]) == 1


def test_line_one_names_the_legacy_bundle_of_this_tag():
    # Clients in the wild download line 1 verbatim and expect the legacy layout inside it.
    url = bundle.legacy_url("v2026.09.1")
    assert url.endswith("/v2026.09.1/vdjdb-legacy-2026.09.1.zip")


def test_a_missing_required_member_fails_loudly(build_dir):
    (build_dir / "legacy/vdjdb.txt").unlink()
    with pytest.raises(FileNotFoundError, match=r"vdjdb\.txt"):
        manifest.resolve(manifest.LEGACY, build_dir, "2026.09.1")


def test_a_missing_optional_member_does_not(build_dir):
    members = manifest.resolve(manifest.LEGACY, build_dir, "2026.09.1")
    names = {Path(n).name for _, n in members}
    assert "cluster_members_tcremp.txt" not in names
    assert names >= VDJMATCH_NEEDS


def test_bundles_are_byte_reproducible(build_dir, tmp_path):
    a = bundle.build(build_dir, tmp_path / "a", "v2026.09.1")
    b = bundle.build(build_dir, tmp_path / "b", "v2026.09.1")
    assert [e["sha256"] for e in a["bundles"]] == [e["sha256"] for e in b["bundles"]]


def test_manifest_names_every_bundle_by_role(build_dir, tmp_path):
    out = tmp_path / "rel"
    bundle.build(build_dir, out, "v2026.09.1")
    m = json.loads((out / "manifest.json").read_text())
    assert {b["role"] for b in m["bundles"]} == {"primary", "legacy", "airr"}
    assert m["tag"] == "v2026.09.1" and m["version"] == "2026.09.1"
    for entry in m["bundles"]:
        assert (out / entry["file"]).stat().st_size == entry["bytes"]


def test_nothing_undeclared_reaches_a_bundle(build_dir, tmp_path):
    # The `cp *.txt` failure mode, asserted directly.
    (build_dir / "legacy/vdjdb_full_filtered.txt").write_text("a side table\n")
    (build_dir / "legacy/vdjdb_broken.txt").write_text("another\n")
    out = tmp_path / "rel"
    bundle.build(build_dir, out, "v2026.09.1")
    with zipfile.ZipFile(out / "vdjdb-legacy-2026.09.1.zip") as zf:
        names = {Path(n).name for n in zf.namelist()}
    assert "vdjdb_full_filtered.txt" not in names
    assert "vdjdb_broken.txt" not in names


def test_checksums_are_the_format_sha256sum_reads(build_dir, tmp_path):
    out = tmp_path / "rel"
    bundle.build(build_dir, out, "v2026.09.1")
    for line in (out / "SHA256SUMS").read_text().splitlines():
        digest, _, name = line.partition("  ")
        assert len(digest) == 64 and (out / name).exists()


def _slim(rows: list[tuple[str, str, str]]) -> str:
    """A minimal `vdjdb.slim.txt`: (cdr3, epitope, comma-joined reference.id)."""
    head = "cdr3\tantigen.epitope\treference.id\n"
    return head + "".join(f"{c}\t{e}\t{r}\n" for c, e, r in rows)


def test_changelog_reports_a_renamed_study_as_a_rename(tmp_path):
    """#347 replaced DOI and preprint URLs with PubMed ids.

    Without content matching the notes read "3 studies withdrawn, 3 added", which is alarming and
    wrong -- measured against the real 2026-06-03 release, that is exactly what it said.
    """
    from vdjdb.release import changelog as cl

    before, after = tmp_path / "before", tmp_path / "after"
    for d in (before, after):
        d.mkdir()
    (before / "vdjdb.slim.txt").write_text(_slim([
        ("CASSA", "GILGFVFTL", "https://doi.org/10.1/x"),
        ("CASSB", "GILGFVFTL", "https://doi.org/10.1/x"),
        ("CASSC", "NLVPMVATV", "PMID:1"),
    ]))
    (after / "vdjdb.slim.txt").write_text(_slim([
        ("CASSA", "GILGFVFTL", "PMID:99"),          # same rows, new identifier
        ("CASSB", "GILGFVFTL", "PMID:99"),
        ("CASSC", "NLVPMVATV", "PMID:1"),
        ("CASSD", "KLGGALQAK", "PMID:2"),           # a genuinely new study
    ]))
    d = cl.diff(before, after)
    assert d["renamed"].to_dicts() == [
        {"from": "https://doi.org/10.1/x", "to": "PMID:99", "rows": 2}]
    assert d["added"]["reference.id"].to_list() == ["PMID:2"]
    assert d["removed"].is_empty()
    notes = cl.render(d, tag="v2026.09.1")
    assert "1 studies re-identified" in notes
    assert "no longer present" not in notes


def test_changelog_splits_comma_joined_references(tmp_path):
    """One slim row can come from several studies; counting the joined string invents a study."""
    from vdjdb.release import changelog as cl

    before, after = tmp_path / "b", tmp_path / "a"
    for d in (before, after):
        d.mkdir()
    (before / "vdjdb.slim.txt").write_text(_slim([("CASSA", "GIL", "PMID:1")]))
    (after / "vdjdb.slim.txt").write_text(_slim([("CASSA", "GIL", "PMID:1,PMID:2")]))
    d = cl.diff(before, after)
    assert d["references_after"] == 2
    assert d["added"]["reference.id"].to_list() == ["PMID:2"]
