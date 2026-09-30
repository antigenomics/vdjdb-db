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

#: Every directory a build writes under `out/`, and the command that writes it, taken from the
#: `build.yml` step list rather than from :data:`manifest.BUNDLES`. The `build_dir` fixture below
#: creates whatever the declaration names, so it can only check that the code agrees with itself:
#: `PRIMARY` named `motifs_new/`, the output directory of the sweeps under `docs/tuning/`, and the
#: fixture dutifully created it, so nothing failed until `vdjdb release` was run against a pipeline
#: build and raised on two required members.
WRITTEN_BY = {
    "tables": "vdjdb build --out out",
    "legacy": "vdjdb build, or vdjdb make legacy",
    "airr": "vdjdb build, or vdjdb convert airr",
    "motifs": "vdjdb motifs",
    "summary": "vdjdb summary, copied by the Dashboard step",
    "": "bundle.stage(): LICENSE and latest-version.txt, from the repository root",
}


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


def test_every_member_comes_from_a_directory_the_pipeline_writes():
    for b in manifest.BUNDLES:
        for m in b.members:
            d = str(Path(m.source).parent).replace(".", "")
            assert d in WRITTEN_BY, (
                f"bundle {b.role!r} member {m.source!r} is under {d!r}/, which no command writes. "
                f"The directories a build produces are {sorted(k for k in WRITTEN_BY if k)}.")


def test_the_motif_members_name_the_directory_the_motifs_command_writes_to():
    """One default, two readers. `vdjdb motifs --out` decides where the files are; the bundles and
    `vdjdb motif-metrics` decide where they are looked for, and all three have to agree."""
    import inspect

    from vdjdb.cli import app

    defaults = {}
    for c in app.registered_commands:
        fn = c.callback
        name = c.name or fn.__name__.replace("_", "-")
        defaults[name] = {k: v.default.default for k, v in inspect.signature(fn).parameters.items()
                          if hasattr(v.default, "default")}

    written = Path(defaults["motifs"]["out"])
    assert written == Path("out/motifs"), (
        f"`vdjdb motifs --out` defaults to {written}, and the bundles read out/motifs")
    assert Path(defaults["motif-metrics"]["motifs"]) == written, (
        "`vdjdb motif-metrics --motifs` must default to the directory `vdjdb motifs` writes")
    for b in manifest.BUNDLES:
        for m in b.members:
            if Path(m.source).name.startswith(("cluster_members", "motif_pwms")):
                assert Path(m.source).parent == written.relative_to("out"), (
                    f"bundle {b.role!r} reads {m.source!r}; `vdjdb motifs` writes to {written}")


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


def test_a_dry_run_leaves_the_tracked_latest_version_alone(tmp_path) -> None:
    """#707. `latest-version.txt` is tracked, so a made-up tag must not reach the repository.

    That is the defect `ROADMAP.md` §3.2 records as having already shipped once, from the other
    direction: for several releases line 1 named the *previous* release. A dry run leaving the tree
    naming a release that does not exist is one `git commit -a` away from the same outcome.
    """
    from vdjdb.release.bundle import legacy_url, prepare_latest

    path = tmp_path / "latest-version.txt"
    before = "https://example.invalid/old.zip\n"
    path.write_text(before)

    content = prepare_latest("v2026.09.9", path, write=False)
    assert path.read_text() == before, "a dry run must not touch the file"
    assert content.splitlines()[0] == legacy_url("v2026.09.9"), (
        "and must still return the line the bundle ships, so the two have one source")

    prepare_latest("v2026.09.9", path)
    assert path.read_text() == content, "a real run writes exactly what the dry run computed"
