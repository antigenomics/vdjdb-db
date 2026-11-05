"""Offline prepared web/structure identity and archive-content checks."""
import gzip
import json
from pathlib import Path

import pytest

from vdjdb.release import bundle, manifest
from vdjdb.release.prepared import validate


def seal(web):
    structure = web / "structure"
    h = "a" * 64
    assets = dict(schema_version=1, web_table_sha256=manifest.sha256(web / "vdjdb.txt.gz"),
                  source_metadata_sha256="b" * 64, hashes=[h], files=[
        dict(path=p.name, bytes=p.stat().st_size, sha256=manifest.sha256(p))
        for p in sorted(structure.iterdir()) if p.name != "asset-manifest.json"])
    (structure / "asset-manifest.json").write_text(json.dumps(assets))
    descriptor = dict(schema_version=1, web_table_sha256=assets["web_table_sha256"],
        source_metadata_sha256=assets["source_metadata_sha256"],
        structures_metadata_sha256=manifest.sha256(web / "structures_metadata.tsv.gz"),
        structure_manifest_sha256=manifest.sha256(structure / "asset-manifest.json"))
    (web / "prepared-assets.json").write_text(json.dumps(descriptor))


@pytest.fixture
def prepared(tmp_path):
    web = tmp_path / "web"
    structure = web / "structure"
    structure.mkdir(parents=True)
    (web / "vdjdb.txt.gz").write_bytes(gzip.compress(("TCR_hash\n" + "a" * 64 + "\n").encode(), mtime=0))
    (web / "structures_metadata.tsv.gz").write_bytes(gzip.compress(("hash\n" + "a" * 64 + "\n").encode(), mtime=0))
    h = "a" * 64
    for name in (f"{h}.html", f"{h}_simplified.html", f"aligned_aligned_{h}.pdb",
                 f"{h}_contacts_aa.txt.gz", f"{h}_aa_coordinates.tsv.gz", f"{h}_all.zip"):
        (structure / name).write_bytes(b"served file")
    import zipfile
    with zipfile.ZipFile(structure / f"{h}_all.zip", "w") as archive:
        for name in (f"aligned_aligned_{h}.pdb", f"{h}_contacts_aa.txt.gz", f"{h}_aa_coordinates.tsv.gz"):
            archive.write(structure / name, name)
    seal(web)
    return web


def test_paired_bundles_include_seals_and_exact_structure_files(prepared, tmp_path):
    import zipfile
    web_bundle = manifest.Bundle("web", "web-{version}.zip", "db-{version}", (
        manifest.Member("web/vdjdb.txt"), manifest.Member("web/structures_metadata.tsv"),
        manifest.Member("web/prepared-assets.json")))
    result = bundle.build(tmp_path, tmp_path / "release", "v2026.11.0", bundles=(web_bundle,))
    assert [b["role"] for b in result["bundles"]] == ["web", "web-structures"]
    for entry in result["bundles"]:
        with zipfile.ZipFile(tmp_path / "release" / entry["file"]) as archive:
            assert set(archive.namelist()) == set(entry["members"])
            if entry["role"] == "web-structures":
                assert set(archive.namelist()) == {f"structure/{p.name}" for p in (prepared / "structure").iterdir()}
            else:
                assert "db-2026.11.0/prepared-assets.json" in archive.namelist()


@pytest.mark.parametrize("changed", ["vdjdb.txt.gz", "structures_metadata.tsv.gz",
                                     "structure/" + "a" * 64 + ".html"])
def test_mismatched_table_metadata_and_maps_fail(prepared, changed):
    (prepared / changed).write_bytes(b"stale")
    with pytest.raises(ValueError, match="digest mismatch"):
        validate(prepared)


def test_missing_required_asset_fails_even_after_resealing(prepared):
    (prepared / "structure" / ("a" * 64 + "_all.zip")).unlink()
    seal(prepared)
    with pytest.raises(ValueError, match="missing required"):
        validate(prepared)


def test_mismatched_source_metadata_identity_fails(prepared):
    descriptor = prepared / "prepared-assets.json"
    data = json.loads(descriptor.read_text())
    data["source_metadata_sha256"] = "c" * 64
    descriptor.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="input identity"):
        validate(prepared)


def test_table_only_allowed_but_unsealed_metrics_fail(tmp_path):
    assert not validate(tmp_path)
    (tmp_path / "structures_metadata.tsv").write_text("metrics")
    with pytest.raises(ValueError, match="missing or unsafe"):
        validate(tmp_path)


def test_prepared_input_archive_stages_safely_and_rejects_overlays(tmp_path):
    import zipfile
    from vdjdb.release.prepared import stage_archive

    source = tmp_path / "input.zip"
    with zipfile.ZipFile(source, "w") as zipped:
        zipped.writestr("web/prepared-assets.json", "{}")
    build = tmp_path / "out"
    stage_archive(source, build)
    assert (build / "web/prepared-assets.json").read_text() == "{}"
    with pytest.raises(ValueError, match="fresh"):
        stage_archive(source, build)
    with zipfile.ZipFile(source, "w") as zipped:
        zipped.writestr("web/structure/../../escape", "bad")
    with pytest.raises(ValueError, match="unsafe"):
        stage_archive(source, tmp_path / "other")
    assert not (tmp_path / "other").exists()


@pytest.mark.parametrize("missing", ["prepared-assets.json", "structure/asset-manifest.json",
                                      "structures_metadata.tsv.gz"])
def test_missing_prepared_sidecars_fail(prepared, missing):
    (prepared / missing).unlink()
    with pytest.raises(ValueError, match="missing"):
        validate(prepared)


def test_changed_structure_manifest_and_unsealed_maps_fail(prepared):
    asset_manifest = prepared / "structure/asset-manifest.json"
    asset_manifest.write_text(asset_manifest.read_text() + "\n")
    with pytest.raises(ValueError, match="digest mismatch"):
        validate(prepared)
    seal(prepared)
    (prepared / "structure/old-map.html").write_text("stale")
    with pytest.raises(ValueError, match="unsealed"):
        validate(prepared)


@pytest.mark.parametrize("counterexample", ["empty_hashes", "absent_table_hash", "absent_metrics_hash", "corrupt_zip"])
def test_resealed_invalid_preparation_is_rejected(prepared, counterexample):
    if counterexample.startswith("absent"):
        name, column = (("vdjdb.txt.gz", "TCR_hash") if counterexample == "absent_table_hash"
                        else ("structures_metadata.tsv.gz", "hash"))
        (prepared / name).write_bytes(gzip.compress((column + "\n" + "c" * 64 + "\n").encode(), mtime=0))
    elif counterexample == "corrupt_zip":
        (prepared / "structure" / ("a" * 64 + "_all.zip")).write_bytes(b"not a zip")
    seal(prepared)
    if counterexample == "empty_hashes":
        path = prepared / "structure/asset-manifest.json"
        data = json.loads(path.read_text())
        data["hashes"] = []
        path.write_text(json.dumps(data))
        path = prepared / "prepared-assets.json"
        data = json.loads(path.read_text())
        data["structure_manifest_sha256"] = manifest.sha256(prepared / "structure/asset-manifest.json")
        path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        validate(prepared)


def test_resealed_download_member_must_equal_served_file(prepared):
    (prepared / "structure" / ("a" * 64 + "_contacts_aa.txt.gz")).write_bytes(b"changed contacts")
    seal(prepared)
    with pytest.raises(ValueError, match="stale structure ZIP member"):
        validate(prepared)
