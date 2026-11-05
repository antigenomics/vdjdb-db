"""Generated text tables obey the byte limit without changing their contents."""
import gzip
import zipfile
from pathlib import Path

import polars as pl
import pytest

from vdjdb.io import text


@pytest.mark.parametrize("size,compressed", [(999_999, False), (1_000_000, False),
                                           (1_000_001, True)])
def test_exact_byte_limit_and_roundtrip(tmp_path, size, compressed):
    path = tmp_path / "table.tsv"
    # Two header bytes, value bytes, and the terminating newline.
    frame = pl.DataFrame({"x": ["a" * (size - 3)]})
    result = text.write_csv(frame, path, separator="\t")
    assert result.name == ("table.tsv.gz" if compressed else "table.tsv")
    assert text.read_csv(path, separator="\t").equals(frame)
    assert not (path if compressed else path.with_suffix(".tsv.gz")).exists()


def test_rebuild_switches_formats_and_gzip_is_deterministic(tmp_path):
    path = tmp_path / "table.txt"
    large = pl.DataFrame({"x": ["é" * 600_000]})
    compressed = text.write_csv(large, path, separator="\t")
    first = compressed.read_bytes()
    assert gzip.decompress(first).decode() == large.write_csv(separator="\t")
    assert text.write_csv(large, path, separator="\t").read_bytes() == first
    small = pl.DataFrame({"x": ["a"]})
    assert text.write_csv(small, path, separator="\t") == path
    assert not compressed.exists()
    text.write_csv(large, path, separator="\t")
    assert not path.exists()


def test_failed_write_preserves_previous_output(tmp_path):
    path = tmp_path / "table.tsv"
    path.write_text("previous\n")
    with pytest.raises(RuntimeError), text.table_writer(path) as stream:
        stream.write("partial")
        raise RuntimeError("failed")
    assert path.read_text() == "previous\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["table.tsv"]


def test_airr_and_motif_exports_and_readers(tmp_path, monkeypatch):
    from vdjdb.emit import airr
    from vdjdb.motifs import emit
    from vdjdb.validate.motif_bench import read_members

    monkeypatch.setattr(text, "MAX_TEXT_BYTES", 1)
    frame = pl.DataFrame({"x": ["value"]})
    written = airr.write_all({"receptor": frame}, tmp_path / "airr")
    assert "vdjdb.receptor.tsv.gz" in written
    assert text.read_csv(tmp_path / "airr/vdjdb.receptor.tsv", separator="\t").equals(frame)
    members = pl.DataFrame({c: ["value"] for c in emit.MEMBER_COLUMNS})
    pwms = pl.DataFrame({c: ["value"] for c in emit.PWM_COLUMNS})
    emit.write(members, pwms, tmp_path / "motifs")
    assert read_members(tmp_path / "motifs/cluster_members.txt").equals(members)
    assert (tmp_path / "motifs/motif_pwms.txt.gz").is_file()
    assert not list(tmp_path.rglob("*.txt"))


def test_corpus_text_compresses_and_parquet_queries_still_roundtrip(tmp_path, monkeypatch):
    from vdjdb.corpus import build

    monkeypatch.setattr(text, "MAX_TEXT_BYTES", 1)
    frame = pl.DataFrame({"x": ["value"]})
    corpus = {name: frame for name in build.TABLES}
    build.write(corpus, tmp_path)
    assert set(build.read(tmp_path)) == set(corpus)
    for name in corpus:
        assert text.read_csv(tmp_path / f"{name}.tsv", separator="\t").equals(frame)
        assert not (tmp_path / f"{name}.tsv").exists()


def test_release_streams_gzip_sources_and_keeps_member_names(tmp_path):
    from vdjdb.compare.diff import Bundle
    from vdjdb.release import manifest

    frame = pl.DataFrame({"x": ["a" * 1_000_001]})
    text.write_csv(frame, tmp_path / "legacy/table.txt", separator="\t")
    text.write_csv(frame, tmp_path / "tables/table.tsv.gz", separator="\t")
    bundle = manifest.Bundle("test", "test-{version}.zip", "test-{version}", (
        manifest.Member("legacy/table.txt"), manifest.Member("tables/table.tsv.gz")))
    first = manifest.write_zip(bundle, tmp_path, tmp_path, "one")
    repeated = manifest.write_zip(bundle, tmp_path, tmp_path / "repeat", "one")
    assert first.read_bytes() == repeated.read_bytes()
    expected = frame.write_csv(separator="\t").encode()
    with zipfile.ZipFile(first) as archive:
        assert set(archive.namelist()) == {"test-one/table.txt.gz", "test-one/table.tsv.gz"}
        for relative in ("legacy/table.txt.gz", "tables/table.tsv.gz"):
            member = "test-one/" + Path(relative).name
            assert archive.read(member) == (tmp_path / relative).read_bytes()
            assert gzip.decompress(archive.read(member)) == expected
    assert Bundle(tmp_path / "legacy").read_bytes("table.txt") == expected
    assert Bundle(first).read_bytes("table.txt") == expected
    assert Bundle(first).read_bytes("table.tsv") == expected
    assert not list(tmp_path.rglob("*.txt"))
