"""Verify the paired web table, prepared metrics and served structure files before release."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path, PurePosixPath

from ..io.text import MAX_TEXT_BYTES, open_text, resolve_path
from .manifest import sha256


def validate(web: Path) -> bool:
    """Return whether a complete sealed structure directory accompanies the fresh web table.

    Table/chart-only releases are allowed. Any supplied structure metrics or assets require
    both descriptors and all files; stale or incomplete preparation fails before packaging.
    """
    def table(name: str) -> Path:
        paths = [p for p in (web / name, web / (name + ".gz")) if p.is_file()]
        if len(paths) != 1:
            raise ValueError(f"missing or ambiguous prepared table: {name}[.gz]")
        return paths[0]

    metrics = resolve_path(web / "structures_metadata.tsv")
    descriptor = web / "prepared-assets.json"
    structure = web / "structure"
    if not (metrics.exists() or descriptor.exists() or structure.exists()):
        return False

    def read(path: Path) -> dict:
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"missing or unsafe prepared asset: {path}")
        try:
            data = json.loads(path.read_text())
        except (ValueError, OSError) as exc:
            raise ValueError(f"invalid prepared descriptor: {path}") from exc
        if not isinstance(data, dict) or data.get("schema_version") != 1:
            raise ValueError(f"unsupported prepared schema: {path}")
        return data

    def digest(path: Path, expected: str | None) -> None:
        if (not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected)
                or path.is_symlink() or not path.is_file() or sha256(path) != expected):
            raise ValueError(f"prepared asset digest mismatch or missing file: {path}")

    prepared = read(descriptor)
    asset_path = structure / "asset-manifest.json"
    assets = read(asset_path)
    web_table = table("vdjdb.txt")
    metrics = table("structures_metadata.tsv")
    digest(web_table, prepared.get("web_table_sha256"))
    digest(metrics, prepared.get("structures_metadata_sha256"))
    digest(asset_path, prepared.get("structure_manifest_sha256"))
    for key in ("web_table_sha256", "source_metadata_sha256"):
        value = prepared.get(key)
        if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value) or assets.get(key) != value:
            raise ValueError(f"prepared input identity mismatch: {key}")
    hashes = assets.get("hashes")
    if (not isinstance(hashes, list) or not hashes or any(not isinstance(h, str) or
            not re.fullmatch(r"[0-9a-f]{64}", h) for h in hashes)
            or hashes != sorted(set(hashes))):
        raise ValueError("invalid prepared structure hashes")
    for path, column, lower in ((web_table, "TCR_hash", True), (metrics, "hash", False)):
        with open_text(path) as stream:
            rows = csv.DictReader(stream, delimiter="\t")
            if column not in (rows.fieldnames or []):
                raise ValueError(f"missing prepared table column: {column}")
            present = {row[column].lower() if lower else row[column] for row in rows}
        if not set(hashes) <= present:
            raise ValueError(f"structure hashes absent from table/metrics: {path}")
    files = assets.get("files")
    if not isinstance(files, list):
        raise ValueError("invalid prepared structure files")
    names = set()
    for entry in files:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise ValueError("invalid prepared structure file entry")
        name = entry["path"]
        relative = PurePosixPath(name)
        if (relative.is_absolute() or len(relative.parts) != 1 or ".." in relative.parts or "\\" in name or
                relative.as_posix() != name or name in names or name == "asset-manifest.json"):
            raise ValueError(f"unsafe or duplicate prepared path: {name}")
        names.add(name)
        path = structure / name
        if structure.is_symlink() or any((structure / PurePosixPath(*relative.parts[:i])).is_symlink()
                                        for i in range(1, len(relative.parts) + 1)):
            raise ValueError(f"unsafe prepared path: {name}")
        digest(path, entry.get("sha256"))
        if path.suffix in {".txt", ".tsv"} and path.stat().st_size > MAX_TEXT_BYTES:
            raise ValueError(f"oversized plain prepared table: {name}")
        if type(entry.get("bytes")) is not int or path.stat().st_size != entry["bytes"]:
            raise ValueError(f"prepared asset size mismatch: {name}")
    actual = {p.relative_to(structure).as_posix() for p in structure.rglob("*") if p.is_file()}
    if actual != names | {"asset-manifest.json"}:
        raise ValueError("prepared structure directory contains missing or unsealed files")
    for h in hashes:
        required = {f"{h}.html", f"{h}_simplified.html", f"aligned_aligned_{h}.pdb", f"{h}_all.zip"}
        for suffix in ("_contacts_aa.txt", "_aa_coordinates.tsv"):
            alternatives = {f"{h}{suffix}", f"{h}{suffix}.gz"} & names
            if len(alternatives) != 1:
                raise ValueError(f"missing required or ambiguous prepared structure assets: {h}")
            required |= alternatives
        if not required <= names:
            raise ValueError(f"missing required prepared structure assets: {h}")
        try:
            with zipfile.ZipFile(structure / f"{h}_all.zip") as archive:
                members = archive.namelist()
                expected = {n for n in required if not n.endswith((".html", ".zip"))}
                if len(members) != len(set(members)) or not expected <= set(members):
                    raise ValueError(f"incomplete/duplicate structure ZIP: {h}")
                for name in members:
                    if PurePosixPath(name).name != name or "\\" in name or name not in names:
                        raise ValueError(f"unsafe/unbound structure ZIP member: {name}")
                    with archive.open(name) as stream:
                        digest_value = hashlib.sha256()
                        for block in iter(lambda: stream.read(1 << 20), b""):
                            digest_value.update(block)
                    if digest_value.hexdigest() != sha256(structure / name):
                        raise ValueError(f"stale structure ZIP member: {name}")
        except (zipfile.BadZipFile, OSError) as exc:
            raise ValueError(f"invalid structure ZIP: {h}") from exc
    map_hashes = {n[:-5] for n in names if re.fullmatch(r"[0-9a-f]{64}\.html", n)}
    if map_hashes != set(hashes):
        raise ValueError("structure map files differ from declared hashes")
    return True


def stage_archive(archive: Path, build: Path) -> None:
    """Extract a verified prepared-input ZIP into a fresh web directory, without overlays.

    Inputs contain only the two prepared sidecars and ``web/structure/*``; the release command
    emits the web table itself and checks its identity against these descriptors.
    """
    import stat
    import zipfile

    with zipfile.ZipFile(archive) as zipped:
        seen = set()
        for entry in zipped.infolist():
            name = entry.filename
            path = PurePosixPath(name)
            allowed = name in {"web/prepared-assets.json", "web/structures_metadata.tsv",
                               "web/structures_metadata.tsv.gz"} or name.startswith("web/structure/")
            if entry.is_dir():
                allowed = name in {"web/", "web/structure/"} or name.startswith("web/structure/")
            if (not allowed or path.is_absolute() or ".." in path.parts or "\\" in name
                    or path.as_posix() != name.rstrip("/") or name in seen
                    or stat.S_ISLNK(entry.external_attr >> 16)):
                raise ValueError(f"unsafe prepared archive member: {name}")
            seen.add(name)
        if (build / "web").exists():
            raise ValueError("prepared input requires a fresh build/web directory")
        zipped.extractall(build)
