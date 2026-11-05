"""Verify the paired web table, prepared metrics and served structure files before release."""
from __future__ import annotations

import json
import re
from pathlib import Path, PurePosixPath

from ..io.text import MAX_TEXT_BYTES, resolve_path
from .manifest import sha256


def validate(web: Path) -> bool:
    """Return whether a complete sealed structure directory accompanies the fresh web table.

    Table/chart-only releases are allowed. Any supplied structure metrics or assets require
    both descriptors and all files; stale or incomplete preparation fails before packaging.
    """
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
    digest(resolve_path(web / "vdjdb.txt"), prepared.get("web_table_sha256"))
    digest(metrics, prepared.get("structures_metadata_sha256"))
    digest(asset_path, prepared.get("structure_manifest_sha256"))
    for key in ("web_table_sha256", "source_metadata_sha256"):
        value = prepared.get(key)
        if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value) or assets.get(key) != value:
            raise ValueError(f"prepared input identity mismatch: {key}")
    hashes = assets.get("hashes")
    if (not isinstance(hashes, list) or any(not isinstance(h, str) or
            not re.fullmatch(r"[0-9a-f]{64}", h) for h in hashes)
            or hashes != sorted(set(hashes))):
        raise ValueError("invalid prepared structure hashes")
    files = assets.get("files")
    if not isinstance(files, list):
        raise ValueError("invalid prepared structure files")
    names = set()
    for entry in files:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise ValueError("invalid prepared structure file entry")
        name = entry["path"]
        relative = PurePosixPath(name)
        if (relative.is_absolute() or ".." in relative.parts or "\\" in name or
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
        if not required <= names or any(not ({f"{h}{suffix}", f"{h}{suffix}.gz"} & names)
                for suffix in ("_contacts_aa.txt", "_aa_coordinates.tsv")):
            raise ValueError(f"missing required prepared structure assets: {h}")
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
