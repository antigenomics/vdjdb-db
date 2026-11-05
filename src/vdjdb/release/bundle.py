"""Assemble a release: the zips, the manifest, the checksums, and ``latest-version.txt``.

Six steps:

1. **plan** -- derive the version from the tag, assert the tag does not already exist;
2. **prepare** -- rewrite ``latest-version.txt`` in the working tree, idempotently;
3. **build** -- write every bundle, all embedding that same ``latest-version.txt``;
4. **verify** -- the comparison against the last release, and the dashboard checks;
5. **publish** -- ``gh release create``, behind an environment with a required reviewer;
6. **finalize** -- commit ``latest-version.txt``, then fetch line 1 and fail on anything but 200.

Steps 1-3 are here; 4-6 belong to the workflow, which has the repository write access and network
access they need.

``latest-version.txt`` needs its own step because the file ships inside the zip it names. The URL is
deterministic from the tag and the tag is chosen before the build, so the self-reference resolves.

The 2026-06-03 prepend was performed on the day and committed, as ``66fd10f`` on the aldan3 GitLab,
which it never left. ``origin`` has two push URLs, so a push can reach one remote and not the other,
and the public GitHub repository went on serving the previous release's URL while the shipped zip
had the correct one. The fix is to commit the line after the release exists and then verify it,
which is what ``release.yml`` does.

The verification compares the tag, not the URL. Measured while this was written, line 1 named
``2026-05-16`` and the published latest was ``2026-06-03-ZENODO``, and both URLs returned 200.
"""
from __future__ import annotations

import os
import re
import tempfile
from datetime import date
from pathlib import Path

from ..io.text import resolve_path, write_text
from . import manifest

# Date tags preserve the public release convention; retain CalVer compatibility.
TAG = re.compile(r"v([0-9]{4})\.(0[1-9]|1[0-2])\.(0|[1-9][0-9]*)")
DATE_TAG = re.compile(r"[0-9]{4}-[0-9]{2}-[0-9]{2}(?:-[A-Za-z0-9][A-Za-z0-9._-]*)?")
RELEASES = "https://github.com/antigenomics/vdjdb-db/releases/download"
LATEST = Path("latest-version.txt")


def version_of(tag: str) -> str:
    """Keep date[-NOTE] tags intact; strip the leading v from CalVer tags."""
    if DATE_TAG.fullmatch(tag):
        try:
            date.fromisoformat(tag[:10])
        except ValueError:
            pass
        else:
            return tag
    if not TAG.fullmatch(tag):
        raise ValueError(
            f"tag {tag!r} is not v<YYYY>.<MM>.<PATCH> or a valid YYYY-MM-DD[-NOTE] tag")
    return tag[1:]


def validate_new_tag(tag: str, existing: list[str]) -> str:
    """Validate an immutable release tag; the first release in a month has patch zero."""
    version = version_of(tag)
    if tag in existing:
        raise ValueError(f"tag {tag!r} already exists; release tags are immutable")
    if DATE_TAG.fullmatch(tag):
        return version
    month, patch = version.rsplit(".", 1)
    monthly = [t for t in existing if TAG.fullmatch(t) and t[1:].rsplit(".", 1)[0] == month]
    if not monthly and patch != "0":
        raise ValueError("the first release in a month must use patch 0")
    return version


def legacy_url(tag: str) -> str:
    """The URL ``latest-version.txt`` line 1 must hold for this release.

    The legacy bundle, while legacy exists: clients in the wild download line 1 verbatim and
    expect that layout inside.
    """
    return f"{RELEASES}/{tag}/{manifest.LEGACY.filename.format(version=version_of(tag))}"


def prepare_latest(tag: str, path: Path = LATEST, *, write: bool = True) -> str:
    """Put this release's URL on line 1, idempotently. Returns the file's new content.

    Idempotent because a re-run of a failed release must not push the same URL twice, and written
    through a temporary file plus :func:`os.replace` because a half-written
    ``latest-version.txt`` is a bundle member that would ship truncated.

    ``write=False`` computes the content and leaves the file alone (#707). :data:`LATEST` is
    **tracked**, and a dry run with a tag no release carries left the repository naming a download
    that 404s - which is the defect `ROADMAP.md` §3.2 records as having already shipped once, from
    the other direction. The content is still returned, so the bundle member and the tree copy have
    one source either way, which is the property §3.2 asks for.
    """
    url = legacy_url(tag)
    lines = path.read_text().splitlines() if path.exists() else []
    if not lines or lines[0].strip() != url:
        lines = [url, *(n for n in lines if n.strip() != url)]
    content = "\n".join(lines) + "\n"
    if not write:
        return content
    with tempfile.NamedTemporaryFile("w", dir=path.parent or Path(), delete=False) as fh:
        fh.write(content)
        tmp = Path(fh.name)
    os.replace(tmp, path)
    return content


def build(build_dir: Path, out: Path, tag: str, *,
          bundles: tuple[manifest.Bundle, ...] = manifest.BUNDLES,
          extras: tuple[Path, ...] = (), source_commit: str | None = None) -> dict:
    """Write every bundle, the manifest and the checksums. Returns the manifest as a dict.

    ``build_dir`` is what ``vdjdb build`` / ``make legacy`` / ``motifs`` / ``summary`` produced,
    plus ``LICENSE`` and ``latest-version.txt`` linked in from the repository root -- every member
    is named in :data:`vdjdb.release.manifest.BUNDLES`, never globbed.

    ``extras`` are assets published beside the zips rather than inside one, so far only the identity
    lifecycle table. They are checksummed but carry no manifest role: nothing selects a release by
    them, and :data:`manifest.BUNDLES` stays the list of downloadable databases.
    """
    version = version_of(tag)
    if any(b.role in {"web", "web-structures"} for b in bundles):
        from .prepared import validate
        complete = validate(build_dir / "web")
        if complete and not any(b.role == "web-structures" for b in bundles):
            bundles += (manifest.WEB_STRUCTURES,)
        if not complete and any(b.role == "web-structures" for b in bundles):
            raise ValueError("web-structures requires complete prepared assets")
    out.mkdir(parents=True, exist_ok=True)
    entries, paths = [], []
    for b in bundles:
        zip_path = manifest.write_zip(b, build_dir, out, version)
        entries.append(manifest.describe(zip_path, b.role))
        paths.append(zip_path)
    metadata = manifest.write_manifest(entries, out, version, tag, source_commit=source_commit)
    manifest.write_checksums(
        paths + [metadata] + [resolve_path(p) for p in extras if resolve_path(p).exists()], out)
    return {"schema_version": 1, "version": version, "tag": tag,
            "source_commit": source_commit, "bundles": entries}


def stage(build_dir: Path, repo: Path = Path()) -> None:
    """Copy the two repository files every bundle includes into the build directory.

    ``LICENSE`` and ``latest-version.txt`` sit in the repository, not in ``out/``; copying them in
    rather than reaching out of the build directory keeps :func:`manifest.resolve` able to say that
    every member of every bundle is under one root.
    """
    for name in ("LICENSE", "latest-version.txt"):
        src, dst = repo / name, build_dir / name
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.suffix == ".txt":
                write_text(dst, src.read_text())
            else:
                dst.write_bytes(src.read_bytes())
