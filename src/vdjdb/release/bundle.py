"""Assemble a release: the zips, the manifest, the checksums, and ``latest-version.txt``.

Six steps, and the sixth is the one that never happened:

1. **plan** -- derive the version from the tag, assert the tag does not already exist;
2. **prepare** -- rewrite ``latest-version.txt`` in the working tree, idempotently;
3. **build** -- write every bundle, all embedding that same ``latest-version.txt``;
4. **verify** -- the difference ledger and the dashboard checks;
5. **publish** -- ``gh release create``, behind an environment with a required reviewer;
6. **finalize** -- *commit* ``latest-version.txt``, then fetch line 1 and fail on anything but 200.

Steps 1-3 are here; 4-6 belong to the workflow, because they touch the repository and the world.

**Why ``latest-version.txt`` needs its own step at all.** The file ships *inside* the zip it names,
so it is self-referential -- but the URL is deterministic from the tag and the tag is chosen before
the build, which is what makes it solvable. The legacy script prepended the line only when a shell
comparison failed, and the prepend either never ran or was never committed: measured today, line 1
names ``2026-05-16`` while the published latest is ``2026-06-03-ZENODO``. Both URLs return 200, so
"does line 1 resolve" would not have caught it -- the check has to compare the *tag*.
"""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

from . import manifest

#: ``v<YYYY>.<MM>.<PATCH>`` -- CalVer, because this is a dataset with no API to break, and the patch
#: component is exactly what the date scheme lacked: ``2024-11-27`` needed an asset called
#: ``vdjdb-2024-11-27-fixed.zip`` because there was nowhere to put a re-release.
TAG = re.compile(r"^v(\d{4})\.(\d{2})\.(\d+)$")
RELEASES = "https://github.com/antigenomics/vdjdb-db/releases/download"
LATEST = Path("latest-version.txt")


def version_of(tag: str) -> str:
    """``v2026.09.1`` -> ``2026.09.1``. Raises on anything that is not the declared scheme."""
    if not TAG.match(tag):
        raise ValueError(
            f"tag {tag!r} is not v<YYYY>.<MM>.<PATCH>. The date scheme it replaces had nowhere to "
            "put a re-release, which is how 2026-06-03, 2026-06-03-ZENODO and "
            "2026-06-03-SNAPSHOT all came to name the same build.")
    return tag[1:]


def legacy_url(tag: str) -> str:
    """The URL ``latest-version.txt`` line 1 must carry for this release.

    The **legacy** bundle, while legacy exists: clients in the wild download line 1 verbatim and
    expect that layout inside.
    """
    return f"{RELEASES}/{tag}/{manifest.LEGACY.filename.format(version=version_of(tag))}"


def prepare_latest(tag: str, path: Path = LATEST) -> str:
    """Put this release's URL on line 1, idempotently. Returns the file's new content.

    Idempotent because a re-run of a failed release must not push the same URL twice, and written
    through a temporary file plus :func:`os.replace` because a half-written
    ``latest-version.txt`` is a bundle member that would ship truncated.
    """
    url = legacy_url(tag)
    lines = path.read_text().splitlines() if path.exists() else []
    if not lines or lines[0].strip() != url:
        lines = [url, *(n for n in lines if n.strip() != url)]
    content = "\n".join(lines) + "\n"
    with tempfile.NamedTemporaryFile("w", dir=path.parent or Path(), delete=False) as fh:
        fh.write(content)
        tmp = Path(fh.name)
    os.replace(tmp, path)
    return content


def build(build_dir: Path, out: Path, tag: str, *,
          bundles: tuple[manifest.Bundle, ...] = manifest.BUNDLES) -> dict:
    """Write every bundle, the manifest and the checksums. Returns the manifest as a dict.

    ``build_dir`` is what ``vdjdb build`` / ``make legacy`` / ``motifs`` / ``summary`` produced,
    plus ``LICENSE`` and ``latest-version.txt`` linked in from the repository root -- every member
    is named in :data:`vdjdb.release.manifest.BUNDLES`, never globbed.
    """
    version = version_of(tag)
    out.mkdir(parents=True, exist_ok=True)
    entries, paths = [], []
    for b in bundles:
        zip_path = manifest.write_zip(b, build_dir, out, version)
        entries.append(manifest.describe(zip_path, b.role))
        paths.append(zip_path)
    manifest.write_manifest(entries, out, version, tag)
    manifest.write_checksums(paths, out)
    return {"version": version, "tag": tag, "bundles": entries}


def stage(build_dir: Path, repo: Path = Path()) -> None:
    """Link the two repository files every bundle carries into the build directory.

    ``LICENSE`` and ``latest-version.txt`` live in the repository, not in ``out/``; copying them in
    rather than reaching out of the build directory keeps :func:`manifest.resolve` able to say that
    every member of every bundle is under one root.
    """
    for name in ("LICENSE", "latest-version.txt"):
        src, dst = repo / name, build_dir / name
        if src.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(src.read_bytes())
