"""What goes in each release bundle, declared once and never globbed.

``release.sh`` gathered a bundle with ``cp *.txt``, which would have swept in
``vdjdb_full_filtered.txt``, three ``*_broken.txt`` and three ``*_scored.txt`` side tables. The
shipped zip contains none of them, so a manual step, not the script, defined the release.
:data:`BUNDLES` is that definition, written down.

Three bundles per release, each named by role rather than by position. ``vdjmatch`` currently takes
the first ``.zip`` asset it finds, so once a release has more than one zip it picks an arbitrary
bundle; ``manifest.json`` is the durable fix, and it ships from the first multi-zip release whether
or not anything reads it yet (``ROADMAP.md`` section 3.1).

Member basenames are a contract. ``vdjmatch`` looks inside the zip for ``vdjdb.txt``,
``vdjdb.slim.txt`` and ``vdjdb_full.txt`` by basename, and ``vdjdb-web`` resolves the rest as
``<database.path>/<name>``. The directory prefix has never mattered to either; the basenames always
have (``docs/outputs.md`` section 2).
"""
from __future__ import annotations

import hashlib
import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

#: Read in blocks rather than all at once: the legacy `vdjdb.txt` is 240 MB.
_CHUNK = 1 << 20


@dataclass(frozen=True)
class Member:
    """One file in a bundle. ``source`` is relative to the build directory."""

    source: str
    #: The name inside the zip. Defaults to the source's basename, which is the contract.
    name: str = ""
    #: A member the bundle is still valid without -- the TCREMP motif tables, which only ship when
    #: the motif stage ran. A required member that is missing fails the release.
    optional: bool = False

    def member_name(self) -> str:
        return self.name or Path(self.source).name


@dataclass(frozen=True)
class Bundle:
    role: str
    #: ``{version}`` is substituted; the result is the asset name.
    filename: str
    #: The single directory every member sits under inside the zip.
    prefix: str
    members: tuple[Member, ...] = field(default_factory=tuple)


LEGACY = Bundle(
    role="legacy",
    filename="vdjdb-legacy-{version}.zip",
    prefix="vdjdb-{version}",
    members=(
        Member("legacy/vdjdb.txt"),
        Member("legacy/vdjdb.meta.txt"),
        Member("legacy/vdjdb.slim.txt"),
        Member("legacy/vdjdb.slim.meta.txt"),
        Member("legacy/vdjdb_full.txt"),
        Member("motifs/cluster_members.txt"),
        Member("motifs/motif_pwms.txt"),
        Member("motifs/cluster_members_tcremp.txt", optional=True),
        Member("motifs/motif_pwms_tcremp.txt", optional=True),
        Member("summary/vdjdb_summary_embed.html"),
        Member("LICENSE"),
        Member("latest-version.txt"),
    ),
)

PRIMARY = Bundle(
    role="primary",
    filename="vdjdb-{version}.zip",
    prefix="vdjdb-{version}",
    members=(
        Member("tables/records.parquet"),
        Member("tables/chains.parquet"),
        Member("tables/evidence.parquet"),
        Member("tables/epitopes.parquet"),
        Member("tables/restriction.parquet"),
        Member("tables/epitope_assessment.parquet"),
        Member("tables/vdjdb.parquet"),
        Member("tables/records.tsv"),
        Member("tables/chains.tsv"),
        Member("tables/evidence.tsv"),
        Member("tables/epitopes.tsv"),
        Member("tables/restriction.tsv"),
        Member("tables/epitope_assessment.tsv"),
        Member("tables/vdjdb.schema.json"),
        # `motifs/`, the directory `vdjdb motifs` writes. It said `motifs_new/` until 2026-09-28,
        # which is the output directory of the clustering sweeps under `docs/tuning/`, so resolving
        # this bundle against a pipeline build raised on two required members and the primary zip
        # could not be built at all. The test fixture was generated from this declaration, so it
        # created the directory it was asserting about.
        Member("motifs/cluster_members.txt"),
        Member("motifs/motif_pwms.txt"),
        Member("motifs/cluster_members_tcremp.txt", optional=True),
        Member("motifs/motif_pwms_tcremp.txt", optional=True),
        Member("summary/vdjdb_summary_embed.html"),
        Member("LICENSE"),
    ),
)

AIRR = Bundle(
    role="airr",
    filename="vdjdb-airr-{version}.zip",
    prefix="vdjdb-airr-{version}",
    members=(
        Member("airr/vdjdb.rearrangement.tsv"),
        Member("airr/vdjdb.receptor.tsv"),
        Member("airr/vdjdb.reactivity.tsv"),
        Member("airr/airr.yaml", optional=True),
        Member("LICENSE"),
    ),
)

BUNDLES: tuple[Bundle, ...] = (PRIMARY, LEGACY, AIRR)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(_CHUNK):
            h.update(block)
    return h.hexdigest()


def resolve(bundle: Bundle, build: Path, version: str) -> list[tuple[Path, str]]:
    """``(source path, member name)`` for every member present. Raises on a missing required one."""
    out: list[tuple[Path, str]] = []
    missing: list[str] = []
    for m in bundle.members:
        src = build / m.source
        if not src.exists():
            if not m.optional:
                missing.append(m.source)
            continue
        out.append((src, f"{bundle.prefix.format(version=version)}/{m.member_name()}"))
    if missing:
        raise FileNotFoundError(
            f"bundle {bundle.role!r} is missing required member(s): {', '.join(missing)} "
            f"(looked under {build})")
    return out


def write_zip(bundle: Bundle, build: Path, out: Path, version: str) -> Path:
    """Write the bundle's zip. Deterministic: members in declared order, fixed timestamps.

    A zip records an mtime per member, so archiving the same bytes twice produces different files
    unless the timestamp is pinned, and "same inputs, same bytes" is hard rule 7.
    """
    target = out / bundle.filename.format(version=version)
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for src, name in resolve(bundle, build, version):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            with src.open("rb") as fh:
                zf.writestr(info, fh.read())
    return target


def describe(path: Path, role: str) -> dict:
    """One manifest entry: role, file, digest, size and the member list."""
    with zipfile.ZipFile(path) as zf:
        members = [n.filename for n in zf.infolist()]
    return {"role": role, "file": path.name, "sha256": sha256(path),
            "bytes": path.stat().st_size, "members": members}


def write_manifest(entries: list[dict], out: Path, version: str, tag: str) -> Path:
    """``manifest.json``: what a consumer reads to pick a bundle by role rather than by order."""
    path = out / "manifest.json"
    path.write_text(json.dumps(
        {"version": version, "tag": tag, "bundles": entries}, indent=2) + "\n")
    return path


def write_checksums(paths: list[Path], out: Path) -> Path:
    """``SHA256SUMS`` in the format ``sha256sum -c`` reads."""
    path = out / "SHA256SUMS"
    path.write_text("".join(f"{sha256(p)}  {p.name}\n" for p in sorted(paths, key=lambda p: p.name)))
    return path
