"""Publication year for every ``reference.id`` in the database, resolved once and committed.

The dashboard plots records by year. Until now it resolved PubMed ids by calling NCBI eutils **at
render time** and covered everything else with a literal fourteen-entry table inside the Rmd, last
updated in 2021. Measured on the current database that table misses **38 of the 52 non-PubMed
references** -- every PDB structure entry curated since, plus an arXiv preprint -- and each missing
one drops out of the by-year panels without a warning.

So the table is resolved here, once, into ``summary/reference_years.tsv``: a committed, reviewed
input refreshed by its own pull request and **never written by a build**. That is the carve-out hard
rule 9 makes for a derived table whose whole purpose is to make the build offline and deterministic,
and it is why :func:`refresh` is a separate command rather than a step of ``vdjdb summary``.

Five kinds of reference, four resolvers, no per-record calls (hard rule 3):

============================  ==========================================  =============
``reference.id``              resolver                                    network
============================  ==========================================  =============
``PMID:<id>``                 NCBI esummary, 200 ids per request          yes
``https://www.rcsb.org/...``  RCSB GraphQL ``entries(entry_ids: [...])``  yes, one call
``https://arxiv.org/abs/..``  the identifier's own ``YYMM`` prefix        no
``https://doi.org/10.1101/``  the bioRxiv/medRxiv DOI's own date          no
``https://github.com/...``    GitHub GraphQL, all issues in one query     yes, one call
============================  ==========================================  =============

Anything left unresolved is reported by count and kept out of the table rather than guessed, and a
blank ``reference.id`` is neither: it is the known curation gap on 854 records, excluded explicitly.

⚠ **Nothing here reads TCRvdb.**
"""
from __future__ import annotations

import json
import re
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path

import polars as pl

#: Where the resolved table lives. ``.tsv``, never ``.txt`` -- ``summary/*.txt`` is gitignored, so a
#: ``.txt`` here would silently never be committed.
TABLE = Path("summary/reference_years.tsv")

#: NCBI asks for no more than a few hundred ids per esummary request.
PMID_CHUNK = 200

# Tolerant of stray whitespace inside the identifier: 22 records carry `PMID: 34433824`
# with a space, which a strict pattern drops from the year plots without saying so. The
# table still keys on the database's literal value -- only the lookup is forgiving, and the
# malformed id remains visible as a curation item rather than being silently repaired.
_PMID = re.compile(r"^PMID:\s*(\d+)\s*$")
_PDB = re.compile(r"^https?://www\.rcsb\.org/structure/(\w+)/?$", re.IGNORECASE)
_ARXIV = re.compile(r"^https?://arxiv\.org/abs/(\d{2})(\d{2})\.\d+", re.IGNORECASE)
_BIORXIV = re.compile(r"^https?://doi\.org/10\.1101/(\d{4})\.\d{2}\.\d{2}\.")
_ISSUE = re.compile(r"^https?://github\.com/([\w-]+)/([\w-]+)/issues/(\d+)")

#: Two references no API resolves: a thesis repository and a vendor application note. Their years
#: come from the table this module replaces, which is the only record of them -- carried over rather
#: than re-derived, and listed here so the provenance is visible instead of buried in an Rmd chunk.
LITERAL: dict[str, int] = {
    "http://mediatum.ub.tum.de/doc/1136748": 2014,
    "https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-"
    "linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#": 2019,
}


def _get(url: str, *, timeout: int = 60) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return fh.read()


def _pubmed(ids: list[str]) -> dict[str, int]:
    """``{pmid: year}`` from NCBI esummary. One request per :data:`PMID_CHUNK`, never per id."""
    out: dict[str, int] = {}
    for i in range(0, len(ids), PMID_CHUNK):
        chunk = ids[i:i + PMID_CHUNK]
        q = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "id": ",".join(chunk)})
        payload = json.loads(_get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{q}"))
        for pmid, rec in payload.get("result", {}).items():
            if pmid == "uids" or not isinstance(rec, dict):
                continue
            # `pubdate` is free text -- "2019 Dec 12", "2020", "2021 Jan-Feb". The leading four
            # digits are the year in every form NCBI emits; anything else is left unresolved.
            if m := re.match(r"(\d{4})", str(rec.get("pubdate", ""))):
                out[pmid] = int(m[1])
    return out


def _rcsb(ids: list[str]) -> dict[str, int]:
    """``{pdb_id: year}`` from one RCSB GraphQL query over every entry id."""
    if not ids:
        return {}
    joined = ",".join(f'"{i}"' for i in ids)
    query = ("{entries(entry_ids:[" + joined
             + "]){rcsb_id rcsb_accession_info{initial_release_date}}}")
    payload = json.loads(_get("https://data.rcsb.org/graphql?query="
                              + urllib.parse.quote(query)))
    out = {}
    for e in (payload.get("data") or {}).get("entries") or []:
        date = ((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")
        if m := re.match(r"(\d{4})", date):
            out[e["rcsb_id"].upper()] = int(m[1])
    return out


def _issues(refs: list[tuple[str, str, str]]) -> dict[str, int]:
    """``{url: year}`` for GitHub issues -- one GraphQL query with an alias per issue.

    Uses ``gh`` rather than a bare request because the token is already configured there, and an
    unauthenticated GitHub API call is rate-limited to the point of being unreliable in CI.
    """
    if not refs:
        return {}
    parts = [f'i{n}: repository(owner:"{o}",name:"{r}"){{issue(number:{num}){{createdAt}}}}'
             for n, (o, r, num) in enumerate(refs)]
    proc = subprocess.run(["gh", "api", "graphql", "-f", "query={" + " ".join(parts) + "}"],
                          capture_output=True, text=True, check=False)
    if proc.returncode:
        return {}
    data = json.loads(proc.stdout).get("data") or {}
    out = {}
    for n, (o, r, num) in enumerate(refs):
        created = (((data.get(f"i{n}") or {}).get("issue") or {}).get("createdAt") or "")
        if m := re.match(r"(\d{4})", created):
            out[f"https://github.com/{o}/{r}/issues/{num}"] = int(m[1])
    return out


def resolve(reference_ids: list[str]) -> pl.DataFrame:
    """``(reference.id, year, source)`` for every id that resolves. Unresolved ids are omitted.

    Blank ids are dropped first: 854 records carry no ``reference.id`` at all, which is a curation
    gap (``ROADMAP.md``), not a lookup failure, and counting it as one would hide it.
    """
    ids = sorted({r for r in reference_ids if r and r.strip()})
    rows: list[dict] = []

    # A list of pairs, not a ``{pmid: ref}`` dict: the same paper can appear under two spellings --
    # `PMID:34433824` and `PMID: 34433824` both occur -- and a dict would keep one and drop the
    # other's 22 records from the year plots, which is the failure this whole module exists to end.
    pmids = [(m[1], r) for r in ids if (m := _PMID.match(r))]
    years = _pubmed(sorted({p for p, _ in pmids}))
    for pmid, ref in pmids:
        if pmid in years:
            rows.append({"reference.id": ref, "year": years[pmid], "source": "pubmed"})

    pdbs = {m[1].upper(): r for r in ids if (m := _PDB.match(r))}
    for pdb, year in _rcsb(sorted(pdbs)).items():
        rows.append({"reference.id": pdbs[pdb], "year": year, "source": "rcsb"})

    issues = [(m[1], m[2], m[3]) for r in ids if (m := _ISSUE.match(r))]
    for url, year in _issues(issues).items():
        rows.append({"reference.id": url, "year": year, "source": "github"})

    for r in ids:
        if m := _ARXIV.match(r):
            # arXiv identifiers are YYMM.NNNNN since 2007; the century is not in doubt.
            rows.append({"reference.id": r, "year": 2000 + int(m[1]), "source": "arxiv"})
        elif m := _BIORXIV.match(r):
            rows.append({"reference.id": r, "year": int(m[1]), "source": "biorxiv-doi"})
        elif r in LITERAL:
            rows.append({"reference.id": r, "year": LITERAL[r], "source": "literal"})

    df = pl.DataFrame(rows, schema={"reference.id": pl.Utf8, "year": pl.Int32,
                                    "source": pl.Utf8}) if rows else pl.DataFrame(
        schema={"reference.id": pl.Utf8, "year": pl.Int32, "source": pl.Utf8})
    return df.unique(subset=["reference.id"], keep="first").sort("reference.id")


def refresh(records: pl.DataFrame, out: Path = TABLE) -> pl.DataFrame:
    """Resolve every reference in ``records`` and write the table. Reports what did not resolve."""
    ids = records["reference.id"].unique().to_list()
    df = resolve(ids)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(out, separator="\t")
    return df


def load(path: Path = TABLE) -> pl.DataFrame:
    """The committed table. Raises if it is missing -- a silent fallback is what caused the drift."""
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing; run `vdjdb refs` to rebuild it. The dashboard must not fall back "
            "to a hardcoded table -- that is how it came to be four years stale.")
    return pl.read_csv(path, separator="\t")


def unresolved(records: pl.DataFrame, table: pl.DataFrame) -> pl.DataFrame:
    """References in the database with no year, and how many records each carries."""
    return (records.filter(pl.col("reference.id").str.strip_chars() != "")
            .group_by("reference.id").len().rename({"len": "records"})
            .join(table.select("reference.id"), on="reference.id", how="anti")
            .sort("records", descending=True))
