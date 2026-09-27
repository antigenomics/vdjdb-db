"""PubMed records and abstracts: one client, one tokeniser, two committed inputs.

The corpus needs a title and an abstract per reference, and the text is an input to the build rather
than an output of it (hard rule 5). So this module fetches the records, tokenises them, and writes two
files that are **committed, reviewed inputs** refreshed by their own pull request, exactly as
``summary/reference_years.tsv`` already is. Every build is then offline and deterministic, and no
third-party abstract text is stored or shipped:

``corpus/pubmed.tsv``
    ``reference.id, pmid, year, journal, title, doi, abstract_words, abstract_sha256``. The digest is
    what says whether an abstract changed between refreshes; the word count is its length. Titles are
    kept because a title is a fact any bibliography carries.

``corpus/text_terms.tsv``
    ``reference.id, term, tf``. Counts only, no prose. The tokeniser is in this module, so the
    transform is reproducible from the text even though the text is not kept.

One transport, shared with :mod:`vdjdb.summary.references`, which asks NCBI a different question
(publication year, via ``esummary``) and should not open a second connection style to do it.
"""
from __future__ import annotations

import hashlib
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree

import polars as pl

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: Ids per request. NCBI's documented ceiling for a POST-sized id list is higher, but 200 matches
#: :data:`vdjdb.summary.references.PMID_CHUNK` so the two callers behave alike, and 610 references
#: are four requests either way. Never one request per id.
CHUNK = 200

#: Requests per second. NCBI allows 3 unauthenticated and 10 with a key, and asks callers to stay
#: under it rather than to discover the limit. Four requests make this politeness rather than
#: throughput, which is the reason to get it right: a tool that ignores the published limit is the
#: reason limits get tightened.
RATE_UNAUTHENTICATED = 3.0
RATE_WITH_KEY = 10.0

#: Retried statuses, and how long to wait. NCBI returns 429 under load and 500 transiently; both are
#: worth one more attempt and neither is worth a loop.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})
BACKOFF = (1.0, 4.0, 10.0)

PUBMED_TABLE = Path("corpus/pubmed.tsv")
TEXT_TERMS_TABLE = Path("corpus/text_terms.tsv")

RECORD_COLUMNS: tuple[str, ...] = (
    "reference.id", "pmid", "year", "journal", "title", "doi", "abstract_words", "abstract_sha256",
)

#: ``PMID:1234``, tolerant of the stray space 22 records carry. The same tolerance
#: :mod:`vdjdb.summary.references` applies, and for the same reason: a malformed id should stay a
#: visible curation item rather than silently drop out of the corpus.
PMID = re.compile(r"^PMID:\s*(\d+)\s*$")

#: A word: letters, digits and inner hyphens. Hyphens are kept because biomedical text carries
#: `T-cell`, `HLA-A2` and `SARS-CoV-2`, and splitting them would scatter one concept across three
#: tokens.
WORD = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*")

#: Shortest token kept. Two-letter words in this corpus are articles, prepositions and single-letter
#: gene suffixes, none of which distinguishes one paper from another.
MIN_WORD = 3

#: English function words, the classic short list used by every text-retrieval baseline. Applied at
#: **query** time and not when the counts are written, so the committed corpus is complete and the
#: filter is reversible: the ``filter_stop_words`` flag the `refsearch` client sends is exactly this
#: choice, and a corpus that had already dropped them could not offer it.
# A readable block rather than a 120-element list literal. Named first so the `.split()` is not a
# call on a literal, which is what ruff's SIM905 objects to.
_STOP_WORD_TEXT = """
about above after again against all also among and any are because been before being between both
but can did does doing during each few for from further had has have having her here hers him his
how into its itself more most not now off once only other our out over own same she should some such
than that the their them then there these they this those through too under until very was were what
when where which while who whom why will with would you your
"""
STOP_WORDS: frozenset[str] = frozenset(_STOP_WORD_TEXT.split())


def _rate() -> float:
    return RATE_WITH_KEY if os.environ.get("NCBI_API_KEY") else RATE_UNAUTHENTICATED


_last_request = 0.0


def http_get(url: str, *, timeout: int = 60) -> bytes:
    """One GET, rate-limited and retried. The only place this package touches the network.

    Sleeps to keep the process under :func:`_rate` requests a second, then retries
    :data:`RETRY_STATUS` on the :data:`BACKOFF` schedule. A status outside that set is not retried: a
    400 means the query is wrong and repeating it wastes someone else's capacity.
    """
    global _last_request
    interval = 1.0 / _rate()
    for wait in (*BACKOFF, None):
        gap = time.monotonic() - _last_request
        if gap < interval:
            time.sleep(interval - gap)
        _last_request = time.monotonic()
        try:
            with urllib.request.urlopen(url, timeout=timeout) as fh:
                return bytes(fh.read())
        except urllib.error.HTTPError as exc:
            if exc.code not in RETRY_STATUS or wait is None:
                raise
            time.sleep(wait)
        except urllib.error.URLError:
            if wait is None:
                raise
            time.sleep(wait)
    raise RuntimeError(f"unreachable: retries exhausted without raising for {url}")


def _params(extra: dict[str, str]) -> str:
    """Query string with NCBI's requested identification.

    ``tool`` always; ``email`` and ``api_key`` only from the environment, because an address committed
    to a public repository becomes someone else's contact for someone else's traffic.
    """
    params = {"tool": "vdjdb-db", **extra}
    if email := os.environ.get("NCBI_EMAIL"):
        params["email"] = email
    if key := os.environ.get("NCBI_API_KEY"):
        params["api_key"] = key
    return urllib.parse.urlencode(params)


def _text(node: ElementTree.Element | None) -> str:
    """All text under a node, including the tails of inline markup.

    An abstract carries ``<i>``, ``<sup>`` and ``<math>`` children, so ``node.text`` alone silently
    truncates at the first italic word.
    """
    return "" if node is None else "".join(node.itertext()).strip()


def _abstract(article: ElementTree.Element) -> str:
    """The abstract, sections in document order, structured labels kept as text.

    ``BACKGROUND``, ``METHODS`` and the rest are words a query can match, and dropping them would lose
    the only signal distinguishing a methods paper from a review.
    """
    pieces = []
    for section in article.iter("AbstractText"):
        label = section.attrib.get("Label", "")
        body = _text(section)
        if body:
            pieces.append(f"{label} {body}".strip() if label else body)
    return "\n".join(pieces)


def parse(xml: bytes) -> list[dict[str, object]]:
    """Every ``PubmedArticle`` in an efetch response, as plain dicts.

    A record missing a field gets an empty string rather than a null, which is the repository's one
    missing marker (hard rule 6).
    """
    root = ElementTree.fromstring(xml)
    out: list[dict[str, object]] = []
    for article in root.iter("PubmedArticle"):
        pmid = _text(article.find(".//MedlineCitation/PMID"))
        if not pmid:
            continue
        year = _text(article.find(".//JournalIssue/PubDate/Year"))
        if not year:
            # Some records carry only `MedlineDate`, free text like "2019 Nov-Dec".
            medline = _text(article.find(".//JournalIssue/PubDate/MedlineDate"))
            year = m[1] if (m := re.match(r"(\d{4})", medline)) else ""
        doi = ""
        for ident in article.iter("ArticleId"):
            if ident.attrib.get("IdType") == "doi":
                doi = _text(ident)
                break
        out.append({
            "pmid": pmid,
            "year": year,
            "journal": _text(article.find(".//Journal/ISOAbbreviation"))
            or _text(article.find(".//Journal/Title")),
            "title": _text(article.find(".//ArticleTitle")),
            "doi": doi,
            "abstract": _abstract(article),
        })
    return out


def fetch(pmids: list[str]) -> list[dict[str, object]]:
    """efetch every id, :data:`CHUNK` at a time. Title, abstract, journal, year and DOI in one call.

    ``esummary`` would give the year alone in less XML, but a second query for the abstract would
    double the requests for no gain: one efetch answers everything the corpus needs.
    """
    out: list[dict[str, object]] = []
    for i in range(0, len(pmids), CHUNK):
        query = _params({"db": "pubmed", "retmode": "xml", "rettype": "abstract",
                         "id": ",".join(pmids[i:i + CHUNK])})
        out.extend(parse(http_get(f"{EUTILS}/efetch.fcgi?{query}")))
    return out


def tokenise(text: str, *, min_word: int = MIN_WORD) -> Counter[str]:
    """Word counts for one document. Lowercased, hyphens kept, no stemming.

    No stemming on purpose. A stemmer is a dependency and it makes a token unauditable: `restrict`
    could have come from `restriction`, `restricted` or `restricting`, and a reader checking why a
    paper ranked where it did cannot tell which. Inverse document frequency already handles the
    inflections that matter, because they co-occur.

    Purely numeric tokens are dropped: a year or a figure number is not what distinguishes one paper
    from another, and `2019` appears in every paper published that year.
    """
    words = (w for w in WORD.findall(text.lower()) if len(w) >= min_word and not w.isdigit())
    return Counter(words)


def pmid_of(reference_id: str) -> str | None:
    """The bare PMID of a ``reference.id``, or ``None`` when the reference is not a PubMed one.

    51 of the 662 references in the current corpus are not: bioRxiv and medRxiv DOIs, an arXiv
    preprint, PDB entries, a thesis and GitHub issues. They become documents with no text family
    rather than being dropped, because their records are still curated evidence.
    """
    return m[1] if (m := PMID.match(reference_id)) else None


def build_tables(reference_ids: list[str]) -> tuple[pl.DataFrame, pl.DataFrame, list[str]]:
    """Fetch, then return ``(records, term_counts, missing)``.

    ``missing`` is the **PubMed** references NCBI returned no record for. Reported and never silently
    dropped: a ``reference.id`` is curated, so a PMID that resolves to nothing is a curation finding.

    A reference that is not a PubMed one is not missing. 51 of the 662 in the corpus are bioRxiv and
    medRxiv DOIs, an arXiv preprint, PDB entries, a thesis and GitHub issues; they become documents
    with no word family, which ``documents.kind`` records, and listing them here every run would bury
    a genuine dead PMID under fifty expected ones.

    Keyed one PMID to **many** references, not one to one. The corpus spells PMID 34433824 two ways,
    ``PMID:34433824`` and ``PMID: 34433824`` with a space, and a one-to-one mapping kept whichever the
    iteration reached last and dropped the other with nothing raised. Both are documents, because the
    document is the ``reference.id`` the records cite and the records cite both spellings.
    """
    wanted: dict[str, list[str]] = {}
    for ref in reference_ids:
        if pid := pmid_of(ref):
            wanted.setdefault(pid, []).append(ref)
    fetched = fetch(sorted(wanted))
    by_pmid = {str(rec["pmid"]): rec for rec in fetched}
    missing = sorted(ref for pid, refs in wanted.items() if pid not in by_pmid for ref in refs)

    record_rows, term_rows = [], []
    for pmid in sorted(by_pmid):
        rec = by_pmid[pmid]
        abstract = str(rec["abstract"])
        counts = tokenise(f"{rec['title']}\n{abstract}")
        digest = hashlib.sha256(abstract.encode("utf-8")).hexdigest()
        for ref in sorted(wanted[pmid]):
            record_rows.append({
                "reference.id": ref, "pmid": pmid, "year": str(rec["year"]),
                "journal": str(rec["journal"]), "title": str(rec["title"]), "doi": str(rec["doi"]),
                "abstract_words": sum(counts.values()), "abstract_sha256": digest,
            })
            term_rows.extend({"reference.id": ref, "term": term, "tf": tf}
                             for term, tf in sorted(counts.items()))
    records = (pl.DataFrame(record_rows, schema={c: pl.Utf8 for c in RECORD_COLUMNS}
                            | {"abstract_words": pl.UInt32})
               if record_rows else
               pl.DataFrame(schema={c: pl.Utf8 for c in RECORD_COLUMNS} | {"abstract_words": pl.UInt32}))
    terms = (pl.DataFrame(term_rows,
                          schema={"reference.id": pl.Utf8, "term": pl.Utf8, "tf": pl.UInt32})
             if term_rows else
             pl.DataFrame(schema={"reference.id": pl.Utf8, "term": pl.Utf8, "tf": pl.UInt32}))
    return (records.select(list(RECORD_COLUMNS)).sort("reference.id"),
            terms.sort("reference.id", "term"), missing)


def write(records: pl.DataFrame, terms: pl.DataFrame, *, root: Path | None = None) -> list[Path]:
    """Write both committed inputs, sorted, so a refresh shows as a reviewable diff."""
    base = root or Path()
    paths = []
    for frame, rel in ((records, PUBMED_TABLE), (terms, TEXT_TERMS_TABLE)):
        path = base / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.write_csv(path, separator="\t", quote_style="necessary")
        paths.append(path)
    return paths


def load_terms(path: Path | None = None) -> pl.DataFrame:
    """``corpus/text_terms.tsv``, or an empty frame when it has not been fetched yet.

    Empty is a usable state: the corpus builds from the receptor, antigen and MHC families alone, and
    says the text family is absent rather than failing. That is what a fork and a first build see.
    """
    target = path or TEXT_TERMS_TABLE
    schema = {"reference.id": pl.Utf8, "term": pl.Utf8, "tf": pl.UInt32}
    if not Path(target).exists():
        return pl.DataFrame(schema=schema)
    return pl.read_csv(target, separator="\t", schema=schema)
