# The reference corpus

One document per publication, twelve token families over it, and tf-idf weights. It reproduces what
[vdjdb.com/refsearch](https://vdjdb.com) serves and answers a question the endpoint cannot: whether a
receptor feature goes with an antigen because of itself, or because of something it travels with.

The corpus is an artifact rather than an index inside a service, so a downstream tool reads three
parquet files and does not have to re-derive a vocabulary.

```bash
uv run vdjdb corpus build --tables out/tables --out out/corpus
uv run vdjdb corpus query --epitope GILGFVFTL
uv run vdjdb corpus lift k:IRS --given e:GILGFVFTL
uv run vdjdb corpus refs --tables out/tables        # network; refreshes the committed PubMed input
```

## Documents

A document is a `reference.id`, and the records citing it are its content. 661 of them on the 2026-09
build:

| Kind | Documents | What it is |
|---|---|---|
| `pubmed` | 610 | a `PMID:` reference, the only kind with a text family |
| `pdb` | 39 | an RCSB structure entry |
| `issue` | 8 | a GitHub issue carrying a curated submission |
| `preprint` | 2 | a bioRxiv, medRxiv or arXiv identifier |
| `other` | 2 | a thesis repository and a vendor application note |

A record with no `reference.id` is not a document: 854 records carry none, all from one chunk, and a
document with no identifier cannot be cited or scored.

## Token families

A token is a prefixed string, so a consumer filters by family with a string comparison and needs no
second table. Each family is present because a question needs it and no other token can stand in.

| Family | Prefix | Example | Terms |
|---|---|---|---|
| word | `w:` | `w:influenza` | needs `vdjdb corpus refs` |
| V gene | `v:` | `v:TRBV9` | 227 |
| J gene | `j:` | `j:TRBJ2-7` | 79 |
| CDR3 k-mer | `k:` | `k:IRS` | 7,708 |
| CDR3 k-mer in V | `kv:` | `kv:IRS@TRBV19` | 243,300 |
| epitope | `e:` | `e:GILGFVFTL` | 2,112 |
| epitope k-mer | `ek:` | `ek:GIL` | 4,954 |
| antigen species | `a:` | `a:InfluenzaA` | 72 |
| antigen gene | `g:` | `g:M` | 624 |
| host species | `s:` | `s:HomoSapiens` | 3 |
| MHC allele | `m:` | `m:HLA-A*02:01` | 164 |
| MHC locus | `ml:` | `ml:HLA-A` | 25 |
| MHC class | `mc:` | `mc:MHCI` | 2 |

k is 3 for both k-mer families. Over 2,118 epitopes of 7 to 25 residues and 180,048 distinct CDR3s, a
3-mer has a document frequency worth an inverse document frequency, while a 5-mer is close to an
identifier of its own sequence and a 2-mer is in nearly every document.

**Why `k:` and `kv:` both.** "Is the CAS motif specific to HIV, or to its TRBV?" is a comparison
between the lift of `k:CAS` on HIV documents and its lift on HIV documents already carrying that V
gene. One token cannot express it and neither can a single search ranking.

**Why `ek:`.** Two epitopes sharing a core, or one epitope reported under two source species, are
linked by their k-mers and by nothing else. Measured: searching `GILGFVFTL` puts nine exact reporters
in the top ten and, tenth, `PMID:27036003`, which reports `GILEFVFTL` and `GILGLVFTL`. Those are
single-residue variants of the epitope being asked about, and an exact-match search cannot find the
altered-peptide-ligand study of its own epitope.

**Why `a:` and `s:` are different families.** `a:HomoSapiens` is a self-antigen; `s:HomoSapiens` is a
human donor. Collapsing them would merge two different claims about a record.

**Why three MHC granularities.** Restriction is a hierarchy and a question picks its level. The lift of
`k:CAS` given `mc:MHCI`, given `ml:HLA-A`, and given `m:HLA-A*02:01` are three different claims, and
one allele token makes the broad ones unaskable. Two fields is the resolution VDJdb curates at; deeper
fields are truncated because keeping them would split one restriction across several tokens.

## The MHC dictionary

`corpus/mhc.parquet`, 213 rows on the current build, is one row per distinct MHC call: the call as
curated, its two-field allele, its locus, the chain the curation placed it on, the curated class, and
the IPD-IMGT/HLA verdict. The three MHC token families are projections of it, so they cannot disagree
with one another, and a tool that wants to group VDJdb by locus reads one mapping rather than writing
a parser.

A locus is everything before the `*`. Murine naming has no field separator, so `H2-Kb` resolves
through a declared five-entry table of haplotype series (`H2-K`, `H2-D`, `H2-L`, `H2-IA`, `H2-IE`)
rather than by stripping a trailing letter, which would turn the IMGT gene `H2-Aa` into `H2-A`.
`H2-Aa`, `H2-Ab1`, `H2-Eb1` and `B2M` are their own locus.

The table doubles as a proofreading report. `status` is
[`vdjdb.assemble.epitopes.mhc_status`](cdr3-fixing.md), not a second check against the same IMGT file:
184 calls are `known` in IPD-IMGT/HLA, 27 `declared` in `proofreading/mhc_nonhuman.tsv` (murine H2,
macaque Mamu and the class I light chain are outside the HLA database), and **none** is `unknown`,
because an unknown call fails the build. The two that used to be were `HLA-A*08:01`, for which no
`HLA-A*08` exists at any resolution, and `HLA-B*12`, a serotype that split into B\*44 and B\*45 and is
no longer an allele group; both are corrected in `patches/mhc.dict`.

## Weighting

| Step | Form | Why |
|---|---|---|
| term frequency | `1 + log tf` | a paper reporting ten thousand receptors would otherwise dominate every receptor token by arithmetic rather than by relevance |
| inverse document frequency | `log((N + 1) / (df + 1)) + 1` | smoothed, so a term in every document still weighs 1 rather than 0 and a term in none cannot divide by zero |
| normalisation | L2 per document | a long abstract and a short one become comparable |

These are scikit-learn's `TfidfVectorizer` conventions, and that is the reason for choosing them:
`tests/unit/test_corpus.py` checks every weight against that library to nine decimal places. An
implementation that agrees with a reference is evidence; one that agrees with its own last run only
proves it has not changed.

## Files

```
corpus/documents.parquet   document_id, reference.id, kind, pmid, year, n_records, n_terms
corpus/terms.parquet       term_id, term, family, df, idf
corpus/postings.parquet    document_id, term_id, tf, weight
corpus/mhc.parquet         mhc, allele, locus, chain, mhc.class, status
```

661 documents, 259,270 terms and 1,040,354 postings on the current build, in 2.8 s.

Long postings rather than a sparse-matrix format, because the consumer is polars or duckdb and the
query is a join. `document_id` is assigned by sorted `reference.id` and `term_id` by sorted term, both
total orders, so the tables are reproducible without a hash and a diff between two builds is readable.
Postings are sorted by `(term_id, document_id)`, which makes a term lookup one contiguous slice.

TSV ships beside each parquet, the same pair the definitive tables use.

## The two questions

### `score`

The sum of a document's matched term weights, which is the quantity the `refsearch` endpoint returns as
`tf_idf`. A sum rather than a cosine against a normalised query vector, because every posting is
already L2-normalised and the query has no length to normalise against. Ties break on `reference.id`,
so the ranking is total and reproducible.

`vdjdb corpus query` takes the fields the `refsearch` client sends - `cdr3`, `antigen.epitope`,
`extra_parameters` drawn from `search_by_antigen` and `filter_stop_words`, and `species_to_search`
defaulting to the three species the client sends when the user picks none. A CDR3 becomes its k-mers
rather than one token, because a query is a motif: `CAS` is what a user types, and it is a k-mer.

### `lift`

How much more often a token appears among documents carrying every token in a condition. 1.0 means
the condition tells you nothing; `None` where nothing carries the condition, which is not a lift of
zero and must never be averaged as one.

Two modes, because they answer different questions.

| Mode | Unit | For |
|---|---|---|
| `documents` | a publication | "who reported this", letting one small study weigh as much as one large one |
| `occurrences` | a token instance, **within the term's own family** | "how much of this antigen's repertoire carries the motif" |

The occurrence mode is scoped to the term's family deliberately. Summing every family into one
denominator puts `k:CAS` over a total including epitope k-mers, MHC tokens and V genes, so its rate
moves with how many epitopes a paper studied rather than with the motif: measured, that dilution
flattens the whole comparison to within 2 % of 1.0.

The document mode cannot answer a common token at all. `k:CAS` is in 614 of 661 documents, so its
document-level lift is bounded near 1 however specific it is, and a test asserts that bound rather
than leaving it as a note.

## Validated against a motif nobody told it about

A tf-idf corpus over receptor k-mers either recovers what immunology already documents, or it is a
table nobody should draw a conclusion from. The case is GILGFVFTL, the influenza A M1 epitope, whose
specific TCRs are known for an RS motif in the beta CDR3. Nothing in the build knows that.

Of the 2,342 CDR3 3-mers with 50 or more occurrences among that epitope's documents:

| | |
|---|---|
| highest-lifting 3-mer | **`k:IRS`, 2.663x** |
| median across all 2,342 | 1.170x |
| RS-bearing 3-mers | 29, of which **25 above the median** |
| next four | `k:DGM` 2.603, `k:DLM` 2.566, `k:FMI` 2.519, `k:SIR` 2.473 |

And the control, on the same instrument: `k:CAS`, the germline-encoded start of nearly every beta
CDR3, lifts **0.969** on HIV-1 documents (28,422 of 739,216 CDR3 3-mer occurrences against 137,751 of
3,473,003 overall) and 1.006 once TRBV9 is held. Slightly depleted, not enriched, which is what a
germline motif should look like. Over the first 400 CDR3 3-mers on HIV-1 the range runs 0.05x to
4.70x, so the instrument has room to move and `k:CAS` genuinely sits at 1.

`tests/release/test_corpus_reproduction.py` pins every number above.

## PubMed records and abstracts

`vdjdb.corpus.pubmed` is the one place this package talks to NCBI, and it shares its transport with
`vdjdb.summary.references`, which asks a different question (publication year, via `esummary`).

* `efetch` with `db=pubmed&retmode=xml`, 200 ids per request, never one per id. One call returns title,
  abstract, journal, year and DOI together, so no second query is needed.
* NCBI etiquette is in the transport rather than in each caller: a `tool` parameter always, `email` and
  `api_key` only from `$NCBI_EMAIL` and `$NCBI_API_KEY`, at most three requests a second without a key
  and ten with one, and retry with backoff on 429 and 5xx. Four requests cover the whole corpus, so
  this is about being a good client rather than about throughput.
* An abstract is assembled from its `AbstractText` sections in order with the structured labels kept as
  text, because `BACKGROUND` and `RESULTS` are words a query can match. Inline `<i>` and `<sup>`
  children are read through, since `node.text` alone truncates a title at the first italic word.
* A PMID that returns no record is reported and never silently dropped: a `reference.id` is curated, so
  one that resolves to nothing is a curation finding.

Tokenisation lowercases, keeps inner hyphens so `HLA-A2-restricted`, `T-cell` and `SARS-CoV-2` stay one
token each, drops tokens shorter than three characters and purely numeric ones, and does not stem. No
stemming on purpose: a stemmer is a dependency and it makes a token unauditable, since `restrict` could
have come from `restriction`, `restricted` or `restricting` and a reader checking a ranking cannot
tell which.

### What is fetched and what is committed

Receptor, antigen and MHC tokens come from `chunks/` and need no network. Only the word family does,
and it follows the rule that already governs backgrounds: the running text is an input to the build and
never an output of it.

| File | Contents | Committed |
|---|---|---|
| `corpus/pubmed.tsv` | `reference.id, pmid, year, journal, title, doi, abstract_words, abstract_sha256` | yes |
| `corpus/text_terms.tsv` | `reference.id, term, tf` | yes |

No abstract text is stored or shipped. `abstract_sha256` says whether an abstract changed between
refreshes and `abstract_words` gives its length, so the corpus is auditable without keeping the prose;
the tokeniser is in the repository, so the transform is reproducible from the text. Titles are kept
because a title is a fact any bibliography carries.

Both are committed, reviewed inputs refreshed by their own pull request, exactly as
`summary/reference_years.tsv` is. A build never runs `vdjdb corpus refs`, and a corpus with neither
file builds from the other eleven families and says the word family is absent - which is the correct
state for a fork and for a first build.

Stop words are filtered at **query** time, not when the counts are written, so the committed corpus is
complete and the filter is reversible. The `filter_stop_words` flag the `refsearch` client sends is
exactly that choice, and a corpus that had already dropped them could not offer it.

## Not built, and what each would buy

The artifact is method-agnostic, so an alternative weighting is an additional file rather than a
rewrite. Both of these are gated on tf-idf reproducing the existing endpoint first, because that is
the only baseline in existence.

**PageRank over the term co-occurrence graph** would rank terms by centrality instead of by rarity.
Inverse document frequency handles receptor tokens badly at the common end: a 3-mer in every document
has an idf near zero and drops out, even when it is the hub connecting two antigen families. Measure it
against tf-idf on the same queries before preferring it.

**An embedding index**, sentence embeddings for the word family and `mir.embedding.TCREmp` for the
receptor families, would answer "which references are about something similar" rather than "which share
a token". It is a dense vector per document in a fourth file.
