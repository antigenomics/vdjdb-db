# CLAUDE.md — working in `vdjdb-db`

Authoritative for how to work in this repo. What is in flight lives in `STATUS.md`; the migration plan
and its measured facts live in `ROADMAP.md`. Do not duplicate either here.

## What this repo is

The curation source and build for **VDJdb**, the TCR–antigen specificity database served at
<https://vdjdb.com>. `chunks/` **is the data**; everything else is machinery for validating it,
assembling it and publishing it.

Consumers: `vdjdb-web` (the site), `vdjmatch` (repertoire annotation), and anyone who downloads a
release zip. The release is the product.

## Layout

| Path | Role |
|---|---|
| `chunks/` | **input, the data.** One file per publication, `PMID_<id>.txt` (230 files, ~203k rows) |
| `patches/` | input, **declared corrections**. `antigen_epitope_species_gene.dict` (epitope → gene, species), `mhc.dict` (allele corrections, optionally scoped to one `reference.id`), `nomenclature.conversions`, `IGM_nomenclature_table.tsv` |
| `res/` | input. `segments.txt`, `segments.aaparts.txt` — germline references for the legacy fixer (retired in phase 5, see `ROADMAP.md`) |
| `proofreading/` | input, **the authority tables**. IMGT TCR alleles, IPD-IMGT/HLA alleles, Arden and antigen aliases — read by `curate/nomenclature.py` since phase 9 |
| `withheld/` | input, quarantined. Old-format submissions excluded from the build |
| `chunks_negative/`, `chunks_unformatted/`, `chunks_with_unconventional_aa/` | **excluded from the build.** Only `chunks/` is read |
| `chunks2/` | empty and untracked. Dead — do not add to it |
| `summary/` | the R dashboard. `.Rmd` and the extractor are input; `*.html`/`*.pdf`/`*.txt` are generated and gitignored |
| `database/` | **generated output**, gitignored except `dummy` and the two `*.meta.txt` files |
| `skills/` | curation workflows (`/vdjdb-extract`, `-format`, `-harmonize`, `-proofread`, `-duplicates`, `-publish`) |
| `src/*.groovy` | **retired Groovy.** Reference only — the one correct spec for the metadata, which phase 1 tests against. Moves to `attic/` in phase 14 |
| `src/vdjdb/` | the build. `assemble/` is the database, `emit/` projects it, `compare/` measures it |

## Commands

```bash
uv sync
uv run vdjdb qc                              # chunk validation, fail-fast
uv run vdjdb build --out out/                # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables # the legacy files, from the tables that shipped
uv run vdjdb convert airr --tables out/tables # AIRR Rearrangement + Reactivity
uv run vdjdb diff <reference.zip> out/legacy # the difference ledger
uv run vdjdb schema --table records          # generated metadata, for any declared table
uv run pytest -q
```

Not yet implemented (the ROADMAP phase that delivers each is printed on invocation):
`motifs`, `summary`, `release`, `refs`, `changelog`.

Output goes to `out/`, **not** `build/` -- `build/` is already gitignored as a Python packaging
convention and using it for release artifacts is confusing.

The pandas pipeline in `py_src/` was retired in phase 4. The 2026-06-03 release zip is the
reference the ledger measures against; `git show 2026-06-03:py_src/` still has the old build if it
is ever needed.

## The data model — `README.md` is authoritative

**`README.md` is the specification until `docs/standards/` replaces it** (ROADMAP phase 13). When
the code and the README disagree, the README wins and the code is the bug. Do not infer the model
from the shape of the data — the shape carries defects.

What it says, and what follows:

- **A chunk is one paper.** `chunks/PMID_<id>.txt` is that publication's report. Two rows in two
  different chunks are **independent reports**, never duplicates, even when every field matches —
  independent replication is a signal, and it is what phase 11 tunes motif clustering against.
- **A chunk row is one record**, and it **reports paired chains**: the alpha and the beta of one
  clone are columns of the same row. `chains` is derived from that, never the other way round.
- **Identity** is the complex-information columns plus the id fields, per the README: *"duplicate
  records (with identical complex information columns) are not allowed, but they will not be
  considered as duplicates in case they have distinct id fields"* — plus the chunk, by the rule
  above. Deduplication is therefore **within** a chunk only.
- **`method.*` and `meta.*` describe the record**, not the act of curating it. They are what the
  publication reports about how the specificity was established, so they belong on the record.
  Only `submitter`, `comment` and `chunk.id` are properties of the curation.
- **Record ids are assigned before CDR3 repair.** Two trimmed sequences that repair to the same
  full one are still two observations; assigning after repair merged 215 pairs of records that the
  publications reported separately.

### The assembly line is tidy; the legacy shapes are projections

`vdjdb.assemble.tables` produces flat tables linked by `record_id` — one observational unit each,
one variable per column, no JSON blobs, no paired alpha/beta columns, no comma-joined sets. Those
are the database.

Every shipped file is a **join and a pivot** away from them, in `vdjdb.emit.*`. Nothing outside
`emit/legacy.py` may know about `complex.id`, the `method`/`meta`/`cdr3fix` blobs, or the positional
column orders. If a legacy quirk leaks into `assemble/`, that is the bug.

## Domain conventions

**`cdr3` in VDJdb is junction space.** Cys104 through Phe/Trp118, **both anchors included**. That is
*not* AIRR's or arda's `cdr3_aa`, which excludes both — `junction_aa` is two residues longer.
Conflating them silently corrupts every coordinate. All `v.end` / `j.start` values are 0-based amino
acid indices in that space.

Four coordinate spaces meet in this codebase. Every conversion between them belongs in one module with
round-trip tests, and nowhere else:

| Space | Basis |
|---|---|
| VDJdb `v.end`/`j.start` | 0-based aa, junction space |
| `vdjtools` `Scenario.v_end`/`j_start` | 0-based half-open, CDR3-**nt** space |
| `arda.annotate.dmap` | 1-based **closed**, junction-nt space |
| AIRR | 1-based closed, **sequence** space |

IMGT nomenclature for V/D/J and MHC. Species vocabulary is `HomoSapiens`, `MusMusculus`,
`RattusNorvegicus`, `MacacaMulatta` (CamelCase, no spaces).

## Hard rules

1. **The legacy motif files are parsed positionally.** `vdjdb-web`'s `Motifs.scala` hands Tablesaw a
   fixed `Array[ColumnType]` with no header check — 27 entries for `motif_pwms.txt`, 19 for
   `cluster_members.txt`. Any inserted, removed or reordered column silently mistypes or shifts the
   whole table. Column order is a contract, not a convention.
2. **`vdjdb.meta.txt` must match `vdjdb.txt` column-for-column, in order.** `vdjdb-web` builds its
   entire schema from the meta file. It has already drifted once.
3. **Never loop a batch API.** `arda.cdr3fix.markup_batch`, `arda.annotate.annotate_records`,
   `TCREmp.embed`, `vdjtools.overlap.tcrnet` — one call with N rows, never N calls with one. A process
   pool around anything that spawns mmseqs2 or BLAS threads deadlocks or thrashes.
4. **Dedupe before an expensive per-record call.** VDJdb has ~192k rows but ~191k distinct
   `(cdr3, v, j, species)` keys and far fewer distinct `(cdr3_aa, v, j)` for junction inference. Run
   on the unique set and join back.
5. **Backgrounds never ship.** `isalgo/airr_control` is CC-BY-NC-ND-4.0 and VDJdb is AGPL-3.0-only.
   Stream at build time; ship derived statistics only.
6. **Empty string is the only missing marker** in the pipeline. The pandas `None`/`NaN`/`""` three-way
   ambiguity is the source of more than one shipped bug.
7. **Every output is reproducible, not merely repeatable.** Same inputs, same bytes — in another
   process, on another host, at another core count. Concretely:
   - **One seed.** `vdjdb.config.SEED`. No `random.seed()` at module scope, no unseeded default, no
     per-call literal. Pass it explicitly to every sampler, shuffler, clustering init and seeded hash.
   - **Sort after anything unordered.** `group_by` without `maintain_order=True`, a `set`, a `dict`
     keyed on strings under `PYTHONHASHSEED`, a thread pool, `os.listdir` — all of them return an
     order you did not choose. Sort, and give ties an explicit tiebreak: an unstable sort over
     equal-labelled groups made the ledger report 158/152/158 changed rows for one comparison.
   - **Never let worker count change the answer.** Split into as many big contiguous slices as there
     are workers, reassemble in slice order, never a pool of small tasks.
   - A determinism test — run it N times, assert one digest — belongs with any stage that samples,
     clusters or parallelises.
8. **Vectorize the hot path before parallelising it.** Reach order, and say which rung you stopped
   at: one polars/numpy expression → one batched call into existing C++ (`seqtree`, `vdjtools`) →
   new C++. Materialising rows into Python is the usual mistake: the ledger's row comparison built
   6.26 million tuples per table until it was replaced with a hashed group-count join in polars,
   which cut the run from 7.6 s to 3.5 s and left 14 of 208,447 keys to inspect in Python.

9. **Nothing computed is ever cached between builds.** Every output is recomputed from `chunks/` on
   every run. A stored intermediate is the first place a build goes wrong: it is authoritative until
   the code that produced it changes, and then it is silently wrong in a way no test sees, because
   the test reads the same stale file. If a stage is slow, make it faster or accept the minutes.
   - Two things this does **not** forbid. **Deduplicating before an expensive per-record call** is
     not caching — it runs the call on the distinct key set *within one build* and joins back, and
     the functions are deterministic in their arguments, so it cannot change an answer (rule 4).
     **Fetching an input** — a release zip, a germline reference, an HF background — is a download,
     not a cache; it is data arriving, not a result being remembered.
   - A derived table that exists to make the build *offline and deterministic*, such as the
     publication-year table the dashboard reads instead of calling NCBI at render time, is a
     **committed, reviewed input**. It is refreshed by its own pull request, never written by a
     build. Do not call it a cache and do not let a build update it in place.

## The issue tracker is mostly a submission queue

440 issues, 130 open, and **103 of the open ones (79 %) are pending papers, preprints and datasets
waiting to be curated** -- not defects. 22 are curation quality, 13 are build infrastructure
(`ROADMAP.md` section 4a).

So: a numbered issue in `ROADMAP.md` is one of the 35 that the build work can close. Do not read the
tracker as a defect list, do not assume an open issue is a bug, and check the label before treating
one as in scope. The queue is worked through the curation skills, not the build.

## Commit conventions

- Every commit that resolves a tracker issue ends with `Closes #N` — the same rule the
  `/vdjdb-publish` skill already applies to chunk commits (`Fixes #issue_id`).
- Gitflow: `master` → `dev` → `feature/*` → `dev` → `master`. `master` only ever receives merges from
  `dev` or `hotfix/*`.
- Chunk submissions are one commit per chunk, referencing that chunk's PMID issue.

## Developer setup gotchas

- `~/hf/airr_control` has every background checked out **except
  `human.trb.aa.vdjtools.tsv.gz`**, which is still an LFS pointer. That is the most important TCRNET
  background. Run `git lfs pull --include=human.trb.aa.vdjtools.tsv.gz` before local motif work.
- `.gitignore` traps: `build/` is ignored (hence `out/`); `summary/*.txt` is ignored (so committed
  tables under `summary/` must be
  `.tsv`); `database/` is ignored but the two `*.meta.txt` files in it are force-added and will shadow
  generated copies until `git rm --cached`-ed.
- The `Dockerfile` installs `colorama` but the code imports `termcolor`. `Dockerfile_2` has it right.

## Open loops

Tracked in `STATUS.md`. Keep it current after each meaningful chunk of work.
