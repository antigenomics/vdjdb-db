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
| `patches/` | input. `antigen_epitope_species_gene.dict`, `nomenclature.conversions`, `IGM_nomenclature_table.tsv` |
| `res/` | input. `segments.txt`, `segments.aaparts.txt` — germline references for the legacy fixer (retired in phase 5, see `ROADMAP.md`) |
| `proofreading/` | input, **skills-only today** — no build code reads it. Alias and IMGT/HLA reference tables |
| `withheld/` | input, quarantined. Old-format submissions excluded from the build |
| `chunks_negative/`, `chunks_unformatted/`, `chunks_with_unconventional_aa/` | **excluded from the build.** Only `chunks/` is read |
| `chunks2/` | empty and untracked. Dead — do not add to it |
| `summary/` | the R dashboard. `.Rmd` and the extractor are input; `*.html`/`*.pdf`/`*.txt` are generated and gitignored |
| `database/` | **generated output**, gitignored except `dummy` and the two `*.meta.txt` files |
| `skills/` | curation workflows (`/vdjdb-extract`, `-format`, `-harmonize`, `-proofread`, `-duplicates`, `-publish`) |
| `src/` | **retired Groovy.** Reference only — it is the one correct spec for the meta files |
| `py_src/` | the current pandas pipeline, being replaced |

## Commands

The rewrite is in progress. Until phase 4 lands, the old path is still the build:

```bash
cd py_src && python runBuidDatabase.py      # assembly (needs 64 GB RAM — see ROADMAP §1)
bash release.sh                              # full release; requires a sibling ../vdjdb-motifs clone
```

The new package (see `ROADMAP.md` for which phases are live):

```bash
uv sync
uv run vdjdb qc chunks/                      # fail-fast chunk validation
uv run vdjdb build --out out/                # the three formats
uv run vdjdb motifs --out out/               # TCRNET + TCREMP
uv run vdjdb summary --out out/              # both dashboards
uv run vdjdb diff --against 2026-06-03       # the difference ledger
uv run pytest -q
```

Output goes to `out/`, **not** `build/` — `build/` is already gitignored as a Python packaging
convention and using it for release artifacts is confusing.

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
- `.gitignore` traps: `build/` is ignored (hence `out/`); `summary/*.txt` is ignored (so caches must be
  `.tsv`); `database/` is ignored but the two `*.meta.txt` files in it are force-added and will shadow
  generated copies until `git rm --cached`-ed.
- The `Dockerfile` installs `colorama` but the code imports `termcolor`. `Dockerfile_2` has it right.

## Open loops

Tracked in `STATUS.md`. Keep it current after each meaningful chunk of work.
