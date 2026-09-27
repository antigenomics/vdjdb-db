# VDJdb build roadmap

Migration of database proofreading, assembly and release from the Docker + GitLab + manual-`scp`
pipeline to a `uv`/polars Python package driven by GitHub Actions.

Deliverables with acceptance criteria, not dates. Sections 0-12 are the plan: the rules a change is
measured against, what ships, the gates, and the phase list.

The execution record is not in this file. Every "phase N result", every measurement and every
parameter sweep is in `ROADMAP_local.md`, with `CHANGELOG_local.md` as the running log. Both are
untracked, because they change every working session. Section numbering is continuous across the
two, so a reference to "section 30" resolves in the local file.

---

## 0. Hard rules

Referenced throughout as "rule N". Breaking one is a defect regardless of what the tests say.


1. The legacy motif files are parsed positionally. `vdjdb-web`'s `Motifs.scala` hands Tablesaw a
   fixed `Array[ColumnType]` with no header check: 27 entries for `motif_pwms.txt`, 19 for
   `cluster_members.txt`. Any inserted, removed or reordered column mistypes or shifts the table
   with no error. Column order is a contract, not a convention.
2. `vdjdb.meta.txt` must match `vdjdb.txt` column-for-column, in order. `vdjdb-web` builds its
   schema from the meta file. It has already drifted once.
3. Never loop a batch API. `arda.cdr3fix.markup_batch`, `arda.annotate.annotate_records`,
   `TCREmp.embed`, `vdjtools.overlap.tcrnet`: one call with N rows, never N calls with one. A process
   pool around anything that spawns mmseqs2 or BLAS threads deadlocks or thrashes.
4. Dedupe before an expensive per-record call. VDJdb has ~192k rows but ~191k distinct
   `(cdr3, v, j, species)` keys and far fewer distinct `(cdr3_aa, v, j)` for junction inference. Run
   on the unique set and join back.
5. Backgrounds never ship. Stream at build time; ship derived statistics only (`count.bg`,
   `total.bg`). A background is an input to the build, never an output of it.
6. Empty string is the only missing marker in the pipeline. The pandas `None`/`NaN`/`""` three-way
   ambiguity is the source of more than one shipped bug.
7. Every output is reproducible, not merely repeatable. Same inputs, same bytes, in another
   process, on another host, at another core count. Concretely:
   - One seed: `vdjdb.config.SEED`. No `random.seed()` at module scope, no unseeded default, no
     per-call literal. Pass it explicitly to every sampler, shuffler, clustering init and seeded hash.
   - Sort after anything unordered. `group_by` without `maintain_order=True`, a `set`, a `dict`
     keyed on strings under `PYTHONHASHSEED`, a thread pool, `os.listdir` all return an order you
     did not choose. Sort, and give ties an explicit tiebreak: an unstable sort over
     equal-labelled groups made `vdjdb diff` report 158/152/158 changed rows for one comparison.
   - Never let worker count change the answer. Split into as many big contiguous slices as there
     are workers, reassemble in slice order, never a pool of small tasks.
   - A determinism test, run it N times and assert one digest, belongs with any stage that samples,
     clusters or parallelises.
8. Vectorize the hot path before parallelising it. Reach order, and say which rung you stopped
   at: one polars/numpy expression → one batched call into existing C++ (`seqtree`, `vdjtools`) →
   new C++. Materialising rows into Python is the usual mistake: the release comparison's row
   matching built 6.26 million tuples per table until it was replaced with a hashed group-count join
   in polars, which cut the run from 7.6 s to 3.5 s and left 14 of 208,447 keys to inspect in Python.

9. Nothing computed is ever cached between builds. Every output is recomputed from `chunks/` on
   every run. A stored intermediate stays authoritative until the code that produced it changes,
   and then it is wrong in a way no test sees, because the test reads the same stale file. If a
   stage is slow, make it faster or accept the minutes.
   - Two things this does not forbid. Deduplicating before an expensive per-record call is
     not caching: it runs the call on the distinct key set within one build and joins back, and
     the functions are deterministic in their arguments, so it cannot change an answer (rule 4).
     Fetching an input, such as a release zip, a germline reference or an HF background, is a
     download rather than a cache.
   - A derived table that exists to make the build offline and deterministic, such as the
     publication-year table the dashboard reads instead of calling NCBI at render time, is a
     committed, reviewed input. It is refreshed by its own pull request, never written by a
     build. Do not call it a cache and do not let a build update it in place.

---

## 1. Current pipeline and its defects

The current build is a Docker image (ubuntu:18.04, Python 3.10 from source, OpenJDK 8 + Groovy 3.0.9,
R 4.x + ~25 CRAN packages, TeXLive, the VDJtools 1.2.1 jar) driven by `release.sh`, triggered from
GitLab CI by `ssh`-ing into aldan3 and `sbatch`-ing a runner. `release.sh` ends by printing an `scp`
command; a human then runs `gh release`.

Measured defects:

| Defect | Evidence |
|---|---|
| The shipped release is not what `release.sh` produces | `release.sh` does `cp *.txt`, which would sweep in `vdjdb_full_filtered.txt`, three `*_broken.txt` and three `*_scored.txt`. The 2026-06-03 zip has exactly 10 files and none of those. |
| Bundle shape drifts every release | 2024: flat + `LICENSE.txt` · 2025-02-21: `vdjdb-<date>/` + 4 extra files · 2025-09-25: `database/` · 2026-06-03: `vdjdb-<date>/` + `LICENSE` |
| Column metadata no longer describes the data | `vdjdb.meta.txt` is a static git file last written by the retired Groovy code: no `TCR_hash` row, `vdjdb.score` in the wrong position, `v.end`/`j.start` wrong place and order in slim |
| `web.cdr3fix.unmp` is wrong on 7,998 of 284,546 rows | `jStart == -1` is truthy in Python; the Groovy tested `jStart > -1` |
| Column definitions duplicated in nine places, already drifted | `ChunkQC.py`, `DefaultDBGenerator.py`, `SlimDBGenerator.py`, `ScoreFactory.py`, two static meta files, `BuildDatabase.groovy`, the README, `template.xls` |
| Build needs 64 GB RAM for 42 MB of data | seven `master_table.T.apply(...)` calls plus `iterrows()` over ~200k rows |
| Stage-II CDR3 fixing never runs | `AlignBestSegments.py` is only called from `BuildDatabase.groovy`; the Python path never invokes it |
| `latest-version.txt` names the *previous* release | line 1 is `2026-05-16`; the published latest is `2026-06-03-ZENODO` |
| The dashboard is non-deterministic | a live NCBI eutils call at render time, plus a hardcoded fallback year table frozen since 2021 |
| Five chunk columns discarded | `submitter`, `chunk.id`, `comment`, `meta.subset.frequency`, `method.pairing`, four of them documented in the README |

## 2. Release assets

Three builds, every release, until legacy is retired:

| Build | Asset | Role |
|---|---|---|
| New VDJdb | `vdjdb-<version>.zip` | primary product: polished legacy, with valid JSON everywhere, generated metadata, stable column order, and the derived `cdr3nt`/`d.*`/Pgen block |
| Legacy | `vdjdb-legacy-<version>.zip` | byte-layout-compatible with today's zip; what `vdjdb-web` and standalone clients consume |
| AIRR | `vdjdb-airr-<version>.zip` | AIRR Rearrangement + Receptor/Reactivity |

Plus `manifest.json` (`{role, file, sha256, bytes, members[]}` per build) and `SHA256SUMS`.

Legacy is a derived export: `vdjdb release` → `make legacy` → zip, reading the new build
directory, never `chunks/`. It is a projection of the primary build and cannot drift from it.

Conversion helpers run both ways: new → AIRR and legacy → AIRR.

## 3. Gates

### 3.1 Cross-repo gate: `vdjmatch` patch and release

`vdjmatch/src/vdjmatch/db/vdjdb.py` no longer reads `latest-version.txt`. It calls the GitHub Releases
API and takes the first `.zip` asset:

```python
def _zip_asset(rel: dict) -> str:
    for a in rel.get("assets", []):
        if a["name"].endswith(".zip"):
            return a["browser_download_url"]
```

The moment a release contains more than one zip, vdjmatch picks an arbitrary one. So, in order:

1. Patch `_zip_asset` to select by role: `manifest.json` `role == "legacy"` → `vdjdb-legacy-*.zip`
   → `vdjdb-<tag>.zip` → first `.zip`. Bump `_HF_TAG` in the same pass.
2. Cut a `vdjmatch` release.
3. Only then cut the first multi-zip VDJdb release.

Corollary: the legacy zip's member basenames must stay `vdjdb.txt`, `vdjdb.slim.txt`,
`vdjdb_full.txt`, which is what vdjmatch's member lookup keys on.

### 3.2 Phase −1: unblocking fixes

None of the three is part of the rewrite.

| Fix | Blocks |
|---|---|
| The `vdjmatch` patch above | "three zips per release" |
| Correct `latest-version.txt` line 1; delete the stale `database/` copy | every client that still reads it |
| Replace the seven `.T.apply` calls and the `ScoreFactory` `iterrows` | "runs on 16 GB GitHub-hosted runners" |

### 3.3 Backgrounds

Backgrounds are streamed at build time and only derived statistics (`count.bg`, `total.bg`) ship:
never the background itself, never a subsampled copy. A background is an input to the build, not
an output of it, and a stale vendored copy would change a call set with no error.

## 4. Phases

`master` → `dev` → `feature/*` → `dev` → `master`. Every phase is independently mergeable and
`dev` stays green. Commits that resolve a tracker issue carry `Closes #N`.

| # | State | Branch | Delivers | Closes | Acceptance |
|---|---|---|---|---|---|
| 0 | merged | `feature/dev-baseline` | `ROADMAP.md`, `pyproject.toml` + `uv.lock`, package skeleton, `chunk-check.yml`, `branch-policy.yml` | #476 | CI green; `release.sh` still works untouched |
| 1 | merged | `feature/schema` | the field registry; `render_meta` | - | reproduces the Groovy `METADATA_LINES` / `SLIM_METADATA_LINES` byte-for-byte; `header == meta names` for all three tables |
| 2 | merged | `feature/golden-harness` | `vdjdb diff` + `expected_diffs.toml` | - | zero diffs against the current pandas build; nothing downstream starts without this |
| 3 | merged | `feature/io-qc` | polars reader, vectorised QC, `--strict` exit-1, chunk header normalisation, `.tsv` rename | #497 | QC report matches the pandas report row-for-row; harness still zero |
| 4 | merged | `feature/pipeline-core` | the definitive tables (`records`, `chains`) + harmonize + score + pairing; the legacy export as a projection of them; deletes `py_src/` | #424, #399 | every difference against the release is a declared rule firing its measured count; peak RSS < 8 GB |
| 5 | merged | `feature/arda-cdr3fix` | `arda.cdr3fix` replaces `Cdr3Fixer.py`; retires `res/segments*.txt` | - | new `expected_diffs.toml` rule, row count measured then frozen |
| 6 | merged | `feature/new-format` | ships the definitive tables as parquet + TSV, adds `evidence`, `vdjdb.schema.json` | - | `make legacy` from the shipped tables still passes the harness |
| 7 | merged | `feature/airr` | `emit/airr.py` (Rearrangement + Reactivity), `convert/coords.py`, `vdjdb convert` | - | `airr.validate_rearrangement` passes on the full table; the legacy path produces nothing the tables path does not |
| 8 | merged | `feature/junction-nt`, `feature/segment-guess`, `feature/dgene` | one branch each | #461, #462, #463 | generated `cdr3nt` back-translates to `cdr3` |
| 9 | merged | `feature/harmonize-rules` | nomenclature rule tables | #327, #389, #347, #368, #564, #467, #561 | each rule gets an `expected_diffs.toml` entry with a measured row count |
| 10 | merged | `feature/motifs-tcrnet` | TCRNET on `vdjtools`, streaming backgrounds | - | deviation report accepted |
| 11 | merged | `feature/motifs-tcremp` | TCREMP + per-epitope DBSCAN; new motif schema; legacy projections | - | beats the shipped `cluster_members_tcremp.txt` re-scored in our harness, per §8.4 |
| 12 | merged | `feature/summary` | Rmd split, ggplot2 4.x fixes, committed publication-year table, data-driven callouts, interactive dashboard | #460 | renders offline; perceptual + structural checks pass |
| 13 | merged | `feature/docs` | Sphinx site, generated schema tables, dashboard tab, Pages | - | zero-warning build, deploys |
| 14 | merged | `feature/release-tooling` | manifest, three zips, checksums, `latest-version.txt`, tag scheme, changelog; retires the legacy CI | #432 | full release dry-run with no unattributed differences |
| 15 | part | `feature/aldan3-runner` | self-hosted runner + `build.yml` retargeting | - | identical canonical digests on both runners. `build.yml` carries the `fromJSON(inputs.runner)` retargeting; **no self-hosted runner is registered** (`actions/runners` returns 0), so the second half of the criterion is unmet |
| 16 | part | `feature/identity` | the four derived id levels, the lifecycle record, `vdjdb identity`, promiscuity columns, one study count | - | every invariant of §10.5 passes; a permuted chunk order changes no id; the dashboard reports 638 of 638 references |
| 17 | planned | `feature/corpus` | the reference corpus: documents, vocabulary, postings, `score` and `lift` | - | the three files reproducible by digest; `score` reproduces the `refsearch` ranking; `lift` answers a specificity question with an n |

Phases 0 to 14 are merged to `master` as of 2026-09-27, and phase 15 is half landed: the comparison against the last release
reads PASS with every difference declared and measured, and the release dry-run produces three
reproducible bundles. `ROADMAP_local.md` carries the per-phase record.

Phase 2 came first: the harness had to show zero diffs against the then-current build before any
behaviour changed, so that later differences could be attributed.

Phases 5, 8, 9, 10 and 11 each introduce exactly one source of deviation, so every difference in the
output has a single attributable cause.

Every phase has a step-by-step subplan in §12. A phase is not startable until its subplan names
the files it creates, the facts it needs (already measured, in §7/§8), and the check that closes it.

## 4a. Issue tracker composition

Measured 2026-09-25 with `gh`: 440 issues, 130 open. Grouped by label, the open ones are

| Category | Open | What they are |
|---|---|---|
| data intake | 103 (79 %) | pending papers (79), preprints (9), paper-pending (3), meta-papers (4), 10x/Immudex sets (5), associations (9), other databases (1), correspondence (2) |
| curation quality | 22 | formatting & proofreading (18), typos, structural, validation |
| build infrastructure | 13 | the build, the summary, maintenance |

Some issues have more than one label, so the columns overlap slightly.

Four out of five open issues are a submission queue, not a defect list. This migration closes
issues from the bottom two rows only, thirteen of them, and nothing it does shortens the first
row. Three decisions follow from that:

* `chunk-check.yml` has a three-minute budget, because it is the job the submission queue runs
  through. The full build's 185 s is paid on `dev` and nightly, with nobody waiting.
* the curation skills (`/vdjdb-extract`, `-format`, `-proofread`, `-publish`) are the tooling with
  the largest backlog pointed at it, and they read `proofreading/`, which until phase 9 no build code
  touched.
* a phase that closes an issue number is not thereby reducing the tracker. Progress on the queue is
  curation throughput, and it is measured separately.

## 5. Release comparison

`vdjdb diff <reference-zip> <candidate-dir>` compares in three passes: file set → two digests per file
(raw and canonical) → row-level classification keyed on
`gene|cdr3|v.segm|j.segm|species|mhc.a|mhc.b|antigen.epitope|reference.id`.

Every changed cell must be attributed to a declared rule in `rules/expected_diffs.toml`:

```toml
[[rule]] id="web-unmp-jstart-minus1" file="vdjdb.txt" column="web.cdr3fix.unmp"
         from="no" to="yes" predicate="cdr3fix.jStart == -1" rows=7998
[[rule]] id="meta-tcr-hash-row"   file="vdjdb.meta.txt" added=["TCR_hash"]
[[rule]] id="meta-score-position" file="vdjdb.meta.txt" moved=["vdjdb.score"]
[[rule]] id="meta-web-data-type"  file="vdjdb.meta.txt" rows=4
[[rule]] id="slim-meta-tcr-hash"  file="vdjdb.slim.meta.txt" added=["TCR_hash"]
[[rule]] id="slim-meta-geometry"  file="vdjdb.slim.meta.txt" moved=["v.end","j.start"]
```

The metadata is wrong rather than merely old: `TCR_hash` has been a `vdjdb.txt` column for years
with no metadata row; `vdjdb.score` and `TCR_hash` are listed after `cdr3fix` but stored before
`method`; and the four `web.*` rows have one value too many, landing `0` in `data.type`. All four
`web.*` rows are `visible = 0`, so nothing user-facing moves.

Any unmatched difference fails. A rule that fires a different number of times than declared also
fails, which is why every rule declares a measured row count rather than a description.

### Canonical equality

`runBuidDatabase.py:49` iterates `os.listdir("../chunks")`, which is readdir order and
host-dependent. The released `vdjdb_full.txt` opens with `PMID:28629751` while the
alphabetically-first chunk is `10xgenomics-2019-07-09.txt`. The 2026-06-03 release therefore cannot
be reproduced byte-for-byte by anyone, including the current pipeline.

Canonical equality is the gate; raw equality is informational. The new reader uses
`sorted(glob("*.tsv"))`, `vdjdb release` writes `chunk-order.txt` into the bundle, and
`--chunk-order <file>` replays a recorded order, which makes raw equality achievable going forward.

## 6. Legacy deprecation path

Legacy stops shipping only when all of the following hold:

1. `manifest.json` role selection is released in every known client (`vdjmatch`, `vdjdb-web`).
2. No consumer reads `latest-version.txt` for anything but historical URLs.
3. `vdjdb-web` reads the new format's generated metadata rather than the legacy `vdjdb.meta.txt`.
4. Two consecutive releases have shipped both, with the new format downloaded at a comparable rate.

Until then, `latest-version.txt` line 1 points at the legacy zip, because clients in the wild
download it verbatim and expect that layout.

## 7. Measured facts

Established 2026-09-25 against the 2026-06-03 release and the current `chunks/`. These are inputs,
not drafts; do not re-derive them.

| Quantity | Value | How |
|---|---|---|
| Chunk rows, raw | 203,308 | 230 files, `chunks/*.txt` |
| Chunk rows after per-chunk dedup on `SIGNATURE_COLS` | 192,753, exactly the released `vdjdb_full.txt` row count | polars |
| Rows matching field-for-field across two chunks | 19 pairs, independent reports rather than duplicates: a chunk is one paper | deduplication is within a chunk; these 19 are evidence (§11.1), and global dedup would delete them |
| polars read + dedup of all 230 chunks | 0.4 s | the pipeline it replaces is documented as needing 64 GB |
| Chunks passing `ChunkQC` | 230 / 230, zero errors | fail-fast needs no quarantine list |
| Non-empty CDR3 cells | 305,031, zero with characters outside the 20 AAs | the TCREMP pre-filter is a guard, not a current defect |
| `arda.cdr3fix` vs shipped `cdr3fix`, 20,000-row sample (seed 42) | `cdr3` 97.81 % · `vFixType` 97.95 % · `vEnd` 96.53 % · `good` 95.68 % · `jFixType` 95.44 % · `jStart` 91.02 % | arda 2.27.0 |
| …of the 1,797 `jStart` disagreements | 556 VDJdb-unmapped → arda-mapped · 1,241 both mapped, arda smaller in every case (mode −2, range −1…−8) · 0 coverage regressions | NW with free end gaps vs k-mer longest-hit |
| `arda.cdr3fix.markup_batch` throughput | 27 µs/record → 5.2 s for 191,440 distinct keys | not vectorised internally; fast enough at this scale |
| `vdjtools.model.infer_nt` | 3.11 ms/record → ~15 min for 284,546 rows | human TRB, warm, single-threaded |
| `vdjdb.txt` row/field shape | 284,546 rows, all exactly 22 fields, zero empty `cdr3fix` | keep and assert; raggedness is not a current defect |
| Quote characters in the release tables | present in every `method` / `meta` / `cdr3fix` cell, because they are JSON. What is absent is a *quoted field*: no field begins with `"` | the plan's "zero `\"` characters" was wrong. `quote_style="never"` is still exact: a default CSV writer would wrap every JSON cell and double its quotes |
| `vdjdb_full.txt` rebuilt by the current pandas pipeline vs the 2026-06-03 release | 119,169,153 bytes both, 192,755 lines both, canonical digest identical, raw digest differs | the reproduction contract holds on the largest table; the raw difference is exactly the `os.listdir` chunk order |
| All five positional column orders, registry vs shipped release | `vdjdb.txt` 22 · `vdjdb.slim.txt` 17 · `vdjdb_full.txt` 35 · `cluster_members.txt` 19 · `motif_pwms.txt` 27, order-for-order identical | validates the phase-1 registry against the artifact consumers actually parse |
| `vdjdb_full.txt` `cdr3fix.alpha` encoding | 122,930 non-empty cells, 100 % Python dict repr, 0 JSON | current defect |
| `web.cdr3fix.unmp` correctness | `(no,no)` 268,546 · `(yes,yes)` 8,002 · `(no,yes)` 7,998 | `jStart` is only ever −1 (8,536) or > 0 (276,010); 0 never occurs |
| Murine MHC-II spellings | `I-Ab` 3,274 · `H2-IAb` 113 · `H2-Ab1` 9 · `H2-IAg7` 333 vs `H2-Ag7` 3 · `H2-Aa` 25 vs `H-2Aa` 19 · `H2-Eb1` 7 vs `H-2Eb1` 7 | mostly in `mhc.b`; class I is clean |
| #327 TRAJ24 | bare `TRAJ24` 1,302 (978 carry `WGKLQF`) · `TRAJ24*01` 113 (75 = 66 % carry `WGKLQF`) · `TRAJ24*02` 34 · malformed `TRAJ24-1` 1 | reproduces the report on a larger set |
| #368 antigen gene/species | patch dict covers 245 epitopes = 154,373 of 203,308 rows, zero swapped | the ~90 reported records are in the uncovered tail of 48,935 rows |
| Dashboard R deps | 15 of 17 installed; `maps`/`scatterpie` used only past the embed cut (lines 812–835 vs marker at 567); `ggh4x` never used | splitting the Rmd drops three deps from the release path |
| Shipped dashboard PNG sizes | 1344×960, 2304×1920, 1152×1920, 1536×1536 → `dpi=96, fig.retina=2` | pin these, or the visual fingerprint is noise |
| Production's 27-row `vdjdb.meta.txt` (`vdjdb-web/test/resources/database/`) | orders `… reference.id method meta cdr3fix vdjdb.score TCR_hash web.*` while the data is `… reference.id vdjdb.score TCR_hash method meta cdr3fix web.*` | the metadata mis-describes the data in production too, not only in the release zip |
| `width="1152"` occurrences in the shipped embed HTML | 0. knitr emits the attribute; pandoc's `--embed-resources` drops it when it inlines the image, so the rewrite ran on a stage where it no longer existed | this is what #460 is |

### 7.1 The corpus since the 2026-06-03 release

`chunks/` is the data, so a difference in the release comparison is either code or curation, and only
one of those two can be argued about. This is the accounting that says which. Established 2026-09-27
between tag `2026-06-03-ZENODO` and `master`.

**230 chunk files before, 230 after. None added, none removed. Zero records added, removed or
edited.** Every file has the same number of non-empty lines at both points.

103 files differ, and every difference is formatting:

| Change | Files | Records affected |
|---|---|---|
| CRLF to LF, content byte-identical after the rewrite | 92 | 0 |
| final newline added, file exactly one byte larger | 8 | 0 |
| two prose captions removed from the header of `PMID_24512815.txt` | 1 | 0, and the columns were verified empty in all 270 rows before removal |
| leading empty column name removed from `PMID_40694338.txt` | 1 | 0, over 2,353 rows |
| `Comment` renamed `comment` in `vandesandt-etal-2019-11-04.txt` | 1 | 0 |

All of it landed in one commit, `b0a479d`, whose 137,538 insertions and 137,538 deletions are equal
because a line-ending rewrite touches every line and adds none.

So every difference the release comparison reports is attributable to the build, which is what makes
`rules/expected_diffs.toml` meaningful: a rule there describes a code change, and there is no
curation change hiding behind it.

Re-derive it with the two properties that matter, rather than by reading a diff stat:

```bash
# files added or removed
diff <(git ls-tree -r --name-only 2026-06-03-ZENODO chunks/) \
     <(git ls-tree -r --name-only HEAD chunks/)
# per file: does the content differ once line endings are normalised?
for f in $(git ls-tree -r --name-only HEAD chunks/); do
  a=$(git show "2026-06-03-ZENODO:$f" | tr -d '\r' | shasum -a 256 | cut -c1-16)
  b=$(git show "HEAD:$f"              | tr -d '\r' | shasum -a 256 | cut -c1-16)
  [ "$a" = "$b" ] || echo "content differs: $f"
done
```

The second loop reports the three header repairs and the eight trailing-newline files, and nothing
else. `CLAUDE.md` carries the rule this section exists to serve: a commit touching `chunks/` says
which files, how many rows, why, and who decided.

## 8. Motif inference findings

Measured 2026-09-25 against the 2026-06-03 slim dump and the shipped motif files.

### 8.1 `vdjtools.overlap.tcrnet` p-value and the legacy statistic

`overlap/tcrnet.py::_score` computes `E = (n_target / max(m_control,1)) * n_control` with no
pseudocount, then `p_enrichment = poisson.sf(degree - 1, E)`. When `n_control == 0`, `E == 0.0`, and
`scipy.stats.poisson.sf(k-1, 0.0) == 0.0` exactly for every k ≥ 1 (verified). So every clonotype
with at least one within-sample neighbour and no background neighbour receives `p = 0.0`, `q = 0.0`,
and sorts to the top of the result.

Measured on human TRB (111,407 unique CDR3s) against the bundled 250k control: `n_control == 0` for
81.0 % of queries; 44,010 of 111,407 rows (39.5 %) get `p_enrichment == 0.0` exactly.

A bigger background does not fix it: mean `n_control` scales linearly with M, so `E` is M-invariant
in expectation and M only controls the zero-inflation (human TRA: 34.7 % zeros at M = 250k, 15.6 % at
M = 2,266,274).

Call `tcrnet()` for `n_control` only, and recompute the legacy statistic, whose pseudocount removes
the pathology:

```
p        = (n_control + 1) / (M + 1)
p_legacy = binom.sf(d_s - 1, N, p) / (1 - (1 - p)**N)
```

This is bit-faithful to the Groovy `DegreeStatisticsAnnotator.computePValue`. File the `E = 0` bug
upstream against `vdjtools`.

### 8.2 Corrections to the plan

- The legacy TCRNET was not V/VJ/VJL-grouped. `CalcDegreeStats.groovy` defaults to
  `-g dummy`, and the legacy Rmd passes only `-o 1,0,1 -g2 dummy`. The new grouping-free `tcrnet()`
  is an exact scope match, not a deviation.
- The 0.941 / 0.569 / 0.709 figure is not a knee-DBSCAN result. It is the shipped
  Leiden/REDCEA production table re-scored, "nothing re-fit". Measured knee-DBSCAN on mirpy
  embeddings traces a strictly lower frontier (pooled human TRB ≥30: 0.910 / 0.368 / 0.524 at
  coef 0.75). The target is the shipped file re-scored in our own harness, not the write-up.
  `coef = 0.75` was calibrated on a different embedding (standalone `tcremp`, ~3000 OLGA
  prototypes, Smith-Waterman) and does not transfer to mirpy's; it must be re-calibrated.

### 8.3 Kneedle at production scale

Pooled human TRB, 112,983 unique clonotypes, PCA(50): Kneedle returns knee = 1 of 112,983
(frac 0.000). The `floor_frac = 0.40` guard fires always, so the operating rule is
`eps = coef × mean(1st-NN distance)` and Kneedle is only a debug cross-check. Do not describe the
method as Kneedle-based.

The two implementations estimate different quantities and must not be conflated:
`knee_eps.py` sorts each column independently and reads column 1, the sorted 1-NN curve, scaled
by 0.75. `mir.bench.metrics.estimate_dbscan_eps` uses the 4-NN curve, scaled by 0.4, with no
degenerate guard. vdjdb-db owns one parameterised implementation, defaulting to the `knee_eps.py`
variant (the one the benchmark numbers came from), and reports both.

### 8.4 Clustering topology

Measured, human TRB ≥30, same embedding, vendored `metrics_lib` semantics:

| topology | purity | retention | cov-F1 |
|---|---|---|---|
| pooled DBSCAN, coef 0.75 | 0.910 | 0.368 | 0.524 |
| pooled DBSCAN, coef 1.10 | 0.787 | **0.545** | **0.644** |
| per-epitope, chain-global eps, coef 1.10 | **1.000** | 0.406 | 0.578 |
| per-epitope, coef 1.30 | **1.000** | 0.452 | 0.623 |
| TCRNET (shipped) reference | 0.950 | 0.230 | 0.370 |

Per-epitope DBSCAN with a chain-global eps is the selected topology: pooled geometry, per-epitope
scope. The eps comes from the pooled k-distance curve, so it is never re-estimated on a per-epitope
n of 30–300, which is the degenerate regime. Purity is 1.000 by construction, because a cid cannot
span epitopes, which is what the output format requires anyway.

Purity 1.000 therefore means cross-epitope contamination is impossible, not that none exists. The
pooled clustering runs in the debug build as the measurement, and its cross-epitope confusion matrix
is a first-class QC artifact.

### 8.5 Defects in the shipped motif files

- 31 of 620 cids in `cluster_members.txt` have no rows in `motif_pwms.txt`, so the web shows those
  clusters in the tree with no logo. Cause: `filter(total.bg > 0)` after a left join against the
  frozen Zenodo background PWM.
- 1.00 % of letter mass is deleted, with no error, across 253 of 589 clusters. The same filter drops
  any `(pos, aa)` absent from the background at that `(v, j, len)`, which are the rarest and most
  informative residues. Per-position `freq` then sums to < 1 (minimum observed 0.111) and
  `I = 1 + Σ freq·log freq / log 20` is computed on a truncated distribution, inflating I.
  `need.impute` is computed *after* the filter, so it is `FALSE` on all 13,456 rows and the
  imputation machinery is dead code.

Fix: a three-level Laplace-smoothed cascade (`vj_len` → `len` → `uniform`) so nothing is ever dropped,
`freq` sums to 1, and `need.impute` becomes a provenance flag. Assert the sign: at the 253
affected clusters the new `Σ I` must be lower. A rewrite that reproduces the inflated values has
reproduced the defect.

Separately, `py_src/MotifsScoresAssembler.py:44-56` sets `cluster.member = 1` for every record of an
(epitope, species, gene): it indexes on those three columns and never on `cdr3`. So
`vdjdb.slim.scored.txt` and `vdjdb.scored.txt` are wrong.

### 8.6 Variable CDR3 length and cluster ids

`vdjdb-web` already splits every cid by `len` before building a cluster
(`Motifs.scala:85-89`, `splitOn(cid)` then `splitOn(len)`), and `MotifCluster` has `vsegm`/`jsegm`
as sequences, so the legacy schema is not the obstacle. Relying on the relaxed `strict = false`
path, however, produces two display clusters sharing one `clusterId`.

Design: a cluster owns one member set and one or more PWM strata (one per CDR3 length with ≥3
members). The default legacy projection emits one legacy cid per (cluster, stratum), `cid =
H.B.<epitope>.<n>L<len>`, so both files satisfy `strict = true` and no two display clusters share a
`clusterId`. `--legacy-shared-cid` reproduces today's shape if the web maintainer prefers it.

### 8.7 Background choice

Use `{species}.{locus}.aa.vdjtools.tsv.gz`, not `*.ntvj` and not `vdjdb-web-control/*`, which is the
opposite of the choice argued in `vdjdb-web-control/README.md`. There the query is a user's
nucleotide-level repertoire, so the control must be nt-level. Here the query is a VDJdb epitope
sample, deduplicated to unique CDR3aa, and `seqtree.Index` is built over aa strings. An ntvj control
inserts the same `cdr3aa` once per nt variant; that inflation is 1.63× for VDJdb-matching clonotypes
against 1.04× overall, concentrated on the germline-proximal public sequences under test. It would
inflate `E` by ~1.6× for those sequences and under-call public motifs.

`seqtree.control.load_control` already streams from `isalgo/airr_control`, filters to the productive
20, reservoir-samples uniformly over unique clonotypes, and content-addresses its download; do
not reimplement it. What it stores is the fetched table, which is an input; no derived control is ever
written (hard rule 9), and the reservoir sample is deterministic in the seed. Never let `tcrnet()`
resolve its own background: it calls `evalue.background(locus, species)` with no `size`, which
indexes the full table.

### 8.8 Memory and throughput

TCREMP embedding, human TRB, measured: single-shot `embed` → `StandardScaler.fit_transform` → `PCA`
peaks at 9.77 GB (12.41 GB with a full-matrix scaler fit). Fitting the scaler and PCA on a 25k
seeded subsample and then embedding in 20k chunks peaks at 3.26 GB. Chunking is mandatory rather
than an optimisation: unchunked it would OOM a 16 GB runner alongside the rest of the build. Chunks
run sequentially; one internally-threaded `embed()` per chunk, never a worker pool.

Throughput: 75,308 rows/s on 16 cores, 61,605 rows/s at `threads=4`; the 113k-row human TRB
embedding takes 1.4 s. Motif stage overall: ~30–60 min cold, ~10–15 min warm, peak 3.3 GB, against a
six-hour CI limit.

Nothing is stored between builds (hard rule 9): the scaler and the PCA are re-fitted every run
from the seeded 25k subsample, which is deterministic in `config.SEED` and the input, so storing them
would buy minutes and risk shipping a fit that no longer matches the data.

A stored artifact was considered for a different reason: raw igraph component numbers are not
stable across releases, so every bookmarked `vdjdb.com` motif URL breaks each time. Cluster ids are
instead derived from cluster content, the same choice `clonotype_id` already makes (a seeded hash of
the clonotype key, never a counter), so a cluster whose membership is unchanged keeps its id with
nothing stored.

### 8.9 `mir.bench` helpers

- `mir.bench.vdjdb.load_vdjdb` keeps only `junction_aa, v_call, j_call, locus, epitope, mhc_class`;
  it drops `species`, `mhc.a`, `mhc.b`, `antigen.gene`, `antigen.species`, `v.end`, `j.start`, seven
  of the 19 legacy columns. Read the slim table directly.
- `mir.bench.metrics.cluster_metrics` uses `recall = tp / n_true_clustered` (clustered records only);
  the vendored `metrics_lib.precision_recall_fscore` folds unclustered records into FN. The two
  metric families are not comparable. Validation uses the vendored `metrics_lib`, verbatim.

## 9. Decisions

Settled 2026-09-25.

- The new format takes ownership of `evidence.*`. Production `vdjdb-web` serves five
  `evidence.*` columns that nothing in this repo produces. The new format produces them, from the
  evidence table in §10.
- The five "dropped" chunk columns are kept. `submitter`, `chunk.id`, `comment`,
  `meta.subset.frequency` and `method.pairing` are all needed for debugging a curation problem, so
  they go into `records.parquet` rather than being discarded.
- The side outputs are produced but not zipped. `vdjdb_full_filtered.txt`, the three
  `*_broken.txt` and the three `*_scored.txt` tables are written to `out/reports/` and uploaded as
  CI artifacts; whether any of them ever ships is deferred.
- Everything the build produces is specified in `docs/outputs.md`, which becomes
  `docs/standards/database-outputs.rst` when the Sphinx site lands.
- How the motif stage is tuned, and what it is for, is specified in `docs/denoising.md`: normative
  for any change to motif clustering, and the source of the two-stage acceptance rule the shipped
  parameters are chosen by. Published as-is on the docs site by `myst-parser`.
- Which clustering algorithms exist, what their parameters do, and what each measures at is
  `docs/clustering.md`: connected components, CPM Leiden, DBSCAN, HDBSCAN, Lumbermark and a
  TCRNET-gated hybrid, with the scorecard that made two of them the default and the rest measured or
  rejected. Both documents use MathJax, so `docs/conf.py` enables
  `myst_enable_extensions = ["dollarmath"]` and publishes them as Markdown rather than converting:
  they are cited by section number from the ROADMAP, the package docstrings and `docs/tuning/`.
- The measurements behind both are committed, in `docs/tuning/`: `scorecard.tsv` is 252
  configurations scored through one harness, and `sweeps.py` / `report.py` / the two `.gp` files
  regenerate it and its figures. They are committed, reviewed inputs to the documents that cite
  them, refreshed by their own pull request and never written by a build (hard rule 9). The sweeps
  need `uv sync --extra tuning`; Lumbermark is measured, not shipped, so it stays out of the build
  container.

## 10. Record identity and the evidence model

### 10.1 Record identity

A VDJdb record has had no identifier: it is a row in a chunk, named only by its content. A curator
fixing a CDR3 typo therefore deletes one record and creates another, indistinguishably from a
deletion plus an addition, and external references, structure links and accumulated evidence break
with no error.

Identity is two-level, in `src/vdjdb/identity/`:

| | |
|---|---|
| `record_id` | opaque, assigned once, never reused, stable across content changes. `VDJDB` + 10 digits. This is what external references point at |
| `content_hash` | sha256 over the canonical content. Changes whenever anything changes. Detects change; never identifies |

Reconciliation, most specific first: exact natural key → amendment (exactly one differing
key field, unambiguous, same chunk and reference; ambiguity is refused, because a wrong link is
worse than a new id) → allocation. Registry entries the build no longer sees are retired,
not deleted. The registry is a committed TSV sorted by `record_id`, so a curation PR shows added,
amended and retired records as a reviewable diff.

Two earlier keys were wrong in opposite directions, and tests now pin the boundary: a narrower key
stopping at `reference.id` collided on 20,769 of 192,753 records (one paper reporting the same TCR
in several donors is several records), and a key without `chunk.file` merged 19 pairs that are two
papers' independent reports.

The natural key is `CHUNK_DEDUP_KEY` plus `chunk.file`. A chunk is one paper, so two matching
rows in two chunks are two independent reports and must keep separate ids. Ids are assigned
before CDR3 repair: two trimmed sequences that repair to the same full one are still two
observations, and assigning afterwards merged 215 pairs the publications reported separately.

Measured: 192,753 rows → 192,753 records; ids stable across rebuilds; a typo fix reported as
`VDJDB0000000101 cdr3.beta: CASSIRSSYEQYF -> CASSIRSSYEQYFF`.

### 10.2 Evidence model

Four tables, normalised so nothing is stored twice (full column lists in `docs/outputs.md`):

```
records.parquet    PK record_id             one row per curated record, with provenance
chains.parquet     PK (record_id, gene)     one row per TCR chain; carries clonotype_id
evidence.parquet   PK (record_id, evidence_id)   long format, one row per piece of evidence
vdjdb.parquet      the joined, pivoted view -- derived, never authored
```

Evidence is long rather than wide because a record can have any number of pieces of any number of
kinds; a wide table would be mostly null. `evidence_type` covers `motif_tcrnet`, `motif_tcremp`,
`structure_native`, `structure_model` and `independent_study`. Structure evidence is keyed on the
legacy `TCR_hash` until the structure store is re-keyed on `record_id`.

`chains` exists so the schema stays non-redundant: folding chains into records forces either
duplicated record fields (what `vdjdb.txt` does) or paired alpha/beta columns (what `vdjdb_full.txt`
does).

### 10.3 Five levels of identity, and which one is allocated

A record is a receptor against a presented peptide. Both halves recur across records and consumers
reference both halves, so each level needs an identifier of its own. Counts are from the current
build.

| Level | Identifies | Key | Distinct | Id |
|---|---|---|---|---|
| clonotype | one receptor chain | `species, gene, cdr3, v.segm, j.segm` | 187,935 over 286,047 chain rows | `CT` + 16 hex, derived |
| clone | one TRA/TRB pair | the two `clonotype_id`s, sorted | 82,266 over 93,294 paired records | `CX` + 16 hex, derived |
| pMHC | one presented peptide | `antigen.epitope, mhc.a, mhc.b` | 2,364 | `PM` + 16 hex, derived |
| epitope | the peptide alone | `antigen.epitope` | 2,118 | `EP` + 16 hex, derived |
| record | one curated line | `NATURAL_KEY`, §10.1 | 192,753 | `VDJDB` + 10 digits, allocated |

99,459 records carry one chain and 93,294 carry two, so `clone_id` is the empty string on slightly
more than half of them, and having no clone is the curated state rather than missing data. Empty
rather than null, because empty is the only missing marker (hard rule 6). `mhc.class` adds nothing to the
pMHC key, which is 2,364 distinct with or without it, because the alleles determine the class; it
stays a derived column.

**Derived against allocated is what settles the chunk-order question.** A derived id is a hash of its
own key and of nothing else, so:

* the order `chunks/` is read in cannot reach it;
* adding a chunk names that chunk's new clonotypes and changes no existing id;
* removing a chunk retires only the ids nothing else supported;
* two hosts, two core counts and two release dates produce the same ids.

An allocated id is a counter and has none of those properties. `record_id` is allocated anyway,
because its purpose is to survive a content change: a curator fixing a CDR3 typo has to keep the
record's id, and no hash of the content can do that. Every other level keys on content with no
separate existence, so a change there is not an amendment but a different clonotype. One registry is
therefore consulted during a build, the record registry, and the other four levels need no history
in order to be correct.

**The hash has to be ours.** `clonotype_id` is `pl.Expr.hash` today, which is polars' xxhash under
our seed. Polars does not specify that output across versions, so an upgrade could renumber every
clonotype with nothing failing anywhere. `vdjdb.identity.ids._hash` is already sha256 over
`\x1f`-joined fields, and the two have to be one function. Measured: sha256 over the 187,935 distinct
clonotype keys costs 83 ms, joining back onto 286,047 chain rows costs 7 ms, against a 185 s build,
with zero collisions at 16 hex digits.

Truncating to 16 digits is a size decision, not a security one. At 187,935 ids the chance of one
collision is around 1 in 10⁹, and the build asserts distinctness per level, so a collision fails a
run rather than corrupting a join.

### 10.4 Lifecycle: what happened to the id I had

An id that vanishes is the failure a consumer cannot diagnose: a reference that used to resolve
returns nothing, and there is no way to tell a correction from a deletion. Every level therefore
carries a lifecycle row, holding only what a build cannot recompute.

| Field | Meaning |
|---|---|
| `id` | the identifier |
| `level` | `clonotype`, `clone`, `pmhc`, `epitope` or `record` |
| `state` | `active` or `retired` |
| `first_release` | the release tag the id first appeared in |
| `last_release` | the newest release carrying it, frozen at retirement |
| `replaced_by` | the id that took over, when the retirement was an amendment; empty otherwise |

The active set and its keys are recomputed from `chunks/` on every build (hard rule 9), so a
lifecycle row never decides what id anything gets. It answers one question, and
`vdjdb identity resolve CT3f9a1c0e8b2d4a67` is that question: active with its key, or retired with
the release that last carried it and what replaced it.

**Retirement is per release, not per build.** A curation branch can add a clonotype and remove it
again before anything ships, and neither event is a lifecycle event, because only the release job
writes lifecycle rows. That is also what keeps a curation pull request readable: nothing in these
files moves when a chunk changes.

Sizes. The four derived levels hold 274,683 ids between them, and the table measures **15.9 MB**.
The record registry is 72.7 MB, because it carries the natural key and the content hash per record.
Both ship as release assets rather than committed files. The retired rows alone are committed, since
they are the part a consumer needs and the part that is otherwise lost, and they are small.

A build with no registry still runs. It allocates record ids from 1 and reports that the run is not
id-stable, which is correct for a fork and for a first build, and is what the code does today by
accident rather than by decision (`ROADMAP_local.md`, known-not-yet-fixed).

### 10.5 The consistency machinery

Four commands and one report shape.

| Command | Does |
|---|---|
| `vdjdb identity build --tables out/tables` | writes the id columns and the lifecycle rows for this build |
| `vdjdb identity resolve <id>` | one id: state, key, first and last release, replacement |
| `vdjdb identity diff <previous> <current>` | added, retired and amended per level, with counts |
| `vdjdb identity check --tables out/tables` | the invariants below; exit 1 on any failure |

`identity check` is what CI runs, and it asserts:

1. every id is distinct within its level;
2. every derived id equals the hash of the key on its own row, recomputed rather than trusted;
3. the ids do not change when `chunks/` is read in a permuted order, nor when a chunk is added and
   removed again;
4. every `record_id` in the previous registry whose natural key is unchanged is unchanged;
5. no id carries a prefix belonging to no level, and no level that was published empties
   silently. Reuse would be an id resolving to a *different* key than the one it was retired under,
   and invariant 2 recomputes every id from the key beside it, so an id pointing anywhere else fails
   there. An id merely reappearing is not reuse: a derived id is the hash of its key, so the same key
   returning has to give the same id;
6. every `clonotype_id` resolves to at least one chain row, every `clone_id` to exactly two chains
   of different `gene`, every `pmhc_id` to at least one record;
7. `TCR_hash` is byte-identical to the previous release on every record whose key is unchanged.

Invariant 3 is the one nothing covers today and the one the chunk-order question asks for. It runs
the assembly twice over a permuted chunk list and asserts one digest per level. It belongs on a
small fixture rather than a full build, so it costs seconds and can run per pull request.

### 10.6 MHC promiscuity is an annotation, never part of an id

An epitope is often presented by several alleles, and the curated allele is not always the one that
binds it best. Both facts belong in the database; neither belongs in a key.

`pmhc_id` keys on the allele **as curated**. A prediction is a moving target: `mhcmatch` ships new
weights, and if a predicted allele were in the key then a model upgrade would renumber pMHC ids and
break every external reference while the build passed. Promiscuity is therefore measured into columns
on `restriction`, which already carries one row per (epitope, antigen species, MHC pair) and 2,373 of
them:

| Column | Meaning |
|---|---|
| `alleles.reported` | distinct alleles VDJdb reports for this epitope |
| `mhc.a.top` | the highest-scoring allele for this epitope under `mhcmatch` |
| `mhc.a.rank` | the curated allele's rank in that ranking |
| `mhc.a.percentile` | the curated allele's binding percentile |
| `promiscuity` | alleles scoring within a declared percentile of the top |
| `mhcmatch.version` | the model that produced the five columns above |

`mhcmatch.version` sits on the row, so a reader can tell which model a number came from, and
re-running under a new model rewrites those columns and no id. Phase 9e already calls `mhcmatch` for
catalogue validation; this keeps its output instead of discarding it.

### 10.7 The legacy structure id is preserved, not redefined

`TCR_hash` is the identifier linking a record to a structure. It comes from the legacy recipe and
nothing else: sha256 over a fixed field list, empty unless every required field is present
(`assemble.master.add_tcr_hash`). That recipe reproduces the values in the released `vdjdb.txt`
exactly, which is why it is preserved rather than redefined, and it is never recomputed under a new
scheme, renamed, or derived from the levels above. Structure evidence keys on it.

The levels above are additive, so a consumer holding a `TCR_hash` is never asked to migrate, and
invariant 7 fails the build if one moves on a clonotype whose key did not.

## 11. Tuning and validation of motif clustering

The tuning set and the validation set are separate datasets and stay separate.

### 11.1 Tuning objective: independent-study support

The signal is how many distinct studies report the same clonotype against the same epitope. Measured
on human records in `chunks/`: 4,129 of 187,238 clonotype-epitope pairs (2.21 %) are supported by
≥2 distinct `reference.id`, spread over 18 epitopes with ≥20 such pairs: GILGFVFTL 2,136,
YLQPRTFLL 619, NLVPMVATV 170, GLCTLVAML 96, RAKFKQLL 61, and a long tail.

A clustering that finds convergent selection should preferentially recover those
independently-replicated clonotypes. That is the objective the DBSCAN `coef` (and any Leiden
resolution) is fitted against, per chain, on the pooled geometry.

Independent replication is also a per-record evidence type in its own right, so the tuning signal and
a shipped evidence column come from the same computation.

Three qualifications, all measured:

*The denominator excludes display-selected records.* A library panned against one pMHC yields
thousands of receptors one substitution apart by construction, all under a single `reference.id`, so
every display clonotype contributes zero independently-replicated pairs while occupying a denominator
slot. Human TRB: all 2,771 replicated pairs sit outside the display set, which fills 29,692 of 116,053
slots, and the same clustering reads lift 3.585 or a much lower figure depending only on which cohort
is scored. Report a lift figure with the cohort it was scored on.

*Lift alone selects for narrowness.* Lift falls monotonically as a clustering's radius widens, so
maximising it prefers the tightest radius, which finds clusters only where the data is densest, i.e.
in fewer epitopes. Measured: the TRB configuration that maximises pooled lift covers 88 of 178
epitopes against the shipped annotation's 103. An epitope with no motif receives no denoising, so
epitope coverage is an admissibility requirement alongside `Q`, purity and precision, not a
tiebreak that can be traded away. `docs/denoising.md` §7.1 is normative; `docs/clustering.md` gives
the frontier for each algorithm.

*The admissibility axes do not, by themselves, exclude doing nothing.* Score the partition that puts
every clonotype of an epitope in one cluster and excludes none: on human TRA it returns `Q` 0.7940,
purity 0.9204, precision 0.9237 and coverage 118 of 118, clearing all four axes, and on TRB
`Q` 0.8632, purity 0.9343, coverage 178 of 178, failing only the legacy-relative purity bar. Its lift
is 1.000 by construction. So `Q` and coverage are guards against shattering rather than evidence of
quality, and an absolute purity floor has to be measured against that partition, per chain per
build: on TRB that is 0.9343, which is why a floor of 0.93 is not a floor. The floor is 0.94 on both
chains, and it is usable only on TRB: the TRB window is [0.9343, 0.9829] while TRA's is empty by
0.0030 of purity, so TRA is guarded by the legacy-relative bar and by stage 2 instead. Setting it
moves no shipped parameter: `coef` 1.55 on TRB still has the highest admissible lift.
`vdjdb.validate.motif_bench.trivial_members` builds the partition; `docs/clustering.md` §8 is the
audit and both windows.

### 11.2 Validation set: TCRvdb

**TCRvdb / MATCHMAKERS is proprietary (academic, non-commercial, no redistribution in whole or in
part) and is never shipped, committed or redistributed.** Read it only via `$VDJDB_TCRVDB`;
`src/vdjdb/validate/guard.py` fails the build if a copy reaches the repository or a bundle, matching
both filename and column fingerprint so a rename does not defeat it.

It is a held-out set: 614 labelled paired HLA-A2 records over exactly two epitopes, YLQPRTFLL
(421 labelled, 40.9 % MM+) and GLCTLVAML (193, 69.4 % MM+). Label definition taken from the author's
own `validation_analytics.Rmd`: MM+ is `padj < 1e-5`, joined on `cdr3` alone with alpha and beta
melted and `min(padj)` per CDR3.

Tuning on it would invalidate it, and two epitopes are too narrow a basis for a global
hyperparameter. Only aggregate metrics may be reported: recall of MM+ at a fixed operating
point, precision against MM-, AUROC. Per-record verdicts may never be emitted.

## 12. Phase subplans

One subplan per phase. Each names the files it creates, the facts it consumes (§7, §8, already
measured, never re-derived), and the single check that closes it. A phase whose check is green is
merged to `dev` and struck here.

Minor decisions taken while executing a subplan are recorded in the execution log rather than escalated.

### Phase 1 - `feature/schema`

1. `src/vdjdb/schema/fields.py` - one `Field` record per column: `name`, `dtype`, the eight
   `vdjdb.meta.txt` attributes (`type`, `visible`, `searchable`, `autocomplete`, `data.type`,
   `title`, `comment`), and its position in each of the four positional orders (§2 of
   `docs/outputs.md`). This is the single declaration the six duplicated column lists collapse into.
2. `render_meta(table)` / `render_slim_meta(table)` - emits `vdjdb.meta.txt` /
   `vdjdb.slim.meta.txt` as text. There is no legacy mode: the shipped metadata does not describe
   the file it belongs to, so reproducing it would ship a known defect. The three fixes are declared
   in `expected_diffs.toml` instead (§5).
3. `header(table)` - the `vdjdb.txt` header derived from the same declaration, restoring the
   `BuildDatabase.groovy:411` invariant the Python port dropped.
4. Tests: parse `src/BuildDatabase.groovy`'s `METADATA_LINES` / `SLIM_METADATA_LINES` at test
   time (never a copy) and assert every difference from `render_meta` is one of the three declared
   fixes; assert `header(t) == [f.name for f in fields(t)]` for all three tables.
5. `vdjdb schema --table {vdjdb,slim,full} --format {meta,header,json}` on the CLI.

**Closes when:** every difference between `render_meta` and the Groovy constants is one of the three
declared fixes, asserted by a test that parses `BuildDatabase.groovy` rather than copying it. No
pipeline behaviour changes.

### Phase 2 - `feature/golden-harness`

1. `vdjdb fetch-reference --tag 2026-06-03-ZENODO --out ref/` - `gh release download`, unzip. An
   input being fetched, not a result being cached.
   `ref/` is gitignored; the harness takes a path, so CI passes an artifact instead.
2. `src/vdjdb/compare/diff.py` - three passes: file set → raw + canonical sha256 per file → row-level
   classification keyed on `gene|cdr3|v.segm|j.segm|species|mhc.a|mhc.b|antigen.epitope|reference.id`,
   bucketed `only-in-reference` / `only-in-candidate` / `changed` / `identical`.
3. `rules/expected_diffs.toml` - every changed cell must match a declared rule and the rule must
   fire exactly its declared `rows` count. Seed it empty: phase 2 starts at zero rules and zero diffs.
4. `vdjdb diff <reference> <candidate> --report out/reports/diff-report.md`, exit 1 on any unattributed
   difference.
5. Run the current pandas build into a scratch dir and diff it against the reference zip.

**Closes when:** the current pandas build diffs to zero against the 2026-06-03 release under canonical
equality, with an empty rule file. Nothing downstream starts before this is green.

### Phase 3 - `feature/io-qc`

1. `src/vdjdb/io/chunks.py` - `read_chunks(dir)` → polars, `sorted(glob(...))` for determinism,
   per-chunk dedup on `CHUNK_DEDUP_KEY`. Target: 192,753 rows (§7).
2. Header normalisation as a one-shot migration commit: CRLF → LF on the 99 files, the prose column
   name in `PMID_24512815.txt`, the bare leading tab in `PMID_40694338.txt`, `Comment` → `comment`,
   and `.txt` → `.tsv` on all 230 (#497). One commit for the rename, one for the content, so
   `git log --follow` still works.
3. `vdjdb qc --strict` exits 1 on any error-level finding, restoring the Groovy behaviour that
   `warnings.warn` replaced. Chunk-level rules ported from `ChunkQC.py` as polars expressions.
4. `chunk-check.yml` switches to `--strict` once the migration commit lands.

**Closes when:** the QC report matches the pandas report row-for-row, and the phase-2 harness is still
zero after the `.tsv` migration.

### Phase 4 - `feature/pipeline-core`

1. `curate/patch.py` - the antigen patch as one join, replacing the per-chunk `T.apply`.
2. `score/confidence.py` - `ScoreFactory` ported to polars expressions; the signature max becomes
   `.max().over(CHUNK_DEDUP_KEY)`.
3. `assemble/master.py` - pairing, `complex.id` allocation, `samples.found` / `studies.found` as
   `len().over(...)` / `n_unique().over(...)`.
4. `emit/legacy.py` - the three legacy tables. Two byte-compatibility helpers belong here and nowhere
   else, each with a comment saying so: `json_column()` (`map_elements(json.dumps)`, because
   `struct.json_encode()` emits `{"a":"x"}` where the release has `{"a": "x"}`) and
   `py_repr_column()` (`vdjdb_full.txt`'s `cdr3fix.*` are Python `dict` repr).
5. Delete `py_src/` in the same commit that makes it redundant, not before. (Done: the last
   thing needed from it was the CDR3 fixer, vendored into `annotate/_legacy_fixer/` until phase 5.)
6. A peak-RSS test: the full build under 8 GB, asserted with `resource.getrusage`.

**Closes when:** the comparison against the release shows only the 7,998 `web.cdr3fix.unmp` rows (§7)
plus the meta-file fixes, each as a declared rule with its measured count.

### Phase 5 - `feature/arda-cdr3fix`

1. Dedupe to distinct `(species, cdr3, v, j)` - 191,440 keys - and make one `markup_batch` call
   per organism, then join back. 5.2 s total (§7); never loop it.
2. `Cdr3Markup.to_cdr3fix()` emits VDJdb's JSON key-for-key; `v_end` / `j_start` are junction-space,
   which is what the `cdr3` column contains.
3. Delete `src/vdjdb/annotate/_legacy_fixer/`, the verbatim copy phase 4 bridged through,
   together with `res/segments.txt` and `res/segments.aaparts.txt`.
4. Measure the difference against the release, then freeze it as declared rule counts. Expected
   shape from §7: `cdr3` ~2.2 % of rows, `jStart` ~9 %, all of it in the direction arda maps more
   and earlier.
5. Assert the one-directional property as a test: arda never loses coverage a VDJdb mapping had.

**Closes when:** the comparison's only new differences are inside the declared `arda-cdr3fix` rules,
and the coverage-regression count is 0.

### Phase 6 - `feature/new-format`

1. `identity/` is already built (§10.1) - wire it into the build so every record carries `record_id`.
2. `emit/vdjdb3.py` writes `records.parquet`, `chains.parquet`, `evidence.parquet` and the derived
   `vdjdb.parquet`, plus a TSV projection of each (`docs/outputs.md` §3).
3. `evidence.parquet`'s first producer is `independent_study` - the same computation as the §11.1
   tuning signal, so one implementation serves both.
4. `emit/legacy.py` is re-pointed to read the new build directory, never `chunks/`.
5. `vdjdb.schema.json` generated from the phase-1 registry.

**Closes when:** `vdjdb make legacy` from the new build still passes the phase-2 harness, which shows
the legacy export is a projection rather than a parallel implementation.

### Phase 7 - `feature/airr`

1. `convert/coords.py` - the only module that converts between the four coordinate spaces
   (section 0 and `README.md`), with round-trip tests. Junction ↔ CDR3 is the two-anchor offset; everything else is
   0-based/1-based and half-open/closed.
2. `emit/airr.py` - Rearrangement + Receptor + Reactivity, both projections of the phase-1 registry's
   AIRR mapping so the two converters cannot diverge.
3. Validate with the `airr` package's own schema validator in CI.
4. Property test: `legacy_to_airr(legacy) == vdjdb3_to_airr(new)` on the fields AIRR can represent.

**Closes when:** `airr.validate_rearrangement` passes on the full table and the round-trip property
holds.

### Phase 8 - `feature/junction-nt`, `feature/segment-guess`, `feature/dgene`

One branch each; all three write new-format columns only, so the harness stays green by
construction.

1. junction-nt (#461): `vdjtools.model.infer_nt` on the unique `(species, cdr3, v, j)` set, four
   big contiguous slices, never a per-record pool. 3.11 ms/record → ~15 min (§7). No cache: the
   ~15 min is inside the budget and the output is authoritative data, not a derived convenience
   (hard rule 9). Test: the generated `cdr3nt` back-translates to the input `cdr3`.
2. segment-guess (#462): kmer candidates vectorised, ties broken by one `pgen_aa_batch` call.
3. dgene: `arda.dpost.posterior_d` (human IGH/TRB/TRD + mouse TRB only; it returns `None`
   elsewhere rather than guessing, and that `None` must be preserved, not defaulted).
4. TCR_hash (#463): keep the legacy hash as-is so structure evidence keeps resolving; the
   re-keying on `record_id` is §10.2's deferred half.
5. AIRR `Receptor` rides along: `vdjtools.model.stitch_*` gives the complete mature variable
   domain its two required columns need, and `receptor_hash` is a sha256 over those.

**Closes when:** each branch's new columns are populated, the harness is unchanged, and the
back-translation test passes.

### Phase 9 - `feature/harmonize-rules`

One rule table, `curate/rules/`, one entry per defect, each with its own `expected_diffs.toml` rule
and measured row count. The counts are already in §7 - do not re-measure:

| Issue | Rule | Rows |
|---|---|---|
| #327 | `TRAJ24` → `*02` where CDR3 carries `WGKLQF` | 75 explicit `*01` + 978 bare `TRAJ24` |
| #389 | `TRAJ24-1` and other malformed gene names | 1 known |
| #564, #467 | MHC-II allele spelling (`DPA` → `DPA1`, `A*24:01`) | see §7 |
| - | murine MHC-II: `I-Ab`/`H2-IAb`/`H2-Ab1` → one IMGT spelling | 3,396 |
| #368 | `antigen.gene` / `antigen.species` in the uncovered 48,935-row tail | ~90 reported |
| #347 | DOI/GitHub `reference.id` → PMID | to measure |
| #561 | identical alpha and beta CDR3 | to measure |

**Closes when:** every rule fires exactly its declared count and nothing else moves.

### Phase 9e - epitope catalogue validation with `mhcmatch`

`mhcmatch` (`~/vcs/code/mhcmatch`, ours, on PyPI) already models the table phase 9d
ships: which peptide a given MHC presents. It goes further than a name lookup: IPD-IMGT/HLA says
whether an allele exists; `mhcmatch` says whether that allele could present that peptide.

Four checks, in increasing strength, each a column on `restriction` and an advisory QC finding:

1. `pseudoseq.normalize_allele` on every `mhc.a` / `mhc.b`. It maps a call to a
   pseudosequence-FASTA key, so a call that does not normalise has no 34-mer binding groove and
   nothing downstream can reason about it. Stricter than the IPD-IMGT/HLA prefix check, which only
   asks whether the name is in the registry.
2. `pseudoseq.load_pseudo("I")` / `load_pseudo("II")` give the allele's class. Cross-check it
   against the recorded `mhc.class`: a class-I allele on an `MHCII` record is a defect the current
   build cannot see.
3. `store.infer_class(peptide)` infers the class from epitope length (MHC-I at ≤ 11). Disagreement
   with `mhc.class` catches a 15-mer filed as MHCI and a 9-mer as MHCII, independent of the allele,
   so it cross-checks check 2 rather than repeating it.
4. The restriction itself. `mhcmatch.predict` / `ligand.presented_span` score presentation, and
   `store.Restriction` gives `p_present`, `rank` and `band`. An epitope the recorded allele cannot
   plausibly present is either a wrong `mhc.a` or a wrong `antigen.epitope`, and neither is visible
   to any nomenclature rule. Advisory and never a rewrite: a presentation model is evidence about
   a pair, not authority over a publication, and its false-positive rate has to be stated with any
   threshold.

Checks 1–3 are deterministic string and length work and belong in the build. Check 4 needs a model
and its reference data (fetched from `isalgo/pmhc_data` on first use), so it runs as its own CI job
over the 2,381 `(epitope, MHC)` pairs and publishes a report, not inside `vdjdb build`, whose
offline determinism (hard rule 9) must not depend on a download.

What these checks would already have caught, from phase 9d's own findings: `HLA-A*08:01` on 74
records (no HLA-A\*08 locus exists, so no pseudosequence), the four null/nonexistent `HLA-A*24:*`
calls, and `HLA-B*12` (a serological antigen with no molecular groove). The IPD-IMGT/HLA prefix check
found those; `mhcmatch` would also have ranked the 74 `HPVTKYIM` records against `HLA-A*08:01` and
reported the peptide as not presented, which points at the fix rather than only the error.

Adds a dependency on `mhcmatch` for the validation extra only, not for the assembly build.

### Phase 10 - `feature/motifs-tcrnet`

1. `motifs/background.py` - `seqtree.control.load_control` against `isalgo/airr_control`, the
   `*.aa.vdjtools.tsv.gz` builds (§8.7). Never let `tcrnet()` resolve its own background: its
   `evalue.background(locus, species)` call takes no `size` and indexes the full table.
2. `motifs/tcrnet.py` - call `tcrnet()` for `n_control` only, then recompute the legacy statistic with
   its pseudocount (§8.1), which is bit-faithful to `DegreeStatisticsAnnotator.computePValue`:
   `p = (n_control + 1) / (M + 1)`, `p_legacy = binom.sf(d_s - 1, N, p) / (1 - (1 - p)**N)`.
3. File the `E = 0` bug upstream against `vdjtools` in the same pass.
4. `motifs/pwm.py` - the three-level Laplace cascade (`vj_len` → `len` → `uniform`, §8.5) so nothing
   is dropped and `freq` sums to 1. Assert the sign: at the 253 affected clusters the new `Σ I`
   must be lower than the shipped value. Reproducing the shipped numbers means reproducing the defect.
5. Fix `MotifsScoresAssembler`'s `(epitope, species, gene)` indexing, which never keys on `cdr3`
   (§8.5) and makes both `*_scored.txt` tables wrong.

**Closes when:** the deviation report accounts for every difference against the shipped files by a
named cause, and the 31 logo-less cids and the 1.00 % deleted letter mass are both gone.

### Phase 11 - `feature/motifs-tcremp`

1. `motifs/tcremp.py` - chunked embedding is mandatory (§8.8): fit `StandardScaler` + `PCA(50)` on
   a 25k seeded subsample, then embed in 20k chunks. Unchunked peaks at 9.77 GB and OOMs a 16 GB
   runner; chunked peaks at 3.26 GB. Chunks run sequentially, one internally-threaded `embed()` each.
2. Per-epitope DBSCAN with a chain-global eps (§8.4): the eps comes from the pooled k-distance
   curve, so it is never re-estimated on a per-epitope n of 30–300.
3. `eps = coef × mean(1st-NN distance)`. Kneedle is a debug cross-check only - it returns knee 1 of
   112,983 at production scale (§8.3). Do not call the method Kneedle-based.
4. Fit `coef` per chain against the §11.1 independent-study objective. Never against TCRvdb.
5. Re-fit the scaler and PCA every build from the seeded subsample; derive cluster ids from cluster
   content so they are stable without storage. Not stored, not cached (hard rule 9). The earlier plan
   used a 22 MB cache keyed on `(cdr3, v, j)`, partly for that id stability (today's igraph component
   numbers are not stable, so every bookmarked motif URL breaks each release).
6. Legacy projection: one legacy cid per `(cluster, stratum)`, `cid = H.B.<epitope>.<n>L<len>`, so
   both files satisfy `vdjdb-web`'s `strict = true` path (§8.6). `--legacy-shared-cid` reproduces
   today's shape.
7. Validation, once, at the end: `$VDJDB_TCRVDB`, aggregate metrics only (§11.2).

**Closes when** both of the following hold, and not before:

1. A floor and a target, scored on one cohort under the vendored `metrics_lib` (§8.9) against
   the REDCEA production clustering re-scored in the same harness:
   - Floor, reproduce. Recall and precision no worse than legacy's. Phase 11 closes here.
     Falling short means the rewrite lost something, and the cause has to be named before it merges.
   - Target, beat. Higher recall at equal or better precision. A trade of more coverage for
     less purity is not the target met; it is a different operating point, and is reported
     as one. This is what the motif optimisation knobs are for, and it is allowed to land after
     phase 11.

   Reference, from the benchmark's own `results/metrics_full.tsv`:

   | chain | method | purity | retention | precision | recall |
   |---|---|---:|---:|---:|---:|
   | TRA | tcremp-redcea | **0.9141** | **0.4481** | **0.9045** | **0.9141** |
   | TRA | tcrnet | 0.8424 | 0.1675 | 0.8136 | 0.8424 |
   | TRB | tcremp-redcea | 0.9509 | **0.6273** | 0.9559 | 0.9509 |
   | TRB | tcrnet | **0.9843** | 0.2602 | **0.9777** | **0.9843** |

   ⚠ Per-epitope clustering has purity 1.000 by construction, so it clears the purity bar with
   no information about quality. Compare on the pooled clustering, where purity is earned, and
   report the per-epitope result beside it, labelled (§8.4).
   ⚠ Not the 0.941/0.569/0.709 write-up: that is the same REDCEA table under a *coverage-aware* F1,
   a different statistic, so quoting it beside these would compare two things.

2. The result is shown to be stable, to the control and to the parameters, with the spread
   reported alongside the point estimate:
   - control choice - re-score at 5–10 background seeds, at {250k, 1M, full}, and against an
     in-silico background (§30.2); report the spread of retention, purity and the called set's
     Jaccard against the reference draw.
   - parameters - `coef`, `min_samples`, PCA components, `n_prototypes` for TCREMP; scope,
     `q`, `MIN_SAMPLE`, `MIN_CLUSTER` for TCRNET (§30.3, §30.4).
   Tuning stays on the §11.1 independent-study objective and never on TCRvdb (§11.2).

### Phase 12 - `feature/summary`

1. Split `summary/vdjdb_summary.Rmd` at the existing `!summary_embed_end!` marker (line 567) into a
   release dashboard and a paper-figures document. That drops `maps`, `scatterpie` and the never-used
   `ggh4x` from the release path (§7).
2. ggplot2 4.x fixes, all of them current breakages against the installed 4.0.2: `guide=F` → `guide="none"`,
   `size=` → `linewidth=`, `..count..` → `after_stat(count)`, `as.tibble` → `as_tibble`, and the
   `g_legend` grob-name grep → `cowplot::get_legend`.
3. Kill the live NCBI eutils call: publication years become `summary/pubmed_years.tsv` (`.tsv`, not
   `.txt`, because `summary/*.txt` is gitignored), refreshed by a scheduled workflow that opens a PR.
4. Data-drive the callouts into `summary/annotations.tsv` carrying `panel, year, label, hjust, vjust`
   and no coordinates (#460).
5. #460 itself: `grep -o 'width="[0-9]*"'` over the shipped embed HTML returns nothing, so
   `MakeEmbedableHtml.py`'s width rewrite is dead code. Delete it; add
   `style="max-width:100%;height:auto"` to the `<img>` tag instead. Done, and the script went with
   it: a pandoc template plus a Lua filter emit the fragment directly, so the three markup guesses it
   made no longer exist.
6. Pin rasterisation: `dpi=96, fig.retina=2`, `dev.args=list(type="cairo")` - explicitly cairo, not
   ragg, which would change font rendering and break parity on day one.
7. Replace the cumulative-by-year `expand.grid` cartesian join (~7M rows) with a first-appearance-year
   `cumsum`. It is the only part of the render that could plausibly OOM at 16 GB.
8. `summary/palette.py` is the single palette source for both renderers, emitting `palette.json` that
   the Rmd reads with `jsonlite::fromJSON`.
9. `summary/check_summary.py` - structural (ordered `<h4>`s, 5 tables, 8 base64 PNGs, IHDR-decoded
   width×height, the three vdjdb-web contracts), style (ColorBrewer anchors within ΔE₇₆ < 3, plus an
   anti-assertion against viridis), perceptual (SSIM vs the previous release, fail < 0.55, warn < 0.80).
10. `summary/preview/index.html` pulls Semantic UI from a CDN and `fetch()`es the fragment, because
    opening it directly shows unstyled tables: the classes come from vdjdb-web's bundle.

**Closes when:** both dashboards render offline and all three check layers pass.

### Phase 13 - `feature/docs`

1. Sphinx + `pydata_sphinx_theme`, `conf.py` copied from `arda/docs/`.
   `html_baseurl = "https://docs.isalgo.dev/vdjdb-db/"` - Pages is already provisioned, no setup step.
2. `docs/_ext/vdjdb_schema.py` provides `.. vdjdb-schema::`, `.. vdjdb-vocabulary::` and
   `.. vdjdb-score-rules::` as directives importing the phase-1 registry at doc-build time. No
   generated `.rst` in the tree means the tables can never be stale.
3. Structure: `getting-started/`, `standards/`, `submission/`, `builds/`, `dashboard/`, `reference/`.
   `docs/outputs.md` becomes `standards/database-outputs.rst`.
4. The dashboard tab is an `<iframe>`, not inline HTML: the fragment's Semantic UI classes and
   plotly's CSS must not leak into the theme. Its inner document is the phase-12 preview harness.
5. The dashboard artifact downloads from the last successful `build.yml` with `continue-on-error` and
   a committed `placeholder.html`, so a 30-minute database build never blocks a typo fix in the docs.
6. README shrinks to a ~60-line front door.

**Closes when:** `sphinx-build -W --keep-going` is clean and Pages deploys.

### Phase 14 - `feature/release-tooling`

1. `io/manifest.py` decides what goes in each bundle - an explicit list, never `cp *.txt`, which is
   how seven unshipped side tables nearly shipped (§1).
2. The six-step release job: plan (derive tag, assert it does not exist) → prepare (rewrite
   `latest-version.txt` in the working tree via temp-file + `os.replace`) → build (all zips embed that
   same content) → verify (release comparison + dashboard checks) → publish (`environment: release`,
   required reviewer) → finalize (commit `latest-version.txt`, then `curl -fsI` line 1 and fail on
   anything but 200). Step 6 has no counterpart in `release.sh`.
3. Tag scheme `v<YYYY>.<MM>.<PATCH>`; tag creation restricted by ruleset to the release environment.
4. No Zenodo step in the workflow. The webhook archived the source tarball rather than the
   release assets, the `-ZENODO` twin tags are the symptom of someone re-triggering it, and the
   minted DOIs are wrong. `.zenodo.json` is removed and the hook is not reinstated: the author
   deposits to Zenodo by hand, so the workflow must not create, mutate or assume a deposition.
   Concretely: never emit a `-ZENODO` twin tag, and treat the Zenodo DOI as an input the author
   supplies after the fact, not something the release derives.
5. `release/changelog.py` - reference diff between releases (#432), cheap because the phase-12
   publication-year table is a committed input.
6. `verify-latest` scheduled job: line 1 returns 200 and its tag equals `releases/latest`.
7. Retire `.gitlab-ci.yml`, `.travis.yml`, `test.sh`, `release.sh`, `docker.sh`, `release_docker.sh`,
   `gitlab/` and both Dockerfiles. Move `src/*.groovy` to `attic/`: it is the only correct
   specification for the meta files and phase 1 parses it at test time.
   `docker_build.log` needs no retiring: it matches `/*.log` in `.gitignore` and was never
   committed. `.zenodo.json` was already removed in `53c2005`.

**Closes when:** a full release dry-run produces three zips, a manifest and no unattributed
differences.

### Phase 15 - `feature/aldan3-runner`

1. Register aldan3 in a runner group scoped to this repo alone.
2. Retarget with `runs-on: ${{ fromJSON(inputs.motifs-runner) }}` - callers pass `'"ubuntu-latest"'`
   or `'["self-hosted","linux","x64","aldan3"]'`. The naive `runs-on: ${{ inputs.runner }}` cannot
   express a multi-label self-hosted target.
3. Guard every self-hostable job with
   `github.event.pull_request.head.repo.full_name == github.repository`. `chunk-check.yml` stays
   `ubuntu-latest` always, because it is the job forks trigger.
4. The existing GitLab runner cannot be reused: `gitlab/runner.slurm` asks for
   `--time=00:10:00 --mem=4G`, and a GitHub Actions runner is a long-lived daemon, not an sbatch job.

**Closes when:** a full build completes on both runners with identical canonical digests.

### Phase 16 - `feature/identity`

Delivers the four derived id levels of §10.3, the lifecycle record of §10.4, the commands and
invariants of §10.5, the promiscuity columns of §10.6, and one definition of "a study".

1. Move `clonotype_id` from `pl.Expr.hash` to `vdjdb.identity.ids` sha256, computed on the distinct
   key set and joined back (rule 4), formatted `CT` + 16 hex. The column is in `chains` and in no
   legacy file, so no shipped legacy byte moves and the release comparison stays PASS on all five
   members. `rules/expected_diffs.toml` gets no entry: it declares differences against the reference
   release, and the reference release has no `chains` table to differ from.
2. Add `clone_id` to `chains`, `pmhc_id` and `epitope_id` to `records`, each declared once in
   `schema/fields.py` so every projection and the documentation tables follow.
3. `src/vdjdb/identity/levels.py`: one function per level, each taking a frame and returning it with
   the id column added. Four declared key tuples and one shared hash-and-join helper. No abstraction
   over the four, because there are four and there will not be more.
4. `src/vdjdb/identity/lifecycle.py`: read a previous lifecycle file, diff it against the current id
   sets, write the new one. The diff is pure and does no I/O, so it is testable on literals.
5. `vdjdb identity build | resolve | diff | check` in `cli.py`.
6. `tests/unit/test_identity_levels.py`: the seven invariants, including the permuted-chunk-order
   test on a three-chunk fixture, plus one test per level that a key change produces a new id and a
   non-key change does not.
7. `identity check` runs in `build.yml` after the assembly step. Not in `chunk-check.yml`: that job
   QCs the changed chunks and never builds the tables, so there is nothing for it to check. The
   order-independence test (invariant 3) runs on a fixture wherever `pytest` runs, which costs
   milliseconds.
8. The release job writes `identity-lifecycle.tsv` and `record-registry.tsv` as release assets and
   commits `identity/retired.tsv`. `vdjdb build --registry <path>` accepts the previous release's
   asset; fetching it is a download, not a cache (hard rule 9's own exception).
9. `restriction` gains the six promiscuity columns of §10.6 from the `mhcmatch` call phase 9e
   already makes, with `mhcmatch.version` recorded per row.
10. **One definition of a study.** `assemble/evidence.support_counts` counts distinct `reference.id`
    on `records`, where the column holds one value per row. The dashboard recomputed it in R from
    `vdjdb.slim.txt`, where references are comma-joined, keeping field 1 only. That reported **527 of
    the 636 distinct non-blank references** in the slim table, missing 109. Rendering the same build
    both ways: human TRA 348 -> 401, human TRB 431 -> 483, mouse TRA 71 -> 84, mouse TRB 88 -> 151,
    macaque TRB unchanged at 2. Five sites in the Rmd carried it, feeding four tables. Fix all
    five to split the field and count every reference. The regression guard is a lint in
    `tests/unit/test_summary.py` asserting the truncating form is absent and the corrected one appears
    five times: `summary/fingerprint.json` records table *headers* rather than cell values, and
    `summary/panels.py` ports the figures rather than the tables, so neither can hold this number
    today.
11. `docs/standards/identity.md`: the five levels, the two mechanisms and why each level has the one
    it has, the lifecycle fields, the seven invariants, how to resolve an id you hold, and the
    statement that `TCR_hash` is preserved from the legacy recipe and never redefined.

**Closes when:** `identity check` passes every invariant; the permuted-order test passes; a rebuild
after adding and then removing a chunk leaves every other id unchanged; and the four dashboard
`Studies` columns report every reference on the row.

**State, 2026-09-27.** Steps 1 to 8, 10 and 11 landed. Measured on the rebuild: 187,935 clonotypes,
82,266 clones over 186,588 paired chain rows, 2,364 pMHCs, 2,118 epitopes, zero collisions at 16 hex
digits, and `identity check` reporting every invariant holding. The release comparison stays PASS on
all five legacy members, so nothing shipped moved. Removing one epitope from a copy of the tables
retires exactly its `pmhc_id` and `epitope_id` with `last_release` frozen, and leaves the other
274,681 ids untouched.

Two carried items. **Step 9, the promiscuity columns, waits on phase 9e**, which is the branch that
introduces the `mhcmatch` call; adding the call twice would put two model versions in one build.
**`record_id` is still not stable across releases**: the record registry is written at release time by
step 8 and nothing has shipped one yet, so until the first release under this scheme a rebuild
reconciles against an empty registry and allocates from 1. The four derived levels do not have that
dependency and are stable from this commit.

### Phase 17 - `feature/corpus`

A reference corpus as a reusable artifact, reproducing what `vdjdb.com/refsearch/` serves and
answering questions the endpoint cannot.

**The contract**, read from `vdjdb-web`'s client in `app/frontend/src/app/pages/refsearch/`. POST to
`/refsearch/` with `{cdr3, antigen.epitope, extra_parameters, species_to_search}`, every value a
space-joined string, `extra_parameters` drawn from `search_by_antigen` and `filter_stop_words`, and
`species_to_search` defaulting to `HomoSapiens MusMusculus MacacaMulatta`. The response is a JSON
array of `{pmid, tf_idf}` and the client keeps the first ten. A second endpoint, POST
`/refsearch/articles` with `{pmid}`, returns `{title, abstract, authors_list, journal,
publication_year}`. The service answering these is not in the `antigenomics` organisation; only the
client is. This phase builds the corpus and the scorer. A service reading them is separate.

**Documents are references, not records.** 662 distinct `reference.id` values in the current build:
610 PubMed ids, 51 others (bioRxiv and medRxiv DOIs, an arXiv preprint, PDB entries, a thesis, and
GitHub issues), and one blank. A document is a reference; the records citing it are its content.

**Three token families in one vocabulary, each prefixed so a consumer can filter.**

| Family | Prefix | Token | Distinct, current build |
|---|---|---|---|
| text | `w:` | a word of the title or abstract, lowercased | to measure; 610 documents have text |
| receptor | `v:` | V gene without allele | 228 |
| receptor | `j:` | J gene without allele | 80 |
| receptor | `k:` | a CDR3 3-mer | 7,713 over 180,048 distinct CDR3s |
| receptor | `kv:` | a CDR3 3-mer scoped to the V gene carrying it, `kv:CAS@TRBV9` | to measure |
| antigen | `e:` | epitope | 2,118 |
| antigen | `a:` | antigen species and antigen gene | from the epitope catalogue |
| antigen | `m:` | MHC allele at two fields | from `restriction` |

`k:` and `kv:` are both in the vocabulary because the question the corpus exists to answer needs
both. "Is the CAS motif specific to HIV, or to its TRBV?" is a comparison between the lift of `k:CAS`
on HIV documents and its lift on HIV documents already carrying that V gene. One token cannot express
that and neither can a single search ranking.

**Weighting.** Sublinear term frequency `1 + log tf`, because a paper reporting 10,000 receptors
would otherwise dominate every receptor token; smoothed inverse document frequency
`log((N + 1) / (df + 1)) + 1`; L2 normalisation per document, so a long abstract and a short one are
comparable. These are scikit-learn's conventions, which makes the implementation checkable against
`TfidfVectorizer` on a fixture rather than only against itself.

**The artifact.** Three parquet files in the release, plus the TSVs beside them:

```
corpus/documents.parquet   document_id, reference.id, kind, year, n_terms
corpus/terms.parquet       term_id, term, family, df, idf
corpus/postings.parquet    document_id, term_id, tf, weight
```

Long postings rather than a sparse-matrix format, because the consumer is polars or duckdb and the
query is a join. `term_id` is assigned by sorted term and `document_id` by sorted reference, so both
are total orders and the files are reproducible without a hash. Postings are sorted by
`(term_id, document_id)`, which makes a term lookup one contiguous slice.

**What is fetched and what is committed.** Receptor and antigen tokens come from `chunks/` and need
no network. Only the text family does, and it follows hard rule 5's pattern: the running text is an
input to the build and never an output of it. `vdjdb refs --abstracts` fetches titles and abstracts,
tokenises them, and writes `corpus/text_terms.tsv` as `reference.id, term, tf` - term counts, no
running text. That file is a committed, reviewed input refreshed by its own pull request, exactly as
`summary/reference_years.tsv` already is, so the build stays offline and deterministic. Estimated at
610 documents and roughly 120 distinct terms each, about 73,000 rows and 2 MB. The tokeniser is in
the repository, so the transform is auditable even though the text is not kept.

**Two functions, not a framework.**

* `score(query_terms) -> (document_id, score)`, the sum of matched weights, which is the `tf_idf`
  the endpoint returns.
* `lift(term, given) -> ratio`, the frequency of `term` among documents carrying every token in
  `given`, over its frequency across all documents. One `group_by` on the postings table. This is
  what answers the specificity question, and it is the reason the corpus is an artifact rather than
  an index inside a service.

**`vdjdb corpus build | query | lift` in `cli.py`**, with `--stop-words` and `--search-by-antigen`
mapping onto the two `extra_parameters` the client sends, so the command and the endpoint take the
same arguments.

**Tests.** `tf_idf` against `TfidfVectorizer` on a five-document fixture to four decimal places; the
three files reproducible across two runs by digest; `term_id` and `document_id` unchanged when the
chunk list is permuted; `lift` equal to 1.0 for a term independent of the condition by construction;
and one end-to-end query whose top document is the paper the epitope was first reported in.

**Follow-ups, named rather than built.** The artifact is method-agnostic, so an alternative weighting
is an additional file and not a rewrite.

* **PageRank over the term co-occurrence graph** ranks terms by centrality instead of by rarity.
  Inverse document frequency handles receptor tokens badly: a 3-mer in every document has idf near
  zero and drops out, even when it is the hub connecting two antigen families. Measure it against
  tf-idf on the same queries before preferring it.
* **An embedding index**, sentence embeddings for the text family and `mir.embedding.TCREmp` for the
  receptor family, answers "which references are about something similar" rather than "which share a
  token". It is a dense vector per document in a fourth file.

Both are gated on tf-idf reproducing the existing endpoint first, because that is the only baseline
in existence.

**Closes when:** `vdjdb corpus build` produces the three files reproducibly; `score` on the
endpoint's own request shape returns a ranking; `lift` answers the CAS-against-TRBV question with a
number and an n; and `docs/standards/corpus.md` documents the vocabulary, the weighting, the two
functions and the file schemas.

### Bootstrap order - done, 2026-09-27

Completed in this order (`ROADMAP_local.md` §44): land the workflows on `master` → create `dev` → let one full `build.yml` run green on `dev` so the
check names exist → then apply the `dev`, `master` and tag rulesets. A required check that has
never run blocks every PR forever. Linear history is deliberately **not** applied: `master` carries 168 reachable merge commits and
the gitflow puts `--no-ff` feature merges on `dev`, so the rule would reject every future update and
train the bypass habit. The four rulesets that did go on are recorded in `ROADMAP_local.md` §44.2.
