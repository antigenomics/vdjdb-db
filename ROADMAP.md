# VDJdb build roadmap

Migration of database proofreading, assembly and release from the Docker + GitLab + manual-`scp`
pipeline to a `uv`/polars Python package driven by GitHub Actions.

Deliverables with acceptance criteria, not dates. `STATUS.md` says what is in flight right now.

---

## 1. Why

The current build is a Docker image (ubuntu:18.04, Python 3.10 from source, OpenJDK 8 + Groovy 3.0.9,
R 4.x + ~25 CRAN packages, TeXLive, the VDJtools 1.2.1 jar) driven by `release.sh`, triggered from
GitLab CI by `ssh`-ing into aldan3 and `sbatch`-ing a runner. `release.sh` ends by printing an `scp`
command; a human then runs `gh release`.

That has produced concrete, measurable defects:

| Defect | Evidence |
|---|---|
| The shipped release is not what `release.sh` produces | `release.sh` does `cp *.txt`, which would sweep in `vdjdb_full_filtered.txt`, three `*_broken.txt` and three `*_scored.txt`. The 2026-06-03 zip has exactly 10 files and none of those. |
| Bundle shape drifts every release | 2024: flat + `LICENSE.txt` · 2025-02-21: `vdjdb-<date>/` + 4 extra files · 2025-09-25: `database/` · 2026-06-03: `vdjdb-<date>/` + `LICENSE` |
| Column metadata no longer describes the data | `vdjdb.meta.txt` is a static git file last written by the retired Groovy code: no `TCR_hash` row, `vdjdb.score` in the wrong position, `v.end`/`j.start` wrong place and order in slim |
| `web.cdr3fix.unmp` is wrong on **7,998 of 284,546 rows** | `jStart == -1` is truthy in Python; the Groovy tested `jStart > -1` |
| Column definitions duplicated in nine places, already drifted | `ChunkQC.py`, `DefaultDBGenerator.py`, `SlimDBGenerator.py`, `ScoreFactory.py`, two static meta files, `BuildDatabase.groovy`, the README, `template.xls` |
| Build needs 64 GB RAM for 42 MB of data | seven `master_table.T.apply(...)` calls plus `iterrows()` over ~200k rows |
| Stage-II CDR3 fixing silently absent | `AlignBestSegments.py` is only called from `BuildDatabase.groovy`; the Python path never invokes it |
| `latest-version.txt` names the *previous* release | line 1 is `2026-05-16`; the published latest is `2026-06-03-ZENODO` |
| The dashboard is non-deterministic | a live NCBI eutils call at render time, plus a hardcoded fallback year table frozen since 2021 |
| Five chunk columns silently discarded | `submitter`, `chunk.id`, `comment`, `meta.subset.frequency`, `method.pairing` — four of them documented in the README |

## 2. What ships

Three builds, every release, until legacy is retired:

| Build | Asset | Role |
|---|---|---|
| **New VDJdb** | `vdjdb-<version>.zip` | primary product: polished legacy — valid JSON everywhere, generated metadata, stable column order, the derived `cdr3nt`/`d.*`/Pgen block |
| **Legacy** | `vdjdb-legacy-<version>.zip` | byte-layout-compatible with today's zip; what `vdjdb-web` and standalone clients consume |
| **AIRR** | `vdjdb-airr-<version>.zip` | AIRR Rearrangement + Receptor/Reactivity |

Plus `manifest.json` (`{role, file, sha256, bytes, members[]}` per build) and `SHA256SUMS`.

Legacy is a **derived export** — `vdjdb release` → `make legacy` → zip — reading the new build
directory, never `chunks/`. It cannot drift from the primary build because it is a projection of it.

Conversion helpers run both ways: new → AIRR and legacy → AIRR.

## 3. Gates

### 3.1 Cross-repo: `vdjmatch` must be patched and released first

`vdjmatch/src/vdjmatch/db/vdjdb.py` no longer reads `latest-version.txt`. It calls the GitHub Releases
API and takes the **first `.zip` asset**:

```python
def _zip_asset(rel: dict) -> str:
    for a in rel.get("assets", []):
        if a["name"].endswith(".zip"):
            return a["browser_download_url"]
```

The moment a release carries more than one zip, vdjmatch picks an arbitrary one.

1. Patch `_zip_asset` to select by role: `manifest.json` `role == "legacy"` → `vdjdb-legacy-*.zip`
   → `vdjdb-<tag>.zip` → first `.zip`. Bump `_HF_TAG` in the same pass.
2. Cut a `vdjmatch` release.
3. **Only then** cut the first multi-zip VDJdb release.

Corollary: the legacy zip's member basenames must stay `vdjdb.txt`, `vdjdb.slim.txt`,
`vdjdb_full.txt` — that is what vdjmatch's member lookup keys on.

### 3.2 Phase −1: three unblocking fixes, none of them part of the rewrite

| Fix | Blocks |
|---|---|
| The `vdjmatch` patch above | "three zips per release" |
| Correct `latest-version.txt` line 1; delete the stale `database/` copy | every client that still reads it |
| Replace the seven `.T.apply` calls and the `ScoreFactory` `iterrows` | "runs on 16 GB GitHub-hosted runners" |

### 3.3 Licence

`isalgo/airr_control` is **CC-BY-NC-ND-4.0**; VDJdb ships **AGPL-3.0-only**. Backgrounds are streamed
at build time and only **derived statistics** (`count.bg`, `total.bg`) ship. Never the background
itself, never a subsampled copy.

## 4. Phases

`master` → `dev` → `feature/*` → `dev` → `master`. Every phase is independently mergeable and `dev`
stays green. Commits that resolve a tracker issue carry `Closes #N`.

| # | Branch | Delivers | Closes | Acceptance |
|---|---|---|---|---|
| 0 | `feature/dev-baseline` | `CLAUDE.md`, `ROADMAP.md`, `STATUS.md`, `pyproject.toml` + `uv.lock`, package skeleton, `chunk-check.yml`, `branch-policy.yml` | #476 | CI green; `release.sh` still works untouched |
| 1 | `feature/schema` | the field registry; `render_meta` | — | reproduces the Groovy `METADATA_LINES` / `SLIM_METADATA_LINES` byte-for-byte; `header == meta names` for all three tables |
| 2 | **`feature/golden-harness`** | `vdjdb diff` + `expected_diffs.toml` | — | **zero diffs against the current pandas build.** Nothing downstream starts without this |
| 3 | `feature/io-qc` | polars reader, vectorised QC, `--strict` exit-1, chunk header normalisation, `.tsv` rename | #497 | QC report matches the pandas report row-for-row; harness still zero |
| 4 | `feature/pipeline-core` | **the definitive tables** (`records`, `chains`) + harmonize + score + pairing; the legacy export as a projection of them; deletes `py_src/` | #424, #399 | every ledger difference is a declared rule firing its measured count; peak RSS < 8 GB |
| 5 | `feature/arda-cdr3fix` | `arda.cdr3fix` replaces `Cdr3Fixer.py`; retires `res/segments*.txt` | — | new ledger rule, row count measured then frozen |
| 6 | `feature/new-format` | ships the definitive tables as parquet + TSV, adds `evidence`, `vdjdb.schema.json` | — | `make legacy` from the shipped tables still passes the harness |
| 7 | `feature/airr` | `emit/airr.py` (Rearrangement + Reactivity), `convert/coords.py`, `vdjdb convert` | — | `airr.validate_rearrangement` passes on the full table; the legacy path produces nothing the tables path does not |
| 8 | `feature/junction-nt`, `feature/segment-guess`, `feature/dgene` | one branch each | #461, #462, #463 | generated `cdr3nt` back-translates to `cdr3` |
| 9 | `feature/harmonize-rules` | nomenclature rule tables | #327, #389, #347, #368, #564, #467, #561 | each rule gets a ledger entry with a measured row count |
| 10 | `feature/motifs-tcrnet` | TCRNET on `vdjtools`, streaming backgrounds | — | deviation report accepted |
| 11 | `feature/motifs-tcremp` | TCREMP + per-epitope DBSCAN; new motif schema; legacy projections | — | beats the **shipped** `cluster_members_tcremp.txt` re-scored in our harness, per §8.4 |
| 12 | `feature/summary` | Rmd split, ggplot2 4.x fixes, committed publication-year table, data-driven callouts, interactive dashboard | #460 | renders offline; perceptual + structural checks pass |
| 13 | `feature/docs` | Sphinx site, generated schema tables, dashboard tab, Pages | — | zero-warning build, deploys |
| 14 | `feature/release-tooling` | manifest, three zips, checksums, `latest-version.txt`, tag scheme, Zenodo, changelog; retires the legacy CI | #432 | full release dry-run with a clean ledger |
| 15 | `feature/aldan3-runner` | self-hosted runner + `build.yml` retargeting | — | identical canonical digests on both runners |

**Phase 2 is load-bearing.** The harness must show zero diffs against the *current* build before any
behaviour changes, or there is no instrument to attribute later differences with.

Phases 5, 8, 9, 10 and 11 each introduce exactly one source of deviation, so every difference in the
output has a single attributable cause.

**Every phase has a step-by-step subplan in §12.** A phase is not startable until its subplan names
the files it creates, the facts it needs (already measured, in §7/§8), and the check that closes it.

## 4a. What the issue tracker actually is

Measured 2026-09-25 with `gh`: **440 issues, 130 open.** Grouped by label, the open ones are

| Category | Open | What they are |
|---|---|---|
| **data intake** | **103 (79 %)** | pending papers (79), preprints (9), paper-pending (3), meta-papers (4), 10x/Immudex sets (5), associations (9), other databases (1), correspondence (2) |
| curation quality | 22 | formatting & proofreading (18), typos, structural, validation |
| build infrastructure | 13 | the build, the summary, maintenance |

Some issues carry more than one label, so the columns overlap slightly.

**Four out of five open issues are a submission queue, not a defect list.** This migration closes
issues from the bottom two rows only -- thirteen of them -- and nothing it does shortens the first
row. That matters for three decisions already taken:

* `chunk-check.yml`'s **three-minute budget is the one that matters**, because it is the job the
  submission queue runs through. The full build's 185 s is paid on `dev` and nightly, by nobody
  waiting.
* the **curation skills** (`/vdjdb-extract`, `-format`, `-proofread`, `-publish`) are the tooling with
  the largest backlog pointed at it, and they read `proofreading/`, which until phase 9 no build code
  touched.
* a phase that closes an issue number is not thereby reducing the tracker. Progress on the queue is
  curation throughput, and it is measured separately.

## 5. The difference ledger

`vdjdb diff <reference-zip> <candidate-dir>` compares in three passes: file set → two digests per file
(raw and canonical) → row-level classification keyed on
`gene|cdr3|v.segm|j.segm|species|mhc.a|mhc.b|antigen.epitope|reference.id`.

Every changed **cell** must be attributed to a declared rule in `rules/expected_diffs.toml`:

```toml
[[rule]] id="web-unmp-jstart-minus1" file="vdjdb.txt" column="web.cdr3fix.unmp"
         from="no" to="yes" predicate="cdr3fix.jStart == -1" rows=7998
[[rule]] id="meta-tcr-hash-row"   file="vdjdb.meta.txt" added=["TCR_hash"]
[[rule]] id="meta-score-position" file="vdjdb.meta.txt" moved=["vdjdb.score"]
[[rule]] id="meta-web-data-type"  file="vdjdb.meta.txt" rows=4
[[rule]] id="slim-meta-tcr-hash"  file="vdjdb.slim.meta.txt" added=["TCR_hash"]
[[rule]] id="slim-meta-geometry"  file="vdjdb.slim.meta.txt" moved=["v.end","j.start"]
```

The metadata fixes are safe because the metadata is currently *wrong* rather than merely old:
`TCR_hash` has been a `vdjdb.txt` column for years with no metadata row; `vdjdb.score` and
`TCR_hash` are listed after `cdr3fix` but stored before `method`; and the four `web.*` rows carry one
value too many, landing `0` in `data.type`. All four `web.*` rows are `visible = 0`, so nothing
user-facing moves.

Any unmatched difference **fails**. A rule that fires a different number of times than declared **also
fails** — that is what turns "we think it is the same" into a gate, and why every rule carries a
measured row count rather than a description.

### "Byte-for-byte" has to be redefined

`runBuidDatabase.py:49` iterates `os.listdir("../chunks")` — readdir order, host-dependent. The
released `vdjdb_full.txt` opens with `PMID:28629751` while the alphabetically-first chunk is
`10xgenomics-2019-07-09.txt`. The 2026-06-03 release therefore cannot be reproduced byte-for-byte by
anyone, including the current pipeline.

So **canonical equality is the gate**; raw equality is informational. The new reader uses
`sorted(glob("*.tsv"))`, `vdjdb release` writes `chunk-order.txt` into the bundle, and
`--chunk-order <file>` replays a recorded order — making raw equality achievable going forward.

## 6. Legacy deprecation path

Legacy stops shipping only when all of the following hold:

1. `manifest.json` role selection is released in every known client (`vdjmatch`, `vdjdb-web`).
2. No consumer reads `latest-version.txt` for anything but historical URLs.
3. `vdjdb-web` reads the new format's generated metadata rather than the legacy `vdjdb.meta.txt`.
4. Two consecutive releases have shipped both, with the new format downloaded at a comparable rate.

Until then, `latest-version.txt` line 1 points at the **legacy** zip, because clients in the wild
download it verbatim and expect that layout.

## 7. Measured facts

Established 2026-09-25 against the 2026-06-03 release and the current `chunks/`. These are inputs, not
drafts — do not re-derive them.

| Quantity | Value | How |
|---|---|---|
| Chunk rows, raw | 203,308 | 230 files, `chunks/*.txt` |
| Chunk rows after per-chunk dedup on `SIGNATURE_COLS` | **192,753** — exactly the released `vdjdb_full.txt` row count | polars |
| Rows matching field-for-field across two chunks | 19 pairs — **independent reports, not duplicates**: a chunk is one paper | deduplication is within a chunk; these 19 are evidence (§11.1), and global dedup would delete them |
| polars read + dedup of all 230 chunks | **0.4 s** | vs a pipeline documented as needing 64 GB |
| Chunks passing `ChunkQC` | **230 / 230**, zero errors | fail-fast needs no quarantine list |
| Non-empty CDR3 cells | 305,031, **zero** with characters outside the 20 AAs | TCREMP pre-filter is a guard, not a live problem |
| `arda.cdr3fix` vs shipped `cdr3fix`, 20,000-row sample (seed 42) | `cdr3` 97.81 % · `vFixType` 97.95 % · `vEnd` 96.53 % · `good` 95.68 % · `jFixType` 95.44 % · `jStart` 91.02 % | arda 2.27.0 |
| …of the 1,797 `jStart` disagreements | 556 VDJdb-unmapped → arda-mapped · 1,241 both mapped, **arda smaller in every case** (mode −2, range −1…−8) · **0** coverage regressions | NW with free end gaps vs k-mer longest-hit |
| `arda.cdr3fix.markup_batch` throughput | 27 µs/record → **5.2 s** for 191,440 distinct keys | not vectorised internally, but fast enough that it does not matter |
| `vdjtools.model.infer_nt` | **3.11 ms/record** → ~15 min for 284,546 rows | human TRB, warm, single-threaded |
| `vdjdb.txt` row/field shape | 284,546 rows, **all exactly 22 fields**, zero empty `cdr3fix` | keep and assert; raggedness is not a live problem |
| Quote characters in the release tables | present in **every** `method` / `meta` / `cdr3fix` cell — they are JSON. What is absent is a *quoted field*: no field begins with `"` | the plan's "zero `\"` characters" was wrong. `quote_style="never"` is still exact, and for the stronger reason: a default CSV writer would wrap every JSON cell and double its quotes |
| `vdjdb_full.txt` rebuilt by the current pandas pipeline vs the 2026-06-03 release | **119,169,153 bytes both**, 192,755 lines both, **canonical digest identical**, raw digest differs | the reproduction contract holds on the largest table; the raw difference is exactly the `os.listdir` chunk order |
| All five positional column orders, registry vs shipped release | `vdjdb.txt` 22 · `vdjdb.slim.txt` 17 · `vdjdb_full.txt` 35 · `cluster_members.txt` 19 · `motif_pwms.txt` 27, **order-for-order identical** | validates the phase-1 registry against the artifact consumers actually parse |
| `vdjdb_full.txt` `cdr3fix.alpha` encoding | 122,930 non-empty cells, **100 % Python dict repr**, 0 JSON | live bug |
| `web.cdr3fix.unmp` correctness | `(no,no)` 268,546 · `(yes,yes)` 8,002 · **`(no,yes)` 7,998** | `jStart` is only ever −1 (8,536) or > 0 (276,010); **0 never occurs** |
| Murine MHC-II spellings | `I-Ab` **3,274** · `H2-IAb` 113 · `H2-Ab1` 9 · `H2-IAg7` 333 vs `H2-Ag7` 3 · `H2-Aa` 25 vs `H-2Aa` 19 · `H2-Eb1` 7 vs `H-2Eb1` 7 | mostly in `mhc.b`; class I is clean |
| #327 TRAJ24 | bare `TRAJ24` 1,302 (978 carry `WGKLQF`) · `TRAJ24*01` 113 (**75 = 66 % carry `WGKLQF`**) · `TRAJ24*02` 34 · malformed `TRAJ24-1` 1 | reproduces the report on a larger set |
| #368 antigen gene/species | patch dict covers 245 epitopes = 154,373 of 203,308 rows, **zero swapped** | the ~90 reported records are in the uncovered tail of 48,935 rows |
| Dashboard R deps | 15 of 17 installed; `maps`/`scatterpie` used only past the embed cut (lines 812–835 vs marker at 567); `ggh4x` never used | splitting the Rmd drops three deps from the release path |
| Shipped dashboard PNG sizes | 1344×960, 2304×1920, 1152×1920, 1536×1536 → `dpi=96, fig.retina=2` | pin it or the visual fingerprint is noise |
| Production's 27-row `vdjdb.meta.txt` (`vdjdb-web/test/resources/database/`) | orders `… reference.id method meta cdr3fix vdjdb.score TCR_hash web.*` while the data is `… reference.id vdjdb.score TCR_hash method meta cdr3fix web.*` | the metadata mis-describes the data **in production too**, not only in the release zip |
| `width="1152"` occurrences in the shipped embed HTML | **0** — the rewrite in `MakeEmbedableHtml.py` is dead code | this is what #460 actually is |

## 8. Motif inference — findings that change the approach

Measured 2026-09-25 against the 2026-06-03 slim dump and the shipped motif files.

### 8.1 `vdjtools.overlap.tcrnet`'s p-value is unusable as shipped — use it as a neighbour counter

`overlap/tcrnet.py::_score` computes `E = (n_target / max(m_control,1)) * n_control` with **no
pseudocount**, then `p_enrichment = poisson.sf(degree - 1, E)`. When `n_control == 0`, `E == 0.0`, and
`scipy.stats.poisson.sf(k-1, 0.0) == 0.0` **exactly** for every k ≥ 1 — verified. So every clonotype
with at least one within-sample neighbour and no background neighbour receives `p = 0.0`, `q = 0.0`,
and sorts to the top of the result.

Measured on human TRB (111,407 unique CDR3s) against the bundled 250k control: `n_control == 0` for
**81.0 %** of queries; **44,010 of 111,407 rows (39.5 %) get `p_enrichment == 0.0` exactly**.

A bigger background does not fix it — mean `n_control` scales linearly with M, so `E` is M-invariant
in expectation and M only controls the zero-inflation (human TRA: 34.7 % zeros at M = 250k, 15.6 % at
M = 2,266,274).

**Therefore:** call `tcrnet()` for `n_control` only, and recompute the legacy statistic, which has the
pseudocount that removes the pathology:

```
p        = (n_control + 1) / (M + 1)
p_legacy = binom.sf(d_s - 1, N, p) / (1 - (1 - p)**N)
```

This is bit-faithful to the Groovy `DegreeStatisticsAnnotator.computePValue`. File the `E = 0` bug
upstream against `vdjtools`.

### 8.2 Two things in the plan were wrong

- **The legacy TCRNET was *not* V/VJ/VJL-grouped.** `CalcDegreeStats.groovy` defaults to
  `-g dummy`, and the legacy Rmd passes only `-o 1,0,1 -g2 dummy`. The new grouping-free `tcrnet()`
  is an **exact scope match**, not a deviation — one less thing to reconcile.
- **The 0.941 / 0.569 / 0.709 figure is not a knee-DBSCAN result.** It is the shipped
  Leiden/REDCEA production table re-scored, "nothing re-fit". Measured knee-DBSCAN on mirpy
  embeddings traces a strictly lower frontier (pooled human TRB ≥30: 0.910 / 0.368 / 0.524 at
  coef 0.75). **The target is the shipped file re-scored in our own harness — beat the file, not the
  write-up.** `coef = 0.75` was calibrated on a different embedding (standalone `tcremp`, ~3000 OLGA
  prototypes, Smith-Waterman) and does not transfer to mirpy's; it must be re-calibrated.

### 8.3 The Kneedle branch is dead at production scale

Pooled human TRB, 112,983 unique clonotypes, PCA(50): Kneedle returns **knee = 1 of 112,983**
(frac 0.000). The `floor_frac = 0.40` guard fires *always*, so the operating rule is
`eps = coef × mean(1st-NN distance)` and Kneedle is only a debug cross-check. Stop describing the
method as Kneedle-based.

The two implementations estimate **different quantities** and must not be conflated:
`knee_eps.py` sorts each column independently and reads column 1 — the sorted **1-NN** curve, scaled
by 0.75. `mir.bench.metrics.estimate_dbscan_eps` uses the true **4-NN** curve, scaled by 0.4, with no
degenerate guard. vdjdb-db owns one parameterised implementation, defaulting to the `knee_eps.py`
variant (the one the benchmark numbers came from), and reports both.

### 8.4 Per-epitope DBSCAN with a chain-global eps is the right topology

Measured, human TRB ≥30, same embedding, vendored `metrics_lib` semantics:

| topology | purity | retention | cov-F1 |
|---|---|---|---|
| pooled DBSCAN, coef 0.75 | 0.910 | 0.368 | 0.524 |
| pooled DBSCAN, coef 1.10 | 0.787 | 0.545 | 0.644 |
| **per-epitope, chain-global eps, coef 1.10** | **1.000** | **0.406** | **0.578** |
| per-epitope, coef 1.30 | 1.000 | 0.452 | 0.623 |
| TCRNET (shipped) reference | 0.950 | 0.230 | 0.370 |

Pooled *geometry*, per-epitope *scope*: the eps comes from the pooled k-distance curve, so it is never
re-estimated on a per-epitope n of 30–300 — exactly the degenerate regime. Purity is 1.000 by
construction because a cid cannot span epitopes, which is what the output format requires anyway.

Stated plainly: that means "no cross-epitope contamination is *possible*", not "none exists". So the
**pooled clustering runs in the debug build** as the honest measurement, and its cross-epitope
confusion matrix is a first-class QC artifact.

### 8.5 Two bugs in the shipped motif files

- **31 of 620 cids in `cluster_members.txt` have no rows in `motif_pwms.txt`** — the web shows those
  clusters in the tree with no logo. Cause: `filter(total.bg > 0)` after a left join against the
  frozen Zenodo background PWM.
- **1.00 % of letter mass is silently deleted, across 253 of 589 clusters.** The same filter drops any
  `(pos, aa)` absent from the background at that `(v, j, len)` — the rarest, most informative
  residues. Per-position `freq` then sums to < 1 (minimum observed 0.111) and
  `I = 1 + Σ freq·log freq / log 20` is computed on a truncated distribution, **inflating I**.
  `need.impute` is computed *after* the filter, so it is `FALSE` on all 13,456 rows — the imputation
  machinery is dead code.

Fix: a three-level Laplace-smoothed cascade (`vj_len` → `len` → `uniform`) so nothing is ever dropped,
`freq` sums to 1, and `need.impute` becomes a real provenance flag. **Assert the sign**: at the 253
affected clusters the new `Σ I` must be *lower*. A rewrite that reproduces the inflated values has
reproduced the bug.

Separately, `py_src/MotifsScoresAssembler.py:44-56` sets `cluster.member = 1` for **every record of an
(epitope, species, gene)** — it indexes on those three columns and never on `cdr3`. So
`vdjdb.slim.scored.txt` and `vdjdb.scored.txt` are wrong.

### 8.6 The variable-length problem is smaller than assumed

`vdjdb-web` **already splits every cid by `len`** before building a cluster
(`Motifs.scala:85-89`, `splitOn(cid)` then `splitOn(len)`), and `MotifCluster` carries `vsegm`/`jsegm`
as sequences. So the legacy schema is not the obstacle. But relying on the relaxed `strict = false`
path produces two display clusters sharing one `clusterId`.

Design: a cluster owns one member set and one or more **PWM strata** (one per CDR3 length with ≥3
members). The default legacy projection emits **one legacy cid per (cluster, stratum)**, `cid =
H.B.<epitope>.<n>L<len>`, so both files satisfy `strict = true` and no two display clusters share a
`clusterId`. `--legacy-shared-cid` reproduces today's shape if the web maintainer prefers it.

### 8.7 Background choice, with the reason

Use `{species}.{locus}.aa.vdjtools.tsv.gz`, not `*.ntvj` and not `vdjdb-web-control/*` — the inverse
of the argument in `vdjdb-web-control/README.md`. There the *query* is a user's nucleotide-level
repertoire, so the control must be nt-level. Here the query is a VDJdb epitope sample, deduplicated to
unique CDR3aa, and `seqtree.Index` is built over aa strings. An ntvj control inserts the same `cdr3aa`
once per nt variant, and that inflation is measured at **1.63× for VDJdb-matching clonotypes against
1.04× overall** — concentrated on exactly the germline-proximal public sequences under test. It would
inflate `E` by ~1.6× for those sequences and systematically **under**-call public motifs.

`seqtree.control.load_control` already streams from `isalgo/airr_control`, filters to the productive
20, reservoir-samples uniformly over unique clonotypes, and content-addresses its **download** — do
not reimplement it. What it stores is the fetched table, which is an input; no derived control is ever
written (hard rule 9), and the reservoir sample is deterministic in the seed anyway. But **never let `tcrnet()` resolve its own background**: it calls
`evalue.background(locus, species)` with no `size`, which indexes the entire table.

### 8.8 Memory is the only real constraint

TCREMP embedding, human TRB, measured: single-shot `embed` → `StandardScaler.fit_transform` → `PCA`
peaks at **9.77 GB** (12.41 GB with a full-matrix scaler fit). Fitting the scaler and PCA on a 25k
seeded subsample and then embedding in 20k chunks peaks at **3.26 GB**. **Chunking is mandatory, not
an optimisation** — unchunked it would OOM a 16 GB runner alongside the rest of the build. Chunks run
sequentially; one internally-threaded `embed()` per chunk, never a worker pool.

Throughput is a non-issue: 75,308 rows/s on 16 cores, 61,605 rows/s at `threads=4`; the full 113k-row
human TRB embedding takes **1.4 s**. Whole motif stage: ~30–60 min cold, ~10–15 min warm, peak 3.3 GB.
The six-hour limit is never in play.

**Nothing is stored between builds** (hard rule 9): the scaler and the PCA are re-fitted every run
from the seeded 25k subsample, which is deterministic in `config.SEED` and the input, so storing them
would buy minutes and risk shipping a fit that no longer matches the data.

The reason a stored artifact was considered is real and needs a different answer: raw igraph component
numbers are **not** stable across releases, so every bookmarked `vdjdb.com` motif URL breaks each
time. The fix is a content-derived cluster id — the same choice `clonotype_id` already makes (a seeded
hash of the clonotype key, never a counter) — so a cluster whose membership is unchanged keeps its id
without anything being remembered.

### 8.9 `mir.bench` helpers are not usable here

- `mir.bench.vdjdb.load_vdjdb` keeps only `junction_aa, v_call, j_call, locus, epitope, mhc_class` —
  it drops `species`, `mhc.a`, `mhc.b`, `antigen.gene`, `antigen.species`, `v.end`, `j.start`, seven
  of the 19 legacy columns. Read the slim table directly.
- `mir.bench.metrics.cluster_metrics` uses `recall = tp / n_true_clustered` (clustered records only);
  the vendored `metrics_lib.precision_recall_fscore` folds unclustered records into FN. **The two
  metric families are not comparable.** Validation uses the vendored `metrics_lib`, verbatim.

## 9. Decisions

Settled 2026-09-25.

- **The new format takes ownership of `evidence.*`.** Production `vdjdb-web` serves five
  `evidence.*` columns that nothing in this repo produces. The new format produces them, from the
  evidence table in §10.
- **The five "dropped" chunk columns are kept.** `submitter`, `chunk.id`, `comment`,
  `meta.subset.frequency` and `method.pairing` all matter for debugging a curation problem, so they
  are carried through to `records.parquet` rather than discarded.
- **The side outputs are produced but not zipped.** `vdjdb_full_filtered.txt`, the three
  `*_broken.txt` and the three `*_scored.txt` tables are written to `out/reports/` and uploaded as
  CI artifacts; whether any of them ever ships is deferred.
- **Everything the build produces is specified** in `docs/outputs.md`, which becomes
  `docs/standards/database-outputs.rst` when the Sphinx site lands.

## 10. Record identity and the evidence model

### 10.1 Record identity

A VDJdb record has had no identifier: it is a row in a chunk, named only by its content. A curator
fixing a CDR3 typo therefore deletes one record and creates another, indistinguishably from a real
deletion plus a real addition. External references, structure links and accumulated evidence all
break silently.

Identity is two-level, in `src/vdjdb/identity/`:

| | |
|---|---|
| `record_id` | opaque, assigned once, never reused, **stable across content changes**. `VDJDB` + 10 digits. This is what external references point at |
| `content_hash` | sha256 over the canonical content. Changes whenever anything changes. Detects change; never identifies |

Reconciliation, most specific first: **exact** natural key → **amendment** (exactly one differing
key field, unambiguous, same chunk and reference — ambiguity is refused, because a wrong link is
worse than a new id) → **allocation**. Registry entries the build no longer sees are **retired**,
not deleted. The registry is a committed TSV sorted by `record_id`, so a curation PR shows added,
amended and retired records as a reviewable diff.

Two earlier keys were wrong in opposite directions, both caught by tests that now pin the
boundary: a narrower one that stopped at `reference.id` collided on **20,769 of 192,753** records
(one paper reporting the same TCR in several donors is several records), and one without
`chunk.file` merged **19** pairs that are two papers' independent reports.

The natural key is `CHUNK_DEDUP_KEY` **plus `chunk.file`**. A chunk is one paper, so two matching
rows in two chunks are two independent reports and must keep separate ids. Ids are assigned
**before** CDR3 repair: two trimmed sequences that repair to the same full one are still two
observations, and assigning afterwards merged 215 pairs the publications reported separately.

Measured: 192,753 rows → 192,753 records; ids stable across rebuilds; a typo fix reported as
`VDJDB0000000101 cdr3.beta: CASSIRSSYEQYF -> CASSIRSSYEQYFF`.

### 10.2 The evidence model

Four tables, normalised so nothing is stored twice (full column lists in `docs/outputs.md`):

```
records.parquet    PK record_id             one row per curated record, with provenance
chains.parquet     PK (record_id, gene)     one row per TCR chain; carries clonotype_id
evidence.parquet   PK (record_id, evidence_id)   long format, one row per piece of evidence
vdjdb.parquet      the joined, pivoted view -- derived, never authored
```

Evidence is long rather than wide because a record carries any number of pieces of any number of
kinds; a wide table would be mostly null. `evidence_type` covers `motif_tcrnet`, `motif_tcremp`,
`structure_native`, `structure_model` and `independent_study`. Structure evidence is keyed on the
legacy `TCR_hash` until the structure store is re-keyed on `record_id`.

`chains` exists so the schema stays non-redundant: folding chains into records forces either
duplicated record fields (what `vdjdb.txt` does) or paired alpha/beta columns (what `vdjdb_full.txt`
does).

## 11. Tuning and validation of motif clustering

**These are separate datasets and must stay separate.**

### 11.1 Tune on independent-study support

The signal is how many distinct studies report the same clonotype against the same epitope. Measured
on human records in `chunks/`: **4,129 of 187,238 clonotype-epitope pairs (2.21 %) are supported by
≥2 distinct `reference.id`**, spread over 18 epitopes with ≥20 such pairs — GILGFVFTL 2,136,
YLQPRTFLL 619, NLVPMVATV 170, GLCTLVAML 96, RAKFKQLL 61, and a long tail.

A clustering that is finding real convergent selection should preferentially recover exactly those
independently-replicated clonotypes. That is the objective the DBSCAN `coef` (and any Leiden
resolution) is fitted against, per chain, on the pooled geometry.

Independent replication is also a per-record evidence type in its own right, so the tuning signal and
a shipped evidence column come from the same computation.

### 11.2 Validate on TCRvdb, held out, touched once

**TCRvdb / MATCHMAKERS is proprietary — academic, non-commercial, no redistribution in whole or in
part — and it is NEVER shipped, committed or redistributed.** Read only via `$VDJDB_TCRVDB`;
`src/vdjdb/validate/guard.py` fails the build if a copy reaches the repository or a bundle, matching
both filename and column fingerprint so a rename does not defeat it.

It is a **held-out** set: 614 labelled paired HLA-A2 records over exactly two epitopes — YLQPRTFLL
(421 labelled, 40.9 % MM+) and GLCTLVAML (193, 69.4 % MM+). Label definition taken from the author's
own `validation_analytics.Rmd`: MM+ is `padj < 1e-5`, joined on `cdr3` alone with alpha and beta
melted and `min(padj)` per CDR3.

Tuning on it would invalidate it, and two epitopes are too narrow a basis for a global
hyperparameter. Only **aggregate** metrics may be reported — recall of MM+ at a fixed operating
point, precision against MM-, AUROC. Per-record verdicts may never be emitted.

## 12. Phase subplans

One subplan per phase. Each names the files it creates, the facts it consumes (§7, §8 — already
measured, never re-derived), and the single check that closes it. A phase whose check is green is
merged to `dev` and struck here.

Minor decisions taken while executing a subplan are recorded in §13 rather than escalated.

### Phase 1 — `feature/schema`

1. `src/vdjdb/schema/fields.py` — one `Field` record per column: `name`, `dtype`, the eight
   `vdjdb.meta.txt` attributes (`type`, `visible`, `searchable`, `autocomplete`, `data.type`,
   `title`, `comment`), and its position in each of the four positional orders (§2 of
   `docs/outputs.md`). This is the single declaration the six duplicated column lists collapse into.
2. `render_meta(table)` / `render_slim_meta(table)` — emits `vdjdb.meta.txt` /
   `vdjdb.slim.meta.txt` as text. There is **no legacy mode**: the shipped metadata does not describe
   the file it belongs to, so reproducing it would ship a known defect. The three fixes are declared
   ledger rules instead (§5).
3. `header(table)` — the `vdjdb.txt` header derived from the same declaration, restoring the
   `BuildDatabase.groovy:411` invariant the Python port dropped.
4. Tests: parse `src/BuildDatabase.groovy`'s `METADATA_LINES` / `SLIM_METADATA_LINES` **at test
   time** (never a copy) and assert every difference from `render_meta` is one of the three declared
   fixes; assert `header(t) == [f.name for f in fields(t)]` for all three tables.
5. `vdjdb schema --table {vdjdb,slim,full} --format {meta,header,json}` on the CLI.

**Closes when:** every difference between `render_meta` and the Groovy constants is one of the three
declared fixes, asserted by a test that parses `BuildDatabase.groovy` rather than copying it. No
pipeline behaviour changes.

### Phase 2 — `feature/golden-harness`

1. `vdjdb fetch-reference --tag 2026-06-03-ZENODO --out ref/` — `gh release download`, unzip. An
   input being fetched, not a result being cached.
   `ref/` is gitignored; the harness takes a path, so CI passes an artifact instead.
2. `src/vdjdb/compare/diff.py` — three passes: file set → raw + canonical sha256 per file → row-level
   classification keyed on `gene|cdr3|v.segm|j.segm|species|mhc.a|mhc.b|antigen.epitope|reference.id`,
   bucketed `only-in-reference` / `only-in-candidate` / `changed` / `identical`.
3. `rules/expected_diffs.toml` — every changed cell must match a declared rule **and** the rule must
   fire exactly its declared `rows` count. Seed it empty: phase 2's whole point is that it starts at
   zero rules and zero diffs.
4. `vdjdb diff <reference> <candidate> --report out/reports/diff-report.md`, exit 1 on any unattributed
   difference.
5. Run the **current pandas build** into a scratch dir and diff it against the reference zip.

**Closes when:** the current pandas build diffs to zero against the 2026-06-03 release under canonical
equality, with an empty rule file. Nothing downstream starts before this is green.

### Phase 3 — `feature/io-qc`

1. `src/vdjdb/io/chunks.py` — `read_chunks(dir)` → polars, `sorted(glob(...))` for determinism,
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

### Phase 4 — `feature/pipeline-core`

1. `curate/patch.py` — the antigen patch as one join, replacing the per-chunk `T.apply`.
2. `score/confidence.py` — `ScoreFactory` ported to polars expressions; the signature max becomes
   `.max().over(CHUNK_DEDUP_KEY)`.
3. `assemble/master.py` — pairing, `complex.id` allocation, `samples.found` / `studies.found` as
   `len().over(...)` / `n_unique().over(...)`.
4. `emit/legacy.py` — the three legacy tables. Two byte-compatibility helpers live here and nowhere
   else, each carrying a comment saying so: `json_column()` (`map_elements(json.dumps)`, because
   `struct.json_encode()` emits `{"a":"x"}` where the release has `{"a": "x"}`) and
   `py_repr_column()` (`vdjdb_full.txt`'s `cdr3fix.*` are Python `dict` repr).
5. Delete `py_src/` in the same commit that makes it redundant, not before. (Done: the last
   thing needed from it was the CDR3 fixer, vendored into `annotate/_legacy_fixer/` until phase 5.)
6. A peak-RSS test: the whole build under 8 GB, asserted with `resource.getrusage`.

**Closes when:** the ledger shows only the 7,998 `web.cdr3fix.unmp` rows (§7) plus the meta-file fixes,
each as a declared rule with its measured count.

### Phase 5 — `feature/arda-cdr3fix`

1. Dedupe to distinct `(species, cdr3, v, j)` — 191,440 keys — and make **one** `markup_batch` call
   per organism, then join back. 5.2 s total (§7); never loop it.
2. `Cdr3Markup.to_cdr3fix()` emits VDJdb's JSON key-for-key; `v_end` / `j_start` are junction-space,
   which is what the `cdr3` column holds.
3. Delete `src/vdjdb/annotate/_legacy_fixer/` -- the verbatim copy phase 4 bridged through --
   together with `res/segments.txt` and `res/segments.aaparts.txt`.
4. Measure the ledger delta, then **freeze it** as declared rule counts. Expected shape from §7:
   `cdr3` ~2.2 % of rows, `jStart` ~9 %, all of it in the direction arda maps more and earlier.
5. Assert the one-directional property as a test: arda never loses coverage a VDJdb mapping had.

**Closes when:** the ledger's only new differences are inside the declared `arda-cdr3fix` rules, and
the coverage-regression count is 0.

### Phase 6 — `feature/new-format`

1. `identity/` is already built (§10.1) — wire it into the build so every record carries `record_id`.
2. `emit/vdjdb3.py` writes `records.parquet`, `chains.parquet`, `evidence.parquet` and the derived
   `vdjdb.parquet`, plus a TSV projection of each (`docs/outputs.md` §3).
3. `evidence.parquet`'s first producer is `independent_study` — the same computation as the §11.1
   tuning signal, so one implementation serves both.
4. `emit/legacy.py` is re-pointed to read the **new build directory**, never `chunks/`.
5. `vdjdb.schema.json` generated from the phase-1 registry.

**Closes when:** `vdjdb make legacy` from the new build still passes the phase-2 harness, proving the
legacy export is a projection rather than a parallel implementation.

### Phase 7 — `feature/airr`

1. `convert/coords.py` — the only module that converts between the four coordinate spaces
   (`CLAUDE.md`), with round-trip tests. Junction ↔ CDR3 is the two-anchor offset; everything else is
   0-based/1-based and half-open/closed.
2. `emit/airr.py` — Rearrangement + Receptor + Reactivity, both projections of the phase-1 registry's
   AIRR mapping so the two converters cannot diverge.
3. Validate with the `airr` package's own schema validator in CI.
4. Property test: `legacy_to_airr(legacy) == vdjdb3_to_airr(new)` on the fields AIRR can represent.

**Closes when:** `airr.validate_rearrangement` passes on the full table and the round-trip property
holds.

### Phase 8 — `feature/junction-nt`, `feature/segment-guess`, `feature/dgene`

One branch each; all three write **new-format columns only**, so the harness stays green by
construction.

1. **junction-nt** (#461): `vdjtools.model.infer_nt` on the unique `(species, cdr3, v, j)` set, four
   big contiguous slices, never a per-record pool. 3.11 ms/record → ~15 min (§7). **No cache** — the
   ~15 min is inside the budget and the output is authoritative data, not a derived convenience
   (hard rule 9). Test: the generated `cdr3nt` back-translates to the input `cdr3`.
2. **segment-guess** (#462): kmer candidates vectorised, ties broken by one `pgen_aa_batch` call.
3. **dgene**: `arda.dpost.posterior_d` (human IGH/TRB/TRD + mouse TRB only — it returns `None`
   elsewhere rather than guessing, and that `None` must be preserved, not defaulted).
4. **TCR_hash** (#463): keep the legacy hash as-is so structure evidence keeps resolving; the
   re-keying on `record_id` is §10.2's deferred half.
5. **AIRR `Receptor`** rides along: `vdjtools.model.stitch_*` gives the complete mature variable
   domain its two required columns need, and `receptor_hash` is a sha256 over those (§18).

**Closes when:** each branch's new columns are populated, the harness is unchanged, and the
back-translation test passes.

### Phase 9 — `feature/harmonize-rules`

One rule table, `curate/rules/`, one entry per defect, each with its own ledger rule and measured row
count. The counts are already in §7 — do not re-measure:

| Issue | Rule | Rows |
|---|---|---|
| #327 | `TRAJ24` → `*02` where CDR3 carries `WGKLQF` | 75 explicit `*01` + 978 bare `TRAJ24` |
| #389 | `TRAJ24-1` and other malformed gene names | 1 known |
| #564, #467 | MHC-II allele spelling (`DPA` → `DPA1`, `A*24:01`) | see §7 |
| — | murine MHC-II: `I-Ab`/`H2-IAb`/`H2-Ab1` → one IMGT spelling | 3,396 |
| #368 | `antigen.gene` / `antigen.species` in the uncovered 48,935-row tail | ~90 reported |
| #347 | DOI/GitHub `reference.id` → PMID | to measure |
| #561 | identical alpha and beta CDR3 | to measure |

**Closes when:** every rule fires exactly its declared count and nothing else moves.

### Phase 9e — validate the epitope catalogue with `mhcmatch`

`~/vcs/code/mhcmatch` is ours, it is on PyPI, and it already models exactly the table phase 9d
ships: which peptide a given MHC presents. It turns the epitope catalogue from a summary into a
**checked** one, and it goes further than a name lookup can — IPD-IMGT/HLA says whether an allele
*exists*; `mhcmatch` says whether that allele could present that peptide.

Four checks, in increasing strength, each a column on `restriction` and an advisory QC finding:

1. **`pseudoseq.normalize_allele`** on every `mhc.a` / `mhc.b`. It maps a call to a
   pseudosequence-FASTA key, so a call that does not normalise has no 34-mer binding groove and
   nothing downstream can reason about it. Stricter than the IPD-IMGT/HLA prefix check, which only
   asks whether the name is in the registry.
2. **`pseudoseq.load_pseudo("I")` / `load_pseudo("II")`** give the allele's class. Cross-check it
   against the recorded `mhc.class`: a class-I allele on an `MHCII` record is a defect the current
   build cannot see.
3. **`store.infer_class(peptide)`** infers the class from epitope length (MHC-I at ≤ 11). Disagreement
   with `mhc.class` catches a 15-mer filed as MHCI and a 9-mer as MHCII — independent of the allele,
   so it cross-checks check 2 rather than repeating it.
4. **The restriction itself.** `mhcmatch.predict` / `ligand.presented_span` score presentation, and
   `store.Restriction` carries `p_present`, `rank` and `band`. An epitope the recorded allele cannot
   plausibly present is either a wrong `mhc.a` or a wrong `antigen.epitope`, and neither is visible
   to any nomenclature rule. **Advisory and never a rewrite**: a presentation model is evidence about
   a pair, not authority over a publication, and its false-positive rate has to be stated with any
   threshold.

Checks 1–3 are deterministic string and length work and belong in the build. Check 4 needs a model
and its reference data (fetched from `isalgo/pmhc_data` on first use), so it runs as its own CI job
over the 2,381 `(epitope, MHC)` pairs and publishes a report — not inside `vdjdb build`, whose
offline determinism (hard rule 9) must not depend on a download.

**What it would already have caught**, from phase 9d's own findings: `HLA-A*08:01` on 74 records (no
HLA-A\*08 locus exists, so no pseudosequence), the four null/nonexistent `HLA-A*24:*` calls, and
`HLA-B*12` (a serological antigen with no molecular groove). The IPD-IMGT/HLA prefix check found
those; `mhcmatch` would also have ranked the 74 `HPVTKYIM` records against `HLA-A*08:01` and said the
peptide is not presented, which is the part that points at the *fix* rather than only the error.

Adds a dependency on `mhcmatch` for the validation extra only, not for the assembly build.

### Phase 10 — `feature/motifs-tcrnet`

1. `motifs/background.py` — `seqtree.control.load_control` against `isalgo/airr_control`, the
   `*.aa.vdjtools.tsv.gz` builds (§8.7). **Never** let `tcrnet()` resolve its own background: its
   `evalue.background(locus, species)` call takes no `size` and indexes the entire table.
2. `motifs/tcrnet.py` — call `tcrnet()` for `n_control` only, then recompute the legacy statistic with
   its pseudocount (§8.1), which is bit-faithful to `DegreeStatisticsAnnotator.computePValue`:
   `p = (n_control + 1) / (M + 1)`, `p_legacy = binom.sf(d_s - 1, N, p) / (1 - (1 - p)**N)`.
3. File the `E = 0` bug upstream against `vdjtools` in the same pass.
4. `motifs/pwm.py` — the three-level Laplace cascade (`vj_len` → `len` → `uniform`, §8.5) so nothing
   is dropped and `freq` sums to 1. **Assert the sign**: at the 253 affected clusters the new `Σ I`
   must be *lower* than the shipped value. Reproducing the shipped numbers means reproducing the bug.
5. Fix `MotifsScoresAssembler`'s `(epitope, species, gene)` indexing, which never keys on `cdr3`
   (§8.5) and makes both `*_scored.txt` tables wrong.

**Closes when:** the deviation report accounts for every difference against the shipped files by a
named cause, and the 31 logo-less cids and the 1.00 % deleted letter mass are both gone.

### Phase 11 — `feature/motifs-tcremp`

1. `motifs/tcremp.py` — chunked embedding is **mandatory** (§8.8): fit `StandardScaler` + `PCA(50)` on
   a 25k seeded subsample, then embed in 20k chunks. Unchunked peaks at 9.77 GB and OOMs a 16 GB
   runner; chunked peaks at 3.26 GB. Chunks run sequentially, one internally-threaded `embed()` each.
2. Per-epitope DBSCAN with a **chain-global** eps (§8.4): the eps comes from the pooled k-distance
   curve, so it is never re-estimated on a per-epitope n of 30–300.
3. `eps = coef × mean(1st-NN distance)`. Kneedle is a debug cross-check only — it returns knee 1 of
   112,983 at production scale (§8.3). Stop calling the method Kneedle-based.
4. Fit `coef` per chain against the §11.1 independent-study objective. **Never against TCRvdb.**
5. Re-fit the scaler and PCA every build from the seeded subsample; derive cluster ids from cluster
   content so they are stable without storage. Not stored, not cached (hard rule 9). Was: cache by `(cdr3, v, j)` —
   22 MB, which also makes cluster ids stable across releases (today's igraph component numbers are
   not, so every bookmarked motif URL breaks each release).
6. Legacy projection: one legacy cid per `(cluster, stratum)`, `cid = H.B.<epitope>.<n>L<len>`, so
   both files satisfy `vdjdb-web`'s `strict = true` path (§8.6). `--legacy-shared-cid` reproduces
   today's shape.
7. Validation, **once**, at the end: `$VDJDB_TCRVDB`, aggregate metrics only (§11.2).

**Closes when:** the clustering beats the **shipped** `cluster_members_tcremp.txt` re-scored in our own
harness under the vendored `metrics_lib` (§8.9) — not the 0.941/0.569/0.709 write-up, which is a
different method's numbers.

### Phase 12 — `feature/summary`

1. Split `summary/vdjdb_summary.Rmd` at the existing `!summary_embed_end!` marker (line 567) into a
   release dashboard and a paper-figures document. That drops `maps`, `scatterpie` and the never-used
   `ggh4x` from the release path (§7).
2. ggplot2 4.x fixes, all live breakages against the installed 4.0.2: `guide=F` → `guide="none"`,
   `size=` → `linewidth=`, `..count..` → `after_stat(count)`, `as.tibble` → `as_tibble`, and the
   `g_legend` grob-name grep → `cowplot::get_legend`.
3. Kill the live NCBI eutils call: publication years become `summary/pubmed_years.tsv` (**`.tsv`, not
   `.txt` — `summary/*.txt` is gitignored**), refreshed by a scheduled workflow that opens a PR.
4. Data-drive the callouts into `summary/annotations.tsv` carrying `panel, year, label, hjust, vjust`
   and **no coordinates** (#460).
5. #460's actual bug: `grep -o 'width="[0-9]*"'` over the shipped embed HTML returns nothing, so
   `MakeEmbedableHtml.py`'s width rewrite is dead code. Delete it; add
   `style="max-width:100%;height:auto"` to the `<img>` tag instead.
6. Pin rasterisation: `dpi=96, fig.retina=2`, `dev.args=list(type="cairo")` — explicitly cairo, not
   ragg, which would change font rendering and break parity on day one.
7. Replace the cumulative-by-year `expand.grid` cartesian join (~7M rows) with a first-appearance-year
   `cumsum`. It is the only part of the render that could plausibly OOM at 16 GB.
8. `summary/palette.py` is the single palette source for both renderers, emitting `palette.json` that
   the Rmd reads with `jsonlite::fromJSON`.
9. `summary/check_summary.py` — structural (ordered `<h4>`s, 5 tables, 8 base64 PNGs, IHDR-decoded
   width×height, the three vdjdb-web contracts), style (ColorBrewer anchors within ΔE₇₆ < 3, plus an
   anti-assertion against viridis), perceptual (SSIM vs the previous release, fail < 0.55, warn < 0.80).
10. `summary/preview/index.html` pulls Semantic UI from a CDN and `fetch()`es the fragment, because
    opening it directly shows unstyled tables — the classes come from vdjdb-web's bundle.

**Closes when:** both dashboards render offline and all three check layers pass.

### Phase 13 — `feature/docs`

1. Sphinx + `pydata_sphinx_theme`, `conf.py` copied from `arda/docs/`.
   `html_baseurl = "https://docs.isalgo.dev/vdjdb-db/"` — Pages is already provisioned, no setup step.
2. `docs/_ext/vdjdb_schema.py` provides `.. vdjdb-schema::`, `.. vdjdb-vocabulary::` and
   `.. vdjdb-score-rules::` as **directives importing the phase-1 registry at doc-build time**. No
   generated `.rst` in the tree means the tables can never be stale.
3. Structure: `getting-started/`, `standards/`, `submission/`, `builds/`, `dashboard/`, `reference/`.
   `docs/outputs.md` becomes `standards/database-outputs.rst`.
4. The dashboard tab is an `<iframe>`, not inline HTML — the fragment's Semantic UI classes and
   plotly's CSS must not leak into the theme. Its inner document is the phase-12 preview harness.
5. The dashboard artifact downloads from the last successful `build.yml` with `continue-on-error` and
   a committed `placeholder.html`: a 30-minute database build must never block a typo fix in the docs.
6. README shrinks to a ~60-line front door.

**Closes when:** `sphinx-build -W --keep-going` is clean and Pages deploys.

### Phase 14 — `feature/release-tooling`

1. `io/manifest.py` decides what goes in each bundle — an explicit list, never `cp *.txt`, which is
   how seven unshipped side tables nearly shipped (§1).
2. The six-step release job: plan (derive tag, assert it does not exist) → prepare (rewrite
   `latest-version.txt` in the working tree via temp-file + `os.replace`) → build (all zips embed that
   same content) → verify (ledger + dashboard checks) → publish (`environment: release`, required
   reviewer) → **finalize (commit `latest-version.txt`, then `curl -fsI` line 1 and fail on anything
   but 200)**. Step 6 is the one that never happened.
3. Tag scheme `v<YYYY>.<MM>.<PATCH>`; tag creation restricted by ruleset to the release environment.
4. Zenodo via the REST API from the workflow (`newversion` → upload → `PUT` metadata → `publish`),
   replacing the webhook that archives the source tarball rather than the assets. Add the missing
   `version` field to `.zenodo.json`.
5. `release/changelog.py` — reference diff between releases (#432), cheap because the phase-12
   publication-year table is a committed input.
6. `verify-latest` scheduled job: line 1 returns 200 **and** its tag equals `releases/latest`.
7. Retire `.gitlab-ci.yml`, `.travis.yml`, `test.sh`, `release.sh`, `docker.sh`, `release_docker.sh`,
   `gitlab/`, both Dockerfiles and the committed 3.7 MB `docker_build.log`. Move `src/*.groovy` to
   `attic/` — it is the only correct specification for the meta files and phase 1 tests against it.

**Closes when:** a full release dry-run produces three zips, a manifest and a clean ledger.

### Phase 15 — `feature/aldan3-runner`

1. Register aldan3 in a runner group scoped to this repo alone.
2. Retarget with `runs-on: ${{ fromJSON(inputs.motifs-runner) }}` — callers pass `'"ubuntu-latest"'`
   or `'["self-hosted","linux","x64","aldan3"]'`. The naive `runs-on: ${{ inputs.runner }}` cannot
   express a multi-label self-hosted target.
3. Guard every self-hostable job with
   `github.event.pull_request.head.repo.full_name == github.repository`. `chunk-check.yml` stays
   `ubuntu-latest` **always** — it is the job forks trigger.
4. The existing GitLab runner cannot be reused: `gitlab/runner.slurm` asks for
   `--time=00:10:00 --mem=4G`, and a GitHub Actions runner is a long-lived daemon, not an sbatch job.

**Closes when:** a full build completes on both runners with identical canonical digests.

### Bootstrap order — protections last

Land the workflows on `master` → create `dev` → let one full `build.yml` run green on `dev` so the
check names exist → **then** apply the `dev`, `master` and tag rulesets. A required check that has
never run blocks every PR forever. Require linear history on `master`; with 85 accumulated branches
that is the rule that stops it getting worse.

## 13. Minor decisions taken while executing

Recorded rather than escalated. Each is reversible and none changes a shipped contract.

| Date | Decision | Why |
|---|---|---|
| 2026-09-25 | Scratch notes and one-off scripts are gitignored at the repo root (`/NOTES*.md`, `/test_*.py`, `/scratch/`, …), anchored so `tests/` and `docs/` are unaffected | keeps working files out of curation PRs |
| 2026-09-25 | The generated metadata fixes the shipped defects rather than reproducing them; no `legacy=True` mode | the shipped metadata does not describe its own file, in the release *and* in production. Declared as ledger rules instead |
| 2026-09-25 | `vdjdb.score`'s title is `Info`, from production, not `Score` from the release | production is what users see; the release file is the stale one |
| 2026-09-25 | The four `web.*` rows get `data.type = factor`; the surplus `0` is dropped | all four are `visible = 0`, so nothing user-facing moves |
| 2026-09-25 | `database/vdjdb.meta.txt` and `.slim.meta.txt` stay tracked until phase 4 emits them | `database/` is gitignored and they were force-added; `git rm --cached` before the emitter exists breaks a fresh clone's legacy build |
| 2026-09-25 | `TOLERATED_DROPPED` renamed `KEPT_CURATION_COLUMNS` | §9 decided they are kept, so the old name asserted the opposite of the decision |
| 2026-09-25 | The motif tables' own vocabulary (`cdr3aa`, `cid`, `csz`, the PWM columns) is declared in the registry too | they have no `.meta.txt`, so the registry is the only place a rename is caught before it mistypes a positionally-parsed file |
| 2026-09-25 | The definitive tables' column orders live in `schema/fields.py` with every other order, not in `assemble/tables.py` | one registry describes every table, so `vdjdb.schema.json` and `schema --table records` fall out for free |
| 2026-09-25 | `clonotype_id` and `evidence_score` are a `UInt64` hash and a `Float64`; `evidence_id` is the readable `<type>:<gene>` | a hash needs no registry and is stable forever; `evidence_score` must hold a model confidence later, and widening a shipped column is worse than choosing the general type now |
| 2026-09-25 | `vdjdb.parquet` declares all six `evidence.*` columns, `false` where no producer exists yet | the view's shape must not change as phases 8-11 land; absent evidence is an honest `false`, not a missing column |
| 2026-09-25 | Record-level evidence (empty `gene`) raises rather than being dropped | nothing produces it yet, and a silent drop of structure evidence is exactly the class of bug the ledger cannot see |
| 2026-09-25 | The record registry becomes a release asset, not a committed file | 72.7 MB per revision (19.8 MB gzipped) against a 42 MB `chunks/` corpus. See §17 |
| 2026-09-25 | AIRR `Receptor` is deferred to phase 8 | its two required domain columns are the complete mature variable domain, which needs germline stitching. See §18 |
| 2026-09-25 | Legacy->AIRR reuses the tables->AIRR emitter rather than being a second mapping | both source shapes already use VDJdb column names, so a second implementation would only add somewhere to drift |
| 2026-09-25 | mypy's `python_version` is 3.12 while `requires-python` stays 3.11 | numpy's stubs use a 3.12 `type` statement and a stub syntax error stops the check before it reaches our code |
| 2026-09-25 | ruff excludes `src/*.py` and `py_src/`, and ignores `B008` | reformatting code that leaves the tree in phases 5 and 14 would bury the real diff; `B008` is typer's idiom |

## 14. Phase 2 result — what the ledger measures against the 2026-06-03 release

Measured 2026-09-25 by rebuilding with the **unchanged** pandas pipeline and diffing against the
release zip. `git diff 2026-06-03 dev -- py_src/ patches/ res/ chunks/` is empty, so the rebuild ran
identical code on identical data; the only variable is the order `os.listdir` returned.

| File | Raw digest | Canonical digest | Rows | Changed rows |
|---|---|---|---|---|
| `vdjdb_full.txt` | differs | **identical** | 192,753 = 192,753 | **0** |
| `vdjdb.slim.txt` | differs | differs | 197,729 = 197,729 | **0** |
| `vdjdb.meta.txt` | identical | identical | 22 = 22 | 0 |
| `vdjdb.slim.meta.txt` | identical | identical | 17 = 17 | 0 |
| `vdjdb.txt` | differs | differs | 284,546 = 284,546 | **28** (14 groups, 7 symmetric swaps) |

Build cost, measured: **344 s wall, peak RSS 2.16 GB** — not the 64 GB the README claims. The 64 GB
figure is wrong by a factor of 30; what the `T.apply` / `iterrows` hot spots cost is *time*, and the
`.loc` lookups in `generate_default_db` are ~95 % of the 344 s.

### The residue: 7 pairs of records the source data cannot distinguish

28 rows of 284,546 differ, over **14 identity groups**, in **28 of 6.26 million cells (0.00045 %)** -- one cell per row -- and every difference is
**symmetric** — A→B paired with B→A. They are seven swaps of two rows each:

| What swaps | Pairs | Example |
|---|---|---|
| `complex.id` of two complexes | 5 | `12092` ↔ `12093`, `40847` ↔ `40848`, `43824` ↔ `43825`, `84536` ↔ `84537`, `84539` ↔ `84540` |
| `meta.structure.id` letter case | 2 | `5EUO` ↔ `5euo`, `6AVF` ↔ `6avf` |
| `meta.epitope.id` numeric form | 1 | `20354.0` ↔ `20354` |

The cause is that `vdjdb.txt` writes **one row per chain**, so a paired record's TRA row carries
nothing about its beta. Two records sharing an alpha and differing only in a field outside
`CHUNK_DEDUP_KEY` melt to two TRA rows that no single-chain key tells apart, and the legacy build
assigns their annotations in whatever order `os.listdir` returned.

Two real data defects fall out of it, both for phase 9:

- the same PDB entry curated in two letter cases (`5EUO`/`5euo`, `6AVF`/`6avf`);
- a float leaking into `meta.epitope.id` (`20354.0` where every other row has `20354`).

Declared as `legacy-melt-ambiguity-*` in `rules/expected_diffs.toml`. They attribute without a
declared count, because the count is a property of the build host's readdir order rather than of the
build. The new format removes the ambiguity outright: `records.parquet` keys on `record_id` and
`chains.parquet` keys on `(record_id, gene)`, so no annotation is ever assigned by position.

### The ledger had to be made reproducible before any of this could be measured

The first three runs of the same comparison reported **158, 152 and 158** changed rows. Three causes,
all of them "an order nobody chose":

| Cause | Effect |
|---|---|
| `set(ref) \| set(cand)` iteration under `PYTHONHASHSEED` | cell order varied between processes |
| `group_by` without `maintain_order=True`, then a sort whose ties are *exactly* the ambiguous groups | canonical `complex.id` numbering varied |
| an unseeded row hash | group identity varied |

Fixed with a sorted key iteration, an explicit positional tiebreak, `maintain_order=True`, and
`SEED = 20260925` from `vdjdb.config`. The measured difference then fell from ~155 cells to 28 — the
instability was manufacturing five times the real difference. Verified: identical cell digest across
four `PYTHONHASHSEED` values, and a unit test asserts one digest over five runs.

This is why `CLAUDE.md` hard rule 7 exists. A rule in the ledger declares a measured row count, and
a count is meaningless against a measurement that moves.

### Cost

The row comparison first materialised 6.26 million Python tuples per table. Replaced with a hashed
`(key, row)` group-count join in polars, so only the keys whose multisets actually disagree — **14 of
208,447** — are ever pulled into Python. Wall time 7.6 s → **3.5 s** for a 425 MB bundle, at 222 % CPU.


## 15. Reconciliation — where this drifted, 2026-09-25

Recorded so the correction is not re-litigated. Each item is a case of inferring the model from the
*shape of the data* instead of reading `README.md`, which is authoritative until `docs/standards/`
replaces it (phase 13).

| Drifted | Correct |
|---|---|
| Called the 19 cross-chunk matches "duplicate submissions" and gave both rows one `record_id` | **A chunk is one paper.** Two chunks are two independent reports, whatever their fields. They get separate ids, and their agreement is evidence — the signal §11.1 tunes motif clustering against |
| Invented a `submissions` table, then a `curation` table, to hold `method.*` and `meta.*` | The README says those columns record how the *publication* established the specificity. They describe the record. There is no submission unit: **one chunk row is one record**, and it reports both chains |
| Built the legacy tables directly from the master frame | The definitive tidy tables are the database; every shipped file is a join and a pivot off them. `emit/legacy.py` is the only module that may know about `complex.id` or the JSON blobs |
| Assigned record ids after CDR3 repair | Ids go on what the publications reported. Repair merged 215 pairs of separately reported records |
| Let phase 6 own the definitive tables while phase 4 owned legacy | Backwards: phase 4 builds the tables and *derives* legacy from them. Phase 6 ships them |

The three earlier phases are unaffected — the field registry, the ledger and the reader make no
claim about what a record is. The cost was confined to phase 4, and the ledger caught every symptom.


## 16. Phase 5 measurement — `arda.cdr3fix` against the vendored k-mer scanner

Measured 2026-09-25 on all **191,447** distinct `(species, cdr3, v, j)` keys of the current corpus,
both engines run through the same pipeline. Not the 20,000-row sample of §7: that sample was drawn
from the *released* table, whose sequences had already been repaired, and every row in it already
had a V and a J.

| Field | Agreement |
|---|---|
| `cdr3` | 97.74 % (187,112) |
| `jFixType` | 96.09 % |
| `vFixType` | 93.69 % |
| `vEnd` | 91.90 % |
| `good` | 91.17 % |
| `jStart` | 88.62 % |

### Coverage is **not** one-directional — this corrects §7

| | legacy maps, arda does not | arda maps, legacy does not | both map |
|---|---|---|---|
| `vEnd` | **7,807** | 3,849 | 176,573 |
| `jStart` | 248 | **6,987** | 183,334 |

§7 recorded "0 coverage regressions" from the 20,000-row sample. On the full corpus arda **loses
7,807 V-end mappings**. The J side is the reverse and much larger in arda's favour. Where both map,
the shift is small and signed as expected: `jStart` mean −0.17 (14,512 smaller, 168,785 equal, 37
larger), `vEnd` mean +0.04.

### Two things arda does not do

1. **It does not guess a segment.** Given a blank `v`, `markup_records` reports `FailedBadSegment`
   rather than proposing one, and the record then fails the legacy "a CDR3 needs a V and a J"
   filter. Swapping both halves at once dropped **13,844 of 284,546** rows. The k-mer guesser
   therefore stays, isolated in `src/vdjdb/_legacy_guess.py`, until #462 replaces it with OLGA Pgen
   scoring in phase 8.
2. **It returns an empty segment id on failure** where the legacy kept the closest match. Passing
   that through cost another 11,619 rows, so the given call is restored when arda declines to name
   one; the coordinates stay −1, which is the honest half of the answer.

### `max_replace = 0`, not the legacy's 1

At 1, arda rewrites **4,486** curated human CDR3s to conform to the germline it was handed, against
649 at 0. The rewrites are wrong for this database: `CAAADSWGKLQF` with `TRAJ24*01` becomes
`CAAADSWGKLEF`, because `WGKLEF` is what `*01` encodes — but `WGKLQF` is the `*02` signature, so the
sequence is right and the **allele call** is wrong (#327: 66 % of explicit `*01` calls carry the
`*02` motif). Substituting the residue destroys the evidence that would fix the call.

It costs nothing: both ends map on **165,223** of 174,630 distinct human keys at either setting, and
`good` is marginally higher at 0. Trimming and extending are unaffected (2,588 against 2,620),
because those repair a truncated sequence rather than contradicting a reported one.

### Ledger position, and why the swap is not yet the default

With `max_replace = 0` and the guesser retained, against the 2026-06-03 release:

| | |
|---|---|
| `vdjdb.txt` rows | 284,546 = 284,546, **14,778 records change their repaired CDR3** |
| `cdr3fix.*` cells | 27,083 |
| `v.end` / `j.start` cells | ~25,000 |

`fix_cdr3(engine=...)` defaults to **`legacy`**, so `dev` stays green and the shipped build is
unchanged. The swap is ready and characterised; accepting it is a decision about shipped sequence
data, not a refactor:

- it changes the repaired `cdr3` of 14,778 records, which changes their identity key and therefore
  their `record_id`;
- it trades 7,807 V-end mappings for 6,987 J-start mappings;
- **#327 should land first.** Repairing against a wrong allele call is what produces the worst of
  these differences, and phase 9 fixes the calls.

### A cross-repo bug found and fixed on the way

`arda.paths._source_root()` decided it was running from an arda checkout by walking up from its own
`__file__` for a directory with `database/` and a project marker. Installed into this project's
`.venv`, that walk reaches **this repository**, which has both — so arda resolved its reference to
`vdjdb-db/database/vdj`, which does not exist, `load_anchors` returned `{}`, and all 191,447 CDR3s
came back `FailedBadSegment` with `vEnd = -1`. Nothing reported a problem.

Fixed in `antigenomics/arda` on `fix/source-root-marker` (commit `d40095c`, 1,118 unit tests pass):
`_source_root()` now requires `database/vdj`, and `load_anchors` raises instead of returning an
empty dict. **Cross-repo gate:** that needs an arda release before CI can rely on it. Until then
`vdjdb.annotate.cdr3fix.ensure_reference()` detects the condition and repoints `$ARDA_HOME` at the
per-user cache, so the build is correct on either arda version.


## 17. Phase 6 result — the definitive tables as they ship

Measured 2026-09-25 on the full corpus, `engine="legacy"` (the arda swap waits for #327, §16).

| Table | Rows | Cols | parquet | TSV |
|---|---|---|---|---|
| `records` | 192,753 | 33 | 1.2 MB | 45.8 MB |
| `chains` | 286,047 | 17 | 9.2 MB | 54.0 MB |
| `evidence` | 53,913 | 8 | 0.1 MB | 8.7 MB |
| `vdjdb` (the joined view) | 286,047 | 54 | 10.8 MB | 129.4 MB |

`vdjdb.schema.json` is 38.7 kB: every column of every declared table, its `vdjdb.meta.txt`
attributes, its position in each table it appears in, and the dtype **read off the written frame**
rather than declared, so the schema cannot claim a type the files do not have.

### The closing criterion holds

```
vdjdb build --out out/                       # tables, then legacy from them
vdjdb make legacy --tables out/tables        # legacy from the parquet that shipped
vdjdb diff ref/vdjdb-2026-06-03.zip out/legacy-made   -> PASS
```

All five legacy members are **byte-identical** whether projected from the in-memory tables or read
back from parquet, and the ledger verdict is unchanged: every difference is a declared rule firing
its exact measured count. The legacy export is a projection of the database, not a second
implementation of it.

### `independent_study` — the first evidence producer

The same computation as the §11.1 tuning objective, in one place, because "two papers found this
receptor against this epitope" is both the strongest evidence a record carries and the signal a
clustering must recover to be believed.

| Scope | Clonotype-epitope pairs | With >= 2 distinct `reference.id` |
|---|---|---|
| human, on `chunks/` as submitted (§11.1) | 187,238 | 4,129 (2.21 %) |
| human, on `chains` after CDR3 repair | 184,660 | 4,974 (2.69 %) |
| all species, after repair | 201,825 | 5,047 (2.50 %) |

The two human rows are the same definition at two stages: repair merges sequences, so pairs fall by
2,578 and replication rises by 845. **The shipped evidence and the tuning objective both use the
post-repair number**, because that is what the database contains; §11.1's figure stands as the
`chunks/`-level measurement it was.

That yields **53,913 evidence rows over 48,893 of 192,753 records (25.4 %)** — far above the 2.5 %
of *pairs*, because the replicated clonotypes are the popular ones and each carries many records.
`evidence_score` ranges 2 to 41 distinct references; `evidence_value` lists *the other* references,
not the one the reader is already holding.

### `clonotype_id`

A seeded hash (`config.SEED`) of `(species, gene, cdr3, v.segm, j.segm)`, not a counter: a counter
renumbers every clonotype the moment a chunk is added, and this id is what accumulated evidence
joins on. **187,984 distinct ids for 187,984 distinct keys** — no collision, no split. Asserted on
every real build (`tests/release/test_tables_contract.py`).

### One defect the contract tests found

34 `chains` rows were the tables' only nulls: a D-segment call with no CDR3, so the fixer was never
handed anything and left no result. Filled — `-1` for the unmapped coordinates, `""` for the
strings, `false` for the flags — because empty string is the only missing marker (CLAUDE.md rule 6)
and a null in a shipped table is the pandas three-way ambiguity coming back. The legacy export gates
all 34 on `cdr3 != ""`, so nothing shipped moved: verified byte-identical against the build made
before the fix.

### `record_id` is stable across builds, not yet across releases

`registry/records.tsv` is not written by the build and not committed. The registry reconciles against
an empty one every time, so ids are deterministic given the corpus but would shift the moment a chunk
is added — which is the failure mode `identity/` exists to prevent.

Measured: the registry is **72.7 MB** for 192,753 records (19.8 MB gzipped), and most of it is the
packed previous natural key that amendment tracing needs. Committing it would add ~20 MB to the repo
per curation PR, against a `chunks/` corpus of 42 MB. So it becomes a **release asset** that the
build fetches, reviewed in the release diff rather than the PR diff — phase 14, where release assets
already live. Until then, ship the tables without leaning on id stability across releases.


## 18. Phase 7 result — AIRR

Measured 2026-09-25 against the `airr` Python package 2.0.0, AIRR schema version 2.0.

| File | Level | Rows | Size |
|---|---|---|---|
| `vdjdb.rearrangement.tsv` | one per chain | 286,047 | 33.8 MB |
| `vdjdb.reactivity.tsv` | one per record | 192,753 | 33.5 MB |

`airr.validate_rearrangement` **passes on the full 286,047-row table** in 1.0 s with no warnings, and
every emitted column is a real `Rearrangement` property. The nucleotide fields (`sequence`,
`junction`, the alignments, the three cigars) are present and empty, which the schema accepts — it
requires the column, not a value. Phase 8 (#461) fills them.

### One implementation, two source shapes

The tidy `chains` table and legacy `vdjdb.txt` already use the **same column names** — `cdr3`,
`v.segm`, `j.segm` — so `rearrangement()` and `reactivity()` take one frame in VDJdb vocabulary and
the "two converters" are two ways of assembling their input, not two mappings that can drift. The
property is therefore stronger than the planned equality, and it holds exactly on the corpus:

| | Rearrangement | Reactivity |
|---|---|---|
| rows the legacy path produces that the tables path does not | **0** | **0** |
| combinations where the legacy path has *more* | **0** | **0** |
| rows the tables path has and legacy does not | 1,501 | 1,141 |

The excess is exactly what the legacy build discards: the **1,467 chains of 1,141 records** whose
chain carried a CDR3 with no V or J (it drops the record whole, both chains) **plus 34 D-only
chains**. So the new format keeps 1,141 records and 1,501 chains that the legacy release loses.

`d_call` is excluded from the comparison: legacy `vdjdb.txt` has no D column at all, so the 42,574
beta chains with a D call are information the file cannot carry, not a disagreement.

### The Reactivity mapping is the spec's own suggestion

AIRR `Reactivity` models a measurement with a value and a unit; VDJdb records a curated assertion.
The spec anticipates exactly this: `reactivity_method` is *"delineated as `annotated` if annotated
from an external source"*, and `reactivity_readout` *"for inferred and annotated methods should
indicate a confidence/quality level"*. So `reactivity_readout = confidence` and
`reactivity_value = vdjdb.score`, which is precisely a confidence in the specificity annotation.

`reactivity_method` is classified from the 61 distinct `method.identification` values in the corpus:

| Keyword class | `reactivity_method` | Records |
|---|---|---|
| tetramer / dextramer / pentamer / multimer / streptamer / monomer / MHC-peptide-beads | `MHC_peptide_multimer` | 147,967 |
| antigen-expressing or antigen-loaded targets, T-Scan | `native_protein` | 12,597 |
| everything else | `annotated` | 32,189 |

This is **not** `emit.legacy._web_method`, which answers a different question (a coarse web filter
class, `sort` / `culture` / `other`). A CD137-expression sort is a `sort` there and is not a multimer
assay here; reusing it would mislabel 462 records.

The two files link by `cell_id = record_id`: a VDJdb record is one publication's report on one T-cell
clone, and a clone is what AIRR's `Cell` names. Both are standard AIRR fields, so nothing invents a
join key.

### `Receptor` is deferred to phase 8, deliberately

AIRR `Receptor` requires `receptor_variable_domain_{1,2}_aa` — the **complete mature variable
domain**, non-nullable. VDJdb has the junction and the allele calls, so producing it means stitching
germline V and J around the junction (`vdjtools.model.stitch_*`), which is the same tooling phase 8
already brings in for `cdr3nt`. Emitting the file now with its two required columns empty would be
worse than not emitting it. Note also that `receptor_hash` is a sha256 over the stitched domains and
is **not** VDJdb's `TCR_hash`, which hashes CDR3s, segments, MHC and epitope.

Separately: `airr` 2.0.0 exposes validators for `Rearrangement` and `Repertoire` only. `Receptor` and
`Reactivity` are in its shipped `airr-schema.yaml` but not in `airr.AIRRSchema`, so the Reactivity
file is gated against our own registry and against the spec's field list read from that YAML, not by
the package's validator. Stated rather than implied.

### One curation defect this surfaced

**854 records carry no `reference.id` at all**, every one from `luciani-samir-etal-hcv-14-09-2018`
(822 at `vdjdb.score` 0, 32 at 1). The QC rule permits a blank one — `_blank("reference.id") | ...`
reads as "blank is acceptable" — which is defensible under the README's *"submitter details in case
unpublished"* but is not what a blank means. Phase 9 / #347.


## 19. Phase 8a result — inferred junction nucleotides (#461)

Measured 2026-09-25 on the full corpus. `vdjtools.model.infer_nt` on the distinct
`(species, gene, cdr3, v.segm, j.segm)` set, joined back. **Nothing cached** (hard rule 9).

| | Chains | With `cdr3nt` |
|---|---|---|
| HomoSapiens | 262,385 | 250,057 (95.3 %) |
| MusMusculus | 21,891 | 11,040 (50.4 %) |
| MacacaMulatta | 1,771 | 0 |
| **Total** | **286,047** | **261,097 (91.3 %)** |

**Zero back-translation mismatches**: every one of the 261,097 inferred sequences translates back to
the junction it was inferred from. That is #461's acceptance criterion and it holds exactly.

`cdr3nt.pgen` spans 3.87 × 10⁻⁶⁶ to 1.62 × 10⁻⁵. `cdr3nt.margin` — the winner's Pgen over the
runner-up's — has median 2.32, but **9.4 % (24,574) fall below 1.1**, where the choice among
synonymous recombination histories was near-arbitrary, and 9,536 had a single candidate (reported as
infinity). The margin column exists so a consumer can filter on exactly that.

### `cdr3nt` is inferred, not observed, and the models disagree

On 600 distinct human TRB keys the OLGA and arda models agree on only **7.2 %** of the nucleotide
sequences they both return (293 both-resolved). They disagree about which synonymous nucleotide
history is most likely, never about the protein. So `cdr3nt` is a plausible representative and must
never be treated as evidence — the field comment in the registry says so.

### Mouse coverage is limited by the model, not by the data

OLGA is human-only (`load_bundled` raises and names arda as the alternative), so mouse must use the
arda source, which is systematically stricter: on the same 600 human TRB keys arda declines **302**
where OLGA declines **30**. Mouse's 50.4 % is that strictness, not a property of murine records.
A mouse model of OLGA's permissiveness would lift ~11,000 chains; that is an upstream `vdjtools`
question, recorded here rather than worked around.

### V/J calls are resolved before the run, and never invented

The model is keyed by allele and raises on anything else, so every distinct call is resolved up
front, three ways: a known allele is used as given; a **gene name whose model carries exactly one
allele** becomes that allele, because there is no choice to make (VDJdb has 2,371 V and 1,784 J calls
with no allele at all, #389); anything else marginalises over that segment. Picking an allele for a
multi-allele gene would be precisely the #327 mistake.

The "anything else" bucket is **3,500 of 284,764 chains (1.2 %)** and reads as a work list for
phase 9:

| Call | Chains | What is wrong |
|---|---|---|
| `TRAV14` | 1,032 | no allele, and the model has several |
| `TRBV21-1*01` / `TRBV21-1` | 303 | pseudogene, carried by no model |
| `TRAV21-DV12` / `TRAV21/DV12*01` | 196 | dash against the model's slash |
| `TRBV13-1*02` (mouse) | 235 | allele absent from the model |
| `TRBJ1-6*02` | 99 | allele absent from the model |
| `TRBJ1.2`, `TRBJ 2-7`, `TRAJ16.5`, `TRAJ01-1*01` | tens | a dot, a space, or zero-padding where a dash belongs |

### Parallelism

Four contiguous slices of the sorted key set, one thread each, reassembled in **slice order** — the
worker count cannot change the answer, and a test asserts that. Measured speedup 2.11× on 4 threads
(the native call releases the GIL only partly), against 2.70 ms per human TRB key single-threaded.
Never a pool of per-record tasks: dispatch on 114k one-row tasks would cost more than it saves.

The whole build goes from **16 s to 170 s** — inference is now 90 % of it. Still half the 344 s the
pandas pipeline took to produce three files and no nucleotides, and nowhere near the runner budget,
which is what makes "recompute every time" (hard rule 9) an easy rule to keep.


## 20. Phase 8b result — the D segment, and how much to believe it

Measured 2026-09-25 on the full corpus. Two sources, each used for the thing it is right about.

**Geometry comes from the junction scenario.** `d.inferred`, `d.start` and `d.end` are taken from the
same `infer_nt` result that produced `cdr3nt`, so the coordinates index the nucleotide sequence
shipped beside them. Taking them from a second model — the original plan's `arda.dpost` — would ship
positions that point at a different hypothetical sequence. 0-based half-open, `cdr3nt` space.

**Confidence comes from `arda.dpost`.** `d.posterior` is the posterior for the gene `d.inferred`
names — not for arda's own winner. The two models name the same D gene on only **78.5 %** of chains,
and the number printed beside a call must be the probability of *that* call, so a disagreement shows
up as a low posterior rather than as a confident wrong answer.

| | Beta chains |
|---|---|
| with an inferred D | 155,019 of 163,117 (95.0 %) |
| `d.posterior` median | 0.728 |
| `d.posterior` below 0.6 | **53,599 (34.6 %)** |
| `d.entropy` median | 0.796 |
| `d.entropy` above 0.9 | **51,950 (33.5 %)** |
| inferred D gene == curated `d.segm` gene | 32,052 of 40,892 (78.4 %) |

**A third of beta D calls are close to undecidable**, which is what a short, heavily trimmed segment
with two similar candidates actually looks like. `d.posterior` is therefore not a decoration: a
consumer that filters on it is doing the only correct thing with a D call. The curated `d.segm` is
left untouched — it is what the publication reported, and it is not this pipeline's to overwrite.

Verified on every build: no `d.end` exceeds its `cdr3nt` length, no `d.start >= d.end`, and no TRA
chain carries a D. arda returns no posterior on 180 of 155,019.

Cost: 9,131 junctions/s, so 13 s on top of the junction inference — the full build is **183 s**.

### The arda reference trap, again

`arda.dpost.posterior_d` returned `None` for **all 200** junctions of the first measurement, because
arda mistakes this repository for its own checkout and loads no anchors — the same defect
`fix/source-root-marker` fixes upstream and `annotate/cdr3fix.ensure_reference()` works around here.
It fails silently, by returning a legitimate-looking "no answer". `add_d_posterior` calls
`ensure_reference()` for exactly this reason. **Any new arda entry point must do the same** until
that release lands (§3 gates).


## 21. Phase 8c result — the V/J guesser (#462), and a bug it uncovered

**VDJdb's V guesser has never worked.** `Cdr3Fixer.guess_id` puts `return ""` **inside** the
five-prime loop, so it tries exactly one prefix length and gives up: measured, **3 non-empty V
guesses in 4,000 sequences**, against 3,797 for J, whose branch has the same statement correctly in
a `for...else`. One level of indentation, and it has been shipping since the file was written.

The consequence is not cosmetic. 711 chains carry a CDR3 with no V. The guesser never supplies one,
so they fail the legacy build's "a CDR3 needs a V and a J" test, and **their records are dropped
from `vdjdb.txt` entirely** — part of the 1,141 records §18 measured the new format keeping.

### Pgen against the k-mer scan, on hidden ground truth

1,200 human chains per locus, curated call hidden, seeded sample:

| Locus | k-mer scan | Pgen (`infer_nt`) |
|---|---|---|
| TRB V | **0 of 1,189 (0.0 %)** | 283 (23.8 %) |
| TRA V | **1 of 1,111 (0.1 %)** | 557 (50.1 %) |
| TRB J | 1,133 (95.9 %) | **1,153 (97.5 %)** |
| TRA J | 882 (79.3 %) | **1,065 (95.8 %)** |

The J rows are a fair comparison and Pgen wins on both. The V rows are not a comparison at all —
they are the bug above. Note what Pgen's V numbers say on their own terms: recovering a V from the
junction alone is right about a quarter of the time for TRB and half for TRA, because TRBV
contributes only a few junction residues. That is 12x chance for TRB, and it is still a guess. The
registry comment carries these numbers so nobody reads `v.inferred` as a call.

### What ships

`v.inferred` and `j.inferred`, filled **only where the curator named no segment** — 686 of the 711
chains with no V, 298 of the 596 with no J — and never beside a curated call, which a release test
asserts. The curated columns are untouched and the ledger still reads PASS.

**Whether the legacy build should start keeping those 711 chains' records is not decided here.** It
is a change to shipped data of exactly the kind #327 is, so it belongs with the nomenclature work in
phase 9, where the allele calls those records carry are being fixed anyway.

Build 183 s -> 185 s.


## 22. Phase 8d result — AIRR `Receptor`, and a stitch that factorises

`Receptor` needs `receptor_variable_domain_{1,2}_aa`: the **complete mature variable domain**, *"from
and including the first AA after the signal peptide to and including the last AA that is completely
encoded by the J gene"*, non-nullable. VDJdb has a junction and two allele calls, so the domain is
rebuilt — V framework 5' of Cys104, the nucleotide junction from phase 8a, J framework 3' of
[FW]118 — and translated. A spot check against IMGT: `TRBV6-1*01` yields
`NAGVTQTPKFQVLKTGQSMTLQC…`, which is the mature sequence, leader correctly absent.

| | Records |
|---|---|
| total | 192,753 |
| paired (TRA **and** TRB) | 93,294 (48.4 %) |
| with both domains rebuilt -> a `Receptor` row | **81,003** |
| distinct `receptor_hash` among them | 71,690 |

Domain lengths run 106–127 aa. Only paired records appear, because a receptor is a two-domain object
and both columns are non-nullable; an unpaired record is not dropped, it simply lives in the
Rearrangement file, which is where AIRR puts a single rearranged sequence. The 9,313-row gap between
81,003 receptors and 71,690 hashes is the same receptor reported by more than one record — the
independent-replication signal of §11.1, at receptor granularity.

### The stitch factorises, so it is one expression instead of 200k calls

`stitch_contig` is a Python function over one sequence, measured at 1,962/s — 2.5 minutes for the
corpus. But its V part depends only on the V allele and its J part only on the J allele, so probing
each allele **once** with a sentinel junction yields two small lookup tables and the contig becomes a
`concat_str`. Verified against `stitch_contig` itself on 3,987 real chains: **identical on all of
them, 0 mismatches**. The stage runs in **2.7 s** for all 286,047 chains.

That is CLAUDE.md rule 8 doing its job — vectorise before parallelising. A thread pool around the
per-sequence call would have been the obvious move, would have been GIL-bound, and would have been
roughly 50x slower than noticing that the function is separable.

`receptor_hash` is AIRR's own: sha256 over the two concatenated domains. It is **not** VDJdb's
`TCR_hash`, which hashes CDR3s, segments, MHC and epitope and is what the structure store is keyed
on. Two hashes, two purposes, both kept (#463 keeps the legacy one as-is).

**Phase 8 is complete**: `cdr3nt` + Pgen + margin, D geometry and confidence, V/J inference for
curation gaps, and the AIRR Receptor. Build 185 s, ledger PASS.


## 23. Phase 9a result — IMGT segment nomenclature (#389)

`proofreading/imgt_alleles.tsv.gz` has been in the repository, unread by any build code, since the
`proofreading/` directory was created. This is what makes it the authority it was collected to be.

**Nothing is invented and nothing is guessed between two candidates.** A call IMGT already knows is
left alone; otherwise a small set of mechanical respellings is generated and the call is rewritten
**only if exactly one of them is an IMGT name for that species**. Of 2,478 non-IMGT calls, 1,928 are
resolved that way — 2,325 records, 105 distinct rewrites.

| Respelling | Example | Records |
|---|---|---|
| restore the `/DV` name IMGT uses for both loci | `TRAV14` → `TRAV14/DV4` | 1,377 |
| drop a D gene's `-1` | `TRBD2-1*01` → `TRBD2*01` | 224 |
| `-DV` → `/DV` | `TRAV21-DV12` → `TRAV21/DV12` | 181 |
| sort a multi-call | `TRBD2,TRBD1` → `TRBD1,TRBD2` | 89 |
| `.` → `-` | `TRBJ1.2` → `TRBJ1-2` | 21 |
| insert the missing slash | `TRAV29DV5` → `TRAV29/DV5` | 10 |
| `TCR` → `TR` | `TCRBD2*02` → `TRBD2*02` | 5 |
| strip a space, including a non-breaking one | `TRAJ12*01 ` → `TRAJ12*01` | 6 |

The 550 calls left are **not spelling problems** and are reported rather than forced: 438 macaque V
calls (1,333 of 1,771 macaque V calls are valid rhesus IMGT, 206 are valid *human* names applied to
macaque records, 232 are neither — rewriting a human gene name to a rhesus one asserts an orthology
this build has no basis for), and names with several IMGT candidates (`TRBV8`, `TRBV7`, `TRAV15`) or
none at all (`TRAJ16.5`).

### What it cascades into, and why that is the point

Correcting the name lets the CDR3 fixer find the germline it never could:

| | Chains |
|---|---|
| V call changed | 1,950 |
| J call changed | 32 |
| **V-end mappings gained** (was `-1`) | **1,511** |
| V-end mappings lost | **0** |
| repaired CDR3 changed | 2 |

`get_closest_id` tries the name, then `<name>*01`, then `<name>-1*01` … `<name>-100*01`. `TRAV14/DV4*01`
is on none of those paths, so 1,032 chains were shipping with no V germline at all. One-directional
gain, 1,511 to 0.

It also moves `TCR_hash` (1,502 cells across the three files, since the hash includes `v.alpha`),
`web.cdr3fix.unmp` (1,340 rows now correctly mapped), `samples.found` (267, the sample signature
includes `v.alpha`), and `vdjdb.score` on 8 records — 2 gaining a score from a neighbour they now
share a signature with, **3 losing one they were being credited with by a record that is not in fact
the same clonotype**. The losses are the more interesting half.

### A new ledger primitive: declared renames

A correction to a column that is part of the identity key produces no changed cell — it removes a row
and adds one, and the cell machinery has nothing to attribute. `[[rename]]` declarations are applied
to the **reference** before keying, so the ledger goes on measuring what *else* moved. Three things
were needed to make that sound, and each was a real failure first:

1. **Apply them simultaneously, per column.** Sequentially, `A → B` and `B → C` chain, and a cell
   that was already `B` comes out `C`: that silently moved mouse `TRAV6-1*01` rows and turned a clean
   comparison into 6,334 phantom unmatched rows in the reference.
2. **Declare the post-fixer value, not the intermediate.** The fixer writes the resolved name back,
   so harmonising `TRAV14` puts `TRAV14/DV4*01` in the file, not `TRAV14/DV4`. A rename declaring the
   intermediate rewrites the reference and rescues no row at all — worse than declaring nothing.
3. **Only declare *injective* renames.** When the fixer does not leave the old spelling alone it has
   already mapped it onto a real allele — `TRAV6-7-DV9` simplifies to `TRAV6` and lands on
   `TRAV6-1*01` — so the reference is indistinguishable from records that genuinely carry that
   allele, and the rename would rewrite both. Those become a declared row delta instead: 430 rows in
   `vdjdb.txt`, 373 in `vdjdb_full.txt`, 356/355 in `vdjdb.slim.txt`.

A rename that matches **nothing** fails the run, so a declaration cannot outlive the data it
describes. That caught 53 stale multi-call declarations immediately: `fix_both` splits a multi-call
and keeps the best member, so what ships is a selection, not a rename.

`vdjdb rules` regenerates the block from a build; the reviewable artifact is its diff in a curation
pull request, and the ledger's complementary job is proving nothing else moved. **Verdict: PASS** with
28 declared rules, 32 renames and 3 row deltas.


## 24. Phase 9b result — TRAJ24*01 vs *02 (#327), the arda gate

`TRAJ24*01` encodes `…GGK**FE**F…` and `*02` `…GGK**LQ**F…`: two residues apart, both inside the
junction, so **the sequence is evidence and the submitted call is not.**

Measured on the corpus, human, across the whole TRAJ24 family (1,444 records):

| `j.alpha` as submitted | Records | carrying `WGKLQF` (*02) | carrying `WGKFEF` (*01) |
|---|---|---|---|
| `TRAJ24` (no allele) | 1,298 | 974 | **0** |
| `TRAJ24*01` | 111 | **73** | **0** |
| `TRAJ24*02` | 34 | 33 | **0** |
| `TRAJ24-1` | 1 | 0 | 0 |

**`WGKFEF` appears zero times in the entire corpus.** The original report was that about two thirds of
explicit `*01` calls are probably `*02`; the sequence says it more strongly — not one of the 111 carries
the `*01` signature, and 73 carry the other one. The 364 with neither have a CDR3 trimmed short of the
anchor: no evidence, no correction.

**1,047 records corrected, and every one of them gains a J germline mapping:**

| | Chains |
|---|---|
| `j.segm` changed | 1,047 |
| **`j.start` mappings gained** (was `-1`) | **1,047** |
| `j.start` mappings lost | **0** |
| repaired CDR3 changed | **0** |

Afterwards `TRAJ24*02` has 1,081 chains of which 1,080 map (99.9 %), while `TRAJ24*01` keeps 439 of
which only 66 map — those are the no-signature records, whose CDR3 genuinely does not reach the
anchor. The correction did not need to alter a single sequence; the repair simply could not place them
against the wrong allele.

`res/segments.txt` carries both alleles, so nothing is lost. OLGA's model carries only `TRAJ24*01`, so
the 1,047 now marginalise over J in junction-nucleotide inference rather than pinning it — a small,
recorded cost of being right.

### This is the gate §16 named

The arda swap was held because *"repairing against a wrong allele call is what produces the worst of
these differences, and phase 9 fixes the calls."* It is now fixed, together with #389's 1,511 gained
V-end mappings. The swap can be re-measured against a corpus whose allele calls are correct, which is
what it was always waiting for.

### Conditional renames

An allele correction is **not injective on value alone**: the fixer resolves a bare `TRAJ24` to
`*01`, so the reference ships the same `TRAJ24*01` for the 1,047 records the CDR3 corrects and the 38
it does not. So a rename may now carry the **same predicate the rule used** —
`when_columns = "cdr3,cdr3.alpha"`, `when_contains = "WGKLQF"` — and the ledger applies it to the
reference under the same evidence. Conditional renames are evaluated against a snapshot taken before
any of them apply, for the same reason the unconditional ones share one mapping: otherwise they chain.

**Verdict: PASS**, 0 unattributed cells, with 11 extra row-delta rows per file where correcting the
allele changed the repair enough to break injectivity.


## 25. Phase 9c result — MHC (#467, #564, and the fragmentation)

`proofreading/mhc.md` states the convention and `proofreading/mhc_alleles.tsv.gz` (46,005 alleles) is
the authority. Neither had been read by build code. Three corrections, 372 records, each with its own
evidence.

### #564 is already fixed — the issue is stale

`HLA-DPA*01:03` against `HLA-DPA1*…`, and `HLA-DRA1*…` against `HLA-DRA*…`. Measured on `chunks/`:
**zero rows carry a malformed class-II gene symbol.** The records use `HLA-DQA1` (8,596), `HLA-DRA`
(3,125) and `HLA-DPA1` (1,519), which are exactly the authority's symbols — note `HLA-DRA` has no
digit and `HLA-DPA1` does. Fixed in the corpus in June 2026, 1,417 rows across 5 chunks. **#564 can be
closed.**

### #467 — an allele that does not exist

All **80** `HLA-A*24:01` records come from one reference, `doi:10.1016/j.xcrm.2023.101017`, which
reports testing in `A*24:02`. IPD-IMGT/HLA lists **no `A*24:01` at any resolution** — 0 rows against
342 for `A*24:02`. Corrected. (The same 80 records are also #347's, since that reference is a DOI.)

### Murine class-II fragmentation

`vdjdb-web` groups motifs by the MHC string, so several spellings of one molecule split its records
and cost the smaller groups their motif badge. Collapsed onto the spelling that both `mhc.md`'s
convention and the data already prefer, so nothing new is introduced:

| From | To | Records | Why |
|---|---|---|---|
| `H2-IAb` | `I-Ab` | 113 | `mhc.md` names class II `I-<locus><haplotype>`, and `I-Ab` dominates 1,368 : 113 |
| `H-2Aa` | `H2-Aa` | 18 | hyphen placement only |
| `H-2Eb1` | `H2-Eb1` | 7 | hyphen placement only |
| `H2-Ag7` | `H2-IAg7` | 3 | the `I` dropped; `H2-IAg7` dominates 333 : 3 |
| `H2-Ed` | `H2-IEd` | 2 | likewise, 30 : 2 |

`H2-Ab1` (9 records in `mhc.b`) is **not** touched: it is the IMGT *gene* symbol for the I-A beta
chain and carries no haplotype, so mapping it to a molecule would need the paired `mhc.a`. Reported,
not guessed.

### The class-II chain order — a finding, not a listed issue

**149 records carry a beta-chain gene in `mhc.a` and an alpha-chain gene in `mhc.b`** — the pair the
wrong way round, led by `(HLA-DRB1*01:01, HLA-DRA*01:01)` on 48. The gene symbol says which chain it
is, so the correction needs no judgement and is applied. A donor typed on the alpha chain would never
have matched those records.

### Two things still open, for the author

1. **Murine class I is `H2-Db` in the data and `H-2Db` in `proofreading/mhc.md`** — 2,451 `H2-Db`,
   2,334 `H2-Kb`, 1,422 `H2-Kd`, and `H-2Db` appears **zero** times. `vdjdb-web` carries a spelling
   repair for exactly this pair. Two authorities disagree about ~6,200 records and the answer changes
   what a user searches for, so it is not taken here.
2. **7,543 mouse records carry `HLA-DQA1*03:01` / `HLA-DQB1*03:02`.** Those are HLA-transgenic mice
   and the combination is correct; any future "species must match the MHC" check has to allow it.

### Swapping two columns at once

A rename declares one column, so the chain-order fix cannot be one. It could have been two
conditional renames — each testing the other column — which is why conditionals now read **both** the
target and the evidence from the pre-rename snapshot: otherwise the second would test a column the
first had already rewritten and the declaration order would decide the answer. The swap is declared
as a row delta regardless, because it is cleaner to read.

**Verdict: PASS**, 0 unattributed cells, 40 renames.


## 26. Phase 9d result — references, the identical-chain report, and input hygiene

### #347 is mostly not a PMID problem

30,977 records carry a non-PMID `reference.id`. Counted by kind:

| Kind | Distinct | Records | Has a PMID? |
|---|---|---|---|
| 10x Genomics application note | 1 | 20,358 | no — a vendor note |
| `github.com/antigenomics/vdjdb-db/issues/*` | 8 | 4,366 | no — **direct submissions**, where the issue *is* the reference |
| DOIs | 4 | 787 | 3 of 4 |
| preprint URLs (bioRxiv, arXiv) | 2 | 322 | 1 of 2 |
| `rcsb.org/structure/*` | 42 | 42 | no — a PDB entry |
| a TUM thesis | 1 | 3 | no |

So the issue's real scope is **668 records across 3 references**, now resolved:

| Reference | PMID |
|---|---|
| `https://doi.org/10.1016/j.xcrm.2023.101017` | PMID:37030296 |
| `https://www.biorxiv.org/content/10.1101/2025.11.05.686789v1.full` | PMID:41279151 |
| `doi:10.1172/jci.insight.174776` | PMID:39024572 |

The bioRxiv preprint had **acquired a PMID since it was submitted**, which is exactly the drift a
committed table catches. `https://doi.org/10.1101/2020.05.04.20085779` is a medRxiv preprint that was
never indexed, and it is recorded as checked-and-unmapped so nobody looks again.

`proofreading/reference_ids.tsv` is a **committed, reviewed input** — the build is offline and
deterministic, so no lookup happens at build time (hard rule 9). Refreshing it is re-running the
resolver and reviewing the diff.

### #561 reports, it does not repair

**99 records carry the same CDR3 on both chains** — the beta sequence copied into the alpha field
with the V and J calls left correct. 98 of them come from two references (PMID:34811538 with 71,
PMID:41610844 with 27). Which chain is wrong cannot be known from the row, so this is an **advisory
QC rule**: it is reported on every `vdjdb qc` run and does not fail the build, because a defect only
a curator can fix must not block a submission.

### #368 is already fixed, for the half that was mechanical

Zero records have `antigen.gene` holding a species name. The other half of the issue — *"sometimes a
human protein name is used instead of a gene symbol"* — is still visible on 9 values, led by
`Trans-sialidase` (284 records), `Nucleocapsid` (171) and `Neuraminidase` (39), and
`proofreading/gene_aliases.tsv` covers only the last. Choosing a gene symbol for the other eight is
curation, not a mechanical rule, so they are reported rather than invented.

Separately, **15 `antigen.species` values sit outside `proofreading/species_aliases.tsv`**, led by
`SIV` (1,771 records), `RotavirusA` (80) and `Synthetic` (62). The vocabulary is incomplete, not the
data wrong; the file should gain them.

### Input hygiene: whitespace forks a value in two

**758 record-cells across 8 columns carried leading or trailing whitespace** — `tetramer-sort `
beside `tetramer-sort` (103 records), `Nucleocapsid ` (171), `HLA-DRB1*15 ` in `meta.donor.MHC`
(336), and one J-gene call with a **non-breaking space**. The reader stripped only the `\r` a CRLF
file leaves; it now strips all surrounding whitespace, which is never meaningful in a TSV cell.

That is a reader change, so it is global and it removed 52 of the segment renames that had existed
only to undo it.

**Phase 9 verdict: PASS**, 0 unattributed cells, 43 renames, 3 row deltas, 301 tests.


## 27. Phase 9d(ii) result — the epitope catalogue, and the patch's own defects

VDJdb now ships **its own list of epitopes and the MHCs that present them**, as two tidy tables:

| Table | Key | Rows |
|---|---|---|
| `epitopes` | `(antigen.epitope, antigen.species)` | **2,132** |
| `restriction` | `(antigen.epitope, antigen.species, mhc.a, mhc.b)` | **2,373** |

`epitopes` carries the antigen gene, the peptide length, the MHC class, and what supports it —
records, chains, distinct clonotypes and **distinct references**. 379 epitopes are reported by two or
more publications.

**The key is the epitope *and* the species**, and that is the point rather than a detail: 13 epitopes
in the corpus are reported under two organisms and none of them is an error — `VEALYLVCG` is insulin
B in both `HomoSapiens`/`INS` and `MusMusculus`/`Ins2`, `LPRWYFYYL` is shared between HCoV-HKU1 and
HCoV-OC43, `KLPDDFMGC` between SARS-CoV and SARS-CoV-2.
`patches/antigen_epitope_species_gene.dict` is keyed on the peptide alone and **cannot express any of
them**; forcing one species through it silently rewrote 79 HCoV-OC43 records to HCoV-HKU1 on the
first attempt. The catalogue can, which is why it is a table and not a view over the patch.

Each allele is checked against IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>, mirrored at
`proofreading/mhc_alleles.tsv.gz`) by **prefix**, because a VDJdb call is two-field and the authority
stores four: 2,236 pairs `known`, 135 `unchecked` (murine and macaque names, which have no such
authority), **2 `unknown`** — down from 10 before the corrections below.

Two checks the catalogue makes possible immediately, ahead of phase 9e's `mhcmatch` work: **20
epitopes are recorded as MHCI but are longer than 11 residues, and 42 as MHCII but shorter than 12.**
Peptide length and MHC class are independent evidence about each other, and they disagree on 62 of
2,132 epitopes.

### Gene symbols, from the source (#368's other half)

According to PubMed, four lookups settled the protein-name-instead-of-gene-symbol cases:

* influenza A nucleoprotein epitopes are **NP** — "HLA-B\*37:01-restricted NP … and HLA-A\*01:01-restricted NP"
  ([DOI](https://doi.org/10.1038/s41467-018-07815-5));
* influenza A neuraminidase epitopes are **NA** — "These neuraminidase-derived peptides, NA(SGPDNGAVAV)…"
  ([DOI](https://doi.org/10.4049/jimmunol.2000689));
* the LCMV epitope is nucleoprotein ([DOI](https://doi.org/10.3389/fimmu.2023.1199064)), and VDJdb
  already writes `NP` on 477 LCMV records;
* the chicken antigen of the D10 structure is **conalbumin**, i.e. ovotransferrin — its MeSH term
  ([DOI](https://doi.org/10.1126/science.286.5446.1913)). One record, no unambiguous short symbol, so
  it is reported rather than renamed.

Separately, the coronavirus genes are normalised to the NCBI RefSeq (NC_045512.2) symbols —
`Spike` → `S`, `Nucleocapsid` → `N`, `Matrix` → `M`, `Envelope` → `E`, `ORF3` → `ORF3a` — **262 patch
entries covering 10,563 records**. That was not a tidy-up: `YLQPRTFLL` was labelled both `S` and
`Spike` on 2,398 records, `LLLDRLNQL` both `N` and `Nucleocapsid` on 1,750, and **39 epitopes carried
two gene symbols under one species**, so a query on either name returned half the data.

### The patch file had four dead entries and eight contradictions

Repairing it was not optional bookkeeping:

* **Four entries used a space where the tab belongs**, so the whole mangled string became the epitope
  key and the entry had **never matched a single record**. One of them is `ALAGIGILTV` — MLANA, among
  the most-studied human epitopes, 180 records — which is why `MART1` survived in the data.
* **The file has no trailing newline**, so the first appended entry was glued onto the last existing
  one. Caught by a parse error; a test now asserts both patch files are well-formed.
* **Eight epitopes are listed twice with different answers**, and `keep="last"` lets file order
  decide, silently. Two are pure nomenclature (`EBNA3B` = `EBNA4`, `EBNA3C` = `EBNA6`); three sit in
  the HIV-1 Gag-Pol frameshift, where both readings are defensible; the rest need a curator. They are
  named in `curate.patch.CONFLICTING_EPITOPES` and asserted, so the set can shrink but not grow.

### MHC corrections became data (`patches/mhc.dict`)

Patch semantics, with a `reference.id` scope so a correction verified against one publication does not
leak to another:

| From | To | Records | Source |
|---|---|---|---|
| `HLA-A*24:01` | `HLA-A*24:02` | 80 | no `A*24:01` in IPD-IMGT/HLA (0 rows against 342) |
| `HLA-A*24:09/:11/:12/:16` | `HLA-A*24:02` | 8 | the paper is *titled* "HLA A\*24:02-restricted T cell receptors…" ([DOI](https://doi.org/10.1172/JCI164535)) |
| the five murine class-II spellings | — | 141 | `proofreading/mhc.md`, and the dominant spelling |

**Checked and deliberately not corrected**, recorded as comments in the file: `HLA-A*08:01` on 74
records — there is no HLA-A\*08 locus at any resolution, so the call is certainly wrong, but the chunk
carries **no `reference.id`** to resolve it against and no source states the intended allele; and
`HLA-B*12` on 1 record, a serological antigen that splits into B44 and B45, which the paper's abstract
does not pin ([DOI](https://doi.org/10.1128/JVI.73.3.2099-2108.1999)).

### The ledger needed an exact-match predicate

`antigen.gene` and `antigen.species` are part of `vdjdb.slim.txt`'s **grouping** key, so correcting a
gene symbol splits and merges slim rows — 13,290 of them undeclared. Two things fixed that:

1. **Renames are derived from the reference bundle, not from `chunks/`.** The release was built with
   the patch as it stood then, so only the entries added since reach the ledger; generating against
   the chunks declared hundreds of corrections the release already carried, and all 72 were stale.
2. **`when_equals`, not `when_contains`.** A 9-mer epitope really is a substring of a 10-mer one
   (`SPRWYFYYL` inside `LSPRWYFYYL`), so a substring predicate fires on the wrong rows.

Slim's undeclared row delta fell from 13,290 to **744**. **Verdict: PASS**, 0 unattributed cells,
323 renames.


## 28. The arda re-measurement, after #327 and #389

§16 measured `arda.cdr3fix` against the vendored k-mer scanner on a corpus whose allele calls were
wrong, and the swap was held for that reason. Re-measured 2026-09-25 on the corrected corpus, both
engines run through the same pipeline on all 192,753 records.

### Agreement rose sharply

| Field | §16, before | After #327 + #389 |
|---|---|---|
| `cdr3` (alpha) | 97.74 % | **99.91 %** |
| `jFixType` | 96.09 % | 95.66 % (alpha) |
| `vFixType` | 93.69 % | 96.47 % (alpha) |
| `good` | 91.17 % | 93.41 % (beta) |

The repaired sequence itself is now **99.91 %** identical between the two engines on alpha chains.
§16's headline "arda changes 14,778 repaired CDR3s" was largely arda and the legacy disagreeing about
sequences the legacy could not place at all.

### Coverage: the alpha side closed, the beta V side did not

| Chain | Field | both | legacy only | arda only | neither |
|---|---|---|---|---|---|
| alpha | `vEnd` | 118,217 | **370** | 3,776 | 567 |
| alpha | `jStart` | 116,982 | 137 | **5,164** | 647 |
| beta | `vEnd` | 152,139 | **7,969** | 1,629 | 1,346 |
| beta | `jStart` | 160,682 | 138 | **1,904** | 359 |

arda now **gains 7,068 J-start mappings and 5,405 V-end mappings**, and loses 275 J-starts. The alpha
V-end loss that §16 reported has essentially closed — 370 chains. **The beta V-end loss has not: 7,969
chains**, and that is the whole of the remaining objection.

### What the 7,969 are

Not a nomenclature problem, and not spread thin: **7,962 of 7,969 are `FailedBadSegment`** — arda
declines the segment outright — while the legacy reports `NoFixNeeded` on 7,919 of them, i.e. it
placed the CDR3 against the germline without changing a residue. They concentrate on common V genes
(`TRBV20-1*01` 1,673, `TRBV3-1*01` 1,058, `TRBV6-1*01` 1,004, `TRBV11-2*01` 823) and on the species
split 6,553 human / 826 mouse / 590 macaque.

The tell is arda's own V call on those rows: `TRBV20`, `TRBV3`, `TRBV6`, `TRBV12-2+TRBV13-2` — **gene
names and ambiguity groups, not alleles**. arda is collapsing the allele and then failing to find an
anchor for the collapsed name. So this is an **arda reference-loading problem on TRBV, not a
disagreement about the sequence** — the same class of defect as `fix/source-root-marker`, which also
presented as a silent "no answer".

### Recommendation

**Do not swap yet, and do not re-litigate it against this corpus either.** The gate §16 named has
been cleared — the allele calls are right, agreement is 99.91 %, and arda is ahead on three of the
four coverage measures. What remains is one upstream defect with a clear signature, and the right
next step is to reproduce those 7,962 `FailedBadSegment` beta V calls in `arda` directly and fix them
there. Then the swap is a one-line default change with nothing left to weigh.
