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
| 4 | `feature/pipeline-core` | harmonize + score + pairing + the three legacy emitters; deletes `py_src/` | #424, #399 | harness shows only the 7,998 `unmp` rows + the meta fixes; peak RSS < 8 GB |
| 5 | `feature/arda-cdr3fix` | `arda.cdr3fix` replaces `Cdr3Fixer.py`; retires `res/segments*.txt` | — | new ledger rule, row count measured then frozen |
| 6 | `feature/new-format` | `emit/vdjdb3.py`, generated meta, `vdjdb.schema.json`, `make legacy` | — | `make legacy` from the new build still passes the harness |
| 7 | `feature/airr` | `emit/airr.py`, `convert/coords.py`, both converters | — | the `legacy_to_airr ≡ vdjdb3_to_airr` round-trip property test |
| 8 | `feature/junction-nt`, `feature/segment-guess`, `feature/dgene` | one branch each | #461, #462, #463 | generated `cdr3nt` back-translates to `cdr3` |
| 9 | `feature/harmonize-rules` | nomenclature rule tables | #327, #389, #347, #368, #564, #467, #561 | each rule gets a ledger entry with a measured row count |
| 10 | `feature/motifs-tcrnet` | TCRNET on `vdjtools`, streaming backgrounds | — | deviation report accepted |
| 11 | `feature/motifs-tcremp` | TCREMP + per-epitope DBSCAN; new motif schema; legacy projections | — | beats the **shipped** `cluster_members_tcremp.txt` re-scored in our harness, per §8.4 |
| 12 | `feature/summary` | Rmd split, ggplot2 4.x fixes, PubMed cache, data-driven callouts, interactive dashboard | #460 | renders offline; perceptual + structural checks pass |
| 13 | `feature/docs` | Sphinx site, generated schema tables, dashboard tab, Pages | — | zero-warning build, deploys |
| 14 | `feature/release-tooling` | manifest, three zips, checksums, `latest-version.txt`, tag scheme, Zenodo, changelog; retires the legacy CI | #432 | full release dry-run with a clean ledger |
| 15 | `feature/aldan3-runner` | self-hosted runner + `build.yml` retargeting | — | identical canonical digests on both runners |

**Phase 2 is load-bearing.** The harness must show zero diffs against the *current* build before any
behaviour changes, or there is no instrument to attribute later differences with.

Phases 5, 8, 9, 10 and 11 each introduce exactly one source of deviation, so every difference in the
output has a single attributable cause.

## 5. The difference ledger

`vdjdb diff <reference-zip> <candidate-dir>` compares in three passes: file set → two digests per file
(raw and canonical) → row-level classification keyed on
`gene|cdr3|v.segm|j.segm|species|mhc.a|mhc.b|antigen.epitope|reference.id`.

Every changed **cell** must be attributed to a declared rule in `rules/expected_diffs.toml`:

```toml
[[rule]] id="web-unmp-jstart-minus1" file="vdjdb.txt" column="web.cdr3fix.unmp"
         from="no" to="yes" predicate="cdr3fix.jStart == -1" rows=7998
```

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
| After global dedup instead | 192,734 (19 cross-chunk duplicates the per-chunk pass keeps) | bears on #390 |
| polars read + dedup of all 230 chunks | **0.4 s** | vs a pipeline documented as needing 64 GB |
| Chunks passing `ChunkQC` | **230 / 230**, zero errors | fail-fast needs no quarantine list |
| Non-empty CDR3 cells | 305,031, **zero** with characters outside the 20 AAs | TCREMP pre-filter is a guard, not a live problem |
| `arda.cdr3fix` vs shipped `cdr3fix`, 20,000-row sample (seed 42) | `cdr3` 97.81 % · `vFixType` 97.95 % · `vEnd` 96.53 % · `good` 95.68 % · `jFixType` 95.44 % · `jStart` 91.02 % | arda 2.27.0 |
| …of the 1,797 `jStart` disagreements | 556 VDJdb-unmapped → arda-mapped · 1,241 both mapped, **arda smaller in every case** (mode −2, range −1…−8) · **0** coverage regressions | NW with free end gaps vs k-mer longest-hit |
| `arda.cdr3fix.markup_batch` throughput | 27 µs/record → **5.2 s** for 191,440 distinct keys | not vectorised internally, but fast enough that it does not matter |
| `vdjtools.model.infer_nt` | **3.11 ms/record** → ~15 min for 284,546 rows | human TRB, warm, single-threaded |
| `vdjdb.txt` row/field shape | 284,546 rows, **all exactly 22 fields**, zero empty `cdr3fix` | keep and assert; raggedness is not a live problem |
| `vdjdb_full.txt` `cdr3fix.alpha` encoding | 122,930 non-empty cells, **100 % Python dict repr**, 0 JSON | live bug |
| `web.cdr3fix.unmp` correctness | `(no,no)` 268,546 · `(yes,yes)` 8,002 · **`(no,yes)` 7,998** | `jStart` is only ever −1 (8,536) or > 0 (276,010); **0 never occurs** |
| Murine MHC-II spellings | `I-Ab` **3,274** · `H2-IAb` 113 · `H2-Ab1` 9 · `H2-IAg7` 333 vs `H2-Ag7` 3 · `H2-Aa` 25 vs `H-2Aa` 19 · `H2-Eb1` 7 vs `H-2Eb1` 7 | mostly in `mhc.b`; class I is clean |
| #327 TRAJ24 | bare `TRAJ24` 1,302 (978 carry `WGKLQF`) · `TRAJ24*01` 113 (**75 = 66 % carry `WGKLQF`**) · `TRAJ24*02` 34 · malformed `TRAJ24-1` 1 | reproduces the report on a larger set |
| #368 antigen gene/species | patch dict covers 245 epitopes = 154,373 of 203,308 rows, **zero swapped** | the ~90 reported records are in the uncovered tail of 48,935 rows |
| Dashboard R deps | 15 of 17 installed; `maps`/`scatterpie` used only past the embed cut (lines 812–835 vs marker at 567); `ggh4x` never used | splitting the Rmd drops three deps from the release path |
| Shipped dashboard PNG sizes | 1344×960, 2304×1920, 1152×1920, 1536×1536 → `dpi=96, fig.retina=2` | pin it or the visual fingerprint is noise |
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
20, reservoir-samples uniformly over unique clonotypes, and content-addresses the cache — do not
reimplement it. But **never let `tcrnet()` resolve its own background**: it calls
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

Freeze the `StandardScaler` + `PCA(50)` as a version-pinned artifact and cache the 50-D vectors by
`(cdr3, v, j)` — 113k × 50 × float32 = 22 MB, well inside the Actions cache limit, where the raw
6000-D matrix at 2.7 GB is not. That also makes cluster ids stable across releases, which today's raw
igraph component numbers are not (every bookmarked vdjdb.com motif URL breaks each release).

### 8.9 `mir.bench` helpers are not usable here

- `mir.bench.vdjdb.load_vdjdb` keeps only `junction_aa, v_call, j_call, locus, epitope, mhc_class` —
  it drops `species`, `mhc.a`, `mhc.b`, `antigen.gene`, `antigen.species`, `v.end`, `j.start`, seven
  of the 19 legacy columns. Read the slim table directly.
- `mir.bench.metrics.cluster_metrics` uses `recall = tp / n_true_clustered` (clustered records only);
  the vendored `metrics_lib.precision_recall_fscore` folds unclustered records into FN. **The two
  metric families are not comparable.** Validation uses the vendored `metrics_lib`, verbatim.

## 9. Open questions

- **Where do the five `evidence.*` columns come from?** Nothing in this repo produces them, yet
  production serves a 27-column `vdjdb.txt` containing them. Three come from
  `vdjdb-web/tools/reconcile_structures.py`; the two `evidence.validation.*` have no located
  generator. Decide whether the new format takes ownership. *(blocks phase 6)*
- **The seven side outputs** (`vdjdb_full_filtered.txt`, three `*_broken.txt`, three `*_scored.txt`)
  are written but never shipped. Recommendation: `build/reports/*`. *(blocks phase 6)*
- **The five discarded chunk columns** — promote, or add to a named `TOLERATED_DROPPED` set so the
  drop is deliberate and testable. `submitter` is attribution data. *(blocks phase 1)*
