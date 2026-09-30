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

**State, checked 2026-09-29: unmet.** `_zip_asset` is still the four lines above. `manifest.json`
ships a `role` per bundle from this side, so step 1 has everything it needs to key on, and the
release dry-run produces the three zips - but until the patch and the release land, a VDJdb release
carrying more than one zip gives vdjmatch an arbitrary one. This gate belongs to the author; nothing
in this repository closes it.

### 3.2 Phase −1: unblocking fixes

None of the three is part of the rewrite.

| Fix | Blocks |
|---|---|
| The `vdjmatch` patch above | "three zips per release" |
| Correct `latest-version.txt` line 1; delete the stale `database/` copy | every client that still reads it. **Done**: line 1 names `2026-06-03/vdjdb-2026-06-03.zip`, which is the published latest release's asset, and `database/latest-version.txt` is untracked - only `dummy` and the two `*.meta.txt` files are. `release.yml` writes line 1 before the build and commits it after the release exists, then fetches it and fails on anything but 200 |
| Replace the seven `.T.apply` calls and the `ScoreFactory` `iterrows` | "runs on 16 GB GitHub-hosted runners" |

### 3.3 Backgrounds

Backgrounds are streamed at build time and only derived statistics (`count.bg`, `total.bg`) ship:
never the background itself, never a subsampled copy. A background is an input to the build, not
an output of it, and a stale vendored copy would change a call set with no error.

### 3.4 Cross-repo gate: `arda` CDR3 repair release - **closed 2026-09-30**

`arda.cdr3fix` declined a class of repair the `Cdr3Fixer` it replaced applied, so a malformed junction
shipped. **arda 2.36.0 and vdjtools 4.8.0 close it**, and the bound here moved with them: the
annotation is one `vdjtools.model.annotate_junctions` call, `annotate/dgene.py` and
`annotate/segments.py` are deleted, and #711 is closed.

What the gate was about, and what it cost: `Anchor.templated_aa` runs Cys104 through [FW]118 inclusive
and stops, so framework past an anchor had nothing to align to. arda's bundled IMGT build already
carried both flanks as named columns in `markup.aa.tsv` - `fwr3` ends at Cys104 on 767 of 775 distinct
human `v_call`, `fwr4` starts at the 118 anchor on 124 distinct `j_call` - and the repair path read
`cdr3_anchors.tsv`, which does not. **No reference data was lost when `res/` was retired.**

Measured on the shipped corpus, 285,794 chains, every junction classified against the germline of the
segment it names:

| | before | after |
|---|--:|--:|
| shipped junctions not canonical against their own germline | 515 | **451** |
| of those, framework kept past an anchor (`under-trimmed`) | 77 | **33** |
| `CASSQSPGGVAFFGQG` shipped | 2 chains | **0** |

The 451 that remain are a proofreading queue in `out/reports/anchors.tsv`, not a repair failure: the
junction disagrees with the germline of the segment the record names, and only a curator can say which
of the two is wrong.

**Two `arda` bottlenecks this build measured are also settled.** `antigenomics/arda#143`
(`markup_records`' `_align` is a pure-Python dynamic program, 35.4 µs/key, 23,556,016 `max()` calls
over 190,624 records) is open and does not gate anything. `antigenomics/arda#142` asked arda for a
batched `posterior_d`; it **dissolved** - see §3.5.

### 3.5 Cross-repo gate: `d.posterior` moves from `arda` to `vdjtools` - **closed 2026-09-30**

The D call and the number beside it used to come from two different models: `d.inferred`, `d.start` and
`d.end` from `vdjtools.model.infer_nt_batch`, and `d.posterior` from `arda.dpost.posterior_d`, which
named a different D gene from the one it annotated on 21.5 % of chains. It was also a
recombination-model computation shipped as a fitted `d_prior.tsv` inside arda's *germline reference*
tree, covering four `(organism, locus)` pairs against a reference covering five organisms.

**`arda.dpost` is gone in arda 2.33.0 and does not reappear in vdjtools.** `d.posterior` now comes from
the same recombination scenario weights that named the call, so the number beside a call is the
probability of that call. `annotate/dgene.py` is deleted rather than re-pointed, and `d.entropy` is
retired: it was the entropy of the second estimator's distribution, and there is no second distribution
any more.

There is no "arda model" that could have kept it instead. arda is the **aligner** and the germline
namespace; every bundled recombination model is `vdjtools`' own fit - `source="olga"` the OLGA
bootstrap, `source="arda"` an EM fit over real non-functional reads in arda's IMGT allele namespace,
which is the only bundled set covering mouse. So the model side is wholly `vdjtools`'.

| | before | after |
|---|--:|--:|
| `add_d_posterior` + `add_junction_nt` | 27.04 s, two stages | **21.43 s, one** |
| `vdjdb build` | 41.50 s | **33.30 s** |
| D gene correct, human TRB against nucleotide truth | 69.67 % | **74.35 %** |
| rows carrying D coordinates | 55.75 % | **99.80 %** |

⚠ Per-row positional precision is the one thing that got worse: `d.start` is exact on 64.71 % of
correctly-called rows against 66.67 % under the retired E-value gate. It is exact on **1,922 rows
rather than 1,278**, because it answers 3,992 rather than 2,230.

`antigenomics/vdjtools#183` and `antigenomics/arda#144` carry the relocation reasoning and
`antigenomics/arda#142` dissolved into it, as predicted: the batched path `#142` asked arda to build
already existed here.

## 4. Phases

`master` → `dev` → `feature/*` → `dev` → `master`. Every phase is independently mergeable and
`dev` stays green. Commits that resolve a tracker issue carry `Closes #N`.

| # | State | Branch | Delivers | Closes | Acceptance |
|---|---|---|---|---|---|
| 0 | merged | `feature/dev-baseline` | `ROADMAP.md`, `pyproject.toml` + `uv.lock`, package skeleton, `chunk-check.yml`, `branch-policy.yml` | #476 | CI green; `release.sh` still works untouched |
| 1 | merged | `feature/schema` | the field registry; `render_meta` | - | reproduces the Groovy `METADATA_LINES` / `SLIM_METADATA_LINES` byte-for-byte; `header == meta names` for all three tables |
| 2 | merged | `feature/golden-harness` | `vdjdb diff` + `expected_diffs.toml` | - | zero diffs against the current pandas build; nothing downstream starts without this |
| 3 | part | `feature/io-qc` | polars reader, vectorised QC, `--strict` exit-1, chunk header normalisation | - | QC report matches the pandas report row-for-row; harness still zero. ⚠ **The `.tsv` rename did not happen and #497 is open**: `chunks/` is 231 files, all `.txt`. Everything the rename was wanted for did land - one canonical 33-column header, the 19 distinct header rows collapsed to one, the 99 CRLF files converted with `*.txt text eol=lf` in `.gitattributes` so it cannot return (#581), and a reader that fails on an unrecognised header. What is left is the extension, which nothing reads to decide the format, against 231 `git mv`s that break every `git log --follow` boundary and every `PMID_<id>.txt` reference in docs, skills, tests and the tracker. Held deliberately: renaming every file in `chunks/` the week curators start opening chunk pull requests is when it costs most. It wants its own branch under the mechanical-repair rule and a quiet period |
| 4 | merged | `feature/pipeline-core` | the definitive tables (`records`, `chains`) + harmonize + score + pairing; the legacy export as a projection of them; deletes `py_src/` | #424, #399 | every difference against the release is a declared rule firing its measured count; peak RSS < 8 GB |
| 5 | done | `feature/arda-cdr3fix`, `feature/retire-res` | `arda.cdr3fix` replaces `Cdr3Fixer.py`; `res/` retired | #658 | new `expected_diffs.toml` rule, row count measured then frozen. `res/segments*.txt` outlived the first branch by three call sites, one of them on every build: `arda.cdr3fix` repairs a junction against a *named* germline and never proposes one, so a blank V or J needed filling first. `feature/retire-res` replaced that with the recombination model falling back to arda's own germline anchor table, deleted `res/` and `annotate/_legacy_fixer/`, and dropped `--engine legacy`. The J proposal gains 347 calls and 469 rows of `vdjdb.txt`; the V proposal is reported as `v.inferred` and does not ship, because a V recovered from a junction alone is right 23.8-50.1 % of the time against a J's 93.6-97.5 %. ⚠ **A junction carrying framework past an anchor is no longer repaired, and that is a live regression** (#711, `antigenomics/arda#141`): `arda.cdr3fix` reads `cdr3_anchors.tsv`, whose `Anchor.templated_aa` runs Cys104 through [FW]118 inclusive and stops, so framework past an anchor has nothing to align it to and the defect is reported without the repair being applied. **No reference data was lost with `res/`** - arda's bundled IMGT build carries both flanks as named columns in `markup.aa.tsv` (`fwr3` ends at Cys104, `fwr4` starts at 118), which is the table arda's own repair path does not read (§3.4). Measured over 190,902 distinct keys: **382** ship a junction that does not run Cys104 to its own segment's anchor, and **2,165** are returned `NoFixNeeded` on a single residue of germline agreement where the retired fixer required a 2-mer hit at offset zero in both sequences. The fix is arda's - it already ships `alleles.fasta` and `anchor_nt` - and reaches this build as a release, per §3.4 |
| 6 | merged | `feature/new-format` | ships the definitive tables as parquet + TSV, adds `evidence`, `vdjdb.schema.json` | - | `make legacy` from the shipped tables still passes the harness |
| 7 | merged | `feature/airr` | `emit/airr.py` (Rearrangement + Reactivity), `convert/coords.py`, `vdjdb convert` | - | `airr.validate_rearrangement` passes on the full table; the legacy path produces nothing the tables path does not |
| 8 | merged | `feature/junction-nt`, `feature/segment-guess`, `feature/dgene` | one branch each | #461, #462, #463 | generated `cdr3nt` back-translates to `cdr3`. The stage was 87.2 % of assembly on `vdjtools` 3.13 and ran as four worker processes over contiguous slices; 4.5 published `infer_nt_batch` (`antigenomics/vdjtools#181`) and it is now one batched call per (species, locus), **114.93 s → 12.44 s**, #656 |
| 9 | merged | `feature/harmonize-rules` | nomenclature rule tables | #327, #389, #347, #368, #564, #467, #561 | each rule gets an `expected_diffs.toml` entry with a measured row count |
| 10 | merged | `feature/motifs-tcrnet` | TCRNET on `vdjtools`, streaming backgrounds | - | deviation report accepted |
| 11 | merged | `feature/motifs-tcremp` | TCREMP + per-epitope DBSCAN; new motif schema; legacy projections | - | beats the shipped `cluster_members_tcremp.txt` re-scored in our harness, per §8.4 |
| 12 | merged | `feature/summary` | Rmd split, ggplot2 4.x fixes, committed publication-year table, data-driven callouts, interactive dashboard | #460 | renders offline; perceptual + structural checks pass |
| 13 | merged | `feature/docs` | Sphinx site, generated schema tables, dashboard tab, Pages | - | zero-warning build, deploys |
| 14 | merged | `feature/release-tooling` | manifest, three zips, checksums, `latest-version.txt`, tag scheme, changelog; retires the legacy CI | #432 | full release dry-run with no unattributed differences |
| 15 | part | `feature/aldan3-runner` | self-hosted runner + `build.yml` retargeting | - | identical canonical digests on both runners. `build.yml` carries the `fromJSON(inputs.runner)` retargeting; **no self-hosted runner is registered** (`actions/runners` returns 0), so the second half of the criterion is unmet |
| 16 | merged | `feature/identity` | the four derived id levels, the lifecycle record, `vdjdb identity`, promiscuity columns, one study count | - | every invariant of §10.5 passes; a permuted chunk order changes no id; the dashboard reports every reference on the row. All three hold, and `record_id` is stable across rebuilds since the registry became a committed input (#674). Step 9, the promiscuity columns, landed 2026-09-29 and completes the phase. It had waited on 9e, and the two never collided once the split was stated: `curate/presentation.py` reads `mhcmatch`'s bundled pseudosequences inside the build, while the model-based scoring stays in `vdjdb promiscuity` and reaches the build only as a join against the committed `proofreading/epitope_promiscuity.tsv` (13,510 rows over 1,729 epitopes and 107 alleles). `restriction` is 15 columns, and 1,793 of its 1,989 class I rows carry a rank |
| 17 | merged | `feature/corpus` | the reference corpus: documents, vocabulary, postings, `score` and `lift` | - | the three files reproducible by digest; `score` reproduces the `refsearch` ranking; `lift` answers a specificity question with an n |

Phases 0 to 14 are merged to `master` as of 2026-09-27, and phase 15 is half landed: the comparison
against the last release reads PASS with every difference declared and measured, and the release dry
run produces three reproducible bundles. `ROADMAP_local.md` carries the per-phase record.

**Promoted to `master` on 2026-09-29, 40 commits, full build green in 12 m 40 s** (#692). It carried
the interactive dashboard, the junction-anchor check, the junction-nt batch call, the profile fix, and
the comparison of the shipped bundle rather than of `out/legacy`. That comparison runs on the assembled
legacy zip over all twelve of its members - it named five with `--only` until then - and reads PASS.
Two declarations make that possible and are part of the release contract: `[members]` for a change of
bundle shape (today the two TCREMP motif tables) and `[measured_elsewhere]` for a member another
instrument gates, with the instrument named. §5 has both.

On `master` with that promotion - #658, #672, #675, #647, #671, #648:

| | What | Measured |
|---|---|--:|
| #658 | `res/` retired; arda's IMGT reference replaces a 2023 import's by-product | J proposal +347 calls |
| #672 | `registry/records.tsv` committed, refreshed and gated | 592 amendments were stale |
| #675 | ten input files given a final newline, two rows padded to 33 fields | 12 files |
| #646 | junctions repaired against their own germline anchor, four commits | **4,838** chains, `anchors.tsv` 5,959 -> 1,132 |
| #647 | the mouse TRAJ47 allele read off the junction | 95 chains, `J allele mismatch` 95 -> **0** |
| #671 | a `[[rename]]` scoped to its organism | 79 mouse rows un-mis-keyed |
| #648 | the spectratype grouped on 500 buckets, not 164,131 CDR3s | dashboard 152 s -> **22.8 s** |

**`dev` ahead of `master` again, four merges as of 2026-09-29.** The build now runs when only
`registry/` changes, which it did not (#691); `out/reports/epitope-sources.tsv` names every peptide
whose source is not single-valued and `RGPGRAFVTI` is re-attributed from `P18-I10` to HIV-1 `GP160`
(#633); and two commits on the B16 chunk (#397), which is the first data in the corpus to record a
clonotype count:

| | What | Measured |
|---|---|--:|
| #397 | 1,556 rows for 485 clonotypes collapse to 485, `method.frequency` as `count/sample total` | records 192,763 -> **192,641** |
| #694 | a spreadsheet autofill had turned the gene `Eef2` into the run `Eef2`..`Eef188` | 187 rows, 65 clonotypes, **122 records that were never real** |
| #397 | `meta.epitope.id` was the constant `p12`, and `meta.subject.cohort` held a culture serial | 459 and 485 rows |
| #397 | `antigen.gene` `Plod1` -> `Plod2`, on `mhcmatch`'s mouse proteome | 49 rows |

QC `duplicate` falls **10,585 -> 9,636** across those, the largest single move the advisory baseline
has recorded.

#646 is not closed by them: 261 chains still carry a germline-supported repair the build proposes and
does not apply, and the reasons are per-case.

**Filed from this work and open:** #693 (the registry retires an id that changed two key fields with
no `replaced_by` forward pointer - the B16 `Plod2` repair produced 49 of them), #696
(`method.frequency` conflates a count with a ratio, and the confidence score needs the count), and
#632 (`antigen.gene` has no authority). #685 is closed and #637's `reference.id` half landed with it;
what is left of #637 is whether the declared `method.identification` vocabulary is enforced at QC
time, which is a decision rather than a measurement.

**Dependency state.** `arda-mapper >= 2.31`, `vdjtools >= 4.7`. The 4.5 bump closed the junction-nt
bottleneck; 4.7 and arda 2.31 carry the two germline-boundary defects this build raised upstream
(`arda#135` `TruncatedGermline`, and `germline_boundary`), and nothing else here reads a 4.x-only API.

Phase 2 came first: the harness had to show zero diffs against the then-current build before any
behaviour changed, so that later differences could be attributed.

Phases 5, 8, 9, 10 and 11 each introduce exactly one source of deviation, so every difference in the
output has a single attributable cause.

Every phase has a step-by-step subplan in §12. A phase is not startable until its subplan names
the files it creates, the facts it needs (already measured, in §7/§8), and the check that closes it.

## 4a. Issue tracker composition

Re-measured 2026-09-29, second pass, with `gh`: 466 issues, **116 open**, against 458 / 123 earlier
the same day and 440 / 130 on 2026-09-25. Grouped by label, one category per issue, intake winning a
tie and maintenance winning over proofreading:

| Category | Open, 2026-09-25 | Open, 2026-09-29 | Open, second pass | What they are |
|---|---:|---:|---:|---|
| data intake | 103 (79 %) | 101 (82 %) | **101 (87 %)** | pending papers, preprints, paper-pending, meta-papers, 10x/Immudex sets, associations, other databases, correspondence |
| curation quality | 22 | 8 | 7 | formatting & proofreading, typos, structural, validation |
| build infrastructure | 13 | 14 | 8 | the build, the summary, maintenance |

The intake row has not moved in eight days. The other two fell because the build work closed what it
had filed: #685, #637's measurable half, #633's epitope-source report, #647, #671, #672, #675, #658.

The build work also files issues from its own measurements - #650 (the profile double-count), #652
(the shipped zips were never compared), #656 (the junction-nt bottleneck), #693, #696 - so the bottom
two rows move in both directions while the intake row stays where it is. The composition is stable
under everything the build does, and that is what this section records.

Six open issues out of seven are a submission queue, not a defect list. This migration closes
issues from the bottom two rows only, and nothing it does shortens the first
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

### What the comparison runs on, and what it is not the instrument for

It runs on the **assembled legacy zip**, so the thing compared is the file a consumer downloads. It
compared `out/legacy` with five members named by `--only` until 2026-09-29, which left seven of the
bundle's twelve members unread; running it on the zip with nothing declared gave **101,877
unattributed cells**, every one in a motif file or in `latest-version.txt` while the five legacy
tables were fully attributed.

Two declarations replace the `--only`, and both are part of the release contract:

```toml
[members]
added = ["cluster_members_tcremp.txt", "motif_pwms_tcremp.txt"]

[measured_elsewhere]
"cluster_members.txt" = "vdjdb motif-metrics, 18 axes; column order by test_reference_contract.py"
```

`[members]` declares a change of bundle shape - an undeclared new or missing member still fails.
`[measured_elsewhere]` names, per member, the instrument that gates it instead; such a member is still
read, digested, row-counted and printed, and is exempt only from cell attribution and the row-bucket
declaration. Four members are listed: the two motif tables, because `cid` carries a position in a
sorted list and one renumbered cluster relabels every cluster after it (32,703 of 55,636 reference
rows key to nothing); `vdjdb_summary_embed.html`, gated by `summary/check_summary.py` on three layers;
and `LICENSE`, which ships verbatim from the repository root.

What replaces the row comparison for the clustering is a per-record measure:
`partition_neighbours_preserved` asks, of every clonotype the last release clustered, what fraction of
its cluster-mates this build still gives it - **0.9224 on TRB and 0.7024 on TRA**. It is per record
rather than per pair because the released TRB clustering puts 19,971 of its 36,906 clonotypes in one
cluster holding 94.7 % of the file's co-clustered pairs, so a pair-weighted score measures that one
blob: the do-nothing partition reads 0.9991 on it (`docs/clustering.md` §8.0).

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
| `vdjtools.model.infer_nt`, per row | 3.11 ms/record on `vdjtools` 3.13, 1.115 ms on 4.5 | human TRB, warm, single-threaded |
| `vdjtools.model.infer_nt_batch` | **0.106 ms/key, 10.5x the per-row loop, and all 3,000 nucleotide sequences identical** | 3,000 distinct human TRB keys, 16 cores, `vdjtools` 4.5. End to end: `add_junction_nt` 114.93 s → 12.44 s, the assemble stage 148.08 s → 37.96 s, `vdjdb build` 45.9 s wall |
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
one of those two can be argued about. This is the accounting that says which. Re-derived 2026-09-29
between tag `2026-06-03-ZENODO` and `dev`; it read *no curation change at all* on 2026-09-27, and
that is no longer true.

**230 chunk files before, 231 after. One added.** `chunks/PMID_18025130.txt`, the KK10 and KK10-L6M
repertoires, recovered from a 2016 branch: 40 records (#161).

**98 files differ once line endings are normalised, and two of them changed row count:**

| File | Rows before | Rows after | Why |
|---|--:|--:|---|
| `Mice_TCRs_..._(Shagina_et_al_2024).txt` | 1,556 | 485 | one row per observation collapsed to one per clonotype, with the count recorded in `method.frequency` (#397) |
| `PMID_18025130.txt` | - | 40 | the file that was added |

The other 96 changed cells and not rows. Fourteen commits produced all 98, one per reason under the
mechanical-repair rule, and the union of the files they touch is exactly 98 - so every file that
differs is accounted for and none differs for a reason nobody wrote down:

| Change | Files | What moved | Issue |
|---|--:|---|---|
| junctions repaired against their own germline anchor, four commits | 78 | 4,838 chains | #646 |
| `meta.structure.id` made a live TCR:pMHC entry or blank | 7 | 2,223 records | #402 |
| final newline added, two `paley-etal` rows padded to 33 fields | 10 | 2 rows | #675 |
| the alpha CDR3 was a copy of the beta: cleared, V and J kept | 3 | 98 records | #561 |
| two class I epitopes recorded under the wrong HLA gene | 2 | 2 epitopes | #597 |
| the space closed in `PMID: 34433824` | 2 | 22 records | #637 |
| `PMID:9971792` Kabat V-beta names converted to IMGT, two J calls corrected | 1 | 47 records | #389 |
| `PMID_15753288` Arden names converted to IMGT | 1 | 35 calls | #302 |
| the two chunks in the row-count table above | 2 | 1,556 -> 485, and 40 landing | #397, #161 |

The 99 CRLF files (#581) are not in the 98 and cannot be: the comparison strips `\r` before hashing,
which is the property that makes a line-ending rewrite invisible to it and a data edit visible.

**Records: 192,753 at the release, 192,641 on `dev`.** Two moves, in opposite directions and both
declared: +40 from `PMID_18025130` landing, and -152 from the B16 collapse, of which 122 are records
that were never real (a spreadsheet autofill had kept 65 clonotypes apart) and 30 are repeated
observations of one clonotype in one mouse.

So a difference the release comparison reports is no longer attributable to the build by default, and
`rules/expected_diffs.toml` says which is which: the curation changes above each carry their own
narrowed rule, declared ahead of the unnarrowed pandas-coercion rules that would otherwise absorb
them, and the three `[[row_delta]]` blocks carry the two chunks that moved rows.

Re-derive it with the two properties that matter, rather than by reading a diff stat:

```bash
# files added or removed
diff <(git ls-tree -r --name-only 2026-06-03-ZENODO chunks/) \
     <(git ls-tree -r --name-only HEAD chunks/)
# per file: does the content differ once line endings are normalised?
for f in $(git ls-tree -r --name-only HEAD chunks/); do
  a=$(git show "2026-06-03-ZENODO:$f" 2>/dev/null | tr -d '\r' | shasum -a 256 | cut -c1-16)
  b=$(git show "HEAD:$f"              | tr -d '\r' | shasum -a 256 | cut -c1-16)
  [ "$a" = "$b" ] || echo "content differs: $f"
done
```

`CLAUDE.md` carries the rule this section exists to serve: a commit touching `chunks/` says which
files, how many rows, why, and who decided. `git log --oneline 2026-06-03-ZENODO..HEAD -- chunks/` is
the list the table above summarises, and it is 18 commits.

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
The record registry is **73.7 MB** over 192,629 rows, because it carries the natural key and the content
hash per record.

⚠ **This section used to say both ship as release assets rather than committed files, because 73.7 MB is
too much to commit per curation PR. That reasoning was wrong and the measurement is in `vdjdb-db#638`:**
it conflates the one-time blob with the per-PR delta. Committed once, the registry costs **18.9 MB** in
the pack; a 40-row amendment to it costs **2.0 KB**, because git deltifies a sorted TSV essentially
perfectly. Committing it uncompressed is therefore affordable and is what `identity/ids.py` was written
for - a curation PR showing the added, amended and retired rows as a reviewable text diff. Gzipping it
would cost more per PR, not less, and lose the diff.

A build with no registry still runs. It allocates record ids from 1 in corpus order and **warns that the
run is not id-stable**, which is correct for a fork, for a first build and for a corpus replayed at an old
tag. That used to be the normal path and the cost was measured: landing one 40-record chunk moved
`record_id` on **168,723 of 192,753 records (87.5 %)**. `registry/records.tsv` has been committed since
2026-09-29 (#674), so it is now the exception, and
`tests/release/test_registry_is_current.py` fails the run on any amendment, retirement or allocation the
committed file does not already hold. Invariant 3 below is the assertion that would have caught the 87.5 %.

**`replaced_by` landed 2026-09-29 (#693), and the amendment rule is unchanged.** It was specified in the
table above and absent from the file: the columns ended at `amended_from_key_hash` and `note`, and
`amended_from_key_hash` points backwards and only for amendments. So when the amendment pass correctly
refused - two key fields moved at once, which is `test_two_field_change_is_a_new_record_not_an_amendment`
- the record retired and a new id was allocated with nothing linking them, which is the failure this
section opens by naming.

The link is the same identity argument the `chunk.row` tie-break already makes, *the same line of the
same file*, with two guards that make it a fact rather than a guess:

* **one in, one out.** A line that retired two ids, or had two allocated against it, links neither.
* **the receptor is unchanged.** The two natural keys must agree on all six of `cdr3.alpha`, `v.alpha`,
  `j.alpha`, `cdr3.beta`, `v.beta`, `j.beta`. A new record at the line an old one left is a coincidence;
  the same TCR at that line with only its annotation moved is not.

Both directions are tested, and the refusals matter as much as the links: a different TCR at that line
writes nothing, and so does a second retirement from it.

`_backfill_successors` runs on every reconciliation and fills retirements that predate the column from
the same evidence, so the field is not empty for everything that already happened. It is idempotent and
it links nothing it cannot prove. **50 of the 242 retired ids now carry a pointer**: the one that
prompted the issue - #633's `RGPGRAFVTI`, where one assertion about a peptide's source moves
`antigen.species` and `antigen.gene` together, `VDJDB0000021117` -> `VDJDB0000192799` - and the 49 from
the B16 `Plod1 -> Plod2` repair, which moves `antigen.gene` and `meta.subject.cohort` in one commit. The
other 192 stay empty and should: 157 of them are the B16 rows that collapsed, which are deletions and
not amendments, and a pointer there would be a lie.

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
   together with `res/segments.txt` and `res/segments.aaparts.txt`. This took a second branch
   (`feature/retire-res`, #658): the scanner was still naming the V and J a record leaves blank,
   because `arda.cdr3fix` repairs against a named germline and never proposes one. What replaced it
   is `vdjdb.annotate.segments.propose` - the recombination model, falling back to arda's germline
   anchor table where the model declines and for species no model covers.
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

1. junction-nt (#461): `vdjtools.model.infer_nt_batch` on the unique `(species, cdr3, v, j)` set,
   one call per (species, locus) and nothing wrapping it - it releases the GIL and partitions the
   batch across its own kernel threads, so a pool of ours would oversubscribe the machine (hard rule
   3). 0.106 ms/key, 12.44 s for the corpus (§7). No cache: the output is authoritative data, not a
   derived convenience (hard rule 9). Test: the generated `cdr3nt` back-translates to the input
   `cdr3`, and the batched result equals the per-row loop field for field.
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

Checks 1-3 are deterministic string and length work and belong in the build. Check 4 needs a model
and its reference data (fetched from `isalgo/pmhc_data` on first use), so it runs as its own CI job
over the 2,381 `(epitope, MHC)` pairs and publishes a report, not inside `vdjdb build`, whose
offline determinism (hard rule 9) must not depend on a download.

**Checks 1-3 landed 2026-09-29** as `vdjdb.curate.presentation`, writing `out/reports/presentation.tsv`
and its summary on every build. Three things the measurement changed about the plan as written:

* **check 3 ships in one direction only.** A class I groove is closed at both ends, so a 14-mer on
  `HLA-A*02:01` is a question worth asking; a class II record carrying a 9-mer is not, because a paper
  reporting the eluted peptide's core rather than the whole 15-to-25-mer is doing something normal.
  Symmetric, it flags 41 pairs of which 18 carry 7,735 records of 9-mer cores. Asymmetric, 23 pairs.
* **a fourth check was added and is the one that found something.** One molecule filed under two
  `mhc.class` values needs no authority at all: `H2-IAb` is `MHCII` on 20 pairs and `MHCI` on the
  77-record `QVYSLIRPNENPAH`, which is also a 14-mer, so three checks agree on it independently.
* **check 1's `nearest` is not a finding.** 87 pairs over 17,836 records resolve by prefix, and
  almost all are an allele *group* the specification allows - `HLA-A*02` completed to `HLA-A*02:01`,
  which is `mhcmatch` guessing rather than the record being wrong. It is carried in `mhc.resolution`
  so a reader can see which scores rest on a guess, and left out of `finding`.

Total: **93 pairs over 1,213 records** of 2,343 and 192,641, against check 2 firing zero times - a
class I allele on an `MHCII` record does not occur, and the check stays because that is what it is
for. 64 of the 71 unreachable calls are murine class II, where `mhcmatch`'s class II pseudosequences
being HLA is a coverage statement rather than a finding against the record.

One upstream defect fell out and is filed as `antigenomics/mhcmatch#3`: `resolve_allele` trims an
allele to two fields and `class2_key` does not, so `HLA-DRB1*11:01:02` builds a key the bundled FASTA
cannot have while `HLA-DRB1*11:01` resolves - the same molecule and the same groove. 29 of the
corpus's 354 class II pairs are reachable only after trimming, which `curate/presentation.py` does as
a workaround with a test that goes when the fix ships.

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

One carried item. **Step 9, the promiscuity columns, waits on phase 9e**, which is the branch that
introduces the `mhcmatch` call; adding the call twice would put two model versions in one build.

**`record_id` became stable on 2026-09-29, and did not need a release to do it** (#638, #672, #674).
This note used to say it could not: the registry was written at release time, nothing had shipped one,
so every rebuild reconciled against an empty registry and allocated from 1. That was the wrong
conclusion from the right observation - what was missing was not a release but a *committed* registry.
`registry/records.tsv` is now an input, 192,763 rows and 76.7 MB in the tree against 14.8 MB in the
pack and kilobytes per amendment, written only by `vdjdb identity update` on a branch that changes the
corpus and read by `vdjdb build`. `tests/release/test_registry_is_current.py` fails the run on any
amendment, retirement or allocation the committed file does not already hold, so a clean clone
reproduces every id. Before it, landing one 40-record chunk moved `record_id` on 168,723 of 192,753
records.

Two defects surfaced from that gate rather than from a release: the registry had been 592 amendments
stale since `52cb4e2` because a harmonisation change moves the natural key exactly as a chunk edit
does (#672), and the amendment pass could not amend a change to `reference.id` at all, retiring the
record silently, because that was the field it bucketed candidates on (#685).

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

**Four token families in one vocabulary, each token prefixed so a consumer can filter by family
without a second table.**

| Family | Prefix | Token | Example |
|---|---|---|---|
| text | `w:` | a word of the title or abstract, lowercased | `w:influenza` |
| receptor | `v:` | V gene, allele dropped | `v:TRBV9` |
| receptor | `j:` | J gene, allele dropped | `j:TRBJ2-7` |
| receptor | `k:` | a CDR3 k-mer | `k:CAS` |
| receptor | `kv:` | a CDR3 k-mer scoped to the V gene carrying it | `kv:CAS@TRBV9` |
| antigen | `e:` | the epitope sequence | `e:GILGFVFTL` |
| antigen | `ek:` | an epitope k-mer | `ek:GIL` |
| antigen | `a:` | antigen source species | `a:InfluenzaA` |
| antigen | `g:` | antigen gene | `g:M` |
| MHC | `m:` | two-field allele, normalised | `m:HLA-A*02:01` |
| MHC | `ml:` | locus | `ml:HLA-A` |
| MHC | `mc:` | class | `mc:MHCI` |

Each family is in the vocabulary because a question needs it and no other token can stand in.

`k:` **and** `kv:` because "is the `CAS` motif specific to HIV-1, or to its TRBV?" is a comparison
between the lift of `k:CAS` on an epitope's documents and its lift on those of
them that already carry that V gene. One token cannot express that and neither can a single search ranking.

⚠ A species condition answers a provenance question and not a specificity one. `a:HIV-1` is a real axis -
"which papers and receptors are about this species" is where most questions start - but the group is a
union over pMHCs, so a lift over it describes the group rather than being a motif for the pathogen, whose
members were shown different antigens. Condition on `e:<epitope>` plus a restriction when the claim is
about recognition. `docs/standards/terminology.md` has the distinction.

`ek:` because two epitopes sharing a core, or one epitope reported under two source species, are
linked by their k-mers and by nothing else. `e:` alone cannot ask "does this motif go with epitopes
containing `GIL`". k is the same k as the CDR3 family, defaulting to 3: 2,118 epitopes over lengths 7
to 25 residues, so a 3-mer is frequent enough to have a document frequency worth an idf.

The three MHC granularities because restriction is a hierarchy and a question picks its level. The lift
of `k:CAS` given `mc:MHCI`, given `ml:HLA-A` and given `m:HLA-A*02:01` are three different claims, and
collapsing them to one allele token makes the broad ones unaskable. A two-field allele is the
granularity the database curates at; deeper fields are noise for this purpose and are truncated.

**The MHC dictionary is its own artifact**, `corpus/mhc.parquet`, one row per distinct MHC call in the
database: the call as curated, the normalised two-field allele, the locus, the chain it sits on, the
curated class, and whether IPD-IMGT/HLA knows the allele (`proofreading/mhc_alleles.tsv.gz`, 46,005
alleles). It exists as a table rather than as three token-generating functions because a downstream
tool that wants to group VDJdb by locus should read one mapping, not re-derive one, and because the
`known` column makes the dictionary its own proofreading report. The three MHC token families are
projections of it, so they cannot disagree with each other.

**Weighting.** Sublinear term frequency `1 + log tf`, because a paper reporting 10,000 receptors
would otherwise dominate every receptor token; smoothed inverse document frequency
`log((N + 1) / (df + 1)) + 1`; L2 normalisation per document, so a long abstract and a short one are
comparable. These are scikit-learn's conventions, which makes the implementation checkable against
`TfidfVectorizer` on a fixture rather than only against itself.

**The artifact.** Three parquet files in the release, plus the TSVs beside them:

```
corpus/documents.parquet   document_id, reference.id, kind, pmid, year, n_terms, n_records
corpus/terms.parquet       term_id, term, family, df, idf
corpus/postings.parquet    document_id, term_id, tf, weight
corpus/mhc.parquet         mhc, allele, locus, chain, mhc.class, known
```

Long postings rather than a sparse-matrix format, because the consumer is polars or duckdb and the
query is a join. `term_id` is assigned by sorted term and `document_id` by sorted reference, so both
are total orders and the files are reproducible without a hash. Postings are sorted by
`(term_id, document_id)`, which makes a term lookup one contiguous slice.

**The PubMed layer.** `src/vdjdb/corpus/pubmed.py` is the one place that talks to NCBI for
bibliographic records, and it reuses the transport `summary/references.py` already has rather than
opening a second one.

* `efetch` with `db=pubmed&retmode=xml`, 200 ids per request, never one per id. One call returns
  title, abstract, journal, year and DOI together, so the corpus needs no second query.
* NCBI etiquette is in the transport, not in each caller: a `tool` and `email` parameter, at most
  three requests a second without an API key, `$NCBI_API_KEY` honoured when present, and retry with
  backoff on 429 and 5xx. Four requests cover the whole corpus, so this is about being a good client
  rather than about throughput.
* An abstract is assembled from its `AbstractText` sections in order, with the structured labels
  (`BACKGROUND`, `METHODS`) kept as text, because they are words a query can match.
* A PMID that returns no record is reported, never silently dropped: `reference.id` values are curated
  and a dead one is a curation finding.

**What is fetched and what is committed.** Receptor, antigen and MHC tokens come from `chunks/` and
need no network at all. Only the text family does, and it follows hard rule 5's pattern: the running
text is an input to the build and never an output of it.

`vdjdb corpus refs` fetches the records and writes two committed, reviewed inputs, refreshed by their
own pull request exactly as `summary/reference_years.tsv` already is, so every build is offline and
deterministic:

| File | Contents |
|---|---|
| `corpus/pubmed.tsv` | `reference.id, pmid, year, journal, title, doi, abstract_words, abstract_sha256` |
| `corpus/text_terms.tsv` | `reference.id, term, tf` - term counts only |

No abstract text is committed or shipped. `abstract_sha256` is what says whether an abstract changed
between refreshes, and `abstract_words` its length, so the corpus is auditable without keeping the
text. Titles are kept because a title is a fact a bibliography carries anyway. Estimated at 610
documents and roughly 120 distinct terms each, about 73,000 rows and 2 MB. The tokeniser is in the
repository, so the transform is reproducible from the text even though the text is not stored.

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

**State, 2026-09-27: done.** 661 documents, 268,160 terms over twelve families, 1,106,074 postings, in
2.8 s. Every weight agrees with `sklearn.feature_extraction.text.TfidfVectorizer` to nine decimal
places, and every document's L2 norm is 1.0.

Four things the work established that the plan did not anticipate.

**The instrument recovers a motif nothing told it about.** Conditioning on `e:GILGFVFTL`, `k:IRS` is
the highest-lifting of the 2,342 CDR3 3-mers with 50 or more occurrences, at **2.663x** against a
median of 1.170x, and 25 of the 29 RS-bearing 3-mers sit above that median. The RS motif of
influenza-M1-specific beta CDR3s is documented immunology and nothing in the build encodes it. The
control holds too: `k:CAS`, the germline start of nearly every beta CDR3, lifts **0.969** on HIV-1
documents and 1.006 with TRBV9 held, so it reads as germline rather than epitope-associated.
`tests/release/test_corpus_reproduction.py` pins all of it.

**A document-level lift cannot answer a common token**, so `lift` has two modes. `k:CAS` is in 614 of
661 documents, which bounds its document-level lift near 1 however specific it is. The occurrence mode
is scoped to the term's own family: pooling every family into one denominator flattens the whole
comparison to within 2 % of 1.0, because the rate then moves with how many epitopes a paper studied.

**The epitope k-mer family earns its place on a measured case.** Searching `GILGFVFTL` puts nine exact
reporters in the top ten and, tenth, `PMID:27036003`, which reports `GILEFVFTL` and `GILGLVFTL`. Those
are single-residue variants, and an exact-epitope search cannot find the altered-peptide-ligand study
of the epitope it is asking about.

**Two spellings of one PMID were silently collapsing.** The corpus carries both `PMID:34433824` and
`PMID: 34433824`, and a one-PMID-to-one-reference map kept whichever the iteration reached last: the
first fetch wrote 609 records where 610 references are PubMed ones, with nothing raised. Fixed to
one-to-many, and `missing` now reports `reference.id` values. 610 records, 65,720 term counts, 8,890
distinct words, 1.6 MB committed.

The MHC dictionary is 213 rows over 25 loci, and doubles as a proofreading report: 184 `known`, 27
`unchecked` (murine H2 and the light chain are outside IPD-IMGT/HLA), and exactly two `unknown` -
`HLA-A*08:01` and `HLA-B*12`, the same two `assemble.epitopes.mhc_status` already finds, because the
dictionary reads that verdict rather than re-deriving one.

### Bootstrap order - done, 2026-09-27

Completed in this order (`ROADMAP_local.md` §44): land the workflows on `master` → create `dev` → let one full `build.yml` run green on `dev` so the
check names exist → then apply the `dev`, `master` and tag rulesets. A required check that has
never run blocks every PR forever. Linear history is deliberately **not** applied: `master` carries 168 reachable merge commits and
the gitflow puts `--no-ff` feature merges on `dev`, so the rule would reject every future update and
train the bypass habit. The four rulesets that did go on are recorded in `ROADMAP_local.md` §44.2.
