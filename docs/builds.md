# Build and release

Use [build integrity](build-integrity.md) for candidate validation, CI artifact checks and
promotion from one release to the next.

The build is a `uv`-managed Python package. It reads `chunks/`, `patches/` and `proofreading/`, and
writes everything else. Nothing computed is stored between builds: every output is recomputed from
`chunks/` on every run.

## Local build

```bash
uv sync --extra motifs --extra summary --extra test
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files, from the tables that shipped
uv run vdjdb convert airr --tables out/tables # AIRR Rearrangement + Reactivity
uv run vdjdb motifs --tables out/tables       # TCRNET + TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, rendered offline and checked
uv run vdjdb identity check --tables out/tables # the identifier invariants
uv run vdjdb corpus build --tables out/tables # the reference corpus: tf-idf over 12 token families
uv run vdjdb diff <reference.zip> out/legacy  # compare against a released zip
uv run vdjdb release --tag vYYYY.MM.P --dry-run  # the three bundles, without touching the tree
VDJDB_REFERENCE_ZIP=reference.zip uv run pytest -q
```

Output goes to `out/`, not `build/`: `build/` is gitignored as a Python packaging convention.

`vdjdb release` writes tracked `latest-version.txt`. Use `--dry-run` to inspect bundles without
changing the repository's download URL. Supply an intended release tag when publishing.

Download the 2026-06-03 release as `reference.zip` before running the comparison or release tests.
A plain `pytest -q` skips the reference checks when `VDJDB_REFERENCE_ZIP` is unset. See
[dashboard requirements](dashboard.md#rendering-it) for R and pandoc installation.

Optional dependency groups: `motifs` for the motif stage, `docs` for this site, `test` for the
suite, `tuning` for the clustering bake-off under `docs/tuning/`.

## Choose the scope of a check

Use the smallest stage that answers the current question. Every invocation recomputes its
results; partial assessment is not a stored intermediate for a later release build.

| Question | Command or workflow | What it establishes |
|---|---|---|
| Is this chunk well formed? | `vdjdb qc <chunk>` | Source validation for the selected files |
| What does it contribute? | `vdjdb submission <chunk>` / `chunk-check.yml` | Corpus-relative scores, novelty and source-review reports, before expensive downstream annotation |
| What does mhcmatch predict for these reported peptides? | `vdjdb assess-epitopes <chunk>` / `assessment.yml` | Selected-pair assessment, independent of junction inference, motifs and dashboard |
| Recompute assessment for a built corpus | `vdjdb assess-epitopes --tables out/tables` | Fresh assessment from the records table, without reading previous predictions |
| Re-emit a format | `vdjdb make legacy` / `vdjdb convert airr` | A projection from the explicitly selected built tables |
| Assess corpus-wide consequences and integrate | `build.yml` | Full assembly, corpus, motifs, dashboard, release comparison and release tests |

For selected-chunk predictions:

```bash
uv run vdjdb epitope-reference --out out/inputs
uv run vdjdb assess-epitopes chunks/PMID_<id>.tsv --out out/assessment \
  --pmhc-reference out/inputs/pmhc/pmhc_full.tsv.gz --jobs 4
```

This writes only `epitope_assessment.parquet`, `epitope_assessment.tsv` and
`assessment-timings.tsv`. Without a reference, reported pairs remain visible with
`reference_not_supplied`. Files are harmonised and deduplicated through the existing master
assembly; support counts describe the selected input only. For corpus-relative scores and
cross-publication checks, use `submission`, which reads the complete source corpus. Neither
partial command writes the identity registry. A subset is not a release bundle.

The opt-in `assessment.yml` workflow takes whitespace-separated tracked `chunks/` paths from
the chosen revision and a `predict` switch. It runs selected QC, the corpus-relative submission
report and the standalone assessment on a hosted runner, uploading one assessment artifact.
It does not run motifs, junction-nucleotide inference, R or the release comparison. Keep ordinary
chunk-check CI fast; request this additional workflow when a peptide/MHC question needs it.
The existing `chunk-check.yml` manual entry point also accepts `assessment_chunks` and
`assessment_predict`, calling the same workflow as an additional job. Ordinary pull requests
leave it disabled. To exercise it on a feature revision:

```bash
gh workflow run chunk-check.yml --ref <branch> \
  -f assessment_chunks=chunks/PMID_<id>.tsv -f assessment_predict=true
```

Presentation calibration still has a fixed per-species/class setup cost for small submissions.
The first full-corpus assessment measured 129.3 seconds on four CI cores, 34.6% of timed
**assembly**, not of the complete workflow. That initial timing used marginal calibration; class-II ranks now use matched-length backgrounds,
so it is historical rather than a current runtime claim. Assessment timing reports also record
`peak_tree_rss_mb`, the aggregate RSS of parent and descendants sampled every 50 milliseconds.
The assembly gate separately budgets this at 8,192 MiB, based on the measured 4,931.4 MiB cold
process-tree peak, while keeping its 4,096 MiB parent budget. Unmeasured tree peaks are blank.
Class-II scorer ownership is bounded to one allele at a time because mhcmatch 1.20.0 retains
frame computations across all queried lengths. This uses fresh public scorer/store objects,
without editing private dependency caches. [mhcmatch #4](https://github.com/antigenomics/mhcmatch/issues/4)
tracks the upstream bounded-memory batch API and removal of this workaround.
Read the standalone timing report for a selected
submission instead of extrapolating that ratio. Full integration must still recompute non-additive
outputs such as clustering, motifs and corpus-wide statistics from the combined corpus.

## Comparing against a release

`vdjdb diff` compares a candidate build against a released bundle in three passes: the file set,
two digests per file (raw and canonical, rows sorted), and a row-level classification that
attributes every changed cell to a declared rule in `rules/expected_diffs.toml`.

Any cell difference not matched by a rule fails, and so does a rule that fires a different number of
times than declared. Every rule therefore declares a measured row count rather than a description.

A correction to an identity column has no cell to attribute: the row leaves the reference bucket and a
different row arrives in the candidate one. Those two counts are declared per file in a `[[row_delta]]`
block, with a `note` giving every reason they are what they are, and the run fails if either moves. **A
new chunk moves all of them**, because its records are rows the reference cannot contain, so landing one
means re-measuring the three declarations and extending each note with the chunk and its record count -
the entry for `PMID_18025130` is the worked example. `vdjdb diff --report` prints declared against
measured per file and flags which one moved.

### What the comparison is, and is not, the instrument for

It runs on the **assembled legacy zip**, so the thing compared is the file a consumer downloads.
Three of that bundle's twelve members cannot be judged by keying their rows against the reference,
and `[measured_elsewhere]` in `rules/expected_diffs.toml` names each one with the instrument that
gates it instead. They are still read, digested, row-counted and printed; what they are exempt from
is the requirement that every changed cell match a rule.

| Member | Why keying it says nothing | What gates it |
|---|---|---|
| `cluster_members.txt` | `cid` is `<species>.<chain>.<epitope>.<n>` and `n` is a position in a sorted list, so one renumbered cluster relabels every cluster after it | `vdjdb motif-metrics`, 18 axes per chain, including `partition_neighbours_preserved` - of every clonotype the release clustered, the fraction of its cluster-mates this build still gives it. Column count and order by `tests/release/test_reference_contract.py` |
| `motif_pwms.txt` | a row is one PWM cell of one cluster, so it has no identity that survives a re-clustering | the same two |
| `vdjdb_summary_embed.html` | a fresh render every build | `summary/check_summary.py`: the ordered headings and tables, PNG dimensions decoded from the IHDR, ColorBrewer anchors, and SSIM against the last release |
| `LICENSE` | not a table | shipped verbatim from the repository root |

Measured before those declarations existed: comparing every member of the legacy zip reported **101,877
unattributed cells**, and every one was in a motif file or in `latest-version.txt` while the five
legacy tables were fully attributed. The workflow's answer had been `--only` naming those five, which
is the same exemption with no reason recorded and no digest taken of the other seven.

A member the reference does not contain is declared in `[members]`. Today that is
`cluster_members_tcremp.txt` and `motif_pwms_tcremp.txt`; an undeclared one still fails.

Canonical equality is the gate; raw equality is informational. The legacy pipeline iterated
`os.listdir("../chunks")`, which is readdir order and therefore filesystem- and host-dependent, so
the 2026-06-03 release cannot be reproduced byte-for-byte on another machine, including by the
pipeline that produced it. The current reader sorts, and a release records its chunk order, so raw
equality is achievable for releases built by this code.

## Testing against older releases

The comparison currently runs one build against one release, 2026-06-03. A single release cannot
distinguish a rule that is correct from a rule that happens to fit that release.

There are 43 published releases, back to 2017-06-13, and each one is a matched pair: the inputs
(`chunks/`, `patches/`, `proofreading/` at that tag, and `res/` for a tag older than the #658
retirement) and the outputs (the zip that shipped).
Each pair is a regression test: take a git worktree at the tag and build it with the current code:

```bash
git worktree add /tmp/vdjdb-2023 2023-06-01      # the corpus as it was
cd /tmp/vdjdb-2023 && uv run --project <repo> vdjdb build --out out/
uv run vdjdb diff <2023-06-01.zip> out/legacy    # does today's code reproduce that release?
```

What a replay across releases adds:

- Code changes separated from corpus changes. Every rule in `expected_diffs.toml` is currently
  declared against one release. A rule that fires the same way across 2021, 2023 and 2026 is a
  property of the code; one that fires only on 2026 is a property of today's data.
- Format drift that the current corpus no longer contains. The 2017–2019 chunks use header
  shapes, MHC spellings and nomenclature that later curation cleaned up, and those are the inputs a
  reader of an old release would hand back to the tool.
- Dates for known defects. `web.cdr3fix.unmp` is wrong on 7,973 rows today; replaying it across
  releases says when it started.

Older releases were produced by older code with different column sets, so the comparison needs
era-scoped rules rather than one global set: `vdjdb.txt` gained columns, the slim table changed
shape, and the motif files did not exist at all before 2018. The first replay therefore produces a
large diff that is mostly format difference, and the work is in classifying it.

Start with one release, 2023-06-01: recent enough to share most of the schema, old enough that a
rule which only fits 2026 will not fit it.

## Protect integration and release branches

The repository has `delete_branch_on_merge` enabled, which is right for feature branches and wrong for
the long-lived ones: a `dev` to `master` pull request has `dev` as its **head**, so merging it used to
delete `dev`. That happened on PR #573 and again on PR #584. Nothing was lost except the ref - every
commit stays reachable from `master` - but a contributor pulling in the window between the merge and the
restore gets a confusing error, and the symptom does not look like a branch-protection question.

It is fixed by a ruleset rather than by turning the setting off, so feature branches still tidy
themselves up:

| Ruleset | Refs | Rules | Bypass actors |
|---|---|---|---|
| `no-deletion-dev-master` | `refs/heads/dev`, the default branch | `deletion` | **none** |
| `protect-dev` | `refs/heads/dev` | `deletion`, `non_fast_forward`, `required_status_checks` | org admin, admin, triage |
| `vdjdb-master` | the default branch | `deletion`, `non_fast_forward`, `pull_request`, `required_status_checks` | org admin |

⚠ **The bypass list is why `protect-dev` alone did not work.** It already carried a `deletion` rule, but
a merge runs as the person merging, and an admin's `always` bypass covers the deletion too. A GitHub
ruleset's bypass list applies to the whole ruleset and cannot be set per rule, so the deletion rule needs
its own ruleset with an empty bypass list. That is what `no-deletion-dev-master` is, and it holds: a
`git push origin --delete dev` from an admin is now rejected with `GH013`.

`protect-dev` keeps its bypass on purpose, so a maintainer can still push a fix to `dev` directly when a
status check is stuck. Deletion is the only thing nobody may do.

If `dev` is ever absent again, the restore is one command, and the commit is the merge's second parent:

```bash
git push origin $(git rev-parse master^2):refs/heads/dev
```

## Reproducibility

Same inputs, same bytes: in another process, on another host, at another core count.

- One seed, `vdjdb.config.SEED`, passed explicitly to every sampler, shuffler and clustering init.
  No module-scope `random.seed()`, no unseeded default.
- Sort after anything unordered: a `group_by` without `maintain_order=True`, a `set`, a `dict`
  keyed on strings under `PYTHONHASHSEED`, a thread pool, `os.listdir`. An unstable sort over
  equal-labelled groups once made the release comparison report three different changed-row counts
  for one comparison.
- Worker count never changes the answer. Split into as many big contiguous slices as there are
  workers, reassemble in slice order, never a pool of small tasks.

## Committed inputs

Three derived tables make the build offline and deterministic. All are committed, reviewed inputs,
refreshed by their own pull request and never written by a build:

- `summary/reference_years.tsv` - the publication year of every reference, so the dashboard does not
  call NCBI while it renders. Rebuilt by `uv run vdjdb refs`.
- `docs/tuning/scorecard.tsv` - the clustering bake-off, 252 configurations scored through one
  harness. Rebuilt by `docs/tuning/sweeps.py`.
- `proofreading/epitope_promiscuity.tsv` - which class I alleles can present each epitope, 13,510
  (epitope, allele) pairs over 1,729 epitopes, of which 1,739 are pairings VDJdb records and 11,771
  are alleles only the predictor reaches. Rebuilt by `uv run vdjdb promiscuity`, which loads
  `mhcmatch`'s reference panel from HuggingFace. It answers a query-side question - filter VDJdb to a
  donor's HLA type and see every record their T cells could have raised - and makes no claim about the
  `mhc.a` a publication reported.

Neither is a cache: a cache is a stored answer to the question the build is currently asking, and
these are inputs that arrive by review.
