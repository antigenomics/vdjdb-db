# Building and releasing

The build is a `uv`-managed Python package. It reads `chunks/`, `patches/` and `proofreading/`, and
writes everything else. Nothing computed is stored between builds: every output is recomputed from
`chunks/` on every run.

## Local build

```bash
uv sync
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files, from the tables that shipped
uv run vdjdb convert airr --tables out/tables # AIRR Rearrangement + Reactivity
uv run vdjdb motifs --tables out/tables       # TCRNET + TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, rendered offline and checked
uv run vdjdb identity check --tables out/tables # the identifier invariants
uv run vdjdb corpus build --tables out/tables # the reference corpus: tf-idf over 12 token families
uv run vdjdb diff <reference.zip> out/legacy  # compare against a released zip
uv run pytest -q
```

Output goes to `out/`, not `build/`: `build/` is gitignored as a Python packaging convention.

Optional dependency groups: `motifs` for the motif stage, `docs` for this site, `test` for the
suite, `tuning` for the clustering bake-off under `docs/tuning/`.

## The junction-nucleotide stage

`annotate.junction.add_junction_nt` was 87.2 % of the assembly stage - 407.64 s of 467.72 s on a
4-vCPU runner over 192,793 records - because `vdjtools.model.infer_nt` wraps a native DP in per-row
Python and the wrapper, not the DP, was the cost. That profile is what
[`antigenomics/vdjtools#181`](https://github.com/antigenomics/vdjtools/issues/181) was opened on, and
`vdjtools` 4.5 answers it with `infer_nt_batch`.

So the stage is **one batched call per (species, locus)** over the distinct
`(species, gene, cdr3, v, j)` keys - 187,055 of them rather than every row, which is rule 4's
deduplication - and nothing wraps it. `infer_nt_batch` releases the GIL and partitions the batch
across its own kernel threads; a pool of our own would oversubscribe the machine and read as
"batching did not help", which is hard rule 3 and section 0e of `CLAUDE.md` both.

Measured on 3,000 distinct human TRB keys from the corpus, 16 cores: **1.115 ms/key serial against
0.106 ms batched, 10.5x, and all 3,000 nucleotide sequences identical.** End to end on the 4-vCPU
runner the stage goes **407.64 s to 108.54 s** and the assembly step 467.72 s to 166.54 s, so it is
65.2 % of that step rather than 87.2 %; on a 16-core laptop it is 12.44 s and `vdjdb build` is 45.9 s
wall. The runner wins less because `infer_nt_batch` defaults to `hardware_concurrency - 2` threads,
which is two there.
`tests/unit/test_junction.py` asserts both halves of that, because each catches a different failure -
identity catches a batch call that is not the same computation, and the ratio catches a regression to
the loop or a batch call that loops internally, neither of which changes an answer.

**Every inferred sequence encodes the junction it came from, and that is by construction rather than
by luck.** The DP enumerates `(V, delV) x (J, delJ) x (D, delD, position)` and picks the best codon
assignment *within* each scenario, so a scenario that cannot spell the given residues has probability
zero and is never a candidate; anything the model cannot encode comes back null rather than wrong.
Probed on human TRB: a stop codon, an `X`, a `Z`, a one- or two-residue junction, an empty string and
a true CDR3 with its anchors stripped are all declined. The one input that survives with a difference
is a lower-case junction, where the nucleotides are right and the comparison is case-sensitive - and
`vdjdb qc` rejects a residue outside the 20 upper-case letters, with zero such chains in the corpus.
Measured on the built corpus: **263,437 of 285,989 chains carry an inferred `cdr3nt`, 0 mismatches, 0
whose length is not exactly three nucleotides per residue**, gated by
`tests/release/test_tables_contract.py`. That check translates the whole column in one threaded native
call, `vdjtools._core.translate_junctions` - 0.023 s against 0.284 s for `vdjtools.model.translate` in
a Python loop, 12.3x, identical on every row.

It previously ran as four worker processes over contiguous parquet slices of the key set, each an
ordinary invocation of a subcommand that existed only to be that worker. The subcommand,
`src/vdjdb/__main__.py`, the slice arithmetic and the worker-count argument are all gone: one batched
call has no worker count, so rule 7's "never let worker count change the answer" holds by
construction rather than by a tiling test.

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

## `dev` and `master` cannot be deleted

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
