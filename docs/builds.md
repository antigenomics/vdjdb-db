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
uv run vdjdb diff <reference.zip> out/legacy  # compare against a released zip
uv run pytest -q
```

Output goes to `out/`, not `build/`: `build/` is gitignored as a Python packaging convention.

Optional dependency groups: `motifs` for the motif stage, `docs` for this site, `test` for the
suite, `tuning` for the clustering bake-off under `docs/tuning/`.

## Comparing against a release

`vdjdb diff` compares a candidate build against a released bundle in three passes: the file set,
two digests per file (raw and canonical, rows sorted), and a row-level classification that
attributes every changed cell to a declared rule in `rules/expected_diffs.toml`.

Any cell difference not matched by a rule fails, and so does a rule that fires a different number of
times than declared. Every rule therefore declares a measured row count rather than a description.

Canonical equality is the gate; raw equality is informational. The legacy pipeline iterated
`os.listdir("../chunks")`, which is readdir order and therefore filesystem- and host-dependent, so
the 2026-06-03 release cannot be reproduced byte-for-byte on another machine, including by the
pipeline that produced it. The current reader sorts, and a release records its chunk order, so raw
equality is achievable for releases built by this code.

## Testing against older releases

The comparison currently runs one build against one release, 2026-06-03. A single release cannot
distinguish a rule that is correct from a rule that happens to fit that release.

There are 43 published releases, back to 2017-06-13, and each one is a matched pair: the inputs
(`chunks/`, `patches/`, `res/`, `proofreading/` at that tag) and the outputs (the zip that shipped).
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

Two derived tables make the build offline and deterministic. Both are committed, reviewed inputs,
refreshed by their own pull request and never written by a build:

- `summary/reference_years.tsv` - the publication year of every reference, so the dashboard does not
  call NCBI while it renders. Rebuilt by `uv run vdjdb refs`.
- `docs/tuning/scorecard.tsv` - the clustering bake-off, 252 configurations scored through one
  harness. Rebuilt by `docs/tuning/sweeps.py`.

Neither is a cache: a cache is a stored answer to the question the build is currently asking, and
these are inputs that arrive by review.
