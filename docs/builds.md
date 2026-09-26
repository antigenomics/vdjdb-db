# Building and releasing

The build is a `uv`-managed Python package. It reads `chunks/`, `patches/` and `proofreading/`, and
writes everything else. Nothing computed is stored between builds: every output is recomputed from
`chunks/` on every run, because a stored intermediate is authoritative until the code that produced
it changes and then silently wrong in a way no test sees.

## Locally

```bash
uv sync
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files, from the tables that shipped
uv run vdjdb convert airr --tables out/tables # AIRR Rearrangement + Reactivity
uv run vdjdb motifs --tables out/tables       # TCRNET + TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, rendered offline and checked
uv run vdjdb diff <reference.zip> out/legacy  # the difference ledger
uv run pytest -q
```

Output goes to `out/`, not `build/` — `build/` is gitignored as a Python packaging convention and
using it for release artifacts is confusing.

Optional dependency groups: `motifs` for the motif stage, `docs` for this site, `test` for the
suite, `tuning` for the clustering bake-off under `docs/tuning/`.

## The difference ledger

`vdjdb diff` compares a candidate build against a released bundle in three passes: the file set,
two digests per file (raw and canonical, rows sorted), and a row-level classification that
attributes **every changed cell** to a declared rule in `rules/expected_diffs.toml`.

Any cell difference not matched by a rule fails. A rule that fires a different number of times than
declared also fails — that is what turns "we think it is the same" into a gate, and it is why every
rule carries a measured row count rather than a description.

**Canonical equality is the gate; raw equality is informational.** The legacy pipeline iterated
`os.listdir("../chunks")`, which is readdir order and therefore filesystem- and host-dependent, so
the 2026-06-03 release cannot be reproduced byte-for-byte on another machine by anyone, including
the pipeline that produced it. The new reader sorts, and a release records its chunk order so raw
equality becomes achievable going forward.

## Reproducibility

Same inputs, same bytes — in another process, on another host, at another core count.

- **One seed**, `vdjdb.config.SEED`, passed explicitly to every sampler, shuffler and clustering
  init. No module-scope `random.seed()`, no unseeded default.
- **Sort after anything unordered**: a `group_by` without `maintain_order=True`, a `set`, a `dict`
  keyed on strings under `PYTHONHASHSEED`, a thread pool, `os.listdir`. An unstable sort over
  equal-labelled groups once made the ledger report three different counts for one comparison.
- **Worker count never changes the answer.** Split into as many big contiguous slices as there are
  workers, reassemble in slice order, never a pool of small tasks.

## Inputs that are committed rather than computed

Two derived tables exist to make the build offline and deterministic. They are committed, reviewed
inputs refreshed by their own pull request, **never written by a build**:

- `summary/reference_years.tsv` — the publication year of every reference, so the dashboard does not
  call NCBI while it renders. Rebuilt by `uv run vdjdb refs`.
- `docs/tuning/scorecard.tsv` — the clustering bake-off, 252 configurations scored through one
  harness. Rebuilt by `docs/tuning/sweeps.py`.

Neither is a cache. A cache is a stored answer to the question the build is currently asking; these
are inputs that arrive by review.
