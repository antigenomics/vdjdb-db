# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 0 — `feature/dev-baseline`** (merged) and **`feature/record-identity`** (on `dev`).

| Item | State |
|---|---|
| `ROADMAP.md` | done |
| `CLAUDE.md` | done |
| `STATUS.md` | done |
| `pyproject.toml` | done (`uv.lock` still to generate — needs a network sync) |
| `src/vdjdb/` skeleton + `vdjdb` CLI entry point | done |
| `vdjdb qc` chunk lint + 15 unit tests | done — 230 chunks, 103 findings (99 CRLF) |
| `.github/workflows/chunk-check.yml` | done |
| `.github/workflows/branch-policy.yml` | done |
| record identity + registry (`src/vdjdb/identity/`) | done — 192,753 rows → 192,734 records in 2.0 s, ids stable |
| `docs/outputs.md` — spec of every produced file | done |
| proprietary-data guard (`src/vdjdb/validate/`) | done — wired into `chunk-check.yml` |
| 65 unit tests | passing |

Nothing in phase 0 changes build behaviour. `release.sh` and the existing pandas pipeline still work
untouched.

## Branches

```
master                   3389001   (origin/master)
dev                      711ede8   phase 0 merged
feature/record-identity  f6dc4d4   ← current; identity + spec + guard
feature/dev-baseline     99779ca   merged into dev
hotfix                   a39fc86   1 commit ahead of master; unmerged
```

Nothing is pushed yet. `dev` and the phase-0 commit are local.

`hotfix` carries `py_src/compute_pdb_cdr3fix.py` plus an untracked `py_src/pdb_cdr3fix.tsv`. Decide
whether that lands on `dev` or is superseded by phase 5 (`arda.cdr3fix` makes the script redundant —
the same lookup falls out of `markup_batch`).

## Decided since

All three questions are settled — see `ROADMAP.md` §9. The new format owns `evidence.*`; the five
debug columns are kept; the side outputs are produced but not zipped.

## Blocked / needs a decision

| Item | Blocks | Question |
|---|---|---|
| `vdjmatch` `_zip_asset` patch | the first multi-zip release | Needs a patch + release in `antigenomics/vdjmatch` first. See `ROADMAP.md` §3.1 |
| Motif `coef` calibration | phase 11 | Needs the motif pipeline running before it can be fitted against the study-support objective |

## Next

1. `uv lock` (needs network) and confirm the dependency set resolves — `vdjtools` and `arda-mapper`
   are the only base deps; `mirpy-lib[bench]` is behind the `motifs` extra.
2. Push `dev`, let `chunk-check` run once so the check name exists, **then** apply branch protection.
   A required check that has never run blocks every PR forever.
3. Phase 1 (`feature/schema`), then phase 2 (`feature/golden-harness`). **Phase 2 must show zero
   diffs against the current pandas build before any behaviour changes land** — it is the instrument
   every later phase is measured with.

## Known, not yet fixed

- 99 chunk files use CRLF line endings; one has a prose sentence as a column name
  (`PMID_24512815.txt`) and one a bare leading tab (`PMID_40694338.txt`). `chunk-check` reports
  these but does not fail on them yet — the `.tsv` migration (#497, phase 3) normalises them, and
  until then a hard gate would block every unrelated submission.
- `src/` is temporarily mixed: the retired Groovy sits beside the new `src/vdjdb/` package. Phase 14
  moves the Groovy to `attic/`. Hatchling only packages `src/vdjdb`, so nothing breaks meanwhile.
