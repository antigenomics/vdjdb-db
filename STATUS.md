# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 0 — `feature/dev-baseline`.** Repo documents, package skeleton, the two fast CI workflows.

| Item | State |
|---|---|
| `ROADMAP.md` | done |
| `CLAUDE.md` | done |
| `STATUS.md` | done |
| `pyproject.toml` + `uv.lock` | in progress |
| `src/vdjdb/` skeleton + `vdjdb` CLI entry point | in progress |
| `.github/workflows/chunk-check.yml` | in progress |
| `.github/workflows/branch-policy.yml` | in progress |

Nothing in phase 0 changes build behaviour. `release.sh` and the existing pandas pipeline still work
untouched.

## Branches

```
master                 3389001   (origin/master)
dev                    3389001   branched from master, not yet pushed
feature/dev-baseline   3389001   ← current
hotfix                 a39fc86   1 commit ahead of master (compute_pdb_cdr3fix.py); unmerged
```

`hotfix` carries `py_src/compute_pdb_cdr3fix.py` plus an untracked `py_src/pdb_cdr3fix.tsv`. Decide
whether that lands on `dev` or is superseded by phase 5 (`arda.cdr3fix` makes the script redundant —
the same lookup falls out of `markup_batch`).

## Blocked / needs a decision

| Item | Blocks | Question |
|---|---|---|
| `evidence.*` columns | phase 6 | Nothing in this repo produces the five columns production serves. Does the new format take ownership? |
| Seven unshipped side outputs | phase 6 | `vdjdb_full_filtered.txt`, three `*_broken.txt`, three `*_scored.txt` → `build/reports/*`, or retire? |
| Five discarded chunk columns | phase 1 | `submitter`, `chunk.id`, `comment`, `meta.subset.frequency`, `method.pairing` — promote or declare `TOLERATED_DROPPED`? |
| `vdjmatch` `_zip_asset` patch | the first multi-zip release | Needs a patch + release in `antigenomics/vdjmatch` first. See `ROADMAP.md` §3.1 |

## Next

Phase 1 (`feature/schema`) then phase 2 (`feature/golden-harness`). **Phase 2 must show zero diffs
against the current pandas build before any behaviour changes land** — it is the instrument every
later phase is measured with.
