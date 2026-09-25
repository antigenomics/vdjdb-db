# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 0 — `feature/dev-baseline`.** Repo documents, package skeleton, the two fast CI workflows.

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

Nothing in phase 0 changes build behaviour. `release.sh` and the existing pandas pipeline still work
untouched.

## Branches

```
master                 3389001   (origin/master)
dev                    711ede8   ← phase 0 merged; not yet pushed
feature/dev-baseline   99779ca   merged into dev
hotfix                 a39fc86   1 commit ahead of master (compute_pdb_cdr3fix.py); unmerged
```

Nothing is pushed yet. `dev` and the phase-0 commit are local.

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
