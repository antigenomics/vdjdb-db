# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 4 — `feature/pipeline-core`**, merging to `dev`.

| Phase | Branch | State |
|---|---|---|
| 0 | `feature/dev-baseline` | merged — docs, package skeleton, chunk lint, CI, identity, guard |
| 1 | `feature/schema` | merged — the field registry every column order projects from |
| 2 | `feature/golden-harness` | merged — `vdjdb diff`, reproducible, validated against the release |
| 3 | `feature/io-qc` | merged — polars reader, vectorised QC, corpus normalised, CI gates `--strict` |
| 4 | `feature/pipeline-core` | **ledger PASS**; `py_src/` retired |

```
uv run vdjdb build --out out/
uv run vdjdb diff ref/vdjdb-2026-06-03.zip out/legacy \
    --only vdjdb.txt,vdjdb.slim.txt,vdjdb_full.txt      # -> PASS
```

## Where the build stands

| | Legacy pandas | Now |
|---|---|---|
| Wall time | 344 s | **12 s** |
| Peak RSS | 2.16 GB | ~1 GB |
| Output | three files, assembled directly | two definitive tables, three files projected from them |

Every difference from the 2026-06-03 release is a declared rule in `rules/expected_diffs.toml`
firing its exact measured count. The largest are the `web.cdr3fix.unmp` truthiness bug (7,973 rows)
and a family of pandas type coercions the all-string reader undoes (~55k cells).

## Next

1. **Phase 5** (`feature/arda-cdr3fix`) — replace `annotate/_legacy_fixer/` with `arda.cdr3fix`;
   measure the ledger delta, then freeze it as declared rule counts.
2. **Phase 6** (`feature/new-format`) — ship the definitive tables, add `evidence`.
3. `uv lock`, push `dev`, let `chunk-check` run once **before** applying branch protection.

## Blocked / needs a decision

| Item | Blocks | Question |
|---|---|---|
| `vdjmatch` `_zip_asset` patch | the first multi-zip release | needs a patch + release in `antigenomics/vdjmatch`. See `ROADMAP.md` §3.1 |
| Motif `coef` calibration | phase 11 | needs the motif pipeline running before it can be fitted |

## Known, not yet fixed

- 14 records are reported by two chunks with the same PDB id in different letter case, and one
  `meta.epitope.id` carries a float. Both are phase 9 nomenclature work; both are visible in the
  ledger today.
- `hotfix` carries `compute_pdb_cdr3fix.py`, which phase 5 supersedes — the same lookup falls out
  of `arda.markup_batch`.
