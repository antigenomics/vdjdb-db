# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 7 — `feature/airr`**, merging to `dev`.

| Phase | Branch | State |
|---|---|---|
| 0 | `feature/dev-baseline` | merged — docs, package skeleton, chunk lint, CI, identity, guard |
| 1 | `feature/schema` | merged — the field registry every column order projects from |
| 2 | `feature/golden-harness` | merged — `vdjdb diff`, reproducible, validated against the release |
| 3 | `feature/io-qc` | merged — polars reader, vectorised QC, corpus normalised, CI gates `--strict` |
| 4 | `feature/pipeline-core` | merged — the definitive tables, legacy as a projection; `py_src/` retired |
| 5 | `feature/arda-cdr3fix` | merged, **behind `engine="legacy"`** — the swap waits for #327 (agreed 2026-09-25) |
| 6 | `feature/new-format` | merged — tables + evidence + `vdjdb.schema.json`; ledger PASS from the shipped tables |
| 7 | `feature/airr` | **`airr.validate_rearrangement` passes on all 286,047 rows** |

```
uv run vdjdb build --out out/                          # tables + every projection, 15 s
uv run vdjdb make legacy --tables out/tables --out out/legacy-made
uv run vdjdb convert airr --tables out/tables --out out/airr
uv run vdjdb diff ref/vdjdb-2026-06-03.zip out/legacy-made \
    --only vdjdb.txt,vdjdb.slim.txt,vdjdb_full.txt     # -> PASS
VDJDB_REFERENCE_ZIP=ref/vdjdb-2026-06-03.zip VDJDB_TABLES=out/tables \
    uv run pytest -q -m release                        # -> 27 passed
```

## Where the build stands

| | Legacy pandas | Now |
|---|---|---|
| Wall time | 344 s | **15 s** |
| Peak RSS | 2.16 GB | ~1 GB |
| Output | three files, assembled directly | three tidy tables + a joined view, with legacy and AIRR projected from them |

`records` 192,753 × 33 · `chains` 286,047 × 17 · `evidence` 53,913 × 8 · `vdjdb` (view) 286,047 × 54.
All five legacy members are byte-identical whether projected from memory or read back from parquet.
Every difference from the 2026-06-03 release is a declared rule in `rules/expected_diffs.toml` firing
its exact measured count — largest are the `web.cdr3fix.unmp` truthiness bug (7,973 rows) and a family
of pandas type coercions the all-string reader undoes (~55k cells). See `ROADMAP.md` §17.

## Next

1. **Phase 8** — `junction-nt` (#461), `segment-guess` (#462), `dgene`, one branch each. These fill
   the `chains` columns phase 6 declared but left absent (`cdr3nt`, `cdr3nt.pgen`, `cdr3nt.margin`,
   `d.start`, `d.end`) and the AIRR nucleotide fields phase 7 left empty. **AIRR `Receptor` rides
   along**: `vdjtools.model.stitch_*` produces the complete variable domain its two required columns
   need (`ROADMAP.md` §18).
2. **Phase 9** (`feature/harmonize-rules`) — **#327 lands here, and it gates the arda swap.**
3. `uv lock`, push `dev`, let `chunk-check` run once **before** applying branch protection.

## Blocked / needs a decision

| Item | Blocks | Question |
|---|---|---|
| `vdjmatch` `_zip_asset` patch | the first multi-zip release | needs a patch + release in `antigenomics/vdjmatch`. See `ROADMAP.md` §3.1 |
| `arda` release with `fix/source-root-marker` | CI without the `$ARDA_HOME` workaround | `d40095c` is committed on a local branch in `~/vcs/code/arda`, not pushed |
| Motif `coef` calibration | phase 11 | needs the motif pipeline running before it can be fitted |

## Known, not yet fixed

- **`record_id` is stable across builds, not across releases.** The registry is not written or
  committed; it reconciles against an empty one every build, so ids would shift the moment a chunk is
  added. It becomes a release asset in phase 14 — 72.7 MB is too much to commit per curation PR
  (`ROADMAP.md` §17). Do not lean on cross-release id stability until then.
- **854 records carry no `reference.id`**, all from `luciani-samir-etal-hcv-14-09-2018` (822 at score
  0, 32 at 1). The QC rule permits a blank one by construction. Phase 9 / #347.
- 14 records are reported by two chunks with the same PDB id in different letter case, and one
  `meta.epitope.id` carries a float. Both are phase 9 nomenclature work; both are visible in the
  ledger today.
- `hotfix` carries `compute_pdb_cdr3fix.py`, which phase 5 supersedes — the same lookup falls out
  of `arda.markup_batch`.
- `summary/MakeEmbedableHtml.py` and `processing/*.ipynb` carry 7 ruff findings. Phase 12 owns them;
  `src/vdjdb` and `tests/` are clean.
