# STATUS

_Last updated: 2026-09-25_

## In flight

**Phase 11 — `feature/motifs-tcremp`**. Phases 0-10 are merged to `dev`; the ledger reads PASS with
every difference declared and measured.

| Phase | Branch | State |
|---|---|---|
| 0 | `feature/dev-baseline` | merged — docs, package skeleton, chunk lint, CI, identity, guard |
| 1 | `feature/schema` | merged — the field registry every column order projects from |
| 2 | `feature/golden-harness` | merged — `vdjdb diff`, reproducible, validated against the release |
| 3 | `feature/io-qc` | merged — polars reader, vectorised QC, corpus normalised, CI gates `--strict` |
| 4 | `feature/pipeline-core` | merged — the definitive tables, legacy as a projection; `py_src/` retired |
| 5 | `feature/arda-cdr3fix` | merged, **behind `engine="legacy"`** — both gates cleared; the swap is now a decision, not a blocker (`ROADMAP.md` §28) |
| 6 | `feature/new-format` | merged — tables + evidence + `vdjdb.schema.json`; ledger PASS from the shipped tables |
| 7 | `feature/airr` | merged — Rearrangement + Reactivity; `airr.validate_rearrangement` passes on all 286,047 rows |
| 8a | `feature/junction-nt` | merged — 261,097 `cdr3nt`, 0 back-translation mismatches |
| 8b | `feature/dgene` | merged — D geometry from the junction scenario, confidence from `arda.dpost` |
| 8c | `feature/segment-guess` | merged — the legacy V guesser has never worked; Pgen fills 686 of 711 gaps |
| 8d | `feature/receptor` | merged — AIRR `Receptor`, 81,003 rows; stitching vectorised to 2.7 s |
| 9 | `feature/harmonize-rules` | **#389, #327, #467, #347, #368 landed; #564 already fixed; #561 advisory. Ships `epitopes` + `restriction`. Ledger PASS** |
| 10 | `feature/motifs-tcrnet` | **merged** — `vdjdb motifs` ships both legacy files; the 137 logo-less cids and the 1.73 % deleted letter mass are gone (`ROADMAP.md` §29) |

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
| Wall time | 344 s | **185 s** (16 s without model inference) |
| Peak RSS | 2.16 GB | ~1 GB |
| Output | three files, assembled directly | three tidy tables + a joined view, with legacy and AIRR projected from them |

Junction-nucleotide inference (#461) is 90 % of the wall time and is **never cached** (hard rule 9):
261,097 of 286,047 chains get a `cdr3nt`, all of which back-translate to the junction they came from.

`records` 192,753 × 33 · `chains` 286,047 × 17 · `evidence` 53,913 × 8 · `vdjdb` (view) 286,047 × 54.
All five legacy members are byte-identical whether projected from memory or read back from parquet.
Every difference from the 2026-06-03 release is a declared rule in `rules/expected_diffs.toml` firing
its exact measured count — largest are the `web.cdr3fix.unmp` truthiness bug (7,973 rows) and a family
of pandas type coercions the all-string reader undoes (~55k cells). See `ROADMAP.md` §17.

## Next

1. **Phase 11, TCREMP** — chunked embedding, per-epitope DBSCAN with a chain-global eps, `coef`
   fitted against the independent-study signal phase 6 already ships. **Never against TCRvdb.**
2. **Optimise the motifs** (`ROADMAP.md` §30) — the paratope motifs are a critical part of VDJdb and
   neither method is tuned. In priority order: the **giant component** (our largest TCRNET cluster
   holds 19,908 of 45,095 clustered records, 44 %; REDCEA's largest holds 3.9 %, so Leiden over
   connected components is the single biggest win), **background re-sampling and size** (`M` is not
   even uniform across chains today), **in-silico backgrounds** from `vdjtools.model` (which would
   also dissolve the CC-BY-NC-ND licence constraint), and the unswept parameters of both methods.
   Phase 10's open human-TRA coverage item lands here.
3. **Decide whether the legacy build should keep the 711 chains with no V** now that #462 can supply
   one for 686 of them (`ROADMAP.md` §21).
4. **Decide the cdr3fix engine swap.** Both `ROADMAP.md` §16 and §28 gates are cleared: the family-
   call defect is fixed and released as **arda-mapper 2.29.0** (bound in `pyproject.toml`), which
   recovers **4,033 beta V-end mappings** (legacy-only 7,969 → 3,936) and puts arda ahead on all four
   coverage measures. The 3,936 that remain are 2,207 ambiguous curation calls arda rightly refuses
   plus 1,218 alleles with no shipped anchor — a curation policy question, not an engine one.
5. **Phase 9e** — validate the epitope catalogue with `mhcmatch` (`ROADMAP.md` §12).
6. `uv lock`, push `dev`, let `chunk-check` run once **before** applying branch protection.

## Blocked / needs a decision

| Item | Blocks | Question |
|---|---|---|
| Murine class I: `H2-Db` or `H-2Db`? | the murine MHC vocabulary | the data says `H2-Db` (6,207 records, `H-2Db` appears zero times), `proofreading/mhc.md` says `H-2Db`, and `vdjdb-web` carries a repair for the pair. `ROADMAP.md` §25 |
| `vdjmatch` `_zip_asset` patch | the first multi-zip release | needs a patch + release in `antigenomics/vdjmatch`. See `ROADMAP.md` §3.1 |
| `arda` release with `fix/source-root-marker` | CI without the `$ARDA_HOME` workaround | `d40095c` is committed on a local branch in `~/vcs/code/arda`, not pushed. **Unrelated to the 2.29.0 family-call fix, which has shipped** |
| Motif `coef` calibration | phase 11 | needs the motif pipeline running before it can be fitted |

## Known, not yet fixed

- **`record_id` is stable across builds, not across releases.** The registry is not written or
  committed; it reconciles against an empty one every build, so ids would shift the moment a chunk is
  added. It becomes a release asset in phase 14 — 72.7 MB is too much to commit per curation PR
  (`ROADMAP.md` §17). Do not lean on cross-release id stability until then.
- **854 records carry no `reference.id`**, all from `luciani-samir-etal-hcv-14-09-2018` (822 at score
  0, 32 at 1). The QC rule permits a blank one by construction, and #347 cannot help: there is
  nothing to resolve. A curation question.
- **`HLA-A*08:01` on 74 records.** No HLA-A\*08 locus exists at any resolution, so the call is
  certainly wrong — but the chunk carries no `reference.id` to resolve it against. Curation.
- **8 epitopes are listed twice in the antigen patch with different answers**, resolved silently by
  file order (`curate.patch.CONFLICTING_EPITOPES`). Two are nomenclature aliases, three are the HIV-1
  Gag-Pol frameshift, the rest need a curator.
- **62 of 2,132 epitopes disagree with their own MHC class on length** — 20 MHCI longer than 11
  residues, 42 MHCII shorter than 12 (`ROADMAP.md` §27).
- **99 records carry the same CDR3 on both chains** (#561), 98 from two references. Reported as an
  advisory QC finding on every run; only a curator can say which chain is wrong.
- **9 `antigen.gene` values are protein names rather than gene symbols** (`Trans-sialidase`, 284
  records; `Nucleocapsid`, 171) and **15 `antigen.species` values are outside
  `proofreading/species_aliases.tsv`** (`SIV`, 1,771). Both are vocabulary gaps for the curation
  skills, not mechanical rules (`ROADMAP.md` §26).
- 14 records are reported by two chunks with the same PDB id in different letter case, and one
  `meta.epitope.id` carries a float. Both are phase 9 nomenclature work; both are visible in the
  ledger today.
- `hotfix` carries `compute_pdb_cdr3fix.py`, which phase 5 supersedes — the same lookup falls out
  of `arda.markup_batch`.
- `summary/MakeEmbedableHtml.py` and `processing/*.ipynb` carry 7 ruff findings. Phase 12 owns them;
  `src/vdjdb` and `tests/` are clean.
