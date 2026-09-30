# Migrating the build onto arda 2.34.0 + vdjtools 4.8.0

2026-09-30. **Everything the build does to annotate a junction is now one call in a library.** This
note says which call, what it returns, and which of this repository's modules it replaces. The
libraries' own design note is `docs/junction_pipeline.md` in `antigenomics/vdjtools`.

## Why

Five of this build's modules re-implement, coordinate, or second-guess work the libraries do:
`annotate/cdr3fix.py` (199 lines), `annotate/junction.py` (203), `annotate/segments.py` (185),
`annotate/dgene.py` (62) and `curate/anchors.py` (388) — **1,037 lines** that call arda and vdjtools
stage by stage, hold the coordinate conversions between them, and in `curate/anchors.py` compute a
junction repair that is then **reported and never applied** (#711). One library call replaces the
annotation part of all five. What stays here is *curation*: which records to flag, and what a curator
does about a contradicted call.

`annotate/segments.py` exists because "`arda.cdr3fix` repairs a junction against a *named* germline and
never proposes one". **arda 2.34.0 proposes one**, which is why the floor below is 2.34.0 and not
2.33.0.

## The one call

```python
from vdjtools.model import annotate_junctions

out = annotate_junctions(
    keys["cdr3"].to_list(),          # junction space: Cys104 .. Phe/Trp118 inclusive
    keys["v"].to_list(),
    keys["j"].to_list(),
    species=keys["species"].to_list(),   # per row; VDJdb names are accepted as-is
)
```

**One row out per row in, in input order**, so it joins positionally. A record the model cannot
explain is present with nulls — never dropped, never an exception. Deduplicate to distinct
`(species, cdr3, v, j)` keys first, as the build already does.

### What comes back

| column | replaces | notes |
|---|---|---|
| `cdr3_repaired` | `annotate/cdr3fix.py` | the repaired junction; `cdr3_aa` is the submission |
| `v_call`, `j_call` | `annotate/cdr3fix.py` | **confirmed or re-called**; feeds `cdr3fix.vId`/`jId` |
| `v_alts`, `j_alts` | — | every allele the junction cannot separate, chosen one first |
| `proposed` | `annotate/segments.py` | which side the submission left blank and the junction named |
| `v_end`, `j_start` | `annotate/cdr3fix.py` | residues; VDJdb's `vEnd` / `jStart` |
| `v_end_nt`, `j_start_nt` | `annotate/junction.py` | nucleotides — no `ceil(nt/3)` conversion here any more |
| `v_flags`, `j_flags`, `good` | `curate/anchors.py` | `mismatch` is the curator's list; `impossible` is a malformed junction |
| `cdr3_nt`, `pgen` | `annotate/junction.py` | the inferred nucleotide junction and its Pgen |
| `d_call`, `d_posterior` | `annotate/dgene.py` | **use as `d.inferred`** — the D GENE, named by the model, with its posterior |
| `d_start_nt`, `d_end_nt` | `annotate/junction.py` | where that gene was placed, 1-based closed in junction space |
| `d_start_aa`, `d_end_aa` | — | the residues whose codons the D touches, recomputed from the nt bounds |
| `np1`, `np2` | `annotate/junction.py` | the N regions either side of the D, sliced from the same bounds |

## Four things that change in the output, and why

**1. The junction repair is applied, not reported (#711).** `curate/anchors.py` computes a repair and
writes it to `out/reports/anchors.tsv`; 261 of 1,037 flagged chains carry a germline-supported repair
and every one ships unrepaired. `cdr3_repaired` is the repair, already gated: measured over all
187,488 keys of the 2026-06-03 release, arda agrees with that release's shipped junction on
**99.8352 %**, reproduces 4,331 of its 4,499 repairs, and ships **677** non-canonical junctions
against the release's 715. Applying `cdr3_repaired` is not the ungated "apply the proposal" that
#711 warns against — the 2,842 rewrites the release never ships are down to 141, of which 35 touch an
already-canonical junction.

**2. `good` means "the junction is well formed", and a germline disagreement no longer contradicts
it.** A residue differing from germline *inside* the templated run is flagged `mismatch` and left
exactly as submitted, because a curation error and an allele IMGT does not record are
indistinguishable from one junction. Treat `mismatch` as a proofreading queue (it is the #681 class),
and `impossible` as the record whose anchors cannot be satisfied.

**3. `d.inferred` comes from `d_call`, and every row that can be drawn has coordinates.** Naming the
D and placing it are separate questions and one estimator answers each: the recombination model names
the gene, the aligner places it greedily and ungated. Measured on 4,000 real human TRB rearrangements
whose D and `DStart`/`DEnd` come from the nucleotide sequence:

| | gene right, of all rows | has coordinates |
|---|---:|---:|
| E-value-gated alignment chooses and places | 47.93 % | 55.75 % |
| gated alignment, model posterior where it declines | 71.40 % | 55.75 % |
| today's `d.inferred` (the length-and-prior posterior) | 69.67 % | — |
| **model names, greedy alignment places** | **74.30 %** | **99.70 %** |

⛔ **`arda.dpost` does not come back, in either library.** The posterior this build reads today is
dominated by a group-by over a call the pipeline already makes — 69.67 % against 74.33 %, 56.3 µs per
key against 30.1 — and it needs a fitted prior table and a per-locus tempering constant that the
replacement does not. `annotate/dgene.py` is replaced by two columns, not re-pointed at a new module.
So `fields.py`'s comment about `d.posterior` needs rewriting rather than renaming: the number now
comes from the model's own scenario weights, normalised over D genes.

⚠ The one thing that gets *worse* is per-row positional precision: `d_start_nt` is exact on 60.88 %
of correctly-called rows against 66.67 % under the gate. It is exact on **1,807 rows rather than
1,278**, because it answers 3,988 rather than 2,230. For a view drawing V/N/D/N/J that is the trade
to take.

**4. The nucleotides are the authority and the amino-acid bounds are recomputed from them.**
`v_end_nt`, `j_start_nt`, `d_start_nt` and `d_end_nt` are all read off the inferred nucleotide
junction; `v_end`, `j_start`, `d_start_aa` and `d_end_aa` follow from those. A view showing both
alphabets therefore cannot draw them disagreeing, and there is no `ceil(nt/3)` conversion left in
this repository.

## Cost

**~280 µs per distinct key**, so under a minute for a 190,000-key corpus in one process. The
nucleotide inference dominates; naming the D costs ~30 µs and placing it ~1.4 µs. Compare what it
replaces: `posterior_d` alone was 15.96 s over 119,034 keys via a Python row loop (35.8 % of the
build), the nucleotide stage ran as four `vdjdb infer-nt` processes over contiguous slices, and
`fix_cdr3` was 7.71 s.

⛔ Do **not** wrap the call in a pool. Every stage is already batched: one `markup_batch`, then one
native threaded `infer_nt_batch` and one `best_aa_scenarios_batch` per `(organism, locus)` on a model
loaded once, then a thin Python loop over arda's C++ aligner. A pool around it would re-import the
libraries and re-load the models per worker.

## Version floors

```toml
"arda-mapper>=2.34.0",   # cdr3fix repair policy, v_alts/j_alts, blank-call proposal, map_d_junction(v_end=, j_start=)
"vdjtools>=4.8.0",       # annotate_junctions
```

`arda.dpost` is **gone** in 2.33.0 and does **not** reappear in vdjtools, so `annotate/dgene.py`'s
`from arda.dpost import posterior_d` breaks on upgrade — that is the intended failure, not a
surprise, and the fix is to read `d_call` / `d_posterior` rather than to re-point the import.
`arda markup --d-posterior` and `--d-prior` are gone with it.

## Suggested order

1. Bump both floors and **delete** `annotate/dgene.py`, reading `d_call` / `d_posterior` from the one
   call instead. The build runs again at this point.
2. Replace `annotate/junction.py`'s `infer` with the `cdr3_nt` / `pgen` / `d_*` columns. Its
   `SPECIES` map and its Pgen provenance comment are worth keeping as documentation of what the
   numbers mean.
3. Replace `annotate/cdr3fix.py` with the stage-1 columns. Keep whatever maps them onto VDJdb's
   `cdr3fix` JSON key names — `Cdr3Markup.to_cdr3fix()` in arda still emits that object key-for-key.
4. Delete `annotate/segments.py` and read `proposed` instead. arda 2.34.0 resolves **2,532 of the
   3,130 blank-call keys** (2,504 of them `good`, with both boundaries placed); the 598 that stay
   refused name neither side, so no locus exists to propose within and the module could not have
   answered them either. An *unresolvable* call — `TRBVnope*01` — is still refused rather than
   proposed for, deliberately: a submission that names something wrong is a defect for a curator, not
   a gap for the junction to fill.
5. **Then** #711: `curate/anchors.py` keeps its classification (which is curation) and drops its
   repair computation (which is annotation), and the build ships `cdr3_repaired`.
6. Re-run `vdjdb diff` against `reference.zip` and read the junction-column changes against the table
   in §"Four things that change" above.
