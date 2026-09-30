# Migrating the build onto arda 2.33.0 + vdjtools 4.8.0

2026-09-30. **Everything the build does to annotate a junction is now one call in a library.** This
note says which call, what it returns, and which of this repository's modules it replaces. The
libraries' own design note is `docs/junction_pipeline.md` in `antigenomics/vdjtools`.

## Why

Four of this build's modules re-implement, coordinate, or second-guess work the libraries do:
`annotate/cdr3fix.py` (199 lines), `annotate/junction.py` (203), `annotate/dgene.py` (62) and
`curate/anchors.py` (388) — 852 lines that call arda and vdjtools stage by stage, hold the coordinate
conversions between them, and in `curate/anchors.py` compute a junction repair that is then
**reported and never applied** (#711). One library call replaces the annotation part of all four. What
stays here is *curation*: which records to flag, and what a curator does about a contradicted call.

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
| `v_end`, `j_start` | `annotate/cdr3fix.py` | residues; VDJdb's `vEnd` / `jStart` |
| `v_end_nt`, `j_start_nt` | `annotate/junction.py` | nucleotides — no `ceil(nt/3)` conversion here any more |
| `v_flags`, `j_flags`, `good` | `curate/anchors.py` | `mismatch` is the curator's list; `impossible` is a malformed junction |
| `cdr3_nt`, `pgen` | `annotate/junction.py` | the inferred nucleotide junction and its Pgen |
| `d_call`, `d_start_nt`, `d_end_nt` | `annotate/junction.py` | D by alignment on those nucleotides |
| `d_start_aa`, `d_end_aa` | — | the residues whose codons the D touches |
| `d_posterior_call`, `d_posterior`, `d_entropy` | `annotate/dgene.py` | `arda.dpost` moved to `vdjtools.model` |
| `d_best`, `d_best_source` | — | **use this as `d.inferred`**: the alignment where it speaks, the posterior where it declines |

## Three things that change in the output, and why

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

**3. `d.inferred` should come from `d_best`, not from either route alone.** Measured on 4,000 real
human TRB rearrangements whose D is called from the nucleotide sequence:

| | called | correct when called | correct over all rows |
|---|---|---|---|
| alignment on the inferred nt | 56.2 % | 85.98 % | 48.30 % |
| the posterior alone (today's `d.inferred`) | 100 % | 69.67 % | 69.67 % |
| `d_best` | 100 % | — | **71.53 %** |

`d.posterior` keeps its meaning and keeps coming from the same model, so the
`fields.py` comment about it needs only its module name updated: `arda.dpost` →
`vdjtools.model.posterior_d`.

## Cost

307 µs per distinct key, so **~58 s for a 190,000-key corpus** in one process, of which the
nucleotide inference is 69 % and the D alignment 2 %. Compare what it replaces: `posterior_d` alone
was 15.96 s over 119,034 keys via a Python row loop (35.8 % of the build), the nucleotide stage ran as
four `vdjdb infer-nt` processes over contiguous slices, and `fix_cdr3` was 7.71 s. Do **not** wrap the
call in a pool: stage 2 already threads across the batch.

## Version floors

```toml
"arda-mapper>=2.33.0",   # cdr3fix repair policy, v_alts/j_alts, map_d_junction(v_end=, j_start=)
"vdjtools>=4.8.0",       # annotate_junctions, posterior_d_batch
```

`arda.dpost` is **gone** in 2.33.0, so `annotate/dgene.py`'s `from arda.dpost import posterior_d`
breaks on upgrade — that is the intended failure, not a surprise. `arda markup --d-posterior` and
`--d-prior` are gone with it.

## Suggested order

1. Bump both floors, and replace `annotate/dgene.py`'s body with the `d_posterior_*` columns the one
   call already returns. The build runs again at this point.
2. Replace `annotate/junction.py`'s `infer` with the `cdr3_nt` / `pgen` / `d_*` columns. Its
   `SPECIES` map and its Pgen provenance comment are worth keeping as documentation of what the
   numbers mean.
3. Replace `annotate/cdr3fix.py` with the stage-1 columns. Keep whatever maps them onto VDJdb's
   `cdr3fix` JSON key names — `Cdr3Markup.to_cdr3fix()` in arda still emits that object key-for-key.
4. **Then** #711: `curate/anchors.py` keeps its classification (which is curation) and drops its
   repair computation (which is annotation), and the build ships `cdr3_repaired`.
5. Re-run `vdjdb diff` against `reference.zip` and read the junction-column changes against the table
   in §"Three things that change" above.
