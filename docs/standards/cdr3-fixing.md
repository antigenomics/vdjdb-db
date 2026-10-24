# CDR3 markup and repair

VDJdb markup describes the junction sequence, its conserved anchors, the assigned V/D/J
segments and their boundaries. Repair normalises sequences reported with missing anchors
or extra framework residues. Nucleotide reconstruction is a separate annotation step.
This page explains how to read these results together; the [column reference](columns.md)
lists their declared types and the [output reference](../outputs.md) describes the files.

## Junction, CDR3 and anchors

**VDJdb's `cdr3` is an amino-acid junction**, from conserved Cys104 through Phe/Trp118,
including both anchors. AIRR and arda call this `junction_aa`; their `cdr3_aa` excludes
both anchors. Removing the anchors shifts every coordinate and shortens the sequence
by two residues. Do not use a VDJdb boundary to slice an anchor-excluded CDR3.

Names below are the exported tidy column names, as in `chains.parquet` and the joined
view. Internal names use dots, for example `v.end.nt` instead of `v_end_nt`.
The legacy releases retain their existing names and coordinate conventions.

## Submitted sequences and repairs

`cdr3` is the sequence that ships after repair. `cdr3_original` preserves the submitted
sequence; when an older report leaves it empty, use `cdr3` as the unchanged submission.
Record identifiers are assigned before repair, so separately reported observations do
not merge merely because their repaired junctions become identical.

The build calls arda markup on distinct `(species, cdr3, v, j)` combinations and joins
the results back to the records. The repair policy is:

- Restore missing anchors and adjoining germline residues when the alignment supports them.
- Trim framework residues outside the junction, such as the extra `GXG` after its terminal `F`.
- Allow a supported terminal-anchor correction within an aligned germline run.
- Retain internal germline mismatches rather than substituting the reference residue.

The configured `max_replace=0` restricts replacement of junction residues before a
supporting germline alignment. It does not disable restoration of an anchor within that
alignment. Consequently `FixReplace` remains a possible report, but is not permission
to rewrite internal junction residues.

Sequence changes and annotation quality are different quantities:

| Field | Interpretation |
|---|---|
| `fix_needed` | The shipped junction differs from the submitted junction. It records a change already made. |
| `fix_good` | The markup engine found neither side impossible under its repair rules. It is not a switch controlling whether the repaired sequence ships. |
| `v_fix_type`, `j_fix_type` | A summary label for each side, described below. |
| `v_canonical` | The shipped junction starts with literal `C`. |
| `j_canonical` | Engine markup checks a terminal `F` or `W`; after J recall, this checks the selected allele's terminal residue instead. |
| `v_canonical_submitted`, `j_canonical_submitted` | Literal `C` and `F/W` anchor checks on the submitted sequence. |
| `cdr3_one_cysteine` | No additional cysteine occurs after the first residue. |

For example, a J-side extension can ship while the V side cannot be annotated:
`fix_needed=true` and `fix_good=false` can therefore occur together. Conversely an
unchanged sequence can have `fix_good=true`. The anchor checks are not a substitute for `fix_good`: an allele with a
nonstandard anchor can satisfy the engine's allele-specific rules while failing a
literal `C` or `F/W` check. J recall also changes which anchor test `j_canonical` uses.

## Segment calls and their sources

The sequence used for markup and the evidence for a segment call must remain distinguishable:

| Fields | Source and use |
|---|---|
| `v_segm_submitted`, `d_segm_submitted`, `j_segm_submitted` | Calls preserved from the submission. |
| `v_segm_arda`, `j_segm_arda` | Calls resolved by the markup engine. |
| `v_segm`, `j_segm` | Calls that ship after harmonisation, allele resolution and the supported J-recall rule. If resolution fails, a reported call is retained. |
| `d_segm` | Curated D call, separate from the inferred D. |
| `v_inferred`, `j_inferred` | Proposals only where the submission left that segment blank. Both proposals ship in `v_segm` and `j_segm`, with inferred provenance retained. |
| `d_inferred`, `d_posterior` | D gene proposed by the recombination model and its gene posterior probability. |

J recall can replace an unsupported J call using the junction's terminal germline
match. It updates the shipped call, `j_start` and `j_canonical`, without editing the
junction or replacing the engine's `j_fix_type` and `fix_good` verdicts. The
[output reference](../outputs.md)
specifies the recall criteria, proposal policy and their measured validation.

## Explicit nucleotide and amino-acid boundaries

Use the eight fields below for a consistent V/N/D/N/J drawing or sequence slicing.
They all use **0-based, half-open junction coordinates**: the start is included and
the end is excluded. Nucleotide coordinates index `cdr3nt`; amino-acid coordinates
index the shipped `cdr3`, including both anchors.

| Segment | Nucleotide span | Amino-acid span |
|---|---|---|
| V | `[0, v_end_nt)` | `[0, v_end_aa)` |
| D | `[d_start_nt, d_end_nt)` | `[d_start_aa, d_end_aa)` |
| J | `[j_start_nt, len(cdr3nt))` | `[j_start_aa, len(cdr3))` |

V/J nucleotide offsets come directly from `vdjtools.model.annotate_junctions`.
D placement comes from that same batch call; its library interval is converted
from 1-based closed to 0-based half-open before publication.

For a nucleotide span `[s, e)`, the amino-acid span is
`[s // 3, (e + 2) // 3)`. Starts round down and ends round up, including every
codon touched. Adjacent segments can therefore share an amino acid. Multiplying
an amino-acid boundary by three cannot recover a boundary inside a codon.

For an illustrative 39-nt, 13-aa junction with V end 4, D span `[5, 7)` and J start 8:

| Region | Nucleotide slice | Amino-acid span touched |
|---|---|---|
| V | `cdr3nt[0:4]` | `[0, 2)` |
| N1 | `cdr3nt[4:5]` | `[1, 2)` |
| D | `cdr3nt[5:7]` | `[1, 3)` |
| N2 | `cdr3nt[7:8]` | `[2, 3)` |
| J | `cdr3nt[8:39]` | `[2, 13)` |

This is a coordinate example, not a reconstructed receptor. N1 and N2 here are the
gaps between the placed germline spans; their widths do not measure how many
nucleotides were inserted during recombination. When D is unavailable, a V-to-J
gap cannot be partitioned into N1, D and N2 from these fields.

All eight fields are **null when unavailable**, including all four D fields on TRA
chains. They are also null if nucleotide reconstruction is absent or fails to
translate to the shipped junction. Test for null before slicing; zero is a valid
boundary. Nullable integer fields use null rather than the empty-string marker
used for missing text.

## Nucleotide reconstruction and D placement

`cdr3nt` is an inferred nucleotide junction, not an observed sequence from the paper.
The recombination model chooses a plausible nucleotide history conditional on the
amino-acid junction and usable segment calls. The automatic model selection uses
bundled OLGA-derived models first, then arda-derived models; a germline scaffold
supports species without a fitted model.

`cdr3nt_pgen` is the reconstruction-stage model score returned as `pgen`;
`cdr3nt_margin` is its ratio to the runner-up score. A ratio near one indicates a
close choice at that stage; a scored winner without a nonzero runner-up has an
infinite margin. The library can subsequently replace germline flanks without
recalculating these scores, and model preparation can adjust recombination
probabilities. Treat them as model scores, not calibrated probabilities or
uncertainty estimates for the final shipped nucleotide string. Synonymous histories
can translate to exactly the same junction, so a valid translation does not
establish the original nucleotides.

The model names the inferred D gene, while a separate alignment places that D within
the inferred V-to-J interval. `d_posterior` concerns the gene assignment, not the
probability that a boundary is exact. D coordinates describe `d_inferred`, which can
differ from the publication's `d_segm`. See the [output reference](../outputs.md)
for the recorded coverage and accuracy measurements.

## Compatibility fields and fallback boundaries

The existing fields remain unchanged so legacy consumers continue to work:

| Fields | Units, source and missing value |
|---|---|
| `v_end`, `j_start` | Protein-alignment boundaries in 0-based junction-aa space, with `-1` for unavailable. V end is exclusive; J start is included. |
| `v_end_inferred`, `j_start_inferred` | Independent germline-alignment fallback, populated only where the corresponding protein boundary is `-1`; otherwise also `-1`. These remain legacy aa boundaries. |
| `d_start`, `d_end` | Aliases of `d_start_nt`, `d_end_nt`, including their nulls. |

`v_end_aa` and `j_start_aa` are derived from nucleotide spans. They need not equal
`v_end` and `j_start`, which answer the protein-alignment question. The fallback
pair is not silently coalesced into either source. In particular, its legacy
conversion rounds a J nucleotide boundary up, whereas `j_start_aa` rounds down to
include a partially templated codon. Do not compare these as interchangeable values.

Legacy positional motif columns and `cdr3fix` reports keep their historical fields.
The explicit coordinates are available in the definitive tables and joined view.
The legacy `cdr3fix` field map is:

| Report key | Tidy field |
|---|---|
| `cdr3`, `cdr3_old` | `cdr3`, `cdr3_original` |
| `fixNeeded`, `good` | `fix_needed`, `fix_good` |
| `vId`, `jId` | `v_segm`, `j_segm` |
| `vEnd`, `jStart` | `v_end`, `j_start` |
| `vFixType`, `jFixType` | `v_fix_type`, `j_fix_type` |
| `vCanonical`, `jCanonical` | `v_canonical`, `j_canonical` |

`vdjdb.txt` serialises the report as JSON. The paired `cdr3fix.alpha` and
`cdr3fix.beta` cells in `vdjdb_full.txt` retain their historical dictionary-string
representation. Prefer the typed flat columns when implementing a new consumer.

## Repair labels

A side's label summarises the engine result; multiple observations can collapse into
one label. It is not a complete edit transcript or an independent quality score.

| Label | Meaning |
|---|---|
| `NoFixNeeded` | No sequence edit reported for this side. Consult `fix_good` and the anchor flags separately. |
| `FixAdd` | Germline-supported residues added. |
| `FixTrim` | Residues outside the junction removed. |
| `FixReplace` | A supported terminal correction applied under the repair policy. |
| `TruncatedGermline` | Boundary supported by an incomplete germline reference; its extent is a lower bound. |
| `FailedBadSegment` | Segment entry cannot support the requested markup. |
| `FailedNoAlignment` | No usable germline alignment found. |
| `FailedReplace` | Required replacement is not permitted by the repair rules. |

`TruncatedGermline` preserves a supported boundary when the IMGT reference stops
inside the anchor region. Residues beyond that supported extent are unattributed,
not established as non-templated. The repair can still be `good`. In the recorded
comparison introducing this label (arda issue 135), 1,144 chains changed from
`FailedBadSegment` with `vEnd=-1`: 1,132 beta and 12 alpha, including 850 assigned
`TRBV11-2*02`. These are counts from that comparison, not a current-build census.

## Interpretation limits

- Literal anchors and a successful repair do not establish a segment call, an exact
  recombination history or experimental specificity. Internal mismatches can remain.
- Nucleotide-derived boundaries describe the selected inferred sequence. They are
  not measurements from an experimentally sequenced rearrangement.
- A null coordinate, or a legacy `-1`, means unavailable. It does not establish a
  zero-length segment or the absence of a germline contribution.
- Amino-acid spans include partial codons and can overlap. Use nucleotide spans
  when drawing fine boundaries; do not force amino-acid spans into a partition.
- AIRR sequence coordinates also require the junction's position in the full
  sequence. See [AIRR mapping](airr-mapping.md) before converting coordinates.
