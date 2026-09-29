# Junction anchors and their repair

> **Authority:** `src/vdjdb/curate/anchors.py`, reading arda's germline reference.
> **Read it with:** `vdjdb submission <chunk>`, or `out/reports/anchors.tsv` after a build.
> **Prose:** [`docs/submission.md`](../docs/submission.md), "Junctions that contradict their own
> germline".

## The rule

VDJdb's `cdr3` is **junction space**: Cys104 through Phe/Trp118, **both anchors included**. A junction
is consistent when its first and last residues are the residues the germline of the V and J the record
names actually encodes.

**The anchor is read from that germline. It is never assumed to be Phe or Trp.** Human `TRAJ35*01`
templates `IGFGNVLHC` and mouse `TRAJ7*01` templates `DYSNNRLTL`, so a correct junction on either ends
in Cys or Leu. Measured over the corpus, a fixed "starts with C, ends with F or W" test calls **481
correct chains broken** and separately **misses 125 that are not** consistent. Earlier revisions of this
file stated that fixed rule and derived its anchors from `chunks/` itself, which is the data validating
itself; both are gone.

`vdjdb qc` does not check this. Its sequence rules check the residue alphabet and a minimum length,
which is why a submission exported in IMGT CDR3 space - short one residue at each end - passes them.

## The defects, as reported

One per end, named for the end and the reason:

| Defect | What the sequence has | Repair proposed |
|---|---|---|
| `V absent anchor` / `J absent anchor` | the anchor residue is missing | the germline residue, added |
| `V corrupt anchor` / `J corrupt anchor` | an anchor is present but is not the germline's | the germline residue, substituted |
| `V under-trimmed` / `J under-trimmed` | framework retained past the anchor | trimmed to the anchor |
| `V allele mismatch` / `J allele mismatch` | the sequence matches a **functional sibling allele** of the named gene | the **call**, not the sequence |
| `V unanchored` / `J unanchored` | no germline to compare against - no call, or a call the reference lacks | none; the universal anchor is all that is left |
| `V unexplained` / `J unexplained` | disagrees with the germline and no rule accounts for it | none |

Three germline residues must agree before a sequence repair is proposed, so a repair is an alignment
result rather than a guess at one letter.

## An allele mismatch is evidence about the call

An ORF or pseudogene allele has a non-canonical anchor by definition, and arda records that faithfully:
of 383 J entries over four organisms only 14 have a `templated_aa` not ending in Phe or Trp, and 13 of
those are marked `ORF` or `P`. So a junction disagreeing with a non-functional allele is evidence that
the **call** is wrong, not the sequence.

95 mouse chains name `TRAJ47`, which resolves to the ORF `*01` (`HYANKMIC`), and every one reads
`DYANKMIF` - exactly `TRAJ47*02`, the functional allele. No corpus record reads the `*01` signature.
This is the same defect as [#327](https://github.com/antigenomics/vdjdb-db/issues/327), where 66 % of
explicit `TRAJ24*01` calls carry the `*02` motif, and rewriting the sequence there would destroy the
evidence for it.

A missing Cys104 deserves a second look even where the rest of the record is clean: no TCR folds
without it, so a first residue that is not Cys while the body still aligns to the V is a sequencing or
transcription error rather than a variant.

## Why the proposal is against the submitted sequence

`arda.cdr3fix` repairs most of these on the way through the build, so the shipped `cdr3` is usually
right and **the chunk keeps the wrong sequence** - the next export of that data is wrong again. The
repair is therefore computed against the value the chunk holds, not the value that ships:
`YLCSSQEGGYGYTFGSG` ships as `YLCSSQEGGYGYTF`, framework trimmed behind the anchor and kept in front
of it.

## Applying one

A repair is a chunk edit: its own branch, its own issue, and a message naming the files, the rows and
the reason. Nothing applies one automatically, and neither should you in bulk -
[`/vdjdb-proofread`](../skills/vdjdb-proofread/SKILL.md) step 3 and
[`/vdjdb-publish`](../skills/vdjdb-publish/SKILL.md) have the procedure.
