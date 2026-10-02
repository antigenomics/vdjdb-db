---
name: vdjdb-proofread
description: Validate a VDJdb chunk with vdjdb qc and vdjdb submission, explain every finding with a specific suggested fix and the authority behind it, resolve method-field and MHC problems the rules cannot decide, read the junction-versus-germline and score reports, and decide whether the chunk lands in chunks/, pending/ or withheld/. Use as the last gate before a chunk is committed, or on any chunk suspected of having a column shift, bad gene or allele names, blank MHC fields or missing method information.
---

# vdjdb-proofread

Validate a chunk, explain every finding, fix what is mechanical, escalate what is judgement, and
decide where the chunk goes. Last of three stages: [extract](../vdjdb-extract/SKILL.md) →
[format](../vdjdb-format/SKILL.md) → **proofread**.

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first.

## Invocation

```
/vdjdb-proofread [path-to-tsv]
```

## Required completeness and duplicate review

Apply [completeness and observation identity](../../docs/standards/chunk-format.md#completeness-and-observation-identity)
to every chunk. At least one junction and the epitope must be reported; neither is imputable.
Strict QC and direct assembly reject missing sequences or incomplete MHC fields. Resolve missing
MHC partners using paper evidence or a documented supported mhcmatch inference before submission.
Inspect built-in segment proposals for missing V/J calls without inventing paper-reported calls.

Audit both-chain receptor-pMHC matches within and across chunks, then compare all observation
metadata. Different donors, references, methods, subsets or other metadata are independent
observations. Record the reviewed groups and decisions; a matching receptor alone is not a duplicate.

## Step 1 - structure, before anything parses it

`vdjdb qc` lints the text before it reads the table, and reports `bom`, `encoding`, `crlf`, `empty`,
`no-data-rows`, `empty-column-name`, `prose-column-name`, `duplicate-column-name`,
`unknown-column`, `missing-required-column`, `reference-id-form`, `missing-file` and `unreadable`.
A header it cannot map at all arrives as `unreadable` with the parser's own message.

One defect it cannot see is a **column shift with matching counts**: a file whose header lost or
gained a column in a way that leaves the row width consistent, so every field lands under the wrong
name and `antigen.species` holds a submitter name. Check content against position:

```python
import csv, re
AA = re.compile(r'^[ARNDCQEGHILKMFPSTWYV]{4,}$')
SPECIES = {'HomoSapiens', 'MusMusculus', 'RattusNorvegicus', 'MacacaMulatta'}

with open(path) as f:
    header = next(csv.reader(f, delimiter='\t'))
    if header[0] != 'chunk.id':
        print(f'first column is {header[0]!r}, not chunk.id')
    f.seek(0)
    for i, row in enumerate(csv.DictReader(f, delimiter='\t'), start=2):
        if (ep := row.get('antigen.epitope', '')) and not AA.match(ep):
            print(f'row {i}: antigen.epitope {ep!r} is not a peptide')
        if (sp := row.get('species', '')) and sp not in SPECIES:
            print(f'row {i}: species {sp!r}')
        if i > 11:
            break
```

If a shift is confirmed, say which direction and by how many columns, then either drop the spurious
header columns or prepend the missing `chunk.id` to the data rows - whichever makes header and data
agree with a shipping chunk (`head -1 chunks/PMID_28423320.tsv`). Re-run step 1 after.

## Step 2 - run the rules

```bash
uv run vdjdb qc <file> --report out/reports/qc.tsv
```

24 row rules over one vectorised pass. `out/reports/qc.tsv` is one row per finding;
`out/reports/qc-summary.tsv` is one row per rule with an `advisory` flag. Under `--strict` (the
default) any non-advisory finding exits 1.

**Fatal rules and what each one means:**

| Rule | Suggested fix |
|---|---|
| `bad cdr3.alpha`, `bad cdr3.beta`, `bad antigen.epitope` | a character outside the 20 canonical amino acids, or 3 residues or fewer. Drop the row and log it - this is a transcription artefact, not a sequence |
| `bad v.alpha`, `bad j.alpha`, `bad v.beta`, `bad d.beta`, `bad j.beta` | wrong locus prefix for the column. Run `/vdjdb-format`; `proofreading/imgt.md` §9 has the older nomenclatures and §10 the common errors |
| `bad species` | not one of the four. If the organism is genuinely a fifth, the chunk goes to `pending/` with an issue naming the missing germline reference |
| `bad mhc.a`, `bad mhc.b` | an `HLA-` string that is not a well-formed allele name. `proofreading/mhc.md` §2, §9, §11 |
| `bad mhc.class` | exactly `MHCI` or `MHCII`. Derive it from the `mhc.a` gene, per `proofreading/mhc.md` §6 |
| `bad antigen.gene` | blank, and only valid when `antigen.species` is `Synthetic`. Run `/vdjdb-harmonize` step 4 |
| `bad reference.id` | `PMID:`, `doi:`, `http://`, `https://` or `unpublished`. Case matters on the first two |
| `no.cdr3` | neither chain has a sequence. Check the extraction - the row has no TCR in it |
| `no.antigen.seq` | no epitope. Required |
| `mhc class/partner mismatch` | Class I requires `mhc.b=B2M`; class II cannot use B2M, and B2M cannot be `mhc.a`. Check the paper restriction, not the CD4/CD8 subset |
| `no.mhc` | one of `mhc.a`/`mhc.b` is blank. `proofreading/mhc.md` §4 for the pairing, and §4.1 for the precedent fills |

**Advisory rules, reported and never fatal.** Each is advisory for a stated reason, so do not
"fix" one into silence:

| Rule | Why it does not fail |
|---|---|
| `non-functional v.alpha`, `non-functional j.alpha`, `non-functional v.beta`, `non-functional j.beta` | IMGT's `ORF`/`P` verdict on the named segment. A P gene can rearrange, and IMGT reclassifies between releases |
| `internal cysteine in cdr3.alpha`, `internal cysteine in cdr3.beta` | a junction has one cysteine, the Cys104 it opens with. A second is rare and not impossible - the Jurkat receptor has one - so the record is kept and flagged. 1,521 + 2,804 corpus rows over 108 chunks; only the submitter's source settles a given one |
| `alpha and beta cdr3 identical` | only a curator can say which of the two chains is the wrong one |
| `segment call with no cdr3` | the call is information; the chain cannot reach an output. Tell the submitter while they can still send the sequence |
| `structure id is not a PDB id` | 2,765 corpus rows hold a figure or table reference there. Blanking them moves 6,004 scores, which is a curation decision |
| `duplicate` | per-chunk deduplication removes these on the way in, so the count is the gap between the raw chunk rows and the released records - 10,555 findings over 59 chunks, declared in `rules/qc_advisories.tsv` |
| `counter in antigen.gene`, `counter in antigen.species`, `counter in mhc.a`, `counter in mhc.b` | a dense run of `prefix`+integer inside one chunk's one epitope, which is what dragging a cell down a spreadsheet column produces. Three confirmed cases: `Eef2`..`Eef188` on one epitope held 65 clonotypes apart (#694), `HLA-A*24:03`..`:20` ran over 18 clonotypes where the paper typed every donor `A*24:02` (#625), and the B16 chunk froze `meta.epitope.id` at `p12`. The repair is a `patches/` entry or a chunk edit, and which one is a curation decision. The corpus carries none today: the #625 chunk now says `HLA-A*24:02` itself |
| `one epitope under two HLA genes in one chunk` | the same peptide under `HLA-A*..` and `HLA-B*..` inside one paper. Either donor typing written into `mhc.a` (`goncharov-various-2023-05-06` did it for `RAKFKQLL`, #597) or a paper that really reports two restrictions, and only the paper says which. The corpus carries 28 rows over 2 chunks today, both declared: `PMID_34793243` (`SIIAYTMSL`, `SLIYSTAAL`, `RLFARTRSM`, which the paper reports under A\*02:01 and B\*07:02) and `goncharov-taa-2020-10-12` (`VQIISCQY`) |
| `frequency disagrees with its count and total` | `method.frequency.count` and `.total` are chunk columns as of #696, so a submitter can report all three - and where they do, `count / total` has to be the frequency they wrote (1 % tolerance, for a two-significant-figure export). Which of the three the paper supports is a curation question and the repair is a chunk edit. Zero findings today: the columns are new, so no submission has used them yet, which makes this a gate on the next one |
| `undeclared method.identification token` | the cell is a comma-separated **set of tokens** and `proofreading/method_vocabulary.tsv` settles all 46 the corpus carries, 9 of them from the specification page. A submission naming a method nobody has seen is a method nobody has seen: the corpus already carries `T-Scan`, `YAMTAD system` and `phage display`, none of which the page ever named. The 36,867 findings today are four `pending` tokens over four chunks, the largest being `magnetic beads` on 36,795 records - one decision, not 36,795 defects (#637) |
| `crlf`, `unknown-column`, `prose-column-name`, `empty-column-name`, `reference-id-form` | the reader produces the right record anyway, and rewriting a corpus file to silence a lint costs more than the lint (invariant 1) |

## Step 3 - run the submission report

```bash
uv run vdjdb submission <file>
```

This assembles the whole corpus, so every number is relative to the database. Four things it gives
that no hand check should be repeating:

**The score distribution, computed.** Per chunk, records at score 0, 1, 2 and 3 from
[`docs/standards/confidence-score.md`](../../docs/standards/confidence-score.md)'s rules. A chunk that
is nearly all score 0 has a method-field problem - go to step 4. A `meta.structure.id` that is a real
PDB entry reaches score 3 on its own.

**Values new to VDJdb**, per column, with the values listed. A new epitope legitimately appears here;
so does a typo, and the two are indistinguishable from the count alone. For an epitope, species, gene
or allele the corpus already has under another spelling, this is the signal - go back to
`/vdjdb-format` or `/vdjdb-harmonize`.

**Junctions that contradict their own V or J germline.** Reported per chain with the defect named
(`J absent anchor`, `J corrupt anchor`, `J under-trimmed`, `V corrupt anchor`, `J allele mismatch`,
`unexplained`) and, where the germline supports one, a proposed repair - against the **submitted**
sequence, because `arda.cdr3fix` may already have fixed one end on the way through the build. Full
table in `out/reports/anchors.tsv`.

Two things to hold on to here. The anchor is read from the germline of the segment the record names,
so there is no "must end in F or W" test to apply (invariant 4). And `J allele mismatch` is evidence
about the **call**, not the sequence: 95 mouse chains name `TRAJ47`, which resolves to the ORF `*01`
(`HYANKMIC`), and every one reads `DYANKMIF`, which is `TRAJ47*02` exactly. Repair the call there and
leave the sequence alone.

**Records that repeat a clonotype another chunk reports.** That is independent replication, which is
what raises `vdjdb.score`, not duplication. But a chunk where nearly every record is an echo may be a
dataset VDJdb already has under another reference - check before landing it.

Applying any proposed repair is a chunk edit: its own branch, its own issue, its own message
(invariant 1). Nothing applies one automatically.

## Step 3a - the junction definition, and the three flags that carry it

A TCR junction **starts with `C`, ends with `F` or `W`, and carries exactly one cysteine**. That is the
definition. `vdjdb qc` does not check the first two, because its sequence rules read the alphabet and a
minimum length; it does check the third, as the advisory `internal cysteine in cdr3.alpha` / `.beta`.

**A record failing the definition is kept and flagged, never dropped.** It can be the best record of
what a publication reported, and consumers filter on the flag. `chains` ships six booleans for it:

| Column | False means |
|---|---|
| `v.canonical`, `j.canonical` | the **shipped** sequence does not open with Cys104 / close with Phe118 or Trp118. 161 and 864 chains |
| `v.canonical.submitted`, `j.canonical.submitted` | the **submitted** sequence did not. 413 and 4,701 chains - so 4,089 repairs are invisible in the shipped pair alone |
| `cdr3.one.cysteine` | a cysteine after the first residue. 4,026 chains. `vdjdb qc` asks the same of the raw chunk cell |

A chain with no CDR3 at all reads `true` on all three: an absent sequence has no anchor to be
missing, and that is `no.cdr3`.

The pair matters because the two answer different questions. Shipped-false tells a consumer this record
is not a canonical junction. Submitted-false tells a curator **the chunk cell is wrong and the next
export of that data will be wrong again**, which is the finding to act on. `arda.cdr3fix` repairs most
of them on the way through, so the shipped `cdr3` is usually right while the chunk keeps the defect.

Then use the germline to say **which** defect a flagged junction has, never whether it is one:

| Germline of the named J | Reading |
|---|---|
| carries the Phe, the sequence does not | the **sequence** is wrong - repair it, ~200 chains |
| is an `ORF` or `P` allele whose anchor is lost | the **call** is wrong - a functional sibling matches, 212 chains, the [#327](https://github.com/antigenomics/vdjdb-db/issues/327) shape |
| there is no J call | nothing adjudicates it; the definition is the only check, 144 chains |

Germline agreement with a pseudogene has not validated the junction. It is the same finding one level
up. `out/reports/anchors.tsv` and `proofreading/cdr3_repair.md` carry the per-defect breakdown.

## Step 4 - the method fields

`vdjdb qc` does not check these, and they set the score. Vocabulary in
[`docs/standards/chunk-format.md`](../../docs/standards/chunk-format.md); mechanical corrections in
[`/vdjdb-format`](../vdjdb-format/SKILL.md) §4. What is left here is the blank
`method.identification`, which costs the row its identification point.

Resolve in this order, and record which rule answered:

1. **Swapped fields.** `method.verification` filled and `method.identification` blank usually means
   the two were transposed: identification is how the epitope-reactive cells were found, verification
   is how the cloned TCR was re-tested. Confirm against the paper before moving the value.
2. **Structural entries.** A `meta.structure.id` that is a PDB entry, with no other method reported:
   `structural` in both fields.
3. **The rest of the chunk.** If every other row uses one identification method and the blank row has
   no anomaly of its own, that method is the answer. Confirm against the paper.
4. **Pre-tetramer papers.** Tetramers arrived in 1996 and were not routine until about 2000. A pre-2000
   paper reporting a handful of T-cell clones is `antigen-loaded-targets,limiting-dilution-cloning`,
   or `limiting-dilution-cloning` alone where the paper describes only the cloning.
5. **The abstract.** Fetch the reference and read it: "tetramer", "sorted", "FACS" →
   `tetramer-sort`; "stimulated", "ELISpot", "IFN-γ", "killing assay" → `antigen-loaded-targets`;
   "transfected", "expressing" → `antigen-expressing-targets`.

If none of the five answers, leave it blank and say so. A guessed method is a wrong score on a
record forever.

## Step 5 - antigen fields

Scan `antigen.gene` and `antigen.species` for the patterns in
[`/vdjdb-harmonize`](../vdjdb-harmonize/SKILL.md) step 2 - a `[species]` annotation, a trailing
` protein`, a `Chain A,` prefix, a UniProt description, a multi-word organism name, a blank. Report
the count and up to five examples per category, then ask whether to run `/vdjdb-harmonize`.

Invoke it as a skill. Do not import its code into this step: two skills sharing a function is a
dependency neither one declares, and it is how one drifts when the other is edited.

## Step 6 - MHC beyond the rules

`proofreading/mhc.md` §11 holds the scan commands for the three MHC-II gene-name defects this corpus
has had - a missing digit in `HLA-DPA*`/`DPB*`/`DQA*`, a spurious one in `HLA-DRA1*`, and a missing
`HLA-` prefix. Run those, and §10's queries for allele membership. §4.1 has the precedent fills for a
blank class II partner chain. §7 has murine, macaque and rat naming, and
`proofreading/mhc_nonhuman.tsv` is the membership test for those.

Two checks worth running by hand because they are cross-field and the rules are per-field:

```bash
awk -F'\t' 'NR>1 && $12=="MHCI"  && $11!="B2M" && $11!="" {print NR,$10,$11,$12}' <file>
awk -F'\t' 'NR>1 && $12=="MHCII" && $11=="B2M"            {print NR,$10,$11,$12}' <file>
```

A collapsed class II pair - both chains in `mhc.a` separated by `/`, `mhc.b` blank - splits into the
two fields, each keeping the `HLA-` prefix.

For an allele the databases do not carry, check what the corpus already uses for the same epitope
before touching the paper's value: same epitope, same restriction is a strong prior. Where the paper
disagrees with the corpus, report the disagreement and ask. Do not overwrite either silently.

## Step 7 - cross-check the earlier stages

If `<basename>_extraction_log.txt` or `<basename>_format_log.txt` exist, read them and confirm the
values they claim to have verified or changed are the values in the file now. Then re-verify five rows
at random against the original source, if it is still available.

## Step 8 - report and decide

```
=== PROOFREAD: <file> ===
rows: N        fatal findings: N        advisory findings: N
score 3: N   score 2: N   score 1: N   score 0: N
junction/germline conflicts: N   (N with a germline-supported repair)
values new to VDJdb: N over M columns
records echoing another chunk: N

findings by rule:
  <rule>  N rows   [fatal|advisory]   -> <fix, and the authority>

escalated to the user:
  <question, with the evidence and the proposed answer>

DECISION: chunks/ | chunks/ after N fixes | pending/ (<missing reference>) | withheld/ (<re-export needed>)
```

`pending/` is for a file that is fine and a build that is not ready - it parses, its header matches a
shipping chunk, and it fails on a reference the build does not have. `withheld/` is for a file that has
to change first, one that predates the current specification and cannot be read at all. A chunk whose
header matches a shipping chunk goes to `pending/`, however far it is from landing. Either way: move
it unchanged, commit it on a chunk branch, comment the path and the blocker on its issue, and leave the
issue **open**.

## When a check this skill runs should become a rule

A finding that is mechanical, has no judgement in it and would apply to any future submission does not
belong in a skill. It belongs in `src/vdjdb/qc/rules.py` with a test in `tests/unit/`, so every
curator and every pull request gets it. Open an issue, say what it detects and on how many corpus rows
it fires today, and add it there. A rule that fires on zero corpus rows is free to make fatal; one that
fires on thousands is advisory until the data is fixed, and `rules/qc_advisories.tsv` declares it.

Do not keep a list of pending checks in this file. A list dated to the month it was written outlives
the thing it describes, and `out/reports/qc-summary.tsv` already says which rules fire and how often.
