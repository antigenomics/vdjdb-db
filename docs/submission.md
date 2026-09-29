# Submitting and curating

To submit a previously published sequence, follow the steps below.

* Create an issue labelled ``paper`` and named by the paper's PubMed id, ``PMID:XXXXXXX``. If the paper is a meta-study, label it ``meta-paper`` and link the issues for its references in a reply to that issue. For unpublished sequences, choose any appropriate issue name and give the submitter details (name, organization) in the issue comments.

* Branch from `dev`, not from `master`, and add one chunk per paper, named ``PMID_XXXXXXX``. One commit per chunk, and close or reference the corresponding issue in the commit message.

* Open a pull request against `dev`. `chunk-check` runs on it and reports what the submission does, as a comment updated in place on every push: the records each changed chunk contributes and their confidence-score histogram, every QC finding by rule, and **every value the chunk introduces that no other chunk carries** - a new epitope, species, gene, MHC allele or reference. That last section is the one to read twice: a mistyped epitope passes every format check and arrives as a new epitope, indistinguishable from a genuine one except that you know which you meant. It also counts the records that repeat a clonotype and pMHC another chunk already reports, which is independent replication rather than duplication and is what raises the score. Fix or remove entries until it is green. A pull request straight into `master` is rejected by `branch-policy`, which only lets `dev` and `hotfix/*` merge there.

* Run `uv run vdjdb identity update` and commit `registry/records.tsv` in the same branch. That file is the record identity registry: it is what keeps a `record_id` pointing at the same record across releases, so external references and accumulated evidence outlive curation. The command is idempotent, so its diff in your pull request is exactly the records your chunk added, amended or retired - a reviewer reads it the way they read the chunk. Skipping it does not fail anything immediately; it leaves the registry describing a corpus that no longer exists.

* Chunks reach `master` with the next `dev` to `master` merge, once the full build has run green on `dev`. The path is `dev` -> chunk branch -> `dev` -> `master`.

The chunk structure is specified in {doc}`standards/chunk-format`. Two rules apply to every submission:

> **STYLE** Avoid spaces in multi-value fields (``TRBV7,TRBV5``, not ``TRBV7, TRBV5``), and leave a field with no information blank rather than filling it with a placeholder. Use only the listed field values. If a critical part of your submission does not fit the current specification, 1) create an issue tagged ``maintainance``, and 2) provide an example, for instance by opening a pull request. Do not put critical information in the comment field.

> **FORMAT** Variable/Joining and MHC names must follow IMGT nomenclature. This does not apply to the donor MHC typing fields.

The ``BuildDatabase`` routine runs in CI on every submission and before every release. It performs table format checks, CDR3 sequence checks and fixes where possible ({doc}`standards/cdr3-fixing`), and confidence score assignment ({doc}`standards/confidence-score`).

Papers not yet processed are listed under the [`paper` label](https://github.com/antigenomics/vdjdb-db/labels/paper).

An [XLS template](https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/template.xls) is available for preparing a chunk.

> **CAUTION** Check that nothing is corrupted on import from the XLS template: ``x/X`` frequencies turned into dates, bad encoding, and similar. The format of every field is pre-set to *text* to prevent this.

## Every change to `chunks/` happens on a branch that names its reason

`chunks/` is the database. A build can be rewritten and re-verified against the last release; a chunk
edit changes what VDJdb *says*, and the only instrument that notices is the comparison against the
last release, which reports it as rows appearing and disappearing with no reason attached. So the
branch carries the reason, and there are exactly two kinds of branch that may touch `chunks/`:

| Branch | Subject | Names |
|---|---|---|
| `chunk/PMID_<id>` | one publication's records | the chunk |
| `proofread/<issue>-<slug>` | one tracker issue about the data | the issue |

No other branch edits `chunks/`. A branch whose subject is the build, the tests, the documentation or
the CI leaves the data alone, however mechanical the edit looks - and a mechanical repair across many
files is not an exception to this, it is the case that most needs it. **A normalisation pass gets its
own issue** saying what it changes and what it must not, and its own `proofread/` branch.

Line endings are the worked example. 99 of the 230 chunks are CRLF, normalising them touches every
line of every one of those files, and when that was once folded into a build branch the commit message
said the content was unchanged and was wrong: 92 files were line-endings only, but 11 carried real
data-line changes, one of them shifting every field of 2,352 rows by dropping an unnamed leading
column. The fix for a change of that shape is to make it checkable - for line endings,
`git diff --ignore-cr-at-eol` on the branch returns empty - and to put it where a reviewer reading
`git log chunks/` sees one line per reason.

Whichever branch it is, the commit message carries four things: **which files**, and per file how many
rows it adds, removes or changes; **why**, naming the paper, the tracker issue or the `proofreading/`
table the change comes from; **what the build shows**, meaning the row-count delta and the score
histogram if it moved; and **who decided**, where the edit is a curation judgement rather than a
mechanical repair. The last one is the part no diff can reconstruct later.

## A submission that cannot land yet goes to `pending/` or `withheld/`

Not every chunk can land when it arrives. It may be in a format that predates the current
specification, carry a species or a nomenclature the build has no germline or allele reference for, or
raise a question only the submitting author can settle.

**Leaving it on a branch is the one wrong answer.** A branch is invisible: nothing indexes it, no build
reads it, and the submission is lost the moment someone tidies the branch list. Two submissions sat
unlanded on branches for eight and ten years and were recovered only by checking every unmerged branch
against the tracker.

Both directories are inputs that no build reads, so the file stays tracked, greppable and reviewable
and the records are there when whatever blocks them is fixed. Which one depends on where the problem
is - in the file, or in the build:

| Directory | The chunk | What has to change | Who changes it |
|---|---|---|---|
| `pending/` | is in the current format and parses cleanly | the **build** gains a reference it lacks - a species vocabulary, a germline set, an allele database | a maintainer |
| `withheld/` | predates the current specification and cannot be read at all | the **file** is re-exported against the current column set | a curator, or the submitting author |

The difference is checkable rather than a judgement. The five chunks in `withheld/` carry 34- or
36-column headers beginning `cdr3.alpha`; the current header is 33 columns beginning `chunk.id`. A
chunk whose header matches a shipping chunk does not belong there, however far it is from landing:
filing it under `withheld/` tells the next curator to re-export a file that needs no re-export.

`pending/PMID_22058411.txt` is the worked example. Its header is identical to a chunk that ships and
`vdjdb qc` parses all 53 rows; every row then fails one rule, `bad species`, because `species` is
`BosTaurus` and three parts of the build have no bovine input. Nothing about the file is wrong.

Either way:

1. `git mv` or copy the chunk under the same name, **unchanged**. Do not repair it on the way in: the
   file should stay what the submitter sent, so the next curator sees the original.
2. Commit it on a chunk branch through `dev`, with a message naming what blocks it and what would
   unblock it.
3. Comment on the issue with the new path, the blocking reason, and the condition that would let it
   land. Leave the issue **open** - it is still a pending submission, and closing it makes a blocked
   chunk indistinguishable from a rejected one. If no issue exists, open one.

Name the blocker in terms someone can act on. "Bad format" is not one; "33-column header predates the
`.tsv` migration, needs re-export from the source table" is.

This applies to a chunk you merely doubt as much as to one that fails a QC rule. A record you are
unsure of is better quarantined with the doubt written down than silently dropped or silently shipped.
Use `pending/` when you expect it to land and `withheld/` when the file itself has to change first.

The repository includes curation skills in `skills/`, for use with [Claude Code](https://claude.ai/code) (Anthropic's CLI agent) and with GitHub Copilot's agent mode. A skill is an instructional document that guides an AI assistant through a multi-step curation, formatting or quality-control task on VDJdb chunks.

## What `chunk-check` alerts on without failing

Two checks in the pull-request report block nothing and both exist because the thing they catch passes
every QC rule.

**Look-alike values.** A value that differs from one VDJdb already has only in case or in a `-`, `_`,
`.` or space will not join it, so every query filtering on one misses the other. `IE1` and `IE-1` are
two `antigen.gene` values today for the same CMV gene. Two values that look alike can also both be
right - one stain against another, `MBP` in human against `Mbp` in mouse - which is why this is an
alert and not a gate.

**Junctions that contradict their own germline.** VDJdb's `cdr3` is **junction space**: Cys104 through
Phe/Trp118, **both anchors included**. That is two residues longer than AIRR's or arda's `cdr3_aa`.
A submission exported in IMGT CDR3 space is therefore short an anchor at each end, and it passes
`vdjdb qc` - those rules check the residue alphabet and a minimum length, not the ends.

Measured on the 2026-09-28 corpus: **5,842 of 285,989 chains (2.04 %)**, and `PMID_15589168.txt` alone
contributes 935 of them, every one short the trailing Phe. 5,097 get a proposed sequence and 95 get a
proposed **allele** instead:

| Defect | Chains | Example | Repair |
|---|--:|---|---|
| J absent anchor | 3,807 | `CASSNEKLF` on `TRBJ1-4` (`TNEKLFF`) | `CASSNEKLFF` |
| J under-trimmed | 650 | `YLCSSQEGGYGYTFGSG` | `CSSQEGGYGYTF` |
| J corrupt anchor | 423 | terminal residue mis-read | germline residue substituted |
| V + J under-trimmed | 93 | framework at both ends | trimmed at both |
| V corrupt anchor | 77 | `GASSDTMNTKIL`, `WAVRDIYTTAKFIL` | `CASSDTMNTKIL` |
| V under-trimmed | 21 | `YFCASSY...` | `CASSY...` |
| J allele mismatch | 95 | `CAADYANKMIF` on mouse `TRAJ47*01` | the call, `TRAJ47*02` |
| unexplained | 650 | - | none proposed |

A missing Cys104 is worth a second look even when the rest of the record is fine: no TCR folds without
it, so a first residue that is not Cys where the body still aligns to the V is a sequencing or
transcription error rather than a variant.

Two things this check is careful about:

- **The anchor residue is read from the germline of the segment the record names, never assumed to be
  Phe or Trp.** Mouse `TRAJ47*01` is `HYANKMIC` and human `TRAJ35*01` is `IGFGNVLHC`, so a junction on
  either ends in Cys. A fixed "ends with F or W" test calls 481 correct chains broken and separately
  misses 125 that are not.
- **Where the junction matches a functional sibling allele of the gene, the call is repaired and the
  sequence is left alone.** An ORF or pseudogene allele has a non-canonical anchor by definition, and
  arda records that faithfully: of 383 J entries over four organisms only 14 have `templated_aa` not
  ending in Phe or Trp and 13 of those are marked `ORF` or `P`. So a junction disagreeing with a
  non-functional allele is evidence about the *call*. 95 mouse chains name `TRAJ47`, which resolves to
  the ORF `*01` (`HYANKMIC`), and every one reads `DYANKMIF` - exactly `TRAJ47*02`, the functional
  allele. No record in the corpus reads the `*01` signature. This is the same defect as
  [#327](https://github.com/antigenomics/vdjdb-db/issues/327), where 66 % of explicit `TRAJ24*01`
  calls carry the `*02` motif, and rewriting the sequence there would destroy the evidence for it -
  the reasoning `MAX_REPLACE = 0` already applies to CDR3 repair.

`arda.cdr3fix` repairs most of these on the way through the build, which is the reason the alert
matters rather than a reason to skip it: the shipped `cdr3` is usually right and **the chunk keeps the
wrong sequence**, so the next export of that data is wrong again. The repair is therefore proposed
against the submitted value, not the shipped one - arda may have fixed one end already, and
`YLCSSQEGGYGYTFGSG` ships as `YLCSSQEGGYGYTF`, framework trimmed behind the anchor and kept in front
of it.

Applying a repair is a chunk edit, so it follows the rule above: its own branch, its own issue, and a
message saying which files and rows moved and why. Nothing applies one automatically.

## Available skills

| Skill | Invocation | Purpose |
|---|---|---|
| `vdjdb-extract` | `/extract [file]` | Extract TCR:pMHC records from raw source files (PDFs, Excel, AIRR-format, 10x Genomics) into VDJdb-format TSV chunks |
| `vdjdb-format` | `/format [file]` | Normalise V/J gene names (IMGT), MHC allele format, species names, and method vocabulary |
| `vdjdb-harmonize` | `/harmonize [file]` | Canonicalize `antigen.gene` and `antigen.species` using `patches/antigen_epitope_species_gene.dict`, `proofreading/gene_aliases.tsv`, and `proofreading/species_aliases.tsv`; detects spurious values and warns about epitope substrings |
| `vdjdb-proofread` | `/proofread [file]` | Run ChunkQC validation, enhanced IMGT gene checks, CDR3 canonical repair, MHC consistency checks, and confidence score estimation |
| `vdjdb-publish` | `/vdjdb-publish` | For each new or modified chunk, find or create a GitHub issue (`PMID:N`), commit the chunk with `Fixes #N`; processes one chunk at a time with user confirmation |
| `vdjdb-duplicates` | `/vdjdb-duplicates` | Audit duplicate TCR records (beta-only and paired), classify by publication source and author overlap, flag high-frequency records from assay artifacts |

## Using with Claude Code

1. Install the [Claude Code](https://claude.ai/code) CLI, or open the repo in the Claude Code desktop app.
2. Invoke any skill by typing `/skill-name [arguments]` in the chat.
3. Each skill asks for confirmation before making a commit or a GitHub API call.

## Using with GitHub Copilot

Skills are plain Markdown documents and can be referenced directly in Copilot chat:

```
@workspace /skills/vdjdb-extract/SKILL.md - extract data from supplementary table X
```

## Reference files for proofreading

| File | Role |
|---|---|
| `proofreading/gene_aliases.tsv` | Free-text antigen gene names → VDJdb canonical symbols (136+ mappings) |
| `proofreading/species_aliases.tsv` | Source organism substrings → canonical CamelCase species names |
| `proofreading/cdr3_repair.md` | CDR3 canonical repair algorithm using V/J germline context |
| `proofreading/imgt.md` | IMGT V/D/J gene naming rules |
| `proofreading/mhc.md` | HLA/MHC allele naming and validation rules |
| `patches/antigen_epitope_species_gene.dict` | Epitope-keyed authority: epitope → (species, gene) |
