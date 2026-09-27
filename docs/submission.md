# Submitting and curating

To submit a previously published sequence, follow the steps below.

* Create an issue labelled ``paper`` and named by the paper's PubMed id, ``PMID:XXXXXXX``. If the paper is a meta-study, label it ``meta-paper`` and link the issues for its references in a reply to that issue. For unpublished sequences, choose any appropriate issue name and give the submitter details (name, organization) in the issue comments.

* Branch from `dev`, not from `master`, and add one chunk per paper, named ``PMID_XXXXXXX``. One commit per chunk, and close or reference the corresponding issue in the commit message.

* Open a pull request against `dev`. `chunk-check` runs on it and reports what the submission does: records added or removed, the confidence-score histogram, and any QC error by rule. Fix or remove entries until it is green. A pull request straight into `master` is rejected by `branch-policy`, which only lets `dev` and `hotfix/*` merge there.

* Chunks reach `master` with the next `dev` to `master` merge, once the full build has run green on `dev`. The path is `dev` -> chunk branch -> `dev` -> `master`.

The chunk structure is specified in {doc}`standards/chunk-format`. Two rules apply to every submission:

> **STYLE** Avoid spaces in multi-value fields (``TRBV7,TRBV5``, not ``TRBV7, TRBV5``), and leave a field with no information blank rather than filling it with a placeholder. Use only the listed field values. If a critical part of your submission does not fit the current specification, 1) create an issue tagged ``maintainance``, and 2) provide an example, for instance by opening a pull request. Do not put critical information in the comment field.

> **FORMAT** Variable/Joining and MHC names must follow IMGT nomenclature. This does not apply to the donor MHC typing fields.

The ``BuildDatabase`` routine runs in CI on every submission and before every release. It performs table format checks, CDR3 sequence checks and fixes where possible ({doc}`standards/cdr3-fixing`), and confidence score assignment ({doc}`standards/confidence-score`).

Papers not yet processed are listed under the [`paper` label](https://github.com/antigenomics/vdjdb-db/labels/paper).

An [XLS template](https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/template.xls) is available for preparing a chunk.

> **CAUTION** Check that nothing is corrupted on import from the XLS template: ``x/X`` frequencies turned into dates, bad encoding, and similar. The format of every field is pre-set to *text* to prevent this.

## A submission that cannot land yet goes to `withheld/`

Not every chunk can land when it arrives. It may be in a format that predates the current
specification, carry a species or a nomenclature the build has no germline or allele reference for, or
raise a question only the submitting author can settle.

**Leaving it on a branch is the one wrong answer.** A branch is invisible: nothing indexes it, no build
reads it, and the submission is lost the moment someone tidies the branch list. Two submissions sat
unlanded on branches for eight and ten years and were recovered only by checking every unmerged branch
against the tracker.

Move it to `withheld/` instead, which is an input directory no build reads, so the file stays tracked,
greppable and reviewable and the records are there when whatever blocks them is fixed.

1. Copy the chunk to `withheld/` under the same name, **unchanged**. Do not repair it on the way in:
   the file should stay what the submitter sent, so the next curator sees the original.
2. Commit it on a chunk branch through `dev`, with a message naming what blocks it and what would
   unblock it.
3. Comment on the issue with the new path, the blocking reason, and the condition that would let it
   land. Leave the issue **open** - it is still a pending submission, and closing it makes a withheld
   chunk indistinguishable from a rejected one. If no issue exists, open one.

Name the blocker in terms someone can act on. "Bad format" is not one; "33-column header predates the
`.tsv` migration, needs re-export from the source table" is.

This applies to a chunk you merely doubt as much as to one that fails a QC rule. A record you are
unsure of is better in `withheld/` with the doubt written down than silently dropped or silently
shipped.

The repository includes curation skills in `skills/`, for use with [Claude Code](https://claude.ai/code) (Anthropic's CLI agent) and with GitHub Copilot's agent mode. A skill is an instructional document that guides an AI assistant through a multi-step curation, formatting or quality-control task on VDJdb chunks.

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
