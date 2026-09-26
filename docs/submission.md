# Submitting and curating


To submit previously published sequence follow the steps below:

* Create an issue(s) labeled as ``paper`` and named by the paper pubmed id, ``PMID:XXXXXXX``. Note that if paper is a meta-study, you can mark it as ``meta-paper`` and link issues for its references in a reply to this issue. Also note that in case submitting unpublished sequences, choose any appropriate issue name with details on submitter (name, organization, etc) in issue comments.

* Create new branch and add chunk(s) for corresponding papers named as ``PMID_XXXXXXX``. Don't forget to close/reference corresponding issues in the commit message.

* Create a pull request for the branch and check if it passes the CI build. If there are any issues, modify them by fixing/removing entries as necessary.

The structure of submission chunk is provided below, but first a couple of notes:

> **STYLE** Try avoiding spaces (e.g. ``TRBV7,TRBV5``, not ``TRBV7, TRBV5``) and leave fields that have no information as blank (don't use any placeholder). Stick to listed field values at all cost! In case a critical part of your submission doesn't fit in current specification: 1) Create an issue in the issues section (and tag it as ``maintainance``), 2) provide us with an example (e.g. open a pull request). Do not insert critical information into the comment field.

> **FORMAT** Please ensure that Variable/Joining and MHC names in your submission come from IMGT nomenclature (this does not apply to donor MHC typing fields).

The ``BuildDatabase`` routine will be executed during CI tests upon each submission and prior to every database release implements table format checks, CDR3 sequence checks and fixes (if possible), and VDJdb confidence score assignment (see below).

To view the list of papers that were not yet processed follow [here](https://github.com/antigenomics/vdjdb-db/labels/paper).

An XLS template is available [here](https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/template.xls).

> **CAUTION** make sure that nothing is messed up (``x/X`` frequencies are transformed to dates, bad encoding, etc) when importing from XLS template. The format of all fields is pre-set to *text* to prevent this case.



This repository includes a set of **skills** (`skills/`) designed for use with [Claude Code](https://claude.ai/code) (Anthropic's CLI agent) and compatible with GitHub Copilot's agent mode. Skills are instructional documents that guide an AI assistant through multi-step curation, formatting, and quality-control tasks on VDJdb chunks.

## Available skills

| Skill | Invocation | Purpose |
|---|---|---|
| **vdjdb-extract** | `/extract [file]` | Extract TCR:pMHC records from raw source files (PDFs, Excel, AIRR-format, 10x Genomics) into VDJdb-format TSV chunks |
| **vdjdb-format** | `/format [file]` | Normalise V/J gene names (IMGT), MHC allele format, species names, and method vocabulary |
| **vdjdb-harmonize** | `/harmonize [file]` | Canonicalize `antigen.gene` and `antigen.species` using `patches/antigen_epitope_species_gene.dict`, `proofreading/gene_aliases.tsv`, and `proofreading/species_aliases.tsv`; detects spurious values and warns about epitope substrings |
| **vdjdb-proofread** | `/proofread [file]` | Run ChunkQC validation, enhanced IMGT gene checks, CDR3 canonical repair, MHC consistency checks, and confidence score estimation |
| **vdjdb-publish** | `/vdjdb-publish` | For each new or modified chunk, find or create a GitHub issue (`PMID:N`), commit the chunk with `Fixes #N`; processes one chunk at a time with user confirmation |
| **vdjdb-duplicates** | `/vdjdb-duplicates` | Audit duplicate TCR records (beta-only and paired), classify by publication source and author overlap, flag high-frequency records from assay artifacts |

## Using with Claude Code

1. Install [Claude Code](https://claude.ai/code) CLI or open the repo in the Claude Code desktop app.
2. Invoke any skill by typing `/skill-name [arguments]` in the chat.
3. Skills guide the agent step-by-step and always ask for confirmation before making commits or GitHub API calls.

## Using with GitHub Copilot

Skills are plain Markdown instructional documents and can be referenced directly in Copilot chat:

```
@workspace /skills/vdjdb-extract/SKILL.md — extract data from supplementary table X
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
