# What decides, and what reads it

Every curation skill in this directory links here once. It holds the parts they share, so no skill
restates a rule another file already owns - a restated rule is a second copy that drifts, and the
copies in these skills had already drifted apart on murine MHC before this file existed.

## One authority per question

| Question | Authority | Read it with |
|---|---|---|
| Which columns may a chunk have, in what order | [`docs/standards/chunk-format.md`](../docs/standards/chunk-format.md) | `vdjdb schema --table records` |
| Is this V/D/J gene or allele a name IMGT has | `proofreading/imgt_alleles.tsv.gz` | `proofreading/imgt.md` §8 has the queries |
| What does this older gene name convert to | `proofreading/arden.tsv`, `patches/nomenclature.conversions` | `proofreading/imgt.md` §9 |
| Is this human HLA allele a name IPD-IMGT/HLA has | `proofreading/mhc_alleles.tsv.gz` | `proofreading/mhc.md` §10 has the queries |
| Is this murine, macaque or rat MHC name the one VDJdb records | `proofreading/mhc_nonhuman.tsv` | `proofreading/mhc.md` §7 |
| What is the declared correction for this MHC name | `patches/mhc.dict` | the file states its own reasoning per row |
| What are this epitope's `antigen.gene` and `antigen.species` | `patches/antigen_epitope_species_gene.dict` | keyed on the epitope |
| What does this free-text antigen name normalise to | `proofreading/gene_aliases.tsv`, `proofreading/species_aliases.tsv` | species fragments are order-sensitive |
| Which anchor residue should this junction end in | the germline of the segment the record names | `vdjdb submission` prints it per chain |
| What score will these method fields earn | [`docs/standards/confidence-score.md`](../docs/standards/confidence-score.md) | `vdjdb submission` computes it |
| Which method vocabulary is recognised | `docs/standards/chunk-format.md`, method columns | - |

When a value has no entry, the fix is to **add the entry** to the authority and say so, not to repair
the chunk cell in isolation. A chunk repaired without a declared rule behind it repeats itself on the
next submission of the same data.

## Five invariants no skill may relax

1. **`chunks/` is the submitter's data.** It is never edited to make a tool pass. A chunk edit is its
   own commit on its own branch named for the chunk or the data issue it answers, with a message
   saying which files and rows moved, why, and who decided. Never bundle a chunk edit with a code
   change. See `CLAUDE.md`.
2. **Empty string is the only missing marker.** Never `NA`, `N/A`, `null`, `nan`, `-`, `.`, or `?`.
3. **`cdr3` is junction space** - Cys104 through Phe/Trp118, **both anchors included**. That is two
   residues longer than AIRR's and arda's `cdr3_aa`. A submission exported in IMGT CDR3 space is
   short one residue at each end and passes `vdjdb qc`, because those rules check the alphabet and a
   minimum length, not the ends. `vdjdb submission` is what catches it.
4. **A TCR junction starts with Cys104 and ends with Phe118 or Trp118. That is the definition, not a
   heuristic**, and it is checked at submission time, where it is cheap to fix. A sequence failing it
   is not a variant: it is an export in IMGT CDR3 space, a mis-read anchor, framework left in, or a
   mis-called allele.

   The germline of the segment the record names says **which** of those it is, never whether it is
   one. Read it that way round. Measured over the corpus, 864 of 285,989 chains fail the definition:
   212 sit on an ORF allele whose anchor is genuinely lost, so the *call* is what needs repairing
   (mouse `TRAJ47`/`TRAJ7`/`TRAJ44` resolve to an ORF `*01` where a functional sibling matches the
   sequence - the same defect as [#327](https://github.com/antigenomics/vdjdb-db/issues/327)); about
   200 name a germline that does carry the Phe, so the *sequence* is wrong; 144 name no J at all, so
   the definition is the only thing left to check them against.

   Germline agreement with a non-functional allele is not a validation. It is the same finding one
   level up.

5. **Never invent a value.** Every amino acid sequence, gene name, allele, species and reference id
   written into a chunk is confirmed present in the source by a search of the source, and the result
   of that search is logged. A PMID is never guessed. A value that cannot be confirmed is marked
   `[UNVERIFIED]` and withheld until the user approves it.

## The commands that replace hand-written checks

The build validates and measures what these skills used to check in prompt-resident Python. Run the
command; read its report. The command is tested, versioned and the same for every curator.

```bash
uv run vdjdb qc <file> --report out/reports/qc.tsv   # 24 row rules, 12 text lints, per row
uv run vdjdb submission <file>                       # records, score distribution, values new to
                                                     # VDJdb, junction/germline conflicts, replication
uv run vdjdb schema --table records                  # the column contract, from the field registry
uv run vdjdb rules                                   # regenerate declared nomenclature renames
```

`vdjdb qc` writes `qc.tsv` (one row per finding) beside `qc-summary.tsv` (one row per rule, with an
`advisory` flag). A rule marked advisory is reported and never fatal; `--strict` fails on the rest.

`vdjdb submission` assembles the whole corpus, so every number it prints is relative to the database
rather than to the file. It stops before the annotation stages, which cost ten times as much and
change nothing a curator decides on.

## What the CLI cannot decide

These are the judgement calls, and they are why the skills exist at all:

- which of two spellings a paper actually meant;
- whether a method the authors describe maps to an existing vocabulary term or needs a new one;
- whether two rows from one paper are two observations or one entered twice;
- whether an epitope that is a substring of a longer one is a truncation artefact or a real shorter
  peptide;
- whether a record should land in `chunks/`, `pending/`, `withheld/`, or not at all.

Each one is escalated to the user with the evidence attached, never resolved by guessing.
