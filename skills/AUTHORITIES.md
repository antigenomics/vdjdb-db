# Curation authorities

Use the specification and reviewed authority tables below. Keep source evidence and decisions in
the local execution record; retain material questions on the publication issue.

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

## Submission invariants

- Preserve literal source sequences and calls when assembly already resolves them. Inspect the
  submitted and shipped values together; an inferred V is not a reported V.
- Use empty cells for missing values. Do not insert placeholders or generated export identifiers.
- The definition of VDJdb junction space includes Cys104 and Phe/Trp118. AIRR `cdr3_aa` excludes both; use
  `junction_aa` when available. Printed cores, flanks and noncanonical anchors require repair
  inspection, not automatic deletion or manual sequence substitution.
- At least one reported junction and the literal peptide are required. Optional provenance may
  stay blank. Partially reported class-II restriction may retain its missing partner blank.
- Report positive, explicitly negative, unassigned and unresolved observations separately.
  Do not infer peptide identity, pairing or restriction from a protein label or sequence similarity.
- Different experiments and publication reports remain separate observations. Compare both chains,
  pMHC and complete metadata before identifying duplicates.
- Chunk edits use a branch for that chunk or data issue, separately from code changes. Reconcile
  identity, document the per-file reason/counts, and validate against the release before publishing.

## Commands and reports

```bash
uv run vdjdb qc <file> --strict --report out/reports/qc.tsv
uv run vdjdb submission <file>
uv run vdjdb build --out out/
uv run vdjdb schema --table records
```

`qc-summary.tsv` declares which findings are advisory. Build reports identify harmonisation,
unresolved nomenclature, junction/germline conflicts and lookalike terms. Read them before proposing
source edits. Unknown optional metadata does not block a supported observation. Multiple reported
chains in the same identified clonotype follow the
[candidate-pair convention](../docs/standards/chunk-format.md#multiple-chains-in-one-clonotype),
retaining the original clonotype ID and marking the uncertainty. Unsupported source-chain
association, peptide/restriction or outcome contradictions block the affected rows only.

Methods describe the experiment. Culture before sequencing, initial identification and later
receptor verification are distinct. Use source evidence for each rather than paper date, token
frequency or another row's method.

For association queries, join records and chains by `record_id`, select species/chain and condition
on the epitope and restriction. State whether the denominator counts observations, clonotypes or
publications. See [record-level queries](../docs/standards/corpus.md#receptor-level-motif-association).

Use the [submission guide](../docs/submission.md) for blockers and quarantine, and
[publish](vdjdb-publish/SKILL.md) for Gitflow, validation and issue closeout. Request a decision only
when source evidence and existing authorization do not resolve it.

Use [structures](vdjdb-structures/SKILL.md) for the local tcren extraction of deposited complexes.
Use [build integrity](vdjdb-build-integrity/SKILL.md) for code, CI and release lifecycle validation;
it does not replace source curation.
