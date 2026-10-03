---
name: vdjdb-format
description: Normalise the controlled-vocabulary fields of a raw VDJdb chunk TSV to the form the database records - IMGT V/D/J gene and allele names including Adaptive ImmunoSEQ and Arden conversions, IPD-IMGT/HLA and murine H2 MHC allele names, species CamelCase, method vocabulary, reference id prefixes - and write a format log naming the authority behind every change. Use after vdjdb-extract and before vdjdb-proofread, or on any TSV whose gene or allele spellings need bringing to IMGT.
---

# vdjdb-format

Standardise every controlled-vocabulary field in a chunk to the spelling VDJdb records. Second of
three stages: [extract](../vdjdb-extract/SKILL.md) → **format** →
[proofread](../vdjdb-proofread/SKILL.md).

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first. This skill applies the authorities in that
table; it does not restate what they say.

## Invocation

```
/vdjdb-format [path-to-tsv]
```

If the file has structural damage - missing columns, wrong separator, a column shift - stop and report
it. Repairing structure is `/vdjdb-proofread` step 1, not this skill.

## The rule this skill runs on

**Normalise toward what the database already records, and prove it with a query.** Every change is
either a row in an authority file or a lookup that returned a hit. A change made from memory of "the
usual spelling" is how two spellings of one molecule entered the corpus in the first place: the
murine `H-2` / `H2-` split cost 768 records their motif badge on the deployed site.

Where the authority has no entry for a value, add the entry and say so in the log. Never repair the
cell alone.

## Preserve source values already handled by assembly

Before rewriting a source cell, inspect the existing harmonisation and batched CDR3-fixer result.
Supported spelling conversion, terminal-flank trimming and subgroup/allele resolution belong in
assembly. Retain the supplied cell when the build resolves it; do not fetch a paper or choose a
family member manually for that case. The rules below describe the target vocabulary and review
steps, not a requirement to duplicate every build conversion at import.

## 1. Species

`species` is one of `HomoSapiens`, `MusMusculus`, `RattusNorvegicus`, `MacacaMulatta` - CamelCase, no
spaces, case-sensitive. Normalise any binomial, common name or abbreviation onto one of the four.

A fifth species is not a formatting problem. `vdjdb qc` fails it as `bad species` because no part of
the build has a germline reference for it, and the chunk belongs in `pending/` with an issue saying
which reference would unblock it. `pending/PMID_22058411.txt` is the worked case: 53 bovine records,
header identical to a shipping chunk, every row failing that one rule.

## 2. V, D and J gene calls

Authority: `proofreading/imgt_alleles.tsv.gz`. Rules and queries: `proofreading/imgt.md` -
§3 naming structure, §8 how to query, §9 older nomenclatures, §10 common errors.

Apply in order:

1. **Strip internal whitespace.** `TRBV 7` → `TRBV7`.
2. **Convert Adaptive ImmunoSEQ names** per `proofreading/imgt.md` §9.2, which has the three
   differences, the five patterns seen in this corpus and the conversion algorithm. Log each as
   `ADAPTIVE <from> -> <to>`.
3. **Convert older nomenclatures.** Arden names (`BV20S1`, `TRBV1S1`) via `proofreading/arden.tsv`;
   everything else via `patches/nomenclature.conversions`. §9.1 and §9.3 cover the patterns.
4. **Look the result up** in `imgt_alleles.tsv.gz` at gene level, then at allele level if an allele is
   named. `proofreading/imgt.md` §8 has both queries.
5. **A name that resolves at neither level is not a formatting problem** - it is nomenclature debt.
   Report it; do not invent a nearest match. The build lists all of them in
   `out/reports/nomenclature.tsv`, separating family-level calls from names the authority does not contain. First inspect the
   build's sequence-supported resolution; only unresolved source ambiguities need a curator.
6. **Ambiguous multi-calls** stay as a comma-separated list with no spaces (`TRBV7-2,TRBV7-3`), each
   part checked separately. The build reads `,`, `;`, `+` and `or` as separators.

**Bound the allele number.** An allele far above any real count is a spreadsheet artefact, not a
call: someone drag-fills a column and gets `*01`, `*02`, `*03` ... `*112`. Two cheap checks catch it,
and both fire on **zero** corpus rows today, so neither costs anything to enforce:

- **Per value.** The highest allele number any TCR gene has in `imgt_alleles.tsv.gz` is 10, and the
  highest in the whole corpus is `*08`. Anything past the named gene's own allele count is suspect;
  anything past about 12 is not a call at all.
- **Across rows.** The drag-fill signature is a *run*: the same gene with allele numbers stepping by
  one down consecutive rows. Measured over `chunks/`, there is no such run of 5 or more anywhere, so
  one appearing is the artefact and nothing else.

Then adjudicate a flagged value against `imgt_alleles.tsv.gz` rather than deleting it on the bound
alone - the bound is the detector, membership is the verdict. That order matters, because the bound
has false positives on genes with many alleles: `TRAV8-4*07`, `TRBV20-1*07`, `TRBV7-9*07` and mouse
`TRAV14D-3/DV8*08` are 19 corpus calls that IMGT does list. The retired build flagged one of them and
had no membership test to resolve it with; now there is one.

A gene whose IMGT functionality is `ORF` or `P` is reported, not rejected. `vdjdb qc` calls it
`non-functional <column>` and treats it as advisory: a pseudogene can rearrange, and IMGT
reclassifies genes between releases.

## 3. MHC alleles

Authorities: `proofreading/mhc_alleles.tsv.gz` for human, `proofreading/mhc_nonhuman.tsv` for
everything else, `patches/mhc.dict` for declared corrections. Rules: `proofreading/mhc.md` - §2
allele structure, §4 class conventions, §7 non-human, §9 historical forms, §10 how to query, §11 the
MHC-II gene-name digit errors this corpus has had.

**Human.** Target `HLA-<GENE>*<field1>:<field2>`. Prefix, `*` and `:` present; low resolution kept as
low resolution with a note; expression suffixes (`N`, `L`, `Q`) kept. `mhc.b` on class I is the
literal `B2M`.

**Murine.** One prefix, `H2-`, for every murine name. `H2-` is the MGI gene symbol prefix; `H-2` is
the classical immunology spelling and is not a symbol. VDJdb records the MGI form and so does the
site - `vdjdb-web`'s `Motifs.scala` maps `h-2` to `h2-` before joining a record to its motif cluster.

| Normalise from | To | Why |
|---|---|---|
| `H-2Db`, `H-2Kb`, `H-2Aa`, `H-2Eb1` | `H2-Db`, `H2-Kb`, `H2-Aa`, `H2-Eb1` | prefix |
| `I-Ab`, `IAb` | `H2-IAb` | prefix; the one murine string of fifteen that lacked it |
| `H2-Ag7`, `H2-Ed` | `H2-IAg7`, `H2-IEd` | the dropped `I` |
| `H-2D^b` | `H2-Db` | superscript flattened, then prefix |

Every one of those is a declared row in `patches/mhc.dict` with its own reasoning. **Do not convert
in the other direction.** An earlier revision of this skill said to write `H-2Db`, which would undo
3,001 rows of corpus curation and re-open the motif-badge split.

`proofreading/mhc_nonhuman.tsv` is the membership test for murine, macaque and rat names, and it
distinguishes two naming levels that are not two spellings: a **molecule** name carries a haplotype
(`H2-IAb`), a **chain** name is an MGI symbol that does not (`H2-Ab1`), so a chain name cannot be
mapped to a molecule without the paired column. Leave chain names as chain names.

**Class correspondence** follows from the `mhc.a` gene, deterministically, with no lookup -
`proofreading/mhc.md` §6. Correct `mhc.class` to match the gene rather than the reverse.

## 4. Method vocabulary

Vocabulary: [`docs/standards/chunk-format.md`](../../docs/standards/chunk-format.md) method columns.
Effect on the score: [`docs/standards/confidence-score.md`](../../docs/standards/confidence-score.md).

**Record what the source says.** Hyphenate and lowercase to the recognised term, and change nothing
else. "tetramer sort" → `tetramer-sort`. "multimer", or "sorted" with no reagent named → `multimer-sort`,
even though tetramers are commoner in VDJdb - the reagent type is not inferable from prevalence.

Fix the mechanical mistakes:

| Wrong | Right | Reason |
|---|---|---|
| a sequencer name in `method.sequencing` (`illumina`, `miseq`) | `amplicon-seq` | those are platforms |
| `RNA-seq`, `Single cell` | `rna-seq`, and set `method.singlecell=yes` | casing, and a separate field |
| a sort term in `method.verification` | the stain form (`tetramer-sort` → `tetramer-stain`) | identification is how cells were found, verification is how the cloned TCR was re-tested |
| software in `method.verification` (`mixcr`, `cellranger`) | blank | not a verification method |
| `antigen-coated-targets` | `antigen-loaded-targets` | recognised term |
| a percentage in `method.frequency` | see below | the field is a count over a total |

A percentage in `method.frequency` needs a decision, not a conversion. If the same value repeats
across every clone for one epitope it is a group-level figure - move it to `meta.subset.frequency` and
blank `method.frequency`. If it varies per clone it may be a repertoire fraction - convert to `N/M`
only when the paper gives the denominator, otherwise leave it blank and log it. Never invent a
denominator; `method.frequency` feeds the score.

A method with no close term is left as the author's wording, logged under vocabulary gaps, and
proposed as a new term. Do not force it into an existing one.

## 5. Reference ids

`https://doi.org/…` and `http://dx.doi.org/…` → `doi:…`. `doi: 10.…` → `doi:10.…`. `pubmed:` or
`PubMed:` → `PMID:`. A bare number is only a PMID once the user confirms it. Case matters: `vdjdb qc`
accepts `PMID:`, `doi:`, `http://`, `https://` and any casing of `unpublished`, and nothing else.

## 6. Antigen fields

`antigen.gene` and `antigen.species` are `/vdjdb-harmonize`'s subject. Run it rather than duplicating
it here. If the epitope is in `patches/antigen_epitope_species_gene.dict`, that entry wins.

## 7. Renumber and check

Renumber `chunk.id` from 1. Then run the two commands and read them, rather than re-deriving what
they report:

```bash
uv run vdjdb qc <file> --report out/reports/qc.tsv
uv run vdjdb submission <file>
```

`vdjdb submission`'s "values new to VDJdb" table is the consistency check this skill used to do with
`cut` and `sort -u`, computed against the whole corpus instead of a grep: a spelling this file
invented shows up there as a new value, and for an established epitope that is the signal to look
again. A genuinely new epitope also shows up there, which is why it is a list and not a gate.

## Output

Write the chunk as `PMID_<pubmed_id>.tsv` where there is a PMID; check the name is not already in
`chunks/`. It goes to `chunks/` only after `/vdjdb-proofread` passes.

`<basename>_format_log.txt` records, per change: the field, the old value, the new value, and **which
authority row or query decided it**. Plus unresolvable values; vocabulary gaps; alleles that exist
only at low resolution; alleles marked `Unconfirmed` in `mhc_alleles.tsv.gz`; and `ORF`/`P`
functionality calls.

## Next step

`/vdjdb-proofread [file]`.
