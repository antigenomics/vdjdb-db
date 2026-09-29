# Terminology: antigen, epitope, and what a receptor recognises

This page exists because confounding *antigen* with *epitope* is the most common imprecision in this
field, and reviewers are right to flag it. VDJdb's own column names carry the loose usage for
historical reasons, so the distinction has to be written down somewhere rather than inferred from
them.

## The two words

**Antigen** - in a T-cell context, the **peptide-MHC complex**. That is what the receptor engages: a
peptide held in the groove of a particular allele. Neither half alone is the antigen.

**Epitope** - the part of the antigen that is specifically recognised. The parallel with the antibody
field is exact: there, the antigen is the molecule and the epitope is the surface the antibody binds.
Here the antigen is the pMHC and the epitope is the portion of it that determines recognition.

So a receptor **does** recognise an antigen. What it does not recognise is a pathogen, a gene, a
protein, or a species. Those are the *provenance* of the peptide - where it came from.

## Provenance is a real axis, and most questions start there

"Which T cells were found against SARS-CoV-2 epitopes?" and "which against human neoantigens?" are
exactly the questions `antigen.species` and `antigen.gene` exist to answer, and they are good
questions. Grouping epitopes by where they came from is how anyone approaches this data, and the
database would be far less useful without it.

What such a group *is*, precisely, is a **union over pMHCs**: every receptor reported against some
epitope of that origin, under whatever restriction each study used. That is well defined, and it is
not one specificity.

The error is only in collapsing the union:

| Statement | Status |
|---|---|
| "this TCR recognises `GILGFVFTL` presented by `HLA-A*02:01`" | precise |
| "receptors in VDJdb reported against SARS-CoV-2 epitopes" | precise, and a union over many pMHCs |
| "T cells recognising human neoantigens" | precise as a set; the members were shown different antigens |
| "this TCR is specific for the influenza M1 protein" | loose: M1 yields many epitopes under many restrictions, and the receptor was shown one |
| "this TCR is influenza-specific" | a property claimed of the receptor that the data does not carry |
| "this TCR recognises the NY-ESO-1 antigen" | conflates provenance with the antigen. The antigen was a peptide on an allele; `NY-ESO-1` is the gene it came from |

The practical consequence is about **aggregates, not about asking**. A motif, a lift or an enrichment
computed over a species group pools receptors that were shown different antigens, so the number
describes the group and cannot be read as a motif *for* that pathogen. `vdjdb.corpus` documents `a:`
on exactly those terms: reach for it to find the papers and the receptors, and condition on
`e:<epitope>` - or on that plus a restriction - when the claim is about recognition.

## How VDJdb's columns map onto this

The `antigen.` prefix predates the distinction and is a published contract, so the names stay. What
each column actually holds:

| Column | What it is | Under the strict terms |
|---|---|---|
| `antigen.epitope` | the peptide sequence | half of the antigen; the epitope is the recognised part of the whole complex |
| `mhc.a`, `mhc.b`, `mhc.class` | the presenting molecule | the other half |
| `antigen.gene` | the gene the peptide came from | **provenance**, not the antigen |
| `antigen.species` | the organism the peptide came from | **provenance**, not the antigen |

The database is keyed accordingly, and this is the part to trust over the naming:

- `restriction` is keyed on `(antigen.epitope, antigen.species, mhc.a, mhc.b)` and `pmhc_id`
  identifies one peptide as presented. Those are the antigen.
- `antigen.gene` and `antigen.species` are deliberately **absent** from the row identity `vdjdb diff`
  keys on, because they annotate the peptide rather than identify the record. A correction to either
  is a changed cell on the same record, never a record retired and another allocated.
- There is no identifier for a gene, a protein or a pathogen anywhere in the schema, and that is not
  an omission.

## The take-home

Be precise, and the precision costs nothing: ask the provenance question freely, then say what the
answer is a set of. Write "epitope-specific" or name the pMHC when the claim is about recognition;
reserve "antigen" for the complex; call a gene or a species what it is, which is where the peptide
came from. `proofreading/epitope_proteome.tsv`
exists to make that provenance checkable - and, because a reference proteome is one genome and a
patient cohort is not, a peptide differing from it is a description of the peptide and never a verdict
on the record.
