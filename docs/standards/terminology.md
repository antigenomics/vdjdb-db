# Terminology: antigen, epitope, and what a receptor recognises

This page exists because confounding *antigen* with *epitope* is the most common imprecision in this
field, and reviewers are right to flag it. VDJdb's own column names carry the loose usage for
historical reasons, so the distinction has to be written down somewhere rather than inferred from
them.

## The two words

**Antigen** - in a T-cell context, the **peptide-MHC complex**. That is what the receptor engages: a
peptide held in the groove of a particular allele. Neither half alone is the antigen.

The word itself is used three ways in print, so the section below states which one this page means and
what the field's own counts are. What no source disputes is the second half of the distinction, which is
the part the database is keyed on: the presenting allele is not optional, and the peptide's origin is
not the antigen.

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

## What the field's own words are, measured

Checked against PubMed on 2026-09-30, because a page about precision should not assert a convention
without knowing whether the field shares it. It largely does not, on this one word:

| Phrase, title and abstract, singular and plural summed | Records |
|---|---|
| `peptide-MHC complex` | 956 |
| `peptide-MHC ligand` | 63 |
| `peptide-MHC antigen` | 12 |

So the field's name for the composite is **complex**, or **ligand**, and "antigen" most often means the
peptide alone. All three senses are in print: antigen as the peptide (the majority), antigen as the
peptide with the epitope being peptide-plus-allele (the inverse of this page), and antigen as the source
protein.

Two consequences for how to write, and neither weakens the distinction above:

- **Where a sentence has to be unassailable, name the object rather than the word**: "the receptor
  engages a peptide-MHC complex" is supported by every source checked and cannot be read three ways.
  Reserve `cognate` for the pair actually measured, `cross-reactive` for further ligands.
- **Say which sense you mean** the first time "antigen" appears in a document. This page means the
  complex. A paper that means the peptide is not wrong, it is the more common usage.

The source gene and species describe peptide provenance.

One thing to carry from the cross-reactivity literature, because it bounds the claim in both
directions. A receptor recognises many peptides, and they usually share a recognition motif rather than
being unrelated: five mouse and human receptors each selected hundreds of reactive peptides whose motifs
resembled the known antigen, which the authors present as surveillance without degenerate recognition of
non-homologous peptides. So "one clonotype, one specificity" is wrong, and so is treating a receptor as
unboundedly promiscuous. The assumption that a receptor has exactly one cognate epitope is named in
*Nature Reviews Immunology* as still embedded in the pre-processing of many prediction models, which is
the same error as collapsing a provenance union onto a receptor.

### References

Retrieved from PubMed, 2026-09-30.

- Sewell AK. Why must T cells be cross-reactive? *Nat Rev Immunol* 2012;12(9):669-77.
  [10.1038/nri3279](https://doi.org/10.1038/nri3279)
- Hudson D, Fernandes RA, Basham M, Ogg G, Koohy H. Can we predict T cell specificity with digital
  biology and machine learning? *Nat Rev Immunol* 2023;23(8):511-521.
  [10.1038/s41577-023-00835-3](https://doi.org/10.1038/s41577-023-00835-3)
- Birnbaum ME, Mendoza JL, Sethi DK, et al. Deconstructing the peptide-MHC specificity of T cell
  recognition. *Cell* 2014;157(5):1073-87.
  [10.1016/j.cell.2014.03.047](https://doi.org/10.1016/j.cell.2014.03.047)

## Shorthand is fine once the precise form is written down

Prose in this repository asks "is the `CAS` motif specific to HIV-1, or to its TRBV?" That question is
understandable and short, and it is also imprecise on both halves. Stated properly it is:

> Is the `CAS` motif specific to any of HIV-1's epitopes, or to its TRBV germline?

Two corrections, and they are the two this page is about. **HIV-1 is provenance**: no receptor is
specific for a species, so what the group can be specific to is some epitope of that species, which is
why the word is "any of". And the alternative to the motif is the **germline** the V gene templates,
not the gene as a label, because that is the thing that would produce `CAS` without any epitope being
involved.

The short form may stand where it already does. What must not happen is the short form being the only
record, because then nobody can tell a shorthand from a claim.

## The take-home

Be precise, and the precision costs nothing: ask the provenance question freely, then say what the
answer is a set of. Write "epitope-specific" or name the pMHC when the claim is about recognition;
reserve "antigen" for the complex; call a gene or a species what it is, which is where the peptide
came from. `proofreading/epitope_proteome.tsv`
exists to make that provenance checkable - and, because a reference proteome is one genome and a
patient cohort is not, a peptide differing from it is a description of the peptide and never a verdict
on the record.
