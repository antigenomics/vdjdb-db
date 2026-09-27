# The chunk format

A chunk is one publication, stored as `chunks/PMID_<id>.txt` with one record per row. A
record reports paired chains: the alpha and the beta of one clone are columns of the same
row. `chains` is derived from that, never the other way round.

Two rows in two different chunks are independent reports, never duplicates, even when
every field matches. The motif stage is tuned against that replication
({doc}`../denoising`).

## Complex information columns (required)

These columns describe the TCR:peptide:MHC complex and are mandatory in any submission.

column name     | description
----------------|-------------
cdr3.alpha | TCR alpha CDR3 amino acid sequence. Give the complete sequence, starting with C and ending with F/W, where possible. Trimmed sequences are fixed at the build stage when sufficient V/J germline parts are present
v.alpha | TCR alpha Variable (V) segment id, to the best resolution available (``TRAVX*XX``, e.g. ``TRAV7``, ``TRAV7*01``, ``TRAV7*02``...). Strictly IMGT nomenclature. May be left blank if unknown.
j.alpha | TCR alpha Joining (J) segment id
cdr3.beta | TCR beta CDR3 amino acid sequence
v.beta | TCR beta V segment id
j.beta | TCR beta J segment id
species | TCR parent species (``HomoSapiens``, ``MusMusculus``,...)
mhc.a | First MHC chain allele, to the best resolution available, ``HLA-X*XX:XX``, e.g. ``HLA-A*02:01``
mhc.b | Second MHC chain allele (``B2M`` for MHCI)
mhc.class | ``MHCI`` or ``MHCII``
antigen.epitope | Amino acid sequence of the epitope
antigen.gene | Parent gene of the epitope sequence (e.g. ``pp24``)
antigen.species | Parent species of the antigen, to the best clade resolution available (e.g. ``HIV-1``, ``HIV-1*HXB2``)
reference.id | Pubmed id, doi, etc
submitter | Name of submitting person/organization

> **Notes:**

> If a record represents a clonotype whose alpha or beta sequence is unknown, leave the missing CDR3/V/(D)/J fields blank.

> V/(D)/J fields may be left blank, in which case the CDR3 fixing and verification procedure is skipped for that record.

> Every record must have at least one of ``cdr3.alpha`` and ``cdr3.beta`` filled.

## Method information columns (optional)

These columns are optional to fill, but should be present in the table header. They set the confidence ranking of an entry: a single confidence score is computed from factors such as the fraction of a given TCRab sequence among the tetramer+ clones sequenced and the verification experiments performed.

column name     | description
----------------|-------------
method.identification | ``tetramer-sort``, ``dextramer-sort``, ``pelimer-sort``, ``pentamer-sort``, etc. for sorting-based identification. For molecular assays use ``antigen-loaded-targets`` (T cell specificity analysed against cells incubated with antigenic peptide) or ``antigen-expressing-targets`` (T cell specificity analysed against cells transformed with an antigenic organism, protein or peptide, e.g. BCL transformed with EBV). For magnetic cell separation use ``beads``. Add ``cultured-T-cells`` or ``limiting-dilution-cloning`` if T cells were cultured before sequencing, since ``method.frequency`` then has a different meaning. For UMI-tagged multimers use ``tetramer-umi``, etc. Separate phrases with a comma.
method.frequency | Frequency in the isolated antigen-specific population, reported as ``X/X`` where possible, e.g. ``7/30`` if a given V/D/J/CDR3 is encountered in 7 out of 30 tetramer+ clones. The formats ``X%``, ``X.X%`` and ``X.X`` are also supported.
method.singlecell | ``yes`` if single cell sequencing was performed, blank otherwise
method.sequencing | Sequencing method: ``sanger``, ``rna-seq`` or ``amplicon-seq``
method.verification | ``tetramer-stain``, ``dextramer-stain``, ``pelimer-stain``, ``pentamer-stain``, etc. for methods that include TCR cloning and re-staining with multimers. For magnetic cell separation use ``beads``. ``restimulation``, ``co-culture``, ``antigen-loaded-targets``, ``antigen-expressing-targets`` for molecular assays that validate the specificity of cloned T-cell receptors. ``direct`` if the affinity of the TCR of a specific T cell to the pMHC is quantified directly. Several comma-separated verification methods may be given.

> **Notes:**

> If ``method.identification`` is left blank, the record is assigned the lowest confidence score possible.

> For special cases such as CD8-null tetramers, which use HLA with mutated residues that abrogate CD8 binding, specify ``cd8null-tetramer`` in ``method.identification`` rather than using the ``mhc.a`` field.

The build collapses the columns above into a JSON string held in a single ``method`` column, e.g.:
```json
{
   "identification":"tetramer-sort",
   "frequency":"5/13",
   "sequencing":"sanger",
   "verification":"antigen-loaded-targets"
}
```

## Meta-information columns (optional)

column name     | description
----------------|-------------
meta.study.id | Internal study id
meta.cell.subset | T-cell subset, free style, e.g. ``CD8+``, ``CD4+CD25+``
meta.subset.frequency | Frequency of a given TCR sequence in the specified cell subset, e.g. ``5%`` means the TCR sequence represents an expanded clone occupying 5% of CD8+ cells
meta.subject.cohort | Subject cohort, free style, e.g. ``healthy`` or ``HIV+``. Where possible, specify to what extent a healthy donor is healthy, e.g. ``CMV-seronegative``.
meta.subject.id | Subject id (e.g. ``donor1``, ``donor2``,...)
meta.replica.id | Replicate sample coming from the same donor, also used for different time points, etc (e.g. ``5mo``)
meta.clone.id | T-cell clone id
meta.epitope.id | Epitope id (e.g. ``FL10``)
meta.tissue | Tissue used to isolate T-cells: ``PBMC``, ``spleen``, etc. or ``TCL`` (T-cell culture) if isolated from re-stimulated T-cells
meta.donor.MHC | Donor MHC list if available, blank otherwise. IMGT nomenclature (e.g. HLA-A*02:01) is preferable. Allele group names (e.g. ``A02``, ``B18``) are also accepted (do not use an asterisk in such cases). Use a comma to separate alleles.
meta.donor.MHC.method | Donor MHC typing method if available, blank otherwise
meta.structure.id | PDB structure ID if one exists, blank otherwise. A record with associated structural data gets the highest confidence score.
comment | Plain text comment, maximum 140 characters

> **Note:**

> These columns are optional, but the subject identifier, replica identifier and the other id fields above are used when scanning a submission for duplicates. Duplicate records, those with identical **complex information** columns, are not allowed, but they are not treated as duplicates when they have distinct id fields.

The build collapses the columns above into a JSON string held in a single ``meta`` column, e.g.:
```json
{
   "cell.subset":"CD8+",
   "subject.cohort":"HSV-2+",
   "subject.id":12,
   "clone.id":46,
   "tissue":"PBMC"
}
```

## Condition association columns (for extended database, TBA)

Condition metadata:

column name    | description
---------------|------------
condition.name | natural language terms like ``T1D``, ``pollen allergy``, ``BRCA`` or ``YF vaccination``
condition.id   | ``ICD-11:5A10`` for ``T1D`` in [ICD-11](https://icd.who.int/browse11/l-m/en) or ``OMIM:114480`` for ``breast cancer`` in [OMIM](https://www.omim.org/entry/114480)
condition.type | ``infection``, ``vaccination``, ``cancer``, ``allergy`` or ``autoimmune``
condition.subtype | natural language terms like ``acute`` or ``poor prognosis`` or ``grade II``

Association metadata:

column name    | description
---------------|------------
condition.freq | fraction of samples matching the entry
condition.count | number of samples matching the entry (can be blank)
population.freq | fraction of controls matching the entry, or Pgen computed by [OLGA/IgOR](https://github.com/statbiophys/OLGA/tree/master/olga)
population.count | number of controls matching the entry (can be blank)
association.pvalue | Association P-value, e.g. enrichment P-value for Fisher's exact test
association.test | ``Fisher``, ``TCRNET``, ``ALICE`` or another statistical method

## Ambiguous antigens (for extended database, TBA)

Columns for peptide pools, long peptides used in T-cell culture expansion, and non-peptide ligands.

column name    | description
---------------|------------
antigen.epitope.long | encompassing protein sequence containing the epitope
antigen.peptide.pool | e.g. ``MIRA COVID19`` TBD
antigen.nonpeptide | ``α-GalCer`` or ``KRN7000`` TBD

## Non TRAB columns (for extended database, TBA)

Columns for non alpha-beta T-cells, CAR-T, etc.

column name    | description
---------------|------------
v.delta | ID of Variable segment in delta chain
cdr3.delta | CDR3 of delta chain
... | ...
v.heavy.shm | CIGAR string of hypermutations in the heavy chain Variable segment
... | ...
