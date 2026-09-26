# The chunk format

A chunk is one publication: `chunks/PMID_<id>.txt`, one record per row, and a record
reports **paired chains** -- the alpha and the beta of one clone are columns of the same
row. `chains` is derived from that, never the other way round.

Two rows in two different chunks are **independent reports**, never duplicates, even when
every field matches: independent replication is a signal, and it is what the motif stage
is tuned against ({doc}`../denoising`).

## Complex information columns (required)

These columns convey full information about TCR:peptide:MHC complex and are mandatory for any submission.

column name     | description
----------------|-------------
cdr3.alpha | TCR alpha CDR3 amino acid sequence. Complete sequence starting with C and ending with F/W should be provided if possible. Trimmed sequences will be fixed at database building stage in case sufficient V/J germline parts are present
v.alpha | TCR alpha Variable (V) segment id, up to best resolution possible (``TRAVX*XX``, e.g. ``TRAV7``, ``TRAV7*01``, ``TRAV7*02``...). Strictly IMGT nomenclature. Can be left blank if unknown.
j.alpha | TCR alpha Joining (J) segment id
cdr3.beta | TCR beta CDR3 amino acid sequence
v.beta | TCR beta V segment id
j.beta | TCR beta J segment id
species | TCR parent species (``HomoSapiens``, ``MusMusculus``,...)
mhc.a | First MHC chain allele, to best resolution possible, ``HLA-X*XX:XX``, e.g. ``HLA-A*02:01``
mhc.b | Second MHC chain allele (``B2M`` for MHCI)
mhc.class | ``MHCI`` or ``MHCII``
antigen.epitope | Amino acid sequence of the epitope
antigen.gene | Parent gene of the epitope sequence (e.g. ``pp24``)
antigen.species | Parent species of the antigen, to the best clade resolution possible (e.g. ``HIV-1``, ``HIV-1*HXB2``)
reference.id | Pubmed id, doi, etc
submitter | Name of submitting person/organization

> **Notes:**

> In case given record represents a clonotype with either TCR alpha or beta sequence unknown, missing CDR3/V/(D)/J fields should be left blank.

> V/(D)/J fields can be left blank, however this will abrogate CDR3 fixing/verification procedure for a given record.

> Any record should have at least one of CDR3 alpha/beta fields that are not blank.

## Method information columns (optional)

Optional columns (i.e. it is not required to fill them, but they **should** be present in table header) that ensure correct confidence ranking of a given entry. Used to calculate a single confidence score based on various factors, e.g. fraction of a given TCRab sequence among tetramer+ clones sequenced and verification experiments performed.

column name     | description
----------------|-------------
method.identification | ``tetramer-sort``, ``dextramer-sort``, ``pelimer-sort``, ``pentamer-sort``, etc for sorting-based identification. For molecular assays use: ``antigen-loaded-targets`` (if T cells specificity was analysed against cells incubatetd with antigenic peptide), ``antigen-expressing-targets`` (if T cells specificity was analysed against cells tranformed with antigenic organism, protein or peptide, e.g. BCL transformed with EBV). For magnetic cell separation use ``beads`` keyword. Add ``cultured-T-cells`` or ``limiting-dilution-cloning`` if T cells were cultured before sequencing as in this case ``method.frequency`` will have completely different meaning. Use comma to separate phrases. For cases that use UMI-tagged multimers use ``tetramer-umi``, etc.
method.frequency | Frequency in isolated antigen-specific population, reported as ``X/X`` if possible, e.g. ``7/30`` if a given V/D/J/CDR3 is encountered in 7 out of 30 tetramer+ clones. Formats ``X%``, ``X.X%`` and ``X.X`` are also supported.
method.singlecell | ``yes`` if single cell sequencing was performed, blank otherwise
method.sequencing | Sequencing method: ``sanger``, ``rna-seq`` or ``amplicon-seq``
method.verification | ``tetramer-stain``, ``dextramer-stain``, ``pelimer-stain``, ``pentamer-stain``, etc for methods that include TCR cloning and re-staining with multimers. For magnetic cell separation use ``beads`` keyword. ``restimulation``, ``co-culture``, ``antigen-loaded-targets``, ``antigen-expressing-targets`` for molecular assays that validate specificity of **cloned** T-cell receptors. ``direct`` in case the affinity of TCRs of specific T-cells to the pMHC is quantified directly in some way. Several comma-separated verification methods can be specified.

> **Notes:**

> In case ``method.identification`` is left blank, the record is automatically assigned with a lowest confidence score possible.

> For special cases such as CD8-null tetramers that utilize HLA with mutated residues that abrogate CD8 binding, specify ``cd8null-tetramer`` in ``method.identification`` field rather than using ``mhc.a`` field.

During database build phase, the information from columns mentioned above is collapsed to a JSON string and stored in a single ``method`` column, e.g.:
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
meta.subset.frequency | Frequency of a given TCR sequence in specified cell subset, e.g. ``5%`` means that the TCR sequence represents an expanded clone occupying 5% of CD8+ cells
meta.subject.cohort | Subject cohort, free style, e.g. ``healthy`` or ``HIV+``. If possible, specify to what extent a healthy donor is healthy, e.g. ``CMV-seronegative``.
meta.subject.id | Subject id (e.g. ``donor1``, ``donor2``,...)
meta.replica.id | Replicate sample coming from the same donor, also applies for different time points, etc (e.g. ``5mo``)
meta.clone.id | T-cell clone id
meta.epitope.id | Epitope id (e.g. ``FL10``)
meta.tissue | Tissue used to isolate T-cells: ``PBMC``, ``spleen``, etc. or ``TCL`` (T-cell culture) if isolated from re-stimulated T-cells
meta.donor.MHC | Donor MHC list if available, blank otherwise. IMGT nomenclature (e.g. HLA-A*02:01) is preferable. Allele group names (e.g. ``A02``, ``B18``) is also acceptable (don't use asterisk in such cases). Use comma to separate alleles.
meta.donor.MHC.method | Donor MHC typing method if available, blank otherwise
meta.structure.id | PDB structure ID if exists, or blank. Records having a structural data associated with them will automatically get the highest confidence score.
comment | Plain text comment, maximum 140 characters

> **Note:**

> While these columns are optional, subject identifier, replica identifier, etc are used when scanning submission for duplicates. Normally duplicate records (with identical **complex information** columns) are not allowed, but they will not be considered as duplicates in case they have distinct id fields mentioned above.

During database build phase, the information from columns mentioned above is collapsed to a JSON string and stored in a single ``meta`` column, e.g.:
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

Peptide pools, long peptides for T-cell culture expansion, non-peptide ligands

column name    | description
---------------|------------
antigen.epitope.long | encompassing protein sequence containing the epitope
antigen.peptide.pool | e.g. ``MIRA COVID19`` TBD
antigen.nonpeptide | ``α-GalCer`` or ``KRN7000`` TBD

## Non TRAB columns (for extended database, TBA)

Information for non alpha-beta T-cells, CAR-T, etc

column name    | description
---------------|------------
v.delta | ID of Variable segment in delta chain
cdr3.delta | CDR3 of delta chain
... | ...
v.heavy.shm | CIGAR string of hypermutations in the heavy chain Variable segment
... | ...
