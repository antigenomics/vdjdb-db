# The chunk format

A chunk is one publication, stored as `chunks/PMID_<id>.tsv` with one record per row. A
record reports paired chains: the alpha and the beta of one clone are columns of the same
row. `chains` is derived from that, never the other way round.

Reports from different publications remain separate, including repeated receptor assays.
Different reference identifiers alone do not establish independent capture or validation.
Check shared experiments, reused constructs and author overlap before assigning that claim
({doc}`../submission`).

## Publication filenames

Use `PMID_<id>.tsv` when PubMed indexes the paper, including indexed preprints. Otherwise use
`DOI_<doi>.tsv`, replacing the DOI slash with `_`, for example
`DOI_10.1101_2021.09.09.459584.tsv`. Strip `https://doi.org/`, publisher URLs,
`/content/`, version suffixes such as `v1`, and page suffixes such as `.full` before naming.
The submitted `reference.id` uses the verified identifier; a filename is not a reference identifier.

Check that the target filename is absent before renaming. Do not overwrite another chunk or
combine two submissions without reviewing their source observations. Perform mechanical renames
on a data-issue branch, retain the source bytes, run `vdjdb identity update`, and check that
record IDs and content hashes are preserved. Source paths and their natural keys can change.

## The declared column order

The 33 columns below, in this order, are what
[`template.tsv`](https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/template.tsv)
carries, and both are generated from the field registry so neither can describe a column set the
build does not read. Chunk columns keep their **dotted** names: the tidy tables ship
`underscore_case`, but a chunk header is the submission contract.

This is an order, not a gate. A submitted chunk may order its columns however it likes and may omit
the optional ones - what is checked is that every column it *does* name is one of these, plus the
three curation columns a chunk may carry and the template does not: `submitter`, `comment` and
`method.pairing`.

```{vdjdb-schema} chunk
:columns: name, title
```

The rest of this page is how to fill each one, which is the part no generated table can carry.

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
mhc.a | First MHC chain allele, to the best resolution available, ``HLA-X*XX:XX``, e.g. ``HLA-A*02:01``. A peptide is presented by one locus, so within one chunk and one ``antigen.epitope`` the HLA gene (A, B, C) should be constant. QC reports a split under ``one epitope under two HLA genes in one chunk``: it is either donor typing written into this column, which is what turned ``RAKFKQLL`` into ``HLA-A*02``, ``HLA-B*08`` and ``HLA-B*07`` in one study (#597), or a paper that reports two restrictions. Only the paper says which
mhc.b | Second MHC chain allele (``B2M`` for MHCI)
mhc.class | ``MHCI`` or ``MHCII``
antigen.epitope | Amino acid sequence of the epitope
antigen.gene | Parent gene of the epitope sequence (e.g. ``pp24``). A property of the peptide, so within one chunk and one ``antigen.epitope`` it should be constant. QC reports a dense run of ``prefix``+integer under ``counter in antigen.gene``: dragging a cell down a spreadsheet column increments it, and that turned one epitope's ``Eef2`` into ``Eef2``..``Eef188``, holding 65 clonotypes apart because the column is part of the deduplication key. The same rule covers ``antigen.species``, ``mhc.a`` and ``mhc.b``, and deliberately not ``meta.clone.id`` or ``meta.subject.id``, where a counter is the content
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
method.frequency | Frequency in the isolated epitope-reactive population, reported as ``X/X`` where possible, e.g. ``7/30`` if a given V/D/J/CDR3 is encountered in 7 out of 30 tetramer+ clones. The population is defined by the pMHC the sort used, not by the antigen it came from - a tetramer is one epitope on one allele. The formats ``X%``, ``X.X%``, ``X.X`` and ``1e-04`` are also supported. Measured over the corpus: 43,231 records write a ratio, 17,700 a percentage, 2,601 a float.
method.frequency.count | **Optional, and preferred over encoding the count in the string above** (#696). Reads, UMIs or cells supporting this clonotype, as an integer. A clonotype supported by 3 reads is better evidence than one supported by 1, and a frequency erases exactly that difference: at any realistic depth both are ``~1/total``, and ``1/33921`` against ``3/33921`` is 2.9e-5 against 8.8e-5, which a two-significant-figure export makes the same number.
method.frequency.total | **Optional.** The sample total ``method.frequency.count`` is out of. Where this pair is left blank the build parses it out of ``method.frequency`` when that is an unambiguous ``x/X``; where it is given, **the submitted value wins and nothing is parsed**. Where all three are present they must agree - ``vdjdb qc`` reports ``frequency disagrees with its count and total`` and does not repair it, because which of the three the paper supports is a curation question.
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
meta.replica.id | Identifier distinguishing replicates, experiments or sampling time points within a study (e.g. ``exp_exploratory``, ``exp_validation_table_s5``, ``5mo``). Use distinct values for separate measurements of the same TCR–pMHC; keep donor identity in ``meta.subject.id`` and assay details in ``method.*``. Prefer identifiers reported by the paper; if assigning experiment labels, document which experiment each label denotes.
meta.clone.id | T-cell clone id
meta.epitope.id | Epitope id (e.g. ``FL10``)
meta.tissue | Tissue used to isolate T-cells: ``PBMC``, ``spleen``, etc. or ``TCL`` (T-cell culture) if isolated from re-stimulated T-cells
meta.donor.MHC | Donor MHC list if available, blank otherwise. IMGT nomenclature (e.g. HLA-A*02:01) is preferable. Allele group names (e.g. ``A02``, ``B18``) are also accepted (do not use an asterisk in such cases). Use a comma to separate alleles.
meta.donor.MHC.method | Donor MHC typing method if available, blank otherwise
meta.structure.id | PDB structure ID if one exists, blank otherwise. A record with associated structural data gets the highest confidence score, so the field is **not free text**: a figure or table reference here awards that score for evidence the reader cannot check. QC reports anything that is not a four-character PDB entry id under ``structure id is not a PDB id``.
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

## Completeness and observation identity

Every record requires at least one reported junction (`cdr3.alpha` or `cdr3.beta`) and a
reported `antigen.epitope`. Neither sequence can be inferred. Missing both junctions or a missing
epitope fails strict QC and direct assembly. Keep the source submission outside the shipping
chunks until the missing evidence is supplied; never silently drop the row during a build.

Missing V/J calls do not justify dropping a sequence-bearing observation. Use the built-in
arda/vdjtools annotation and inspect its proposals. Preserve the submitted call separately from
inference; an inferred V proposal is not a paper-reported V call.

Resolve every reported MHC name and `mhc.class` before submission. Class I uses `B2M` as its second
chain. Strict QC and direct assembly reject class I without `B2M`, class II with
`B2M`, and `B2M` in the first-chain field. Assembly also checks both harmonised
chain names against their declared MHC class. For class II, distinguish an explicitly reported pair from a single-chain or haplotype
label. Use the installed mhcmatch naming and partner-inference functions where supported;
record the original label, inferred partner and inference basis in the review. Its
`pseudoseq.class2_key` supports eligible DP/DQ beta-only typings through `alpha_prior`;
unsupported or ambiguous typings remain unresolved. A prediction of peptide binding is not
proof of the restriction reported by a paper. Never infer restriction solely from peptide length.

A source reporting only one class-II chain can be submitted with its partner blank.
QC reports a partial restriction; the build preserves the reported chain and does
not invent its partner. Both chains absent, invalid nonempty names and incomplete
class-I pairs still fail. Unknown peptide parent genes stay blank and are reported
as missing provenance. Bovine receptors use `BosTaurus`; unavailable germline
annotations stay empty and do not borrow another species' model.

For duplicate review, compare both chains together: `v.alpha`, `j.alpha`, `cdr3.alpha`,
`v.beta`, `j.beta`, `cdr3.beta`, epitope and both MHC chains, within species. Compare the complete
observation metadata for every matching group, within and across files. Different references,
donors, methods, subsets, tissues, clone IDs or other reported metadata distinguish observations.
File boundaries alone do not establish independence: the same paper can occur in an aggregate
and its own chunk. Merge only confirmed duplicates, preserving complementary information and
recording the decision. Review both submitted and harmonised values without silently replacing
the source values.

Large chunk files may use `.tsv.gz`; validation reads the same TSV contents after decompression.
Negative observations remain in `chunks_negative/` and are excluded from the positive build.
