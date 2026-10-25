# Search publications and trace a specificity claim

Use this workflow when asked for literature, an explanation of a record's experimental
support, or research context for a downstream task. Reading papers does not authorize
editing source chunks or submitting a correction. Route that work to curator skills.

## Search from the biological question

1. Define species, chain or reported pair, junction sequence, peptide and MHC, plus the
   requested question: assay support, cross-reactivity, neighbourhoods or structures.
   Missing restriction stays missing. Start from the record's `reference_id`; retrieve
   its PMID/DOI before expanding to related work.
2. Use PubMed for biomedical publication discovery and metadata. Search exact peptide
   with receptor/MHC terms, then peptide aliases or source protein if needed. Exact
   receptor sequences are often found only in supplements or deposited tables, so a
   zero-hit sequence search does not establish absence. Example queries:

   ```text
   "GILGFVFTL" AND ("T-cell receptor" OR TCR)
   "GILGFVFTL" AND ("HLA-A*02:01" OR "HLA-A2")
   TCREMP AND receptor
   (TCRNET OR ALICE) AND "T-cell receptor"
   TCRen AND structure
   ```

3. Retrieve abstracts and, where available, full text, supplements and linked public
   repositories. Search methods/results for the receptor construct, chain pairing,
   ligand, restriction and actual readout. Separate multimer binding, activation,
   cytotoxicity, affinity measurements and structural observations. Engineered receptor
   or peptide variants are separate constructs, even when named after the same clone.
4. Follow references and later independent reports for that pMHC when relevant. Keep
   assay evidence, repeated reporting and reused datasets separate. A publication that
   imports the same original observations is not a new validation experiment.
5. Return a compact evidence table: record/query ID, source PMID/DOI, receptor/ligand
   construct, assay and reported result, supporting figure/table, and unresolved context.
   Preserve experimental negatives if provided; untested combinations are unlabeled.

Record the search date, exact queries, inspected source identifiers and database snapshot.
State whether this was a targeted search or a systematic review; do not claim exhaustive
coverage from a first page of search results. Page through results needed for the task,
and disclose unavailable or truncated sources. Copy citation metadata from retrieved
records, retaining PMID and DOI. Cite the supporting publication next to the claim.

## Verified starting literature

Selected PubMed records and abstracts checked on 2026-10-25. This is a task-oriented
reading set, not a systematic review. Paper titles retain their published terminology;
use the database's current pMHC terminology when interpreting a record.

| Paper | Why read it | Identifiers |
|---|---|---|
| VDJdb: a curated database of T-cell receptor sequences with known antigen specificity. | Published observations and repertoire annotation | [PMID 28977646](https://pubmed.ncbi.nlm.nih.gov/28977646/), [DOI](https://doi.org/10.1093/nar/gkx760) |
| VDJdb in 2019: database extension, new analysis infrastructure and a T-cell receptor motif compendium. | Batch lookup and the motif compendium | [PMID 31588507](https://pubmed.ncbi.nlm.nih.gov/31588507/), [DOI](https://doi.org/10.1093/nar/gkz874) |
| VDJdb in 2026: boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure prediction. | Embedding and structural evidence in the updated resource | [PMID 42773777](https://pubmed.ncbi.nlm.nih.gov/42773777/), [DOI](https://doi.org/10.1093/nar/gkag904) |
| A Framework for Annotation of Antigen Specificities in High-Throughput T-Cell Repertoire Sequencing Studies. | Combining neighbourhood enrichment, HLA context and database matches; inspect expanded clones separately | [PMID 31616409](https://pubmed.ncbi.nlm.nih.gov/31616409/), [DOI](https://doi.org/10.3389/fimmu.2019.02159) |
| Neighborhood enrichment for the identification of antigen-specific T-cell receptors. | Comparative ALICE/TCRNET/GLIPH2/tcrdist3 workflows; check the studied populations before transferring a result | [PMID 40996146](https://pubmed.ncbi.nlm.nih.gov/40996146/), [DOI](https://doi.org/10.1093/bib/bbaf495) |
| TCREMP: A Bioinformatic Pipeline for Efficient Embedding of T-cell Receptor Sequences from Immune Repertoire and Single-cell Sequencing Data. | Prototype-based receptor representation and single-cell/repertoire examples | [PMID 40368275](https://pubmed.ncbi.nlm.nih.gov/40368275/), [DOI](https://doi.org/10.1016/j.jmb.2025.169205) |
| Structure-based prediction of T cell receptor recognition of unseen epitopes using TCRen. | Structure-based candidate ranking and contact preferences | [PMID 38987378](https://pubmed.ncbi.nlm.nih.gov/38987378/), [DOI](https://doi.org/10.1038/s43588-024-00653-0) |
| VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. | Repertoire statistics, overlap and cohort post-analysis | [PMID 26606115](https://pubmed.ncbi.nlm.nih.gov/26606115/), [DOI](https://doi.org/10.1371/journal.pcbi.1004503) |

These papers motivate the downstream routes; they do not specify the current CLI,
score-model version or release schema. Use the selected software version's documentation
for those contracts. Do not attribute a later implementation to an older paper.
