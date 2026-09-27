"""The field registry -- the single declaration every column order is projected from.

Before this module, the column set was written out in nine places (``ChunkQC.py``,
``DefaultDBGenerator.py``, ``SlimDBGenerator.py``, ``ScoreFactory.py``, the two static
``*.meta.txt`` files, ``BuildDatabase.groovy``, the README and ``template.xls``) and had already
drifted apart. Everything here is a projection of :data:`FIELDS`:

* ``vdjdb.meta.txt`` / ``vdjdb.slim.meta.txt`` -- :func:`render_meta`, :func:`render_slim_meta`
* the ``vdjdb.txt`` header -- :func:`header`, restoring ``BuildDatabase.groovy:411``
  (``METADATA_LINES[1..-1].collect { it.split("\\t")[0] }``), the invariant the Python port dropped
* the four positional column orders the release and ``vdjdb-web`` depend on
* the docs tables and the AIRR mapping, in later phases

Three defects in the shipped metadata are fixed here deliberately, each with a declared rule in
``rules/expected_diffs.toml``:

1. no ``TCR_hash`` row, though ``vdjdb.txt`` has had the column for years;
2. ``vdjdb.score`` and ``TCR_hash`` placed after ``cdr3fix`` in the metadata but before ``method`` in
   the data -- so the metadata does not describe the file it belongs to, in the release and in
   production;
3. the four ``web.*`` rows have one value too many, landing ``0`` in ``data.type``, ``factor`` in
   ``title`` and ``Internal`` in ``comment``. Corrected to ``data.type = factor``. They are
   ``visible = 0``, so nothing user-facing moves.

The ``web.method`` row additionally has a space where a tab belongs. That one is a git-edit
regression rather than a generator defect -- ``BuildDatabase.groovy:367`` emits a proper tab -- and
it disappears the moment the file is generated again.
"""
from __future__ import annotations

from dataclasses import dataclass

TXT, SEQ = "txt", "seq"


@dataclass(frozen=True, slots=True)
class Field:
    """One column, with the seven ``vdjdb.meta.txt`` attributes ``vdjdb-web`` builds its schema from."""

    name: str
    type: str = TXT
    visible: int = 1
    searchable: int = 1
    autocomplete: int = 1
    data_type: str = "factor"
    title: str = ""
    comment: str = ""
    #: The AIRR field this column maps to, or empty. Declared here so the AIRR emitter and the
    #: docs mapping table are two projections of one statement rather than two hand-written lists.
    #: Empty is not "no counterpart exists" but "none that is the same quantity": VDJdb's
    #: `v.end` is an amino-acid offset in junction space and AIRR's `v_sequence_end` a nucleotide
    #: offset in sequence space, so declaring them equal would be wrong (see `convert.coords`).
    airr: str = ""

    def meta_row(self) -> str:
        return "\t".join((self.name, self.type, str(self.visible), str(self.searchable),
                          str(self.autocomplete), self.data_type, self.title, self.comment))


def _f(name: str, **kw: object) -> tuple[str, Field]:
    return name, Field(name, **kw)  # type: ignore[arg-type]


#: Every column the build knows about, keyed by name. Order here is documentation; the output
#: orders are the tuples below.
FIELDS: dict[str, Field] = dict([
    # -- flat-table identity -------------------------------------------------------------------
    _f("complex.id", visible=0, searchable=0, autocomplete=0, data_type="complex.id",
       title="complex.id",
       comment="TCR alpha and beta chain records having the same complex identifier belong to the "
               "same T-cell clone."),
    _f("gene", airr="locus", title="Gene", comment="TCR chain: alpha or beta."),
    _f("cdr3", airr="junction_aa", type=SEQ, autocomplete=0, data_type="cdr3", title="CDR3",
       comment="TCR complementarity determining region 3 (CDR3) amino acid sequence."),
    _f("v.segm", airr="v_call", title="V", comment="TCR Variable segment allele."),
    _f("j.segm", airr="j_call", title="J", comment="TCR Joining segment allele."),
    _f("species", title="Species", comment="TCR parent species."),
    _f("mhc.a", airr="mhc_allele_1", title="MHC A", comment="First MHC chain allele."),
    _f("mhc.b", airr="mhc_allele_2", title="MHC B",
       comment="Second MHC chain allele (defaults to Beta2Microglobulin for MHC class I)."),
    _f("mhc.class", airr="mhc_class", title="MHC class", comment="MHC class (I or II)."),
    _f("antigen.epitope", airr="peptide_sequence_aa", type=SEQ, data_type="peptide", title="Epitope",
       comment="Amino acid sequence of the epitope."),
    _f("antigen.gene", airr="antigen", title="Epitope gene",
       comment="Representative parent gene of the epitope."),
    _f("antigen.species", airr="antigen_source_species", title="Epitope species",
       comment="Representative parent species of the epitope."),
    _f("reference.id", airr="reactivity_refs", data_type="url", title="Reference",
       comment="Pubmed reference / URL / or submitter details in case unpublished."),
    _f("vdjdb.score", airr="reactivity_value", autocomplete=0, data_type="uint", title="Info",
       comment="VDJdb confidence score, the higher is the score the more confidence we have in the "
               "antigen specificity annotation of a given TCR clonotype/clone. Zero score indicates "
               "that there are insufficient method details to draw any conclusion."),
    # `visible=0` is what production serves (`vdjdb-web/test/resources/database/vdjdb.meta.txt`),
    # and this row reproduces that one field for field. `vdjdb-web` overrides the flag anyway --
    # `DatabaseMetadata.ForcedColumns` copies `visible = true` onto it so the structure viewer can
    # read the hash, and the search page then hides the column again
    # (`search-table.service.ts:198`) -- so the value here only tells a standalone reader that the
    # hash is not a column to display.
    _f("TCR_hash", visible=0, searchable=0, autocomplete=0, data_type=TXT, title="TCR hash",
       comment="SHA256 hash of TCR structure used for structure visualization."),
    _f("method", searchable=0, autocomplete=0, data_type="method.json", title="Method",
       comment="Details on method used to assay TCR specificity."),
    _f("meta", searchable=0, autocomplete=0, data_type="meta.json", title="Meta",
       comment="Various meta-information: cell subset, donor status, etc."),
    _f("cdr3fix", searchable=0, autocomplete=0, data_type="fixer.json", title="CDR3fix",
       comment="Details on CDR3 sequence fixing (if applied) and consistency between V, J and "
               "reported CDR3 sequence."),

    # -- internal filtering fields, never displayed --------------------------------------------
    _f("web.method", visible=0, searchable=0, autocomplete=0, title="Internal",
       comment="Internal: coarse identification method for fast filtering."),
    _f("web.method.seq", visible=0, searchable=0, autocomplete=0, title="Internal",
       comment="Internal: coarse sequencing method for fast filtering."),
    _f("web.cdr3fix.nc", visible=0, searchable=0, autocomplete=0, title="Internal",
       comment="Internal: CDR3 has non-canonical V or J anchor residues."),
    _f("web.cdr3fix.unmp", visible=0, searchable=0, autocomplete=0, title="Internal",
       comment="Internal: CDR3 could not be mapped onto V or J germline."),

    # -- evidence, served by production vdjdb-web; owned by the new format (ROADMAP 9) ---------
    # All five carry `0/0/0`, reproducing the metadata production serves. `autocomplete` is the one
    # that has an effect: `DatabaseColumnInfo.scala:37` copies a column's whole value list into the
    # metadata JSON sent to the browser when it is `1`, and unlike the `web.*` columns these five
    # survive the invisible-column filter, because `DatabaseMetadata.ForcedColumns` forces them
    # through for the Evidence badges. They are badge inputs, not search facets.
    _f("evidence.validation.same.study", visible=0, searchable=0, autocomplete=0,
       title="Validation same study",
       comment="Antigen specificity validated within the same study."),
    _f("evidence.validation.independent", visible=0, searchable=0, autocomplete=0,
       title="Validation independent",
       comment="Antigen specificity independently validated in another study."),
    _f("evidence.structure.native", visible=0, searchable=0, autocomplete=0,
       title="Structure native",
       comment="Native (experimental) TCR-pMHC structure available."),
    _f("evidence.structure.contacts", visible=0, searchable=0, autocomplete=0,
       title="Structure model with contacts",
       comment="Structural model with annotated TCR-pMHC contacts available."),
    _f("evidence.structure.quality", visible=0, searchable=0, autocomplete=0,
       title="Structure good quality model",
       comment="Good-quality structural model available."),

    # -- CDR3 geometry, slim and full only ------------------------------------------------------
    _f("v.end", searchable=0, autocomplete=0, data_type="uint", title="V end",
       comment="Last amino acid position of the V germline part of CDR3, 0-based, junction space."),
    _f("j.start", searchable=0, autocomplete=0, data_type="uint", title="J start",
       comment="First amino acid position of the J germline part of CDR3, 0-based, junction space."),

    # -- chunk columns, paired form -------------------------------------------------------------
    _f("cdr3.alpha", type=SEQ, data_type="cdr3", title="CDR3 alpha"),
    _f("v.alpha", title="V alpha"),
    _f("j.alpha", title="J alpha"),
    _f("cdr3.beta", type=SEQ, data_type="cdr3", title="CDR3 beta"),
    _f("v.beta", title="V beta"),
    _f("d.beta", title="D beta"),
    _f("j.beta", title="J beta"),
    _f("cdr3fix.alpha", searchable=0, autocomplete=0, data_type="fixer.json", title="CDR3fix alpha"),
    _f("cdr3fix.beta", searchable=0, autocomplete=0, data_type="fixer.json", title="CDR3fix beta"),

    # -- assay description ----------------------------------------------------------------------
    _f("method.identification", title="Identification method"),
    _f("method.frequency", title="Frequency"),
    _f("method.singlecell", title="Single cell"),
    _f("method.sequencing", title="Sequencing"),
    _f("method.verification", title="Verification"),
    _f("method.pairing", title="Pairing",
       comment="How alpha and beta chains were paired. Kept for debugging; dropped by the legacy build."),

    # -- sample description ---------------------------------------------------------------------
    _f("meta.study.id", title="Study id"),
    _f("meta.cell.subset", title="Cell subset"),
    _f("meta.subset.frequency", title="Subset frequency",
       comment="Kept for debugging; dropped by the legacy build."),
    _f("meta.subject.cohort", title="Subject cohort"),
    _f("meta.subject.id", title="Subject id"),
    _f("meta.replica.id", title="Replica id"),
    _f("meta.clone.id", title="Clone id"),
    _f("meta.epitope.id", title="Epitope id"),
    _f("meta.tissue", title="Tissue"),
    _f("meta.donor.MHC", title="Donor MHC"),
    _f("meta.donor.MHC.method", title="Donor MHC method"),
    _f("meta.structure.id", title="Structure id"),

    # -- motif tables ----------------------------------------------------------------------------
    # These have no ``.meta.txt`` of their own: ``Motifs.scala`` hands Tablesaw a fixed
    # ``Array[ColumnType]`` and never reads a header. Declared here so the docs tables and the
    # positional-width tests have one source, and so a renamed column fails a test rather than
    # mistyping the file with no error.
    _f("cdr3aa", type=SEQ, data_type="cdr3", title="CDR3",
       comment="CDR3 amino acid sequence of the cluster member."),
    _f("x", data_type="float", title="x", comment="Layout x coordinate of the member in the cluster graph."),
    _f("y", data_type="float", title="y", comment="Layout y coordinate of the member in the cluster graph."),
    _f("cid", title="Cluster id",
       comment="Motif cluster identifier, <species>.<chain>.<epitope>.<n>."),
    _f("csz", data_type="uint", title="Cluster size", comment="Number of members in the cluster."),
    _f("v.segm.repr", title="V representative",
       comment="Modal V allele of the cluster; the PWM is pinned to it."),
    _f("j.segm.repr", title="J representative",
       comment="Modal J allele of the cluster; the PWM is pinned to it."),
    _f("aa", title="Residue", comment="Amino acid at this PWM position."),
    _f("pos", data_type="uint", title="Position", comment="0-based CDR3 position."),
    _f("len", data_type="uint", title="Length", comment="CDR3 length the PWM stratum is pinned to."),
    _f("count", data_type="uint", title="Count",
       comment="Observed occurrences of this residue at this position in the cluster."),
    _f("count.bg", data_type="uint", title="Background count",
       comment="Occurrences in the matched background at the same (V, J, length)."),
    _f("total.bg", data_type="uint", title="Background total",
       comment="Background sequences contributing at this position."),
    _f("count.bg.i", data_type="uint", title="Imputed background count",
       comment="Background count after the Laplace cascade; equals count.bg where no imputation was needed."),
    _f("total.bg.i", data_type="uint", title="Imputed background total",
       comment="Background total after the Laplace cascade."),
    _f("need.impute", title="Imputed",
       comment="Whether this (position, residue) was absent from the matched background and imputed."),
    _f("freq", data_type="float", title="Frequency",
       comment="Within-cluster residue frequency. Sums to 1 across residues at each position."),
    _f("freq.bg", data_type="float", title="Background frequency",
       comment="Background residue frequency at the same position."),
    _f("I", data_type="float", title="Information",
       comment="Information content in bits, 1 + sum(freq * log freq) / log 20."),
    _f("I.norm", data_type="float", title="Normalised information",
       comment="Information content relative to the background distribution."),
    _f("height.I", data_type="float", title="Letter height",
       comment="Sequence-logo letter height, freq * I."),
    _f("height.I.norm", data_type="float", title="Normalised letter height",
       comment="Sequence-logo letter height against the background, freq * I.norm."),

    # -- the definitive tables (ROADMAP phase 6) -------------------------------------------------
    # ``records``, ``chains`` and ``evidence`` are the database; every shipped file is a
    # projection of them. Their columns are declared here with everything else so that one registry
    # describes every table, and so a rename fails a test instead of drifting.
    _f("record_id", searchable=0, autocomplete=0, title="Record id",
       comment="Stable VDJdb record identifier. Assigned once, never reused, and it survives a "
               "content change -- a curator fixing a typo amends a record rather than deleting one "
               "and creating another."),
    _f("clonotype_id", searchable=0, autocomplete=0, data_type="uint", title="Clonotype id",
       comment="Identifies a receptor chain: a hash of species, gene, CDR3, V and J. Records "
               "reporting the same chain share it, and motif evidence attaches at this level."),
    _f("d.segm", airr="d_call", title="D", comment="TCR Diversity segment allele."),
    _f("cdr3nt", type=SEQ, searchable=0, autocomplete=0, data_type="cdr3", airr="junction",
       title="CDR3 nucleotide",
       comment="Most plausible nucleotide junction behind the amino-acid one, inferred from the "
               "recombination model -- not observed. Two models agree on only 7.2 % of these, so it "
               "is a representative history, never evidence."),
    _f("cdr3nt.pgen", searchable=0, autocomplete=0, data_type="float", title="CDR3nt Pgen",
       comment="Generation probability of the inferred nucleotide junction."),
    _f("d.inferred", searchable=0, title="D inferred", airr="d_call",
       comment="D allele of the recombination scenario that produced cdr3nt. Not the curated "
               "d.segm, which it matches at gene level on 76.5 % of beta chains."),
    _f("d.start", searchable=0, autocomplete=0, data_type="uint", title="D start",
       comment="First nucleotide of the D segment in cdr3nt, 0-based, half-open with d.end."),
    _f("d.end", searchable=0, autocomplete=0, data_type="uint", title="D end",
       comment="One past the last nucleotide of the D segment in cdr3nt."),
    _f("v.inferred", searchable=0, title="V inferred",
       comment="V call proposed by the recombination model, filled only where the curator named "
               "none. Recovers the curated V on 23.8 % of human TRB and 50.1 % of TRA when it is "
               "hidden -- the junction carries little V. Never overwrites v.segm."),
    _f("j.inferred", searchable=0, title="J inferred",
       comment="J call proposed by the recombination model, filled only where the curator named "
               "none. Recovers the curated J on 97.5 % of human TRB and 95.8 % of TRA."),
    _f("d.posterior", searchable=0, autocomplete=0, data_type="float", title="D posterior",
       comment="Posterior probability of the gene d.inferred names, from arda.dpost. Median 0.791 "
               "and below 0.6 on 21.8 % of beta chains -- filter on it."),
    _f("d.entropy", searchable=0, autocomplete=0, data_type="float", title="D entropy",
       comment="Entropy of the posterior over D genes. Above 0.9 on 28.9 % of beta chains, where "
               "TRBD1 and TRBD2 are essentially undecidable from the junction."),
    _f("cdr3nt.margin", searchable=0, autocomplete=0, data_type="float", title="CDR3nt margin",
       comment="How far the inferred junction beat the runner-up: its Pgen divided by the next "
               "candidate's. Near 1 means the choice among synonymous histories was near-arbitrary."),
    _f("cdr3.original", type=SEQ, autocomplete=0, data_type="cdr3", title="CDR3 as submitted",
       comment="The CDR3 as the reference publication reported it, before repair."),
    _f("fix.needed", searchable=0, autocomplete=0, data_type="bool", title="Fix needed",
       comment="Whether the repaired CDR3 differs from the submitted one."),
    _f("fix.good", searchable=0, autocomplete=0, data_type="bool", title="Fix good",
       comment="Whether the CDR3 could be placed on both germline segments."),
    _f("v.fix.type", searchable=0, title="V fix type",
       comment="How the V side was repaired: NoFixNeeded, FixAdd, FixTrim, FixReplace, or a "
               "Failed* reason."),
    _f("j.fix.type", searchable=0, title="J fix type", comment="How the J side was repaired."),
    _f("v.canonical", searchable=0, autocomplete=0, data_type="bool", title="V anchor canonical",
       comment="Whether the CDR3 begins with the Cys104 the V germline predicts."),
    _f("j.canonical", searchable=0, autocomplete=0, data_type="bool", title="J anchor canonical",
       comment="Whether the CDR3 ends with the Phe/Trp118 the J germline predicts."),
    _f("chunk.file", searchable=0, title="Chunk file",
       comment="The chunk the record was read from. One chunk is one publication."),
    _f("chunk.row", searchable=0, autocomplete=0, data_type="uint", title="Chunk row",
       comment="0-based row within the chunk; with chunk.file it points at the curated line."),
    # -- the epitope catalogue (ROADMAP phase 9d) -----------------------------------------------
    _f("epitope.length", searchable=0, autocomplete=0, data_type="uint", title="Epitope length",
       comment="Residues in the epitope. MHC-I presents 8-11, MHC-II 12-25, so it cross-checks "
               "mhc.class independently of the allele."),
    _f("records", searchable=0, autocomplete=0, data_type="uint", title="Records",
       comment="Curated records supporting this row."),
    _f("chains", searchable=0, autocomplete=0, data_type="uint", title="Chains",
       comment="TCR chains across those records."),
    _f("clonotypes", searchable=0, autocomplete=0, data_type="uint", title="Clonotypes",
       comment="Distinct receptor chains, by clonotype_id -- records minus the replication."),
    _f("references", searchable=0, autocomplete=0, data_type="uint", title="References",
       comment="Distinct publications reporting this row. Two or more is independent replication."),
    _f("mhc.a.status", searchable=0, title="MHC A status",
       comment="known / unknown / unchecked against IPD-IMGT/HLA by prefix. `unchecked` means there "
               "is no authority for the name -- a murine molecule, or B2M."),
    _f("mhc.b.status", searchable=0, title="MHC B status",
       comment="As mhc.a.status, for the second chain."),

    _f("evidence_id", searchable=0, autocomplete=0, title="Evidence id",
       comment="Identifies one piece of evidence within a record: a hash of its type, chain, "
               "source and value, so the same evidence keeps the same id across releases."),
    _f("evidence_type", title="Evidence type",
       comment="independent_study, motif_tcrnet, motif_tcremp, structure_native or "
               "structure_model."),
    _f("evidence_source", title="Evidence source",
       comment="Where the evidence came from: a release tag, PDB, or a model set identifier."),
    _f("evidence_value", searchable=0, title="Evidence value",
       comment="The evidence itself: the other reference ids, a cluster id, or a PDB id."),
    _f("evidence_score", searchable=0, autocomplete=0, data_type="float", title="Evidence score",
       comment="Its strength: distinct supporting references, cluster size, or model confidence."),
    _f("first_seen_release", title="First seen release",
       comment="The release in which this piece of evidence first appeared."),

    # -- curation provenance, kept (ROADMAP 9) ---------------------------------------------------
    _f("submitter", title="Submitter", comment="Kept for debugging; dropped by the legacy build."),
    _f("chunk.id", title="Chunk id", comment="Kept for debugging; dropped by the legacy build."),
    _f("comment", title="Comment", comment="Kept for debugging; dropped by the legacy build."),
])

# ---------------------------------------------------------------------------------------------
# Positional orders. Each of these is a contract with a consumer -- see CLAUDE.md, hard rule 1.
# ---------------------------------------------------------------------------------------------

#: Chunk input, in README order: complex (15) + method (5) + meta (11).
COMPLEX_COLUMNS: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "d.beta", "j.beta",
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species", "reference.id",
)
METHOD_COLUMNS: tuple[str, ...] = (
    "method.identification", "method.frequency", "method.singlecell",
    "method.sequencing", "method.verification",
)
META_COLUMNS: tuple[str, ...] = (
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.epitope.id", "meta.tissue",
    "meta.donor.MHC", "meta.donor.MHC.method", "meta.structure.id",
)
ALL_COLUMNS: tuple[str, ...] = COMPLEX_COLUMNS + METHOD_COLUMNS + META_COLUMNS

#: ``vdjdb.txt`` as shipped -- 22 columns. Verified positionally against the 2026-06-03 release.
VDJDB_COLUMNS: tuple[str, ...] = (
    "complex.id", "gene", "cdr3", "v.segm", "j.segm",
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species", "reference.id",
    "vdjdb.score", "TCR_hash", "method", "meta", "cdr3fix",
    "web.method", "web.method.seq", "web.cdr3fix.nc", "web.cdr3fix.unmp",
)

#: The five columns production ``vdjdb-web`` appends. Owned by the new format from phase 6.
EVIDENCE_COLUMNS: tuple[str, ...] = (
    "evidence.validation.same.study", "evidence.validation.independent",
    "evidence.structure.contacts", "evidence.structure.quality", "evidence.structure.native",
)

#: ``vdjdb.txt`` as production serves it -- 27 columns.
VDJDB_WEB_COLUMNS: tuple[str, ...] = VDJDB_COLUMNS + EVIDENCE_COLUMNS

#: ``vdjdb.slim.txt`` -- 17 columns. Note ``j.start`` precedes ``v.end``, which the shipped
#: ``vdjdb.slim.meta.txt`` gets backwards (and two positions too early).
SLIM_COLUMNS: tuple[str, ...] = (
    "gene", "cdr3", "species", "antigen.epitope", "antigen.gene", "antigen.species",
    "complex.id", "v.segm", "j.segm", "mhc.a", "mhc.b", "mhc.class",
    "reference.id", "vdjdb.score", "TCR_hash", "j.start", "v.end",
)

#: ``vdjdb_full.txt`` -- 35 columns: the 31 chunk columns plus four derived ones.
FULL_COLUMNS: tuple[str, ...] = (*ALL_COLUMNS, "cdr3fix.alpha", "cdr3fix.beta", "vdjdb.score", "TCR_hash")

# ---------------------------------------------------------------------------------------------
# The definitive tables. Tidy: one observational unit per table, one variable per column. These are
# the database -- see `vdjdb.assemble.tables` for what each unit is and why.
# ---------------------------------------------------------------------------------------------

#: The pMHC and its parent antigen -- what the receptor recognises.
RECORD_ANTIGEN: tuple[str, ...] = (
    "species", "mhc.a", "mhc.b", "mhc.class",
    "antigen.epitope", "antigen.gene", "antigen.species",
)

#: The donor and sample a record was observed in. Part of its identity: the same TCR against the
#: same epitope in two donors is two records, not one seen twice.
RECORD_SAMPLE: tuple[str, ...] = (
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.tissue",
)

#: Where the record was written down. Part of the record: one row, one paper, one report.
RECORD_CURATION: tuple[str, ...] = ("chunk.file", "chunk.row", "chunk.id", "submitter", "comment")

#: ``records`` -- one row per curated record, PK ``record_id``.
RECORD_COLUMNS: tuple[str, ...] = (
    "record_id", *RECORD_ANTIGEN, "reference.id", *RECORD_SAMPLE,
    "meta.epitope.id", "meta.donor.MHC", "meta.donor.MHC.method", "meta.structure.id",
    "meta.subset.frequency",
    *METHOD_COLUMNS, "method.pairing",
    "vdjdb.score",
    *RECORD_CURATION,
)

#: ``chains`` -- one row per TCR chain, PK ``(record_id, gene)``. ``cdr3fix`` is flattened: every
#: member of the legacy JSON blob is a column, because a blob is not a variable.
CHAIN_COLUMNS: tuple[str, ...] = (
    "record_id", "gene", "clonotype_id",
    "cdr3", "v.segm", "d.segm", "j.segm",
    "v.end", "j.start",
    "cdr3nt", "cdr3nt.pgen", "cdr3nt.margin",
    "v.inferred", "j.inferred",
    "d.inferred", "d.start", "d.end", "d.posterior", "d.entropy",
    "cdr3.original", "fix.needed", "fix.good",
    "v.fix.type", "j.fix.type", "v.canonical", "j.canonical",
    "TCR_hash",
)

#: ``epitopes`` -- one row per antigen, PK ``(antigen.epitope, antigen.species)``. A peptide is not
#: unique to one organism: 13 epitopes are reported under two species and none of them is an error.
EPITOPE_COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "antigen.gene", "epitope.length",
    "mhc.class", "records", "chains", "references", "clonotypes",
)

#: ``restriction`` -- one row per (antigen, presenting MHC), each allele checked against
#: IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>).
RESTRICTION_COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "mhc.a", "mhc.b", "mhc.class",
    "mhc.a.status", "mhc.b.status", "records", "references",
)

#: ``evidence`` -- one row per piece of evidence, PK ``(record_id, evidence_id)``. Long rather than
#: wide: a record may have any number of pieces of evidence of any number of kinds, and the wide
#: form would be mostly empty.
EVIDENCE_TABLE_COLUMNS: tuple[str, ...] = (
    "record_id", "gene", "evidence_id", "evidence_type",
    "evidence_source", "evidence_value", "evidence_score", "first_seen_release",
)

#: ``cluster_members.txt`` -- 19 columns, parsed positionally by ``Motifs.scala`` with a fixed
#: ``Array[ColumnType]`` and no header check. Order is a contract.
CLUSTER_MEMBERS_COLUMNS: tuple[str, ...] = (
    "species", "antigen.epitope", "antigen.gene", "antigen.species",
    "mhc.a", "mhc.b", "mhc.class", "gene", "cdr3aa", "x", "y", "cid", "csz",
    "v.segm", "j.segm", "v.end", "j.start", "v.segm.repr", "j.segm.repr",
)

#: ``motif_pwms.txt`` -- 27 columns, same positional contract.
MOTIF_PWMS_COLUMNS: tuple[str, ...] = (
    "species", "antigen.epitope", "gene", "aa", "pos", "len",
    "v.segm.repr", "j.segm.repr", "cid", "csz",
    "count", "count.bg", "total.bg", "count.bg.i", "total.bg.i", "need.impute",
    "freq", "freq.bg", "I", "I.norm", "height.I", "height.I.norm",
    "antigen.gene", "antigen.species", "mhc.a", "mhc.b", "mhc.class",
)

#: Present in submitted chunks and discarded by the legacy build with no error. Kept from phase 6.
KEPT_CURATION_COLUMNS: tuple[str, ...] = (
    "chunk.id", "submitter", "comment", "meta.subset.frequency", "method.pairing",
)

#: Per-chunk deduplication key. Not the score signature, which is a different 11-column key
#: sharing the name ``SIGNATURE_COLS`` in the legacy code -- two keys, one name.
CHUNK_DEDUP_KEY: tuple[str, ...] = (
    *COMPLEX_COLUMNS,
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.tissue",
)

#: VDJdb column -> AIRR field, projected from the registry so it cannot be restated anywhere.
AIRR_MAP: dict[str, str] = {}   # populated below, after FIELDS is complete

SPECIES: frozenset[str] = frozenset({
    "HomoSapiens", "MusMusculus", "RattusNorvegicus", "MacacaMulatta",
})

_META_HEADER = "name\ttype\tvisible\tsearchable\tautocomplete\tdata.type\ttitle\tcomment"

TABLES: dict[str, tuple[str, ...]] = {
    "vdjdb": VDJDB_COLUMNS,
    "vdjdb-web": VDJDB_WEB_COLUMNS,
    "slim": SLIM_COLUMNS,
    "full": FULL_COLUMNS,
    "cluster_members": CLUSTER_MEMBERS_COLUMNS,
    "motif_pwms": MOTIF_PWMS_COLUMNS,
    "records": RECORD_COLUMNS,
    "chains": CHAIN_COLUMNS,
    "evidence": EVIDENCE_TABLE_COLUMNS,
    "epitopes": EPITOPE_COLUMNS,
    "restriction": RESTRICTION_COLUMNS,
}


AIRR_MAP.update({f.name: f.airr for f in FIELDS.values() if f.airr})


def fields(table: str) -> tuple[Field, ...]:
    """The :class:`Field` records of ``table``, in positional order."""
    try:
        names = TABLES[table]
    except KeyError:
        raise KeyError(f"unknown table {table!r}; known: {', '.join(sorted(TABLES))}") from None
    return tuple(FIELDS[n] for n in names)


def header(table: str) -> str:
    """The tab-separated header line of ``table``."""
    return "\t".join(TABLES[table])


def render_meta(table: str = "vdjdb") -> str:
    """``vdjdb.meta.txt`` for ``table`` -- eight columns, one row per data column, trailing newline.

    ``header(table)`` is the ``name`` column of this output by construction, which is the invariant
    ``BuildDatabase.groovy:411`` had and the Python port lost.
    """
    return "\n".join([_META_HEADER, *(f.meta_row() for f in fields(table))]) + "\n"


def render_slim_meta(table: str = "slim") -> str:
    """``vdjdb.slim.meta.txt`` -- two columns only, which is all standalone clients read."""
    rows = "\n".join(f"{f.name}\t{f.type}" for f in fields(table))
    return f"name\ttype\n{rows}\n"


def schema_json(dtypes: dict[str, dict[str, str]] | None = None, *, indent: int = 2) -> str:
    """The full registry, machine-readably: every column, its attributes and where it appears.

    Ships in the new-format bundle as ``vdjdb.schema.json``. A consumer reading ``records.parquet``
    can answer "what is this column, and which other tables have it" without scraping the docs.

    ``dtypes`` maps table name -> column name -> the physical dtype the build wrote, so the declared
    schema and the shipped files cannot disagree: it is read off the frames, not asserted.
    """
    import json

    dtypes = dtypes or {}
    out = []
    for name in sorted({c for cols in TABLES.values() for c in cols}):
        f = FIELDS[name]
        where = {t: cols.index(name) for t, cols in TABLES.items() if name in cols}
        entry: dict[str, object] = {
            "name": f.name, "type": f.type, "visible": f.visible, "searchable": f.searchable,
            "autocomplete": f.autocomplete, "data_type": f.data_type, "title": f.title,
            "comment": f.comment, "airr": f.airr, "position": where,
        }
        seen = {dtypes[t][name] for t in where if name in dtypes.get(t, {})}
        if seen:
            entry["dtype"] = sorted(seen)[0] if len(seen) == 1 else sorted(seen)
        out.append(entry)
    return json.dumps({"tables": {t: list(c) for t, c in TABLES.items()}, "fields": out},
                      indent=indent) + "\n"
