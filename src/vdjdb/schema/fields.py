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

Three defects in the shipped metadata are **fixed here, deliberately**, each with a declared rule in
``rules/expected_diffs.toml``:

1. no ``TCR_hash`` row, though ``vdjdb.txt`` has carried the column for years;
2. ``vdjdb.score`` and ``TCR_hash`` placed after ``cdr3fix`` in the metadata but before ``method`` in
   the data -- so the metadata does not describe the file it belongs to, in the release *and* in
   production;
3. the four ``web.*`` rows carry one value too many, landing ``0`` in ``data.type``, ``factor`` in
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
    _f("gene", title="Gene", comment="TCR chain: alpha or beta."),
    _f("cdr3", type=SEQ, autocomplete=0, data_type="cdr3", title="CDR3",
       comment="TCR complementarity determining region 3 (CDR3) amino acid sequence."),
    _f("v.segm", title="V", comment="TCR Variable segment allele."),
    _f("j.segm", title="J", comment="TCR Joining segment allele."),
    _f("species", title="Species", comment="TCR parent species."),
    _f("mhc.a", title="MHC A", comment="First MHC chain allele."),
    _f("mhc.b", title="MHC B",
       comment="Second MHC chain allele (defaults to Beta2Microglobulin for MHC class I)."),
    _f("mhc.class", title="MHC class", comment="MHC class (I or II)."),
    _f("antigen.epitope", type=SEQ, data_type="peptide", title="Epitope",
       comment="Amino acid sequence of the epitope."),
    _f("antigen.gene", title="Epitope gene", comment="Representative parent gene of the epitope."),
    _f("antigen.species", title="Epitope species",
       comment="Representative parent species of the epitope."),
    _f("reference.id", data_type="url", title="Reference",
       comment="Pubmed reference / URL / or submitter details in case unpublished."),
    _f("vdjdb.score", autocomplete=0, data_type="uint", title="Info",
       comment="VDJdb confidence score, the higher is the score the more confidence we have in the "
               "antigen specificity annotation of a given TCR clonotype/clone. Zero score indicates "
               "that there are insufficient method details to draw any conclusion."),
    _f("TCR_hash", searchable=0, autocomplete=0, data_type=TXT, title="TCR hash",
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
    _f("evidence.validation.same.study", title="Validation same study",
       comment="Antigen specificity validated within the same study."),
    _f("evidence.validation.independent", title="Validation independent",
       comment="Antigen specificity independently validated in another study."),
    _f("evidence.structure.native", title="Structure native",
       comment="Native (experimental) TCR-pMHC structure available."),
    _f("evidence.structure.contacts", title="Structure model with contacts",
       comment="Structural model with annotated TCR-pMHC contacts available."),
    _f("evidence.structure.quality", title="Structure good quality model",
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
    # silently mistyping a whole file.
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

#: Present in real chunks and silently discarded by the legacy build. Kept from phase 6 on.
KEPT_CURATION_COLUMNS: tuple[str, ...] = (
    "chunk.id", "submitter", "comment", "meta.subset.frequency", "method.pairing",
)

#: Per-chunk deduplication key. **Not** the score signature, which is a different 11-column key
#: sharing the name ``SIGNATURE_COLS`` in the legacy code -- two keys, one name.
CHUNK_DEDUP_KEY: tuple[str, ...] = (
    *COMPLEX_COLUMNS,
    "meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
    "meta.replica.id", "meta.clone.id", "meta.tissue",
)

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
}


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
