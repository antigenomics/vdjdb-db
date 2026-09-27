"""AIRR export: Rearrangement and Reactivity, projected from the definitive tables.

Two files, at the two levels AIRR models:

===============================  ==============================================================
``vdjdb.rearrangement.tsv``      one row per chain -- AIRR Rearrangement
``vdjdb.reactivity.tsv``         one row per record -- AIRR Reactivity
``airr.yaml``                    the AIRR schema version the files conform to
===============================  ==============================================================

They are linked by ``cell_id``, which is the ``record_id``: a VDJdb record is one publication's
report on one T-cell clone, and a clone is the unit AIRR's ``Cell`` names. Both fields are standard
AIRR, so nothing here invents a join key.

Both source shapes use the same column names. The tidy ``chains`` table uses ``cdr3``, ``v.segm``,
``j.segm`` -- exactly what legacy ``vdjdb.txt`` calls them -- so :func:`rearrangement` and
:func:`reactivity` take one frame in VDJdb vocabulary and there is no second implementation for the
legacy path to drift from. :func:`from_tables` and :func:`from_legacy` differ only in how they
assemble that input.

``Receptor`` joins them once the nucleotide junction exists (phase 8): its two domain columns are
the complete mature variable domain and non-nullable, so they are rebuilt by stitching germline V
and J around the junction (:mod:`vdjdb.annotate.contig`). A receptor is a two-domain object by
definition, so only paired records appear there; a single chain is a Rearrangement, which is already
shipped.

One thing is still deliberately absent:

* **nucleotide fields** (``sequence``, ``junction``, the alignments and cigars) are present as
  columns and empty, which ``airr.validate_rearrangement`` accepts: the schema requires the column,
  not a value. Phase 8 (#461) fills them.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..convert.coords import cdr3_aa

#: The AIRR schema version these files conform to. Asserted against the installed ``airr`` package
#: in the tests: the build itself must not need ``airr`` installed.
AIRR_VERSION = "2.0"

#: AIRR Rearrangement. The first fourteen are the schema's ``required`` set, in its order -- a test
#: asserts that against ``airr.RearrangementSchema.required`` so a schema bump is caught here rather
#: than by a consumer.
REARRANGEMENT_COLUMNS: tuple[str, ...] = (
    "sequence_id", "sequence", "rev_comp", "productive",
    "v_call", "d_call", "j_call",
    "sequence_alignment", "germline_alignment",
    "junction", "junction_aa",
    "v_cigar", "d_cigar", "j_cigar",
    "locus", "cdr3_aa", "cell_id",
)

#: AIRR Reactivity. ``reactivity_id`` through ``reactivity_unit`` are its ``required`` set.
REACTIVITY_COLUMNS: tuple[str, ...] = (
    "reactivity_id", "cell_id", "ligand_type", "antigen_type", "antigen",
    "antigen_source_species", "peptide_sequence_aa",
    "mhc_class", "mhc_allele_1", "mhc_allele_2",
    "reactivity_method", "reactivity_readout", "reactivity_value", "reactivity_unit",
    "reactivity_refs",
)

#: AIRR Receptor. ``receptor_id`` through ``receptor_variable_domain_2_locus`` are its ``required``
#: set, all non-nullable -- which is why unpaired records cannot appear.
RECEPTOR_COLUMNS: tuple[str, ...] = (
    "receptor_id", "receptor_hash", "receptor_type",
    "receptor_variable_domain_1_aa", "receptor_variable_domain_1_locus",
    "receptor_variable_domain_2_aa", "receptor_variable_domain_2_locus",
)

#: Domain 1 is the heavy/beta/delta chain, domain 2 the light/alpha/gamma one -- the schema's
#: controlled vocabularies say so, and swapping them would validate while being wrong.
_DOMAIN = {1: "TRB", 2: "TRA"}

#: VDJdb spells the MHC class without the hyphen AIRR's controlled vocabulary requires.
MHC_CLASS = {"MHCI": "MHC-I", "MHCII": "MHC-II"}

#: An assay naming any of these determined specificity with a peptide-MHC multimer. Taken from the
#: 61 distinct ``method.identification`` values in the corpus, not from a guess: ``tetramer-sort``
#: (94,613 records), ``dextramer-sort`` (34,818), and the pentamer / streptamer / monomer tail.
_MULTIMER = ("tetramer", "dextramer", "pentamer", "multimer", "streptamer", "monomer",
             "pelimer", "mhc-peptide-beads")

#: An assay presenting the antigen as protein on a target cell rather than as a peptide-MHC complex.
_NATIVE = ("expressing-targets", "loaded-targets", "loaded-targed", "t-scan")

_CSV = {"separator": "\t", "line_terminator": "\n", "include_header": True}


def _reactivity_method() -> pl.Expr:
    """``method.identification`` -> one of AIRR's recommended ``reactivity_method`` keywords.

    Not :func:`vdjdb.emit.legacy._web_method`, which answers a different question (a coarse filter
    class for the web front end, ``sort`` / ``culture`` / ``other``). A CD137-expression sort is a
    ``sort`` there and is not a multimer assay here, so reusing it would mislabel 462 records.

    Everything the corpus does not name as a multimer or a target assay becomes ``annotated``, which
    is what the spec asks for: *"delineated as `annotated` if annotated from an external
    source"*. VDJdb curates from publications, so ``annotated`` is the correct default, not a gap.
    """
    m = pl.col("method.identification").str.to_lowercase()
    multimer = pl.any_horizontal(*[m.str.contains(k, literal=True) for k in _MULTIMER])
    native = pl.any_horizontal(*[m.str.contains(k, literal=True) for k in _NATIVE])
    return (pl.when(multimer).then(pl.lit("MHC_peptide_multimer"))
            .when(native).then(pl.lit("native_protein"))
            .otherwise(pl.lit("annotated")))


def rearrangement(chains: pl.DataFrame) -> pl.DataFrame:
    """One AIRR Rearrangement row per chain.

    ``chains`` must have ``record_id``, ``gene``, ``cdr3``, ``v.segm``, ``j.segm`` and ``d.segm``
    under those names, which both the tidy table and legacy ``vdjdb.txt`` already do.
    """
    d = pl.col("d.segm") if "d.segm" in chains.columns else pl.lit("")
    return chains.select(
        # Unique per chain, and it resolves back to the record it came from.
        pl.concat_str("record_id", "gene", separator=":").alias("sequence_id"),
        pl.lit("").alias("sequence"),
        # Every VDJdb record is a productive, forward-strand rearrangement: it is an expressed
        # receptor with a curated specificity, which is what got it into the database.
        pl.lit("F").alias("rev_comp"),
        pl.lit("T").alias("productive"),
        pl.col("v.segm").alias("v_call"), d.alias("d_call"), pl.col("j.segm").alias("j_call"),
        pl.lit("").alias("sequence_alignment"), pl.lit("").alias("germline_alignment"),
        # The inferred nucleotide junction (#461) when this source has one: the tidy `chains`
        # table does, legacy `vdjdb.txt` never did.
        (pl.col("cdr3nt") if "cdr3nt" in chains.columns else pl.lit("")).alias("junction"),
        # VDJdb's `cdr3` is the junction: Cys104..Phe/Trp118 inclusive. The identity here and the
        # two-residue trim below are what `convert.coords` exists for.
        pl.col("cdr3").alias("junction_aa"),
        pl.lit("").alias("v_cigar"), pl.lit("").alias("d_cigar"), pl.lit("").alias("j_cigar"),
        pl.col("gene").alias("locus"),
        cdr3_aa("cdr3"),
        pl.col("record_id").alias("cell_id"),
    ).select(REARRANGEMENT_COLUMNS)


def reactivity(records: pl.DataFrame) -> pl.DataFrame:
    """One AIRR Reactivity row per record.

    ``reactivity_value`` / ``reactivity_unit`` hold ``vdjdb.score``, which is what the spec asks a
    non-physical assay for: *"For inferred and annotated methods this should indicate a
    confidence/quality level"*, recommended keyword ``confidence``. VDJdb's score is a confidence in
    the specificity annotation, 0 to 3.
    """
    return records.select(
        pl.col("record_id").alias("reactivity_id"),
        pl.col("record_id").alias("cell_id"),
        pl.lit("MHC:peptide").alias("ligand_type"),
        pl.lit("peptide").alias("antigen_type"),
        # Non-nullable. The parent gene is the antigen; where a chunk leaves it blank the epitope
        # itself is the most specific true answer, and never empty.
        pl.when(pl.col("antigen.gene") != "").then(pl.col("antigen.gene"))
        .otherwise(pl.col("antigen.epitope")).alias("antigen"),
        pl.col("antigen.species").alias("antigen_source_species"),
        pl.col("antigen.epitope").alias("peptide_sequence_aa"),
        pl.col("mhc.class").replace_strict(MHC_CLASS, default="").alias("mhc_class"),
        pl.col("mhc.a").alias("mhc_allele_1"), pl.col("mhc.b").alias("mhc_allele_2"),
        _reactivity_method().alias("reactivity_method"),
        pl.lit("confidence").alias("reactivity_readout"),
        pl.col("vdjdb.score").cast(pl.Int64).alias("reactivity_value"),
        pl.lit("vdjdb.score").alias("reactivity_unit"),
        # An AIRR array, flattened to the comma-separated form AIRR TSVs use for `*_call` lists. A
        # VDJdb record is one publication's report, so it is always one reference.
        pl.col("reference.id").alias("reactivity_refs"),
    ).select(REACTIVITY_COLUMNS)


def receptor(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """One AIRR Receptor row per paired record whose two variable domains can be rebuilt.

    A receptor is a two-domain object: both ``receptor_variable_domain_*_aa`` are required and
    non-nullable, so a record with one chain has no Receptor row. It is not dropped -- its chain is
    in the Rearrangement file, which is where AIRR puts a single rearranged sequence.

    ``receptor_hash`` is AIRR's: sha256 over the concatenated domain sequences. It is not VDJdb's
    ``TCR_hash``, which hashes CDR3s, segments, MHC and epitope and is what the structure store is
    keyed on. Two hashes, two purposes, both kept.
    """
    import hashlib

    from ..annotate.contig import variable_domains

    domains = variable_domains(chains, records).filter(pl.col("vdomain_aa") != "")
    wide = None
    for n, locus in _DOMAIN.items():
        side = domains.filter(pl.col("gene") == locus).select(
            "record_id", pl.col("vdomain_aa").alias(f"receptor_variable_domain_{n}_aa"))
        wide = side if wide is None else wide.join(side, on="record_id", how="inner")
    assert wide is not None

    return wide.sort("record_id").select(
        pl.col("record_id").alias("receptor_id"),
        pl.concat_str("receptor_variable_domain_1_aa", "receptor_variable_domain_2_aa")
          .map_elements(lambda s: hashlib.sha256(s.encode()).hexdigest(), return_dtype=pl.Utf8)
          .alias("receptor_hash"),
        pl.lit("TCR").alias("receptor_type"),
        "receptor_variable_domain_1_aa", pl.lit(_DOMAIN[1]).alias("receptor_variable_domain_1_locus"),
        "receptor_variable_domain_2_aa", pl.lit(_DOMAIN[2]).alias("receptor_variable_domain_2_locus"),
    ).select(RECEPTOR_COLUMNS)


def from_tables(tables: dict[str, pl.DataFrame]) -> dict[str, pl.DataFrame]:
    """The definitive tables -> the three AIRR frames."""
    return {"rearrangement": rearrangement(tables["chains"]),
            "receptor": receptor(tables["chains"], tables["records"]),
            "reactivity": reactivity(tables["records"])}


def from_legacy(vdjdb_txt: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """Legacy ``vdjdb.txt`` -> the two AIRR frames, for users holding an old release zip.

    ``vdjdb.txt`` is already one row per chain in VDJdb vocabulary, so the only work is minting a
    ``record_id`` it does not have, pulling ``method.identification`` out of the ``method`` blob, and
    collapsing the duplicated record fields back to one row per record.
    """
    # A CSV reader turns an empty field into a null; empty string is the only missing marker here
    # (CLAUDE.md rule 6). 854 records ship with no `reference.id` at all, and without this fill they
    # compare unequal to the same records in the tables.
    df = vdjdb_txt.fill_null("").with_columns(
        pl.col("method").str.json_decode(
            dtype=pl.Struct([pl.Field("identification", pl.Utf8)])
        ).struct.field("identification").fill_null("").alias("method.identification"),
        pl.lit("").alias("d.segm"),
    )
    # The legacy file has no record identity. `complex.id` pairs two chains of one clone and is 0
    # for every unpaired record, so it cannot stand in alone: pair on it where it is non-zero, and
    # fall back to the row's own position where it is not.
    df = df.with_row_index("__i").with_columns(
        pl.when(pl.col("complex.id").cast(pl.Int64) > 0)
        .then(pl.concat_str(pl.lit("c"), pl.col("complex.id")))
        .otherwise(pl.concat_str(pl.lit("r"), pl.col("__i")))
        .alias("record_id")
    )
    return {"rearrangement": rearrangement(df),
            # One row per record: the record fields are duplicated across a clone's two chains.
            "reactivity": reactivity(df.unique(subset="record_id", keep="first",
                                               maintain_order=True))}


def write_all(frames: dict[str, pl.DataFrame], out: Path) -> dict[str, Path]:
    """Both TSVs plus ``airr.yaml``."""
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    for name, frame in frames.items():
        path = out / f"vdjdb.{name}.tsv"
        frame.write_csv(path, **_CSV)
        written[path.name] = path
    y = out / "airr.yaml"
    y.write_text(f'Info:\n  title: VDJdb AIRR export\n  version: {AIRR_VERSION}\n')
    written[y.name] = y
    return written
