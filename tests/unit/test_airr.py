"""The AIRR export: schema conformance, the value mappings, and the two source paths agreeing."""
from __future__ import annotations

from pathlib import Path

import airr as A
import polars as pl
import pytest
from airr.schema import RearrangementSchema as RS

from vdjdb.emit import airr
from vdjdb.schema import AIRR_MAP

CHAINS = pl.DataFrame({
    "record_id": ["VDJDB0000000001", "VDJDB0000000001", "VDJDB0000000002"],
    "gene": ["TRA", "TRB", "TRB"],
    "cdr3": ["CAVSDLEPNSSASKIIF", "CASSIRSSYEQYF", "CASSQQQGGF"],
    "v.segm": ["TRAV12-2*01", "TRBV10-3*01", "TRBV19*01"],
    "d.segm": ["", "TRBD2*01", ""],
    "j.segm": ["TRAJ3*01", "TRBJ2-7*01", "TRBJ2-7*01"],
})

RECORDS = pl.DataFrame({
    "record_id": ["VDJDB0000000001", "VDJDB0000000002"],
    "species": ["HomoSapiens", "HomoSapiens"],
    "antigen.epitope": ["GILGFVFTL", "NLVPMVATV"],
    "antigen.gene": ["M", ""],
    "antigen.species": ["InfluenzaA", "CMV"],
    "mhc.a": ["HLA-A*02:01", "HLA-DRA*01:01"],
    "mhc.b": ["B2M", "HLA-DRB1*04:01"],
    "mhc.class": ["MHCI", "MHCII"],
    "method.identification": ["tetramer-sort", "CD137 expression"],
    "vdjdb.score": [2, 0],
    "reference.id": ["PMID:28629751", ""],
})


# -- schema conformance ------------------------------------------------------------------------

def test_the_declared_required_columns_are_the_schemas_required_set():
    """A schema bump that adds a required field must fail here, not at a consumer."""
    assert set(airr.REARRANGEMENT_COLUMNS[:14]) == set(RS.required)


def test_every_emitted_column_is_a_real_airr_field():
    assert [c for c in airr.REARRANGEMENT_COLUMNS if c not in RS.properties] == []


def test_the_declared_schema_version_matches_the_installed_package():
    import yaml
    spec_path = Path(A.__file__).parent / "specs" / "airr-schema.yaml"
    with spec_path.open() as fh:
        spec = yaml.safe_load(fh)
    assert str(spec["Info"]["version"]) == airr.AIRR_VERSION


def test_the_written_file_passes_the_airr_validator(tmp_path):
    """The acceptance criterion. Nucleotide fields are empty until phase 8, which the schema allows:
    it requires the column, not a value."""
    airr.write_all(airr.from_tables({"records": RECORDS, "chains": CHAINS}), tmp_path)
    assert A.validate_rearrangement(tmp_path / "vdjdb.rearrangement.tsv")


# -- the mappings ------------------------------------------------------------------------------

def test_junction_aa_is_vdjdbs_cdr3_and_cdr3_aa_is_two_residues_shorter():
    """VDJdb's `cdr3` is the junction. Getting this backwards shifts every coordinate by two."""
    r = airr.rearrangement(CHAINS)
    assert r["junction_aa"].to_list() == CHAINS["cdr3"].to_list()
    assert (r["junction_aa"].str.len_chars() - r["cdr3_aa"].str.len_chars()).to_list() == [2, 2, 2]
    assert AIRR_MAP["cdr3"] == "junction_aa"


def test_a_rearrangement_resolves_back_to_its_record():
    r = airr.rearrangement(CHAINS)
    assert r["sequence_id"].to_list() == ["VDJDB0000000001:TRA", "VDJDB0000000001:TRB",
                                         "VDJDB0000000002:TRB"]
    assert r["cell_id"].to_list() == CHAINS["record_id"].to_list()
    assert r["locus"].to_list() == ["TRA", "TRB", "TRB"]


def test_mhc_class_gets_the_hyphen_airrs_vocabulary_requires():
    assert airr.reactivity(RECORDS)["mhc_class"].to_list() == ["MHC-I", "MHC-II"]


def test_antigen_is_never_empty_because_the_schema_forbids_it():
    """`antigen` is non-nullable. Where a chunk leaves the gene blank the epitope is the most
    specific true answer available."""
    assert airr.reactivity(RECORDS)["antigen"].to_list() == ["M", "NLVPMVATV"]


def test_the_score_is_reported_as_a_confidence_readout():
    """The spec asks a non-physical assay for a confidence level, which is what the score is."""
    rx = airr.reactivity(RECORDS)
    assert rx["reactivity_readout"].unique().to_list() == ["confidence"]
    assert rx["reactivity_value"].to_list() == [2, 0]
    assert rx["reactivity_unit"].unique().to_list() == ["vdjdb.score"]


@pytest.mark.parametrize("identification,expected", [
    ("tetramer-sort", "MHC_peptide_multimer"),
    ("dextramer-sort,cultured-T-cells", "MHC_peptide_multimer"),
    ("antigen-loaded-targets,dextramer-sort", "MHC_peptide_multimer"),
    ("antigen-expressing-targets", "native_protein"),
    ("T-Scan", "native_protein"),
    # the case that stops this reusing `emit.legacy._web_method`: a `sort` there, and 462 records
    ("CD137 expression", "annotated"),
    ("phage display, magnetic beads", "annotated"),
    ("cultured-T-cells", "annotated"),
])
def test_the_reactivity_method_classification(identification, expected):
    df = RECORDS.head(1).with_columns(pl.lit(identification).alias("method.identification"))
    assert airr.reactivity(df)["reactivity_method"][0] == expected


# -- the two source paths ----------------------------------------------------------------------

def _legacy_frame() -> pl.DataFrame:
    """A minimal legacy `vdjdb.txt`: one row per chain, record fields duplicated, `method` as JSON."""
    import json
    rows = []
    for i, (rid, gene, cdr3, v, j) in enumerate(zip(
            CHAINS["record_id"], CHAINS["gene"], CHAINS["cdr3"], CHAINS["v.segm"],
            CHAINS["j.segm"], strict=True)):
        rec = RECORDS.filter(pl.col("record_id") == rid).row(0, named=True)
        rows.append({
            "complex.id": "1" if rid == "VDJDB0000000001" else "0",
            "gene": gene, "cdr3": cdr3, "v.segm": v, "j.segm": j,
            "antigen.epitope": rec["antigen.epitope"], "antigen.gene": rec["antigen.gene"],
            "antigen.species": rec["antigen.species"], "mhc.a": rec["mhc.a"],
            "mhc.b": rec["mhc.b"], "mhc.class": rec["mhc.class"],
            "reference.id": rec["reference.id"], "vdjdb.score": str(rec["vdjdb.score"]),
            "method": json.dumps({"identification": rec["method.identification"]}),
            "__i": i,
        })
    return pl.DataFrame(rows).drop("__i")


def test_the_legacy_path_produces_nothing_the_tables_path_does_not():
    """The property that matters: legacy is a projection, so its AIRR is a subset of the tables'.

    `d_call` is excluded because legacy `vdjdb.txt` has no D column at all -- it is information the
    file cannot carry, not a disagreement. On the real corpus this holds exactly, with the tables
    carrying 1,501 extra chains: the 1,467 chains of 1,141 records the legacy build drops for a CDR3
    with no V or J, plus 34 D-only chains.
    """
    keys = ["locus", "v_call", "j_call", "junction_aa", "cdr3_aa"]
    a = airr.from_tables({"records": RECORDS, "chains": CHAINS})["rearrangement"]
    b = airr.from_legacy(_legacy_frame())["rearrangement"]
    assert b.select(keys).join(a.select(keys), on=keys, how="anti").is_empty()
    assert b.height == a.height


def test_the_legacy_path_collapses_the_duplicated_record_fields():
    """`vdjdb.txt` repeats a record's fields on each of its chains; Reactivity is per record."""
    b = airr.from_legacy(_legacy_frame())
    assert b["reactivity"].height == 2
    keys = [c for c in airr.REACTIVITY_COLUMNS if c not in ("reactivity_id", "cell_id")]
    a = airr.from_tables({"records": RECORDS, "chains": CHAINS})["reactivity"]
    assert (b["reactivity"].select(keys).with_columns(pl.col("reactivity_value").cast(pl.Int64))
            .join(a.select(keys), on=keys, how="anti").is_empty())


def test_an_empty_reference_id_is_not_read_back_as_a_null():
    """854 records ship with no `reference.id`; a CSV reader makes that a null, and then the two
    paths compare unequal on a record that is in fact identical."""
    b = airr.from_legacy(_legacy_frame().with_columns(pl.lit(None, pl.String).alias("reference.id")))
    assert b["reactivity"]["reactivity_refs"].to_list() == ["", ""]


# -- Receptor ----------------------------------------------------------------------------------

PAIRED = pl.DataFrame({
    "record_id": ["VDJDB0000000001", "VDJDB0000000001", "VDJDB0000000002"],
    "gene": ["TRB", "TRA", "TRB"],
    "cdr3": ["CASSEGWHSYEQYF", "CADLGSQGNLIF", "CASSIRSSYEQYF"],
    "v.segm": ["TRBV6-1*01", "TRAV21*01", "TRBV10-3*01"],
    "j.segm": ["TRBJ2-7*01", "TRAJ42*01", "TRBJ2-7*01"],
    "d.segm": ["", "", ""],
    "cdr3nt": ["TGTGCCAGCAGTGAAGGGTGGCACTCCTACGAGCAGTACTTC",
               "TGTGCAGACCTAGGAAGCCAAGGAAATCTCATCTTT",
               "TGTGCCAGTTCTATTAGGAGCTCCTACGAGCAGTACTTC"],
})
PAIRED_RECORDS = pl.DataFrame({"record_id": ["VDJDB0000000001", "VDJDB0000000002"],
                               "species": ["HomoSapiens", "HomoSapiens"]})


def test_a_receptor_needs_both_domains_so_unpaired_records_have_none():
    """Not a loss: a single chain is a Rearrangement, which is the file AIRR puts it in."""
    got = airr.receptor(PAIRED, PAIRED_RECORDS)
    assert got["receptor_id"].to_list() == ["VDJDB0000000001"]
    assert tuple(got.columns) == airr.RECEPTOR_COLUMNS


def test_domain_one_is_the_beta_chain_and_domain_two_the_alpha():
    """The schema's controlled vocabularies pin this; swapping them is silently wrong, not rejected."""
    got = airr.receptor(PAIRED, PAIRED_RECORDS).row(0, named=True)
    assert got["receptor_variable_domain_1_locus"] == "TRB"
    assert got["receptor_variable_domain_2_locus"] == "TRA"
    assert "CASSEGWHSYEQY" in got["receptor_variable_domain_1_aa"]
    assert "CADLGSQGNLI" in got["receptor_variable_domain_2_aa"]


def test_the_receptor_hash_is_airrs_and_not_vdjdbs_tcr_hash():
    import hashlib

    got = airr.receptor(PAIRED, PAIRED_RECORDS).row(0, named=True)
    joined = got["receptor_variable_domain_1_aa"] + got["receptor_variable_domain_2_aa"]
    assert got["receptor_hash"] == hashlib.sha256(joined.encode()).hexdigest()


def test_the_domain_is_the_mature_variable_region_not_just_the_junction():
    """AIRR asks for everything from after the signal peptide to the end of the J gene."""
    got = airr.receptor(PAIRED, PAIRED_RECORDS).row(0, named=True)
    d1 = got["receptor_variable_domain_1_aa"]
    assert len(d1) > 100, "a variable domain is ~110 aa; a junction is ~14"
    assert d1.endswith("GPGTRLTVT") or d1.endswith("GPGTRLTV")   # J framework 4


def test_the_legacy_path_cannot_produce_a_receptor():
    """It has no nucleotide junction to stitch around, and the domain columns are non-nullable."""
    assert "receptor" not in airr.from_legacy(_legacy_frame())
