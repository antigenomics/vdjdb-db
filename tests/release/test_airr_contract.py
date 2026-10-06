"""AIRR conformance and the legacy-vs-tables property, on a real build.

Marked ``release``: needs a built directory (``VDJDB_TABLES``, default ``out/tables``) and the legacy
projection beside it. The unit tests cover the mappings; these two properties only mean something at
corpus scale.
"""
from __future__ import annotations

import os
from pathlib import Path

import airr as A
import polars as pl
import pytest

from vdjdb.emit import airr
from vdjdb.emit.vdjdb3 import read_tables

pytestmark = pytest.mark.release

#: Chains the legacy build drops: the ones whose chain carried a CDR3 with no V or J (it drops the
#: record whole, both chains), plus 34 D-only chains with no CDR3 at all.
#:
#: **Both numbers falling is the outcome to want**, and they have fallen twice. The tables path always
#: carried these records; what changes is that legacy carries them too.
#:
#: ====================================  ======  =======
#: after                                 chains  records
#: ====================================  ======  =======
#: `res/` retired (#658)                 1,501    1,141
#: the k-mer J proposal replaced         1,032      845
#: `annotate_junctions` (arda 2.36)        970      755
#: ====================================  ======  =======
#:
#: The last step moved the blank-call proposal from this repository's own `annotate/segments.py` into
#: `arda.cdr3fix`, which proposes the **locus** as well as the call (2.36) and so answers the 461 keys
#: that named neither side and had no locus to look one up under. 62 more chains of 90 more records
#: stop failing the legacy "a CDR3 needs a V and a J" filter.
#: #845 restores nine additional single-chain observations with unresolved segment calls.
#: All earlier observations remain; these nine are retained by the tables and AIRR projection.
#: #893 adds 336 records (542 chains) whose source calls remain incomplete after assembly.
#: All are retained in the definitive tables and AIRR; legacy requires both V and J.
#: #934 adds 174 records (239 chains) with incomplete calls after assembly.
#: These observations remain in the definitive tables and AIRR.
#: PMID34290408 adds eight TRB observations with incomplete calls; tables and AIRR retain them.
#: PMID25681349 adds 19 incomplete-call chains retained in definitive tables and AIRR.
#: PMID36516854 adds four incomplete-call observations retained by tables and AIRR.
#: Final publication additions retain 22 further incomplete-V chains in tables/AIRR.
# Three new source observations retain six chains without resolvable segment calls.
# Primary-table segment reconciliation recovers ten chain rows and nine record rows
# for the legacy projection; the definitive tables retain every observation.
# PMID38100526 full-protein calls recover the two previously unpaired, blank-V observations
# through their reported pair. Two fewer chains and records fail the legacy segment filter.
# Finalization retains 98 further single-chain observations with unresolved segments.
LEGACY_DROPS_CHAINS = 1915
LEGACY_DROPS_RECORDS = 1422


@pytest.fixture(scope="module")
def built() -> tuple[dict[str, pl.DataFrame], pl.DataFrame]:
    d = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    legacy = d.parent / "legacy" / "vdjdb.txt"
    if not (d / "records.parquet").exists() or not legacy.exists():
        pytest.skip(f"no build at {d.parent}; run `vdjdb build --out out/`")
    return read_tables(d), pl.read_csv(legacy, separator="\t", infer_schema=False, quote_char=None)


def test_the_whole_rearrangement_table_passes_the_airr_validator(built, tmp_path):
    """Phase 7's acceptance criterion."""
    tables, _ = built
    airr.write_all(airr.from_tables(tables), tmp_path)
    assert A.validate_rearrangement(tmp_path / "vdjdb.rearrangement.tsv")


def test_the_legacy_path_is_a_strict_subset_of_the_tables_path(built):
    """`d_call` is excluded: legacy `vdjdb.txt` has no D column, so it is information the file cannot
    carry rather than a disagreement. 42,574 beta chains have a D call in the tables."""
    tables, legacy = built
    keys = ["locus", "v_call", "j_call", "junction_aa", "cdr3_aa"]
    a = airr.from_tables(tables)["rearrangement"].select(keys).group_by(keys).len()
    b = airr.from_legacy(legacy)["rearrangement"].select(keys).group_by(keys).len()
    j = a.join(b, on=keys, how="full", suffix="_legacy", coalesce=True).fill_null(0)

    assert j.filter(pl.col("len") == 0)["len_legacy"].sum() == 0, "legacy invented rows"
    assert j.filter(pl.col("len_legacy") > pl.col("len")).is_empty(), "legacy has more of a row"
    excess = (j["len"] - j["len_legacy"])
    assert excess.filter(excess > 0).sum() == LEGACY_DROPS_CHAINS


def test_reactivity_agrees_and_the_tables_keep_the_records_legacy_drops(built):
    tables, legacy = built
    keys = [c for c in airr.REACTIVITY_COLUMNS if c not in ("reactivity_id", "cell_id")]
    a = airr.from_tables(tables)["reactivity"].select(keys).group_by(keys).len()
    b = (airr.from_legacy(legacy)["reactivity"].select(keys)
         .with_columns(pl.col("reactivity_value").cast(pl.Int64)).group_by(keys).len())
    j = a.join(b, on=keys, how="full", suffix="_legacy", coalesce=True).fill_null(0)

    assert j.filter(pl.col("len") == 0)["len_legacy"].sum() == 0
    excess = (j["len"] - j["len_legacy"])
    assert excess.filter(excess > 0).sum() == LEGACY_DROPS_RECORDS


def test_receptors_are_paired_records_with_both_domains_rebuilt(built):
    """A receptor is a two-domain object, so the file is smaller than the record table by design:
    93,294 of 192,753 records are paired, and 81,003 of those have both variable domains rebuilt."""
    tables, _ = built
    rec = airr.from_tables(tables)["receptor"]
    assert rec.filter((pl.col("receptor_variable_domain_1_aa") == "")
                      | (pl.col("receptor_variable_domain_2_aa") == "")).is_empty()
    assert rec["receptor_id"].n_unique() == rec.height
    assert set(rec["receptor_variable_domain_1_locus"]) == {"TRB"}
    assert set(rec["receptor_variable_domain_2_locus"]) == {"TRA"}
    # a mature TCR variable domain is roughly 110 aa; a junction is ~14
    lens = rec["receptor_variable_domain_1_aa"].str.len_chars()
    assert lens.min() > 90 and lens.max() < 160

    paired = (tables["chains"].group_by("record_id")
              .agg(pl.col("gene").n_unique().alias("n")).filter(pl.col("n") == 2))
    assert rec.height <= paired.height


def test_every_chain_has_a_rearrangement_row(built):
    tables, _ = built
    assert airr.from_tables(tables)["rearrangement"].height == tables["chains"].height
