"""Record associations must survive papers reporting several unrelated antigens."""
import polars as pl
import pytest

from vdjdb.corpus.query import receptor_lift


def test_receptor_lift_keeps_epitope_chain_and_species_on_the_same_record():
    records = pl.DataFrame({
        'record_id': ['a', 'b', 'c', 'd'],
        'reference.id': ['PMID:1'] * 4,
        'species': ['HomoSapiens'] * 3 + ['MusMusculus'],
        'antigen.epitope': ['TARGET', 'OTHER', 'TARGET', 'TARGET'],
    })
    chains = pl.DataFrame({
        'record_id': ['a', 'a', 'b', 'c', 'd'],
        'gene': ['TRA', 'TRB', 'TRB', 'TRB', 'TRB'],
        'cdr3': ['CNDDDF', 'CIRSF', 'CNDDDF', 'CIRSF', 'CNDDDF'],
    })
    result = receptor_lift(records, chains, species='HomoSapiens', gene='TRB',
                           epitope='TARGET')
    irs = result.filter(pl.col('term') == 'k:IRS').row(0, named=True)
    ndd = result.filter(pl.col('term') == 'k:NDD').row(0, named=True)
    assert (irs['units'], irs['given_units'], irs['both'], irs['term_units']) == (3, 2, 2, 2)
    assert irs['lift'] == pytest.approx(1.5)
    assert ndd['both'] == 0
    assert ndd['term_units'] == 1
    assert receptor_lift(records, chains, species='HomoSapiens', gene='TRB',
                          epitope='ABSENT')['lift'].null_count() == result.height
    with pytest.raises(ValueError, match='TRA or TRB'):
        receptor_lift(records, chains, species='HomoSapiens', gene='', epitope='TARGET')
