"""Built chains are complete, and incomplete input stops publication."""
import polars as pl
import pytest

from vdjdb.assemble.master import FIX_FIELDS, assert_complete_chains


def master():
    row = {'record_id': 'r1'}
    for gene in ('alpha', 'beta'):
        row.update({f'cdr3.{gene}': 'CASSF' if gene == 'beta' else '',
                    f'v.{gene}': 'TRBV1*01' if gene == 'beta' else '',
                    f'j.{gene}': 'TRBJ1-1*01' if gene == 'beta' else ''})
        for _, name, kind in FIX_FIELDS:
            row[f'{name}.{gene}'] = {str: '', bool: False, int: -1}[kind]
        for name in ('cdr3', 'v', 'j'):
            row[f'__{name}.{gene}'] = row[f'{name}.{gene}']
        row[f'__vfix.{gene}'] = row[f'__jfix.{gene}'] = 'NoFixNeeded'
    return pl.DataFrame([row])


def test_complete_chain_and_absent_partner_are_accepted():
    assert_complete_chains(master())


@pytest.mark.parametrize('column,value', [
    ('v.beta', ''), ('j.beta', ''), ('cdr3.beta', ''),
    ('__vend.beta', None), ('__cdr3.beta', ''), ('__v.beta', ''), ('__j.beta', ''),
    ('__jfix.beta', ''), ('__good.beta', None),
])
def test_partial_chains_stop_the_build(column, value):
    with pytest.raises(ValueError, match=r'r1.*beta'):
        assert_complete_chains(master().with_columns(pl.lit(value).alias(column)))


def test_zero_chain_records_stop_the_build():
    empty = master().with_columns(pl.lit('').alias(c) for c in ('cdr3.beta', 'v.beta', 'j.beta'))
    with pytest.raises(ValueError, match=r'r1.*zero chains'):
        assert_complete_chains(empty)
