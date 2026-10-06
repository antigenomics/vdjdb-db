"""Reported chain typing can be incomplete without inventing its partner."""

import polars as pl
import pytest

from vdjdb.assemble.epitopes import assert_mhc_class, assert_mhc_resolves


def record(a, b, cls):
    return pl.DataFrame({'mhc.a': [a], 'mhc.b': [b], 'mhc.class': [cls],
                         'chunk.file': ['reported.tsv']})


@pytest.mark.parametrize('a,b', [('', 'HLA-DPB1*03'), ('HLA-DQA1*05', '')])
def test_partial_class_II_checks_the_reported_chain(a, b):
    r = record(a, b, 'MHCII')
    assert_mhc_resolves(r)
    assert_mhc_class(r)
    assert r['mhc.a'][0] == a and r['mhc.b'][0] == b


@pytest.mark.parametrize('a,b,cls', [('', '', 'MHCII'), ('HLA-A*02', '', 'MHCI'),
                                    ('', 'HLA-DPB1*99:99', 'MHCII')])
def test_missing_complete_restrictions_and_invalid_named_alleles_still_fail(a, b, cls):
    with pytest.raises(ValueError):
        assert_mhc_resolves(record(a, b, cls))


def test_bovine_source_serotype_is_declared():
    r = record('MHC-A*10', 'B2M', 'MHCI')
    assert_mhc_resolves(r)
    assert_mhc_class(r)
