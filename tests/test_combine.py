import pandas as pd
import pytest
import jq

from pgscatalog_utils.download.score import query_score


def test_combine_scorefiles(combined_scorefile, _n_variants):
    df = pd.read_table(combined_scorefile)
    cols = {'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type', 'accession'}
    assert set(df.columns).issubset(cols)
    assert df.shape[0] == _n_variants


def test_liftover(lifted_scorefiles):
    df = pd.read_table(lifted_scorefiles)
    assert df.shape[0] > 50000  # approx size


@pytest.fixture
def _n_variants(pgs_accessions):
    json = query_score(pgs_accessions)
    n: list[int] = jq.compile("[.results][][].variants_number").input(json).all()
    return sum(n)
