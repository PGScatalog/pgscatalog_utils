from unittest.mock import patch

import pandas as pd
import pytest
import jq

from pgscatalog_utils.download.score import query_score
from pgscatalog_utils.scorefile.combine_scorefiles import combine_scorefiles


def test_combine_scorefiles(combined_scorefile, _n_variants):
    df = pd.read_table(combined_scorefile)
    cols = {'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
            'is_duplicated', 'accession'}
    assert set(df.columns).issubset(cols)
    assert df.shape[0] == _n_variants


def test_liftover(lifted_scorefiles):
    df = pd.read_table(lifted_scorefiles)
    assert df.shape[0] == 832  # approx size


def test_fail_combine(scorefiles, tmp_path_factory):
    # these genomes are in build GRCh37, so combining with -t GRCh38 will raise an exception
    with pytest.raises(Exception):
        out_path = tmp_path_factory.mktemp("scores") / "combined.txt"
        args: list[str] = ['combine_scorefiles', '-t', 'GRCh38', '-s'] + scorefiles + ['-o', str(out_path.resolve())]
        with patch('sys.argv', args):
            combine_scorefiles()


@pytest.fixture
def _n_variants(pgs_accessions):
    json = query_score(pgs_accessions)
    n: list[int] = jq.compile("[.results][][].variants_number").input(json).all()
    return sum(n)

