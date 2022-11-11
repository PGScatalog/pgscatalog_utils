import gzip

import pytest
import os
import polars as pl

from unittest.mock import patch

from pgscatalog_utils.match.combine_matches import combine_matches
from pgscatalog_utils.match.match_variants import match_variants


def test_combine_matches_pass(mini_scorefile, only_matches, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['combine_matches', '-s', mini_scorefile,
                       '-m', only_matches,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--min_overlap', '0.9',
                       '--ignore_strand_flips',
                       '--keep_first_match',
                       '--keep_multiallelic']

    with patch('sys.argv', args):
        combine_matches()

    # and double check the output format of scorefiles
    with gzip.open(os.path.join(tmp_path, 'test_ALL_additive_0.scorefile.gz')) as f:
        scores = pl.read_csv(f, sep='\t')
    # pl.Categorical vs pl.Utf8 doesn't matter for this test
    assert scores.schema == {'ID': pl.Utf8, 'effect_allele': pl.Utf8, 'PGS001229_22': pl.Float64}



def test_combine_matches_fail(mini_scorefile, only_matches, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['combine_matches', '-s', mini_scorefile,
                       '-m', only_matches,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--min_overlap', '1.0',
                       '--ignore_strand_flips',
                       '--keep_first_match',
                       '--keep_multiallelic']

    with pytest.raises(Exception) as excinfo:
        with patch('sys.argv', args):
            combine_matches()

    assert "No valid matches found" in str(excinfo.value)


@pytest.fixture
def only_matches(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--only_match']

    with pytest.raises(SystemExit, match='0'):
        with patch('sys.argv', args):
            match_variants()

    return os.path.join(out_dir, 'matches', 'test_match_0.ipc.zst')

