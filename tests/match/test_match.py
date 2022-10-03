""" Test that match strategies return the expected match results"""

from unittest.mock import patch

import polars as pl
import pytest

from pgscatalog_utils.match.match import get_all_matches
from pgscatalog_utils.match.match_variants import match_variants


def test_match_fail(combined_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', combined_scorefile,
                       '-t', target_path,
                       '-m', 1,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with pytest.raises(Exception):
        with patch('sys.argv', args):
            match_variants()


def test_match_pass(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-m', 0,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with patch('sys.argv', args):
        match_variants()


def _cast_cat(scorefile, target) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    with pl.StringCache():
        scorefile = scorefile.with_columns([
            pl.col("chr_name").cast(pl.Utf8).cast(pl.Categorical),
            pl.col("effect_allele").cast(pl.Categorical),
            pl.col("other_allele").cast(pl.Categorical),
            pl.col("effect_type").cast(pl.Categorical),
            pl.col("effect_allele_FLIP").cast(pl.Categorical),
            pl.col("other_allele_FLIP").cast(pl.Categorical),
            pl.col("accession").cast(pl.Categorical)
        ])
        target = target.with_columns([
            pl.col("#CHROM").cast(pl.Utf8).cast(pl.Categorical),
            pl.col("REF").cast(pl.Categorical),
            pl.col("ALT").cast(pl.Categorical)
        ])
        return scorefile.lazy(), target.lazy()


def test_match_strategies(small_scorefile, small_target):
    scorefile, target = _cast_cat(small_scorefile, small_target)

    # check unambiguous matches
    df = (get_all_matches(scorefile, target, skip_flip=True, remove_ambiguous=False, keep_first_match=False)
          .filter(pl.col('ambiguous') == False)).collect()
    assert set(df['ID'].to_list()).issubset({'3:3:T:G', '1:1:A:C'})
    assert set(df['match_type'].to_list()).issubset(['altref', 'refalt'])

    # when keeping ambiguous and flipping alleles
    flip = (get_all_matches(scorefile, target, skip_flip=False, remove_ambiguous=False, keep_first_match=False)
            .filter(pl.col('ambiguous') == True)).collect()

    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'altref', 'refalt_flip'})


def test_no_oa_match(small_scorefile_no_oa, small_target):
    scorefile, target = _cast_cat(small_scorefile_no_oa, small_target)

    df = (get_all_matches(scorefile, target, skip_flip=True, remove_ambiguous=False, keep_first_match=False)
          .filter(pl.col('ambiguous') == False)).collect()

    assert set(df['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])
    assert set(df['match_type'].to_list()).issubset(['no_oa_alt', 'no_oa_ref'])

    # check ambiguous matches
    flip = (get_all_matches(scorefile, target, skip_flip=False, remove_ambiguous=False, keep_first_match=False)
            .filter(pl.col('ambiguous') == True)).collect()
    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'no_oa_alt', 'no_oa_ref_flip'})


def test_flip_match(small_flipped_scorefile, small_target):
    scorefile, target = _cast_cat(small_flipped_scorefile, small_target)

    df = get_all_matches(scorefile, target, skip_flip=True, remove_ambiguous=False, keep_first_match=False).collect()
    assert set(df['ambiguous']) == {True}
    assert set(df['match_type']) == {'refalt'}

    flip = (get_all_matches(scorefile, target, skip_flip=False, remove_ambiguous=False, keep_first_match=False)
            .filter(pl.col('ambiguous') == False)).collect()
    assert flip['match_type'].str.contains('flip').all()
    assert set(flip['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])



