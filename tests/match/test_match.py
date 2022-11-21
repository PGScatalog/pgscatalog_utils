""" Test that match strategies return the expected match results"""
import os
from unittest.mock import patch

import polars as pl
import pytest

from pgscatalog_utils.match.label import label_matches
from pgscatalog_utils.match.match import get_all_matches
from pgscatalog_utils.match.match_variants import match_variants


def test_only_match_pass(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-d', 'test',
                       # '--min_overlap', '0.5',
                       '--only_match',
                       '--outdir', out_dir]
                       # '--keep_ambiguous', '--keep_multiallelic']

    with patch('sys.argv', args):
        with pytest.raises(SystemExit) as se:
            match_variants()
        assert se.value.code == 0

    assert os.path.exists(os.path.join(out_dir, "matches/test_match_0.ipc.zst"))


def test_match_pass(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-d', 'test',
                       '--min_overlap', '0.95',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with patch('sys.argv', args):
        match_variants()

    assert os.path.exists(os.path.join(out_dir, "test_summary.csv"))
    assert os.path.exists(os.path.join(out_dir, "test_log.csv.gz"))
    assert os.path.exists(os.path.join(out_dir, "test_ALL_additive_0.scorefile.gz"))


def test_match_fail(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-d', 'test',
                       '--min_overlap', '1',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with pytest.raises(Exception) as excinfo:
        with patch('sys.argv', args):
            match_variants()
            assert "No valid matches found" in str(excinfo.value)


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

    params = {'skip_flip': True, 'remove_ambiguous': False, 'keep_first_match': False, 'remove_multiallelic': False}
    # check unambiguous matches
    df: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                        .pipe(label_matches, params=params)
                        .filter(pl.col('ambiguous') == False)
                        .collect())
    assert set(df['ID'].to_list()).issubset({'3:3:T:G', '1:1:A:C'})
    assert set(df['match_type'].to_list()).issubset(['altref', 'refalt'])

    # when keeping ambiguous and flipping alleles
    flip_params = {'skip_flip': False, 'remove_ambiguous': False, 'keep_first_match': False, 'remove_multiallelic': False}
    flip: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                          .pipe(label_matches, params=flip_params)
                          .filter(pl.col('ambiguous') == True)
                          .collect())

    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'altref', 'refalt_flip'})


def test_no_oa_match(small_scorefile_no_oa, small_target):
    scorefile, target = _cast_cat(small_scorefile_no_oa, small_target)

    no_ambig = {'skip_flip': True, 'remove_ambiguous': False, 'keep_first_match': False, 'remove_multiallelic': False}
    df: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                        .pipe(label_matches, params=no_ambig)
                        .filter(pl.col('ambiguous') == False)
                        .collect())

    assert set(df['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])
    assert set(df['match_type'].to_list()).issubset(['no_oa_alt', 'no_oa_ref'])

    # check ambiguous matches
    ambig = {'skip_flip': False, 'remove_ambiguous': False, 'keep_first_match': False, 'remove_multiallelic': False}
    flip: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                          .pipe(label_matches, ambig)
                          .filter(pl.col('ambiguous') == True)
                          .collect())
    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'no_oa_alt', 'no_oa_ref_flip'})


def test_flip_match(small_flipped_scorefile, small_target):
    scorefile, target = _cast_cat(small_flipped_scorefile, small_target)
    params = {'skip_flip': True, 'remove_ambiguous': False, 'keep_first_match': False, 'remove_multiallelic': False}
    df: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                        .pipe(label_matches, params=params)
                        .collect())

    assert df['ambiguous'].to_list() == [True, False, True, False]
    assert df['match_type'].to_list() == ['refalt', 'refalt_flip', 'altref_flip', 'altref_flip']
    assert df['match_status'].to_list() == ['matched', 'excluded', 'not_best', 'excluded']  # flipped -> excluded

    no_flip_params = {'skip_flip': False, 'remove_ambiguous': False, 'keep_first_match': False,
                      'remove_multiallelic': False}
    flip: pl.DataFrame = (pl.concat(get_all_matches(scorefile, target))
                          .pipe(label_matches, params=no_flip_params)
                          .filter(pl.col('ambiguous') == False)
                          .collect())
    assert flip['match_type'].str.contains('flip').all()
    assert set(flip['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])



