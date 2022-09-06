from unittest.mock import patch

import polars as pl
import pytest

from pgscatalog_utils.match.match import get_all_matches, _cast_categorical
from pgscatalog_utils.match.match_variants import match_variants
from pgscatalog_utils.match.preprocess import complement_valid_alleles


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


def _cast_cat(scorefile, target):
    with pl.StringCache():
        return _cast_categorical(scorefile, target)


def test_match_strategies(small_scorefile, small_target):
    scorefile, target = _cast_cat(small_scorefile, small_target)

    # check unambiguous matches
    df = get_all_matches(scorefile, target, skip_flip=True).filter(pl.col('ambiguous') == False)
    assert set(df['ID'].to_list()).issubset({'3:3:T:G', '1:1:A:C'})
    assert set(df['match_type'].to_list()).issubset(['altref', 'refalt'])

    # when keeping ambiguous and flipping alleles
    flip = (get_all_matches(scorefile, target, skip_flip=False).filter(pl.col('ambiguous') == True))

    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'altref', 'refalt_flip'})


def test_no_oa_match(small_scorefile_no_oa, small_target):
    scorefile, target = _cast_cat(small_scorefile_no_oa, small_target)

    df = get_all_matches(scorefile, target, skip_flip=True).filter(pl.col('ambiguous') == False)

    assert set(df['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])
    assert set(df['match_type'].to_list()).issubset(['no_oa_alt', 'no_oa_ref'])

    # check ambiguous matches
    flip = (get_all_matches(scorefile, target, skip_flip=False)
            .filter(pl.col('ambiguous') == True))
    assert set(flip['ID'].to_list()).issubset({'2:2:T:A'})
    assert set(flip['match_type'].to_list()).issubset({'no_oa_alt', 'no_oa_ref_flip'})


def test_flip_match(small_flipped_scorefile, small_target):
    scorefile, target = _cast_cat(small_flipped_scorefile, small_target)

    df = get_all_matches(scorefile, target, skip_flip=True)
    assert set(df['ambiguous']) == {True}
    assert set(df['match_type']) == {'refalt'}

    flip = get_all_matches(scorefile, target, skip_flip=False).filter(pl.col('ambiguous') == False)
    assert flip['match_type'].str.contains('flip').all()
    assert set(flip['ID'].to_list()).issubset(['3:3:T:G', '1:1:A:C'])


@pytest.fixture
def small_scorefile():
    df = pl.DataFrame({"accession": ["test", "test", "test"],
                       "row_nr": [1, 2, 3],
                       "chr_name": [1, 2, 3],
                       "chr_position": [1, 2, 3],
                       "effect_allele": ["A", "A", "G"],
                       "other_allele": ["C", "T", "T"],
                       "effect_weight": [1, 2, 3],
                       "effect_type": ["additive", "additive", "additive"]})

    return complement_valid_alleles(df, ["effect_allele", "other_allele"])


@pytest.fixture
def small_scorefile_no_oa(small_scorefile):
    return small_scorefile.with_column(pl.lit(None).alias('other_allele'))


@pytest.fixture
def small_flipped_scorefile(small_scorefile):
    # simulate a scorefile on the wrong strand
    return (complement_valid_alleles(small_scorefile, ['effect_allele', 'other_allele'])
            .drop(['effect_allele', 'other_allele'])
            .rename({'effect_allele_FLIP': 'effect_allele', 'other_allele_FLIP': 'other_allele'})
            .pipe(complement_valid_alleles, ['effect_allele', 'other_allele']))


@pytest.fixture
def small_target():
    return pl.DataFrame({"#CHROM": [1, 2, 3],
                         "POS": [1, 2, 3],
                         "REF": ["A", "T", "T"],
                         "ALT": ["C", "A", "G"],
                         "ID": ["1:1:A:C", "2:2:T:A", "3:3:T:G"],
                         "is_multiallelic": [False, False, False]})
