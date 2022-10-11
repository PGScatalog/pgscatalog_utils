""" Test that matches have the correct labels, which is important for edge case handling and summary stats """

import logging
import pytest
import polars as pl

from pgscatalog_utils.match.label import label_matches
from pgscatalog_utils.match.match import get_all_matches
from tests.match.test_match import _cast_cat

logger = logging.getLogger(__name__)


def test_label_best_match(multiple_match_types):
    """ Test that multiple match candidates are correctly prioritised """
    # both matches are flagged as candidates
    assert multiple_match_types['match_candidate'].to_list() == [True, True]
    # but the matches have different match types
    assert multiple_match_types['match_type'].to_list() == ["altref", "refalt_flip"]
    # only one match candidate can survive!
    assert multiple_match_types['best_match'].to_list() == [True, False]
    assert multiple_match_types['match_status'].to_list() == ["matched", "not_best"]
    # however, exclude is _only_ for omitting a 'best match' from the final results, e.g. because of duplication
    assert multiple_match_types['exclude'].to_list() == [False, False]
    # match candidates are filtered by best_match == True and exclude == False


def test_label(small_scorefile, small_target):
    """ Test typical labels for match candidates with one match per position """
    scorefile, target = _cast_cat(small_scorefile, small_target)

    # get_all_matches calls label_matches
    params = {'skip_flip': True, 'remove_ambiguous': True, 'remove_multiallelic': False, 'keep_first_match': False}
    labelled: pl.DataFrame = (get_all_matches(scorefile=scorefile, target=target)
                              .pipe(label_matches, params=params)
                              .collect())

    logger.debug(labelled.select(['ID', 'match_type', 'best_match', 'ambiguous', 'match_status', 'exclude']))

    assert labelled['best_match'].to_list() == [True, True, True, False]
    assert labelled['ambiguous'].to_list() == [False, True, False, True]
    assert labelled['exclude'].to_list() == [False, True, False, True]
    assert labelled['match_status'].to_list() == ["matched", "excluded", "matched", "not_best"]


def test_ambiguous_label(small_flipped_scorefile, small_target):
    """ Test ambiguous variant labels change when they're kept for match candidates with one match per position """
    scorefile, target = _cast_cat(small_flipped_scorefile, small_target)
    no_flip = {'skip_flip': True, 'remove_ambiguous': True, 'remove_multiallelic': False, 'keep_first_match': False}
    no_ambiguous: pl.DataFrame = (get_all_matches(scorefile=scorefile, target=target)
                                  .pipe(label_matches, params=no_flip)
                                  .collect())

    # 2:2:T:A -> refalt      -> ambiguous     -> excluded (best match but ambiguous)
    # 1:1:A:C -> refalt_flip -> not ambiguous -> excluded (best match but skip_flip)
    # 2:2:T:A -> refalt_flip -> ambiguous     -> not_best (refalt priority so not best and excluded)
    # 3:3:T:G -> refalt_flip -> not ambiguous -> excluded (best match but skip_flip)
    assert no_ambiguous['best_match'].to_list() == [True, True, False, True]
    assert no_ambiguous['ambiguous'].to_list() == [True, False, True, False]
    assert no_ambiguous['exclude'].to_list() == [True, True, True, True]
    assert no_ambiguous['match_status'].to_list() == ["excluded", "excluded", "not_best", "excluded"]

    # otherwise, ambiguous variants are kept
    flip_params = {'skip_flip': True, 'remove_ambiguous': False, 'remove_multiallelic': False,
                   'keep_first_match': False}
    labelled = (get_all_matches(scorefile=scorefile, target=target)
                .pipe(label_matches, params=flip_params)
                .collect())

    # 2:2:T:A -> refalt      -> ambiguous     -> matched
    # 1:1:A:C -> refalt_flip -> not ambiguous -> excluded (best match but skip_flip)
    # 2:2:T:A -> refalt_flip -> ambiguous     -> not_best (refalt priority so not best and excluded)
    # 3:3:T:G -> refalt_flip -> not ambiguous -> excluded (best match but skip_flip)
    assert labelled['best_match'].to_list() == [True, True, False, True]
    assert labelled['ambiguous'].to_list() == [True, False, True, False]
    assert labelled['exclude'].to_list() ==  [False, True, True, True]
    assert labelled['match_status'].to_list() == ["matched", "excluded", "not_best", "excluded"]


def test_duplicate_ID(duplicated_matches, request):
    # these matches come from different lines in the original scoring file
    assert duplicated_matches["row_nr"].to_list() == [1, 4]
    # but they have the same ID!
    assert duplicated_matches["ID"].to_list() == ["1:1:A:C", "1:1:A:C"]
    # and they're matched with the same match type
    assert duplicated_matches["match_type"].to_list() == ["refalt", "refalt"]
    # oh dear, they're both the best match
    assert duplicated_matches["best_match"].to_list() == [True, True]
    # however, we've flagged them as duplicate IDs
    assert duplicated_matches['duplicate_ID'].to_list() == [True, True]

    if request.node.callspec.id == "keep_first_match":
        # and correctly label _the first occurring match_ as best match
        assert duplicated_matches['exclude'].to_list() == [False, True]
        assert duplicated_matches['match_status'].to_list() == ["matched", "excluded"]
    elif request.node.callspec.id == "delete_both":
        # and correctly labelled all duplicate instances for exclusion (default behaviour)
        assert duplicated_matches['exclude'].to_list() == [True, True]
        assert duplicated_matches['match_status'].to_list() == ["excluded", "excluded"]


def test_duplicate_best_match(duplicate_best_match):
    # all best matches come from the same row number in the original scoring file
    assert duplicate_best_match['row_nr'].to_list() == [1, 1, 1]
    # and the match type is duplicated, so we can't prioritise
    assert duplicate_best_match['match_type'].to_list() == ['no_oa_alt', 'no_oa_alt', 'no_oa_ref_flip']
    # find the duplicate best matches (with the same match type)
    assert duplicate_best_match['duplicate_best_match'].to_list() == [True, True, False]
    # and only keep the first occurring best match. the worse match type is correctly set to not_best too.
    assert duplicate_best_match['match_status'].to_list() == ["matched", "not_best", "not_best"]
    assert duplicate_best_match['best_match'].to_list() == [True, False, False]


@pytest.fixture(params=[True, False], ids=["keep_first_match", "delete_both"])
def duplicated_matches(small_scorefile, small_target, request) -> pl.DataFrame:
    # pgs catalog scorefiles can contain the same variant remapped to multiple rows
    # this happens after liftover to a different genome build
    # row_nrs will be different, but other information may be the same
    dups = (pl.concat([small_scorefile, small_scorefile])
            .with_column(pl.Series(list(range(1, 7)))
                         .alias('row_nr'))
            .filter(pl.col('chr_name') == 1))

    scorefile, target = _cast_cat(dups, small_target)

    params = {'skip_flip': False, 'remove_ambiguous': False, 'remove_multiallelic': False,
              'keep_first_match': request.param}
    return (get_all_matches(scorefile=scorefile, target=target)
            .pipe(label_matches, params=params)
            .collect())


@pytest.fixture
def multiple_match_types(small_target, small_scorefile) -> pl.DataFrame:
    # skip flip will return two candidate matches for one target position: refalt + refalt_flip
    scorefile, target = _cast_cat(small_scorefile, small_target)

    params = {'skip_flip': False, 'remove_ambiguous': False, 'remove_multiallelic': False, 'keep_first_match': False}
    return (get_all_matches(scorefile=scorefile, target=target)
            .pipe(label_matches, params=params)
            .filter(pl.col('chr_name') == '2')
            .collect())


@pytest.fixture
def duplicate_best_match(small_target, small_scorefile_no_oa) -> pl.DataFrame:
    # this type of target genome can sometimes occur when the REF is different at the same position
    odd_target = {'#CHROM': [1, 1], 'POS': [1, 1], 'REF': ['T', 'C'], 'ALT': ['A', 'A'], 'ID': ['1:1:T:C', '1:1:A:A'],
                  'is_multiallelic': [False, False]}
    scorefile, target = _cast_cat(small_scorefile_no_oa, pl.DataFrame(odd_target))

    params = {'skip_flip': False, 'remove_ambiguous': False, 'remove_multiallelic': False, 'keep_first_match': False}
    return (get_all_matches(scorefile=scorefile, target=target)
            .pipe(label_matches, params=params)
            .collect())
