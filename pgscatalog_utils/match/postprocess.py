import logging
from functools import reduce

import polars as pl

from pgscatalog_utils.match.preprocess import complement_valid_alleles

logger = logging.getLogger(__name__)


def postprocess_matches(df: pl.DataFrame) -> pl.DataFrame:
    """ Label match candidates with additional metadata. Column definitions:

    - match_candidate: All input variants that were returned from match.get_all_matches() (always True in this function)
    - best_match: True if row is the best possible match type (refalt > altref > ...)
    - duplicate: True if >1 scoring file line matches to the same variant ID
    - ambiguous: True if ambiguous
    """
    return (df.with_column(pl.lit(True).alias('match_candidate'))
            .pipe(_label_biallelic_ambiguous)
            .pipe(_label_pruned_matches))


def _label_biallelic_ambiguous(df: pl.DataFrame) -> pl.DataFrame:
    logger.debug("Labelling ambiguous variants")
    df = df.with_columns([
        pl.col(["effect_allele", "other_allele", "REF", "ALT", "effect_allele_FLIP", "other_allele_FLIP"]).cast(str),
        pl.lit(True).alias("ambiguous")
    ]).pipe(complement_valid_alleles, ["REF"])

    return (df.with_column(
        pl.when(pl.col("REF_FLIP") == pl.col("ALT"))
        .then(pl.col("ambiguous"))
        .otherwise(False)))


def _label_pruned_matches(df: pl.DataFrame) -> pl.DataFrame:
    best_matches = (df.pipe(_label_best_match)
                    .pipe(_label_duplicates))

    # check that duplicates were correctly labelled
    u_counts = best_matches.filter(pl.col('duplicate') == False).groupby(['accession', 'ID']).count()
    assert (u_counts['count'] == 1).all(), \
        "Duplicate effect weights for a variant: {}".format(list(u_counts['accession'].unique()))

    labelled = (df.join(best_matches, how='left', on=['row_nr', 'accession', 'ID'])
                .select(pl.exclude("^.*_right$")))
    assert labelled.shape[0] == df.shape[0]  # don't want to lose any rows from the input df

    return labelled


def _label_duplicates(df: pl.DataFrame) -> pl.DataFrame:
    """ Label scorefile (accession) matches with only one ID match (singletons) vs. multiple (duplicates)"""
    logger.debug('Labelling multiple accession - ID rows as duplicates')

    join_cols = ['accession', 'ID']
    counted = df.groupby(join_cols).count()
    singletons = (counted.filter(pl.col('count') == 1)[:, join_cols]
                         .join(df, on=join_cols, how='left')
                         .with_column(pl.lit(False).alias('duplicate')))
    duplicates = (counted.filter(pl.col('count') > 1)[:, join_cols]
                         .join(df, on=join_cols, how='left')
                         .with_column(pl.lit(True).alias('duplicate')))

    return pl.concat([singletons, duplicates])


def _label_best_match(df: pl.DataFrame) -> pl.DataFrame:
    match_priority = ['refalt', 'altref', 'refalt_flip', 'altref_flip', 'no_oa_ref', 'no_oa_alt', 'no_oa_ref_flip',
                      'no_oa_alt_flip']
    match: list[pl.DataFrame] = []
    for match_type in match_priority:
        logger.debug(f"Selecting matches with match type {match_type}")
        match.append(df.filter(pl.col("match_type") == match_type))

    logger.debug("Labelling best match type (refalt > altref > ...)")
    best_match: pl.DataFrame = reduce(lambda x, y: _prioritise_best_match(x, y), match)
    return best_match.with_column(pl.lit(True).alias('best_match'))


def _prioritise_best_match(x: pl.DataFrame, y: pl.DataFrame) -> pl.DataFrame:
    # variants in dataframe x have a higher priority than dataframe y
    # when concatenating the two dataframes, use an anti join to first remove variants in y that are in x
    not_in: pl.DataFrame = y.join(x, how='anti', on=['accession', 'ID', 'row_nr'])
    return pl.concat([x, not_in])
