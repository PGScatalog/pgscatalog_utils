import logging

import polars as pl

from pgscatalog_utils.match.preprocess import complement_valid_alleles

logger = logging.getLogger(__name__)


def label_matches(df: pl.LazyFrame, params: dict[str, bool]) -> pl.LazyFrame:
    """ Label match candidates with additional metadata. Column definitions:

    - match_candidate: All input variants that were returned from match.get_all_matches() (always True in this function)
    - best_match: True if row is the best possible match type (refalt > altref > ...)
    - duplicate: True if more than one best match exists for the same accession and ID
    - ambiguous: True if ambiguous
    """
    assert set(params.keys()) == {'keep_first_match', 'remove_ambiguous', 'remove_multiallelic', 'skip_flip'}
    labelled = (df.with_column(pl.lit(False).alias('exclude'))  # set up dummy exclude column for _label_*
                .pipe(_label_best_match)
                .pipe(_label_duplicate_best_match)
                .pipe(_label_duplicate_id, params['keep_first_match'])
                .pipe(_label_biallelic_ambiguous, params['remove_ambiguous'])
                .pipe(_label_multiallelic, params['remove_multiallelic'])
                .pipe(_label_flips, params['skip_flip'])
                .with_column(pl.lit(True).alias('match_candidate')))

    return _encode_match_priority(labelled)


def _encode_match_priority(df: pl.LazyFrame) -> pl.LazyFrame:
    """ Encode a new column called match status containing matched, unmatched, excluded, and not_best """
    return (df.with_columns([
        # set false best match to not_best
        pl.col('best_match').apply(lambda x: {None: 0, True: 1, False: 3}[x]).alias('match_priority'),
        pl.col('exclude').apply(lambda x: {None: 0, True: 2, False: 0}[x]).alias('excluded_match_priority')
    ])
            .with_column(pl.max(["match_priority", "excluded_match_priority"]))
            .with_column(pl.col("max")
                         .apply(lambda x: {0: 'unmatched', 1: 'matched', 2: 'excluded', 3: 'not_best'}[x])
                         .alias('match_status')
                         .cast(pl.Categorical)).drop(["max", "excluded_match_priority", "match_priority"]))


def _label_best_match(df: pl.LazyFrame) -> pl.LazyFrame:
    """ Best matches have the lowest match priority type. Find the best matches and label them.  """
    logger.debug("Labelling best match type (refalt > altref > ...)")
    match_priority = {'refalt': 0, 'altref': 1, 'refalt_flip': 2, 'altref_flip': 3, 'no_oa_ref': 4, 'no_oa_alt': 5,
                      'no_oa_ref_flip': 6, 'no_oa_alt_flip': 7}
    match_priority_rev = {v: k for k, v in match_priority.items()}

    # use a groupby aggregation to guarantee the number of rows stays the same
    # rows were being lost using an anti join + reduce approach
    prioritised: pl.LazyFrame = (df.with_column(pl.col('match_type')
                                                .apply(lambda x: match_priority[x])
                                                .alias('match_priority'))
                                 .with_column(pl.col("match_priority")
                                              .min()
                                              .over(["accession", "row_nr"])
                                              .apply(lambda x: match_priority_rev[x])
                                              .alias('best_match_type'))
                                 .with_column(pl.when(pl.col('best_match_type') == pl.col('match_type'))
                                              .then(pl.lit(True))
                                              .otherwise(pl.lit(False))
                                              .alias('best_match')))

    return prioritised.drop(['match_priority', 'best_match_type'])


def _label_duplicate_best_match(df: pl.LazyFrame) -> pl.LazyFrame:
    """ A scoring file row_nr in an accession group can be duplicated if a target position has different REF, e.g.:

    ┌────────┬────────────────────────┬────────────┬────────────────┬─────┬────────────┐
    │ row_nr ┆ accession              ┆ match_type ┆ ID             ┆ REF ┆ best_match │
    │ ---    ┆ ---                    ┆ ---        ┆ ---            ┆ --- ┆ ---        │
    │ i64    ┆ cat                    ┆ str        ┆ cat            ┆ str ┆ bool       │
    ╞════════╪════════════════════════╪════════════╪════════════════╪═════╪════════════╡
    │ 38557  ┆ PGS000012_hmPOS_GRCh37 ┆ no_oa_alt  ┆ 3:29588979:A:G ┆ A   ┆ true       │
    ├╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┤
    │ 38557  ┆ PGS000012_hmPOS_GRCh37 ┆ no_oa_alt  ┆ 3:29588979:T:G ┆ T   ┆ true       │
    └────────┴────────────────────────┴────────────┴────────────────┴─────┴────────────┘

    Label the first row with best_match = true, and duplicate rows with best_match = false
    """
    logger.debug("Labelling duplicated best match: keeping first instance as best_match = True")
    labelled: pl.LazyFrame = (df.with_column(pl.col('best_match')
                                             .count()
                                             .over(['accession', 'row_nr', 'best_match'])
                                             .alias('count'))
                              .with_column(pl.when((pl.col('count') > 1) & (pl.col('best_match') == True))
                                           .then(pl.lit(True))
                                           .otherwise(pl.lit(False))
                                           .alias('duplicate_best_match'))
                              .drop('count')
                              .rename({'row_nr': 'score_row_nr'})
                              .with_row_count()  # add temporary row count to get first variant
                              .with_column(pl.when((pl.col("best_match") == True) &
                                                   (pl.col("duplicate_best_match") == True) &
                                                   (pl.col("row_nr") > pl.min("row_nr")).over(
                                                       ["accession", "score_row_nr"]))
                                           .then(False)  # reset best match flag for duplicates
                                           .otherwise(pl.col("best_match"))  # just keep value from existing column
                                           .alias('best_match_duplicate_row_nr'))
                              .drop(['row_nr', 'best_match'])
                              .rename({'score_row_nr': 'row_nr', 'best_match_duplicate_row_nr': 'best_match'}))

    return labelled


def _label_duplicate_id(df: pl.LazyFrame, keep_first_match: bool) -> pl.LazyFrame:
    """ Label best match duplicates made when the scoring file is remapped to a different genome build

    ┌─────────┬────────────────────────┬─────────────┬────────────────┬─────┬────────────┐
    │ row_nr  ┆ accession              ┆ match_type  ┆ ID             ┆ REF ┆ best_match │
    │ ---     ┆ ---                    ┆ ---         ┆ ---            ┆ --- ┆ ---        │
    │ i64     ┆ cat                    ┆ str         ┆ cat            ┆ str ┆ bool       │
    ╞═════════╪════════════════════════╪═════════════╪════════════════╪═════╪════════════╡
    │ 1194115 ┆ PGS002244_hmPOS_GRCh37 ┆ altref      ┆ 3:50924580:C:A ┆ C   ┆ true       │
    ├╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌╌╌╌╌╌╌┤
    │ 1194132 ┆ PGS002244_hmPOS_GRCh37 ┆ refalt_flip ┆ 3:50924580:C:A ┆ C   ┆ true       │
    └─────────┴────────────────────────┴─────────────┴────────────────┴─────┴────────────┘

    refalt > altref > ... prioritisation doesn't fix this problem because row_nr is different (duplicated by remapping)
    """
    logger.debug('Labelling multiple scoring file lines (accession/row_nr) that best_match to the same variant')

    # the window in .over() starts with accession + ID
    # best_match is added to not count: same row_nr, different match_type (_label_best_match)
    # duplicate_best_match is added to not count: same row_nr, same match_type (_label_duplicate_best_match)
    duplicates = (df.with_column(pl.count("ID")
                                 .over(['accession', 'ID', 'best_match', 'duplicate_best_match'])
                                 .alias('count'))
                  .with_column(pl.when((pl.col('count') > 1) & (pl.col('best_match') == True))
                               .then(pl.lit(True))
                               .otherwise(pl.lit(False))
                               .alias('duplicate_ID')))

    if keep_first_match:
        logger.debug("Keeping first duplicate, labelling others with exclude flag ")
        # set first duplicate (with the smallest row_nr) to exclude = false
        labelled = duplicates.with_column(pl.when((pl.col("duplicate_ID") == True) &
                                                  (pl.col("row_nr") != pl.min("row_nr")
                                                   .over(["accession", "ID", "duplicate_ID"])))
                                          .then(True)
                                          .otherwise(False)
                                          .alias('exclude_duplicate'))
    else:
        logger.debug("Labelling all duplicates with exclude flag")
        labelled = duplicates.with_column(pl.col('duplicate_ID').alias('exclude_duplicate'))

    # get the horizontal maximum to combine the exclusion columns for each variant
    return (labelled.with_column(pl.max(["exclude", "exclude_duplicate"]))
            .drop(["exclude", "exclude_duplicate"])
            .rename({"max": "exclude"}))


def _label_biallelic_ambiguous(df: pl.LazyFrame, remove_ambiguous) -> pl.LazyFrame:
    logger.debug("Labelling ambiguous variants")
    ambig = ((df.with_columns([
        pl.col(["effect_allele", "other_allele", "REF", "ALT", "effect_allele_FLIP", "other_allele_FLIP"]).cast(str),
        pl.lit(True).alias("ambiguous")])
              .pipe(complement_valid_alleles, ["REF"]))
             .with_column(pl.when(pl.col("REF_FLIP") == pl.col("ALT"))
                          .then(pl.col("ambiguous"))
                          .otherwise(False)))

    if remove_ambiguous:
        logger.debug("Labelling ambiguous variants with exclude flag")
        return (ambig.with_column(pl.when(pl.col('ambiguous') == True)
                                  .then(True)
                                  .otherwise(False)
                                  .alias('exclude_ambiguous'))
                .with_column(pl.max(["exclude", "exclude_ambiguous"]))
                .drop(["exclude", "exclude_ambiguous"])
                .rename({"max": "exclude"}))
    else:
        return (ambig.with_column(pl.lit(False).alias('exclude_ambiguous'))
                .with_column(pl.max(["exclude", "exclude_ambiguous"]))
                .drop(["exclude", "exclude_ambiguous"])
                .rename({"max": "exclude"}))


def _label_multiallelic(df: pl.LazyFrame, remove_multiallelic: bool) -> pl.LazyFrame:
    """ Label multiallelic variants with exclude flag

    (Multiallelic variants are already labelled with the "is_multiallelic" column in match.preprocess)
    """
    if remove_multiallelic:
        logger.debug("Labelling multiallelic matches with exclude flag")
        return df.with_column(pl.when(pl.col('is_multiallelic') == True)
                              .then(True)
                              .otherwise(pl.col('exclude'))  # don't overwrite existing exclude flags
                              .alias('exclude'))
    else:
        logger.debug("Not excluding multiallelic variants")
        return df


def _label_flips(df: pl.LazyFrame, skip_flip: bool) -> pl.LazyFrame:
    df = df.with_column(pl.when(pl.col('match_type').str.contains('_flip'))
                        .then(True)
                        .otherwise(False)
                        .alias('match_flipped'))
    if skip_flip:
        logger.debug("Labelling flipped matches with exclude flag")
        return df.with_column(pl.when(pl.col('match_flipped') == True)
                              .then(True)
                              .otherwise(pl.col('exclude'))  # don't overwrite existing exclude flags
                              .alias('exclude'))
    else:
        logger.debug("Not excluding flipped matches")
        return df