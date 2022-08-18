from functools import reduce
import polars as pl
import logging

from pgscatalog_utils.match.preprocess import complement_valid_alleles

logger = logging.getLogger(__name__)


def postprocess_matches(df: pl.DataFrame, remove_ambiguous: bool) -> pl.DataFrame:
    df = _label_biallelic_ambiguous(df)
    if remove_ambiguous:
        logger.debug("Removing ambiguous matches")
        return df.filter(pl.col("ambiguous") == False)
    else:
        logger.debug("Keeping best possible match from ambiguous matches")
        ambiguous: pl.DataFrame = df.filter((pl.col("ambiguous") == True) & \
                                            (pl.col("match_type").str.contains('flip').is_not()))
        unambiguous: pl.DataFrame = df.filter(pl.col("ambiguous") == False)
        return pl.concat([ambiguous, unambiguous])


def _label_biallelic_ambiguous(df: pl.DataFrame) -> pl.DataFrame:
    logger.debug("Labelling ambiguous variants")
    df = df.with_columns([
        pl.col(["effect_allele", "other_allele", "REF", "ALT", "effect_allele_FLIP", "other_allele_FLIP"]).cast(str),
        pl.lit(True).alias("ambiguous")
    ]).pipe(complement_valid_alleles, ["REF"])

    return (df.with_column(
        pl.when(pl.col("REF_FLIP") == pl.col("ALT"))
        .then(pl.col("ambiguous"))
        .otherwise(False))).pipe(_get_distinct_weights)


def _get_distinct_weights(df: pl.DataFrame) -> pl.DataFrame:
    """ Select single matched variant in target for each variant in the scoring file (e.g. per accession) """
    count: pl.DataFrame = df.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele']).count()
    singletons: pl.DataFrame = (count.filter(pl.col('count') == 1)[:, "accession":"effect_allele"]
                                .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele'], how='left'))

    dups: pl.DataFrame = (count.filter(pl.col('count') > 1)[:, "accession":"effect_allele"]
                          .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele'], how='left'))

    if dups:
        distinct: pl.DataFrame = pl.concat([singletons, _prioritise_match_type(dups)])
    else:
        distinct: pl.DataFrame = singletons

    assert all(distinct.groupby(['accession', 'ID']).count()['count'] == 1), "Duplicate effect weights for a variant"

    return distinct


def _prioritise_match_type(duplicates: pl.DataFrame) -> pl.DataFrame:
    dup_oa: pl.DataFrame = duplicates.filter(pl.col("other_allele") != None)
    dup_no_oa: pl.DataFrame = duplicates.filter(pl.col("other_allele") == None)
    best_matches: list[pl.DataFrame] = []

    if dup_oa:
        match_priority: list[str] = ['refalt', 'altref', 'refalt_flip', 'altref_flip']
        logger.debug(f"Prioritising matches in order {match_priority}")
        best_matches.append(_get_best_match(dup_oa, match_priority))

    if dup_no_oa:
        match_priority: list[str] = ['no_oa_ref', 'no_oa_alt', 'no_oa_ref_flip', 'no_oa_alt_flip']
        logger.debug(f"Prioritising matches in order {match_priority}")
        best_matches.append(_get_best_match(dup_no_oa, match_priority))

    return pl.concat(best_matches)


def _get_best_match(df: pl.DataFrame, match_priority: list[str]) -> pl.DataFrame:
    match: list[pl.DataFrame] = []
    for match_type in match_priority:
        match.append(df.filter(pl.col("match_type") == match_type))
    logger.debug("Filtering best match types")
    return reduce(lambda x, y: _join_best_match(x, y), match)


def _join_best_match(x: pl.DataFrame, y: pl.DataFrame) -> pl.DataFrame:
    # variants in dataframe x have a higher priority than dataframe y
    # when concatenating the two dataframes, use an anti join to first remove variants in y that are in x
    not_in: pl.DataFrame = y.join(x, how='anti',
                                  on=['accession', 'chr_name', 'chr_position', 'effect_allele', 'other_allele'])
    return pl.concat([x, not_in])
