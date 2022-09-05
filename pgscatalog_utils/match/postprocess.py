import logging
from functools import reduce

import polars as pl

from pgscatalog_utils.match.preprocess import complement_valid_alleles

logger = logging.getLogger(__name__)


def postprocess_matches(df: pl.DataFrame, remove_ambiguous: bool, keep_first_match: bool) -> pl.DataFrame:
    """ Clean up match candidates ready for writing out, including:

    - Label ambiguous variants
    - Prune match candidates to select the best match for each variant in the scoring file
    - Optionally remove ambiguous variants
    """
    df = _label_biallelic_ambiguous(df).pipe(_prune_matches, keep_first_match)

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
        .otherwise(False)))


def _prune_matches(df: pl.DataFrame, keep_first_match: bool = True) -> pl.DataFrame:
    """ Select the best match candidate in the target for each variant in the scoring file

    - In a scoring file (accession), each variant ID with the same effect allele and weight *must be unique*
    - The variant matching process normally returns multiple match candidates for each variant ID, e.g.:
        refalt > altref > refalt_flip > altref_flip
    - When multiple match candidates for an ID exist, they must be prioritised and pruned to be unique
    - If it's impossible to prioritise match candidates (i.e. same strategy is used), drop all matches by default

    :param df: A dataframe containing multiple match candidates for each variant
    :param drop_duplicates: If it's impossible to make match candidates unique, drop all candidates?
    :return: A dataframe containing the best match candidate for each variant
    """
    logger.debug("First match pruning: prioritise by match types")
    prioritised = _prioritise_match_type(df)
    singletons: pl.DataFrame = _get_singleton_variants(prioritised)
    dups: pl.DataFrame = _get_duplicate_variants(prioritised)

    if dups:
        if keep_first_match:
            logger.debug("Final match pruning: keeping first match")
            distinct: pl.DataFrame = pl.concat([singletons, dups.unique(maintain_order=True)])
        else:
            logger.debug("Final match pruning: dropping remaining duplicate matches")
            distinct: pl.DataFrame = singletons
    else:
        logger.debug("Final match pruning unnecessary")
        distinct: pl.DataFrame = singletons

    assert all(distinct.groupby(['accession', 'ID']).count()['count'] == 1), "Duplicate effect weights for a variant"
    logger.debug("Match pruning complete")

    return distinct.with_column(pl.lit(True).alias('passes_pruning'))


def _get_singleton_variants(df: pl.DataFrame) -> pl.DataFrame:
    """ Return variants with only one row (match candidate) per variant ID """
    return (df.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele', 'other_allele'])
            .count()
            .filter(pl.col('count') == 1)[:, "accession":"other_allele"]
            .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele', 'other_allele'], how='left'))


def _get_duplicate_variants(df: pl.DataFrame) -> pl.DataFrame:
    """ Return variants with more than one row (match candidate) per variant ID """
    return (df.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele', 'other_allele'])
            .count()
            .filter(pl.col('count') > 1)[:, "accession":"other_allele"]
            .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele', 'other_allele'], how='left'))


def _prioritise_match_type(duplicates: pl.DataFrame) -> pl.DataFrame:
    # first element has the highest priority and last element has the lowest priority
    match_priority = ['refalt', 'altref', 'refalt_flip', 'altref_flip', 'no_oa_ref', 'no_oa_alt', 'no_oa_ref_flip',
                      'no_oa_alt_flip']
    return _get_best_match(duplicates, match_priority)


def _get_best_match(df: pl.DataFrame, match_priority: list[str]) -> pl.DataFrame:
    match: list[pl.DataFrame] = []
    for match_type in match_priority:
        logger.debug(f"Selecting matches with match type {match_type}")
        match.append(df.filter(pl.col("match_type") == match_type))
    logger.debug("Prioritising match types (refalt > altref > ...)")
    return reduce(lambda x, y: _join_best_match(x, y), match)


def _join_best_match(x: pl.DataFrame, y: pl.DataFrame) -> pl.DataFrame:
    # variants in dataframe x have a higher priority than dataframe y
    # when concatenating the two dataframes, use an anti join to first remove variants in y that are in x
    not_in: pl.DataFrame = y.join(x, how='anti', on=['accession', 'ID'])
    return pl.concat([x, not_in])
