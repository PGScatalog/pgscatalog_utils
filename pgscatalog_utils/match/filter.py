import logging

import polars as pl

from pgscatalog_utils.match.log import write_log

logger = logging.getLogger(__name__)


def filter_scores(scorefile: pl.DataFrame, matches: pl.DataFrame, remove_ambiguous: bool, keep_first_match: bool,
                  min_overlap: float, dataset: str) -> pl.DataFrame:
    """ Remove scores that don't match well """
    scorefile: pl.DataFrame = scorefile.with_columns([
        pl.col('effect_type').cast(pl.Categorical),
        pl.col('accession').cast(pl.Categorical)])  # same dtypes for join

    # matches may contain more than one row per variant in the scoring file
    # e.g., one ambiguous match and one clear match, or duplicates may be in the scoring file
    filtered_matches: pl.DataFrame = _filter_matches(matches, remove_ambiguous, keep_first_match)
    match_log: pl.DataFrame = _join_matches(filtered_matches, scorefile, dataset)
    match_log['best_match'] = match_log['best_match'].fill_null(False)

    fail_rates: pl.DataFrame = _calculate_match_rate(match_log)

    scores: list[pl.DataFrame] = []
    for accession, rate in zip(fail_rates['accession'].to_list(), fail_rates['fail_rate'].to_list()):
        if rate < (1 - min_overlap):
            df: pl.DataFrame = pl.DataFrame({'accession': [accession], 'score_pass': [True], 'match_rate': [1 - rate]})
            logger.debug(f"Score {accession} passes minimum matching threshold ({1 - rate:.2%}  variants match)")
            scores.append(df)
        else:
            df: pl.DataFrame = pl.DataFrame({'accession': [accession], 'score_pass': [False], 'match_rate': [1 - rate]})
            logger.error(f"Score {accession} fails minimum matching threshold ({1 - rate:.2%} variants match)")
            scores.append(df)

    (match_log.with_column(pl.col('accession').cast(str))
     .join(pl.concat(scores), on='accession', how='left')).pipe(write_log, dataset)  # write log to gzipped CSV

    return (filtered_matches.with_column(pl.col('accession').cast(str))
            .join(pl.concat(scores), on='accession', how='left'))


def _calculate_match_rate(df: pl.DataFrame) -> pl.DataFrame:
    logger.debug("Calculating overlap between target genome and scoring file")
    return (df.groupby('accession')
            .agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
            .with_column((pl.col('no_match') / pl.col('count')).alias('fail_rate')))


def _filter_matches(df: pl.DataFrame, remove_ambiguous: bool, keep_first_match: bool) -> pl.DataFrame:
    logger.debug("Final match candidate filtering")
    return (df.filter(pl.col('best_match') == True)
            .pipe(_handle_ambiguous, remove_ambiguous)
            .pipe(_handle_duplicates, keep_first_match))


def _handle_ambiguous(df: pl.DataFrame, remove_ambiguous: bool) -> pl.DataFrame:
    if remove_ambiguous:
        logger.debug("Filtering: Removing ambiguous matches")
        return df.filter(pl.col("ambiguous") == False)
    else:
        logger.debug("Filtering: Keeping best possible match from ambiguous matches")
        ambiguous: pl.DataFrame = df.filter((pl.col("ambiguous") == True) & \
                                            (pl.col("match_type").str.contains('flip').is_not()))
        unambiguous: pl.DataFrame = df.filter(pl.col("ambiguous") == False)
        return pl.concat([ambiguous, unambiguous])


def _handle_duplicates(df: pl.DataFrame, keep_first_match: bool) -> pl.DataFrame:
    singletons = df.filter(pl.col('duplicate') == False)
    if keep_first_match:
        logger.debug("Filtering: keeping first match")
        first = (df.filter(pl.col('duplicate') == True)
                 .groupby(["accession", "ID"])
                 .agg([pl.col("row_nr").first()])
                 .join(df, on=['accession', 'row_nr'], how='left'))
        return pl.concat([singletons, first.select(singletons.columns)])
    else:
        logger.debug("Filtering: dropping any duplicate matches")
        return singletons


def _join_matches(matches: pl.DataFrame, scorefile: pl.DataFrame, dataset: str) -> pl.DataFrame:
    return (scorefile.join(matches, on=['accession', 'row_nr'], how='left')
            .with_column(pl.lit(dataset).alias('dataset'))
            .select(pl.exclude("^.*_right$")))


def _match_keys() -> list[str]:
    return ['chr_name', 'chr_position', 'effect_allele', 'other_allele',
            'accession', 'effect_type', 'effect_weight']
