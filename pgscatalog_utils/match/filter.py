import logging

import polars as pl

logger = logging.getLogger(__name__)


def filter_scores(scorefile: pl.LazyFrame, matches: pl.LazyFrame, min_overlap: float,
                  dataset: str) -> tuple[pl.LazyFrame, pl.LazyFrame]:
    """ Check overlap between filtered matches and scorefile, remove scores that don't match well and report stats """
    filtered_matches: pl.LazyFrame = _filter_matches(matches)
    match_log: pl.LazyFrame = (_join_filtered_matches(filtered_matches, scorefile, dataset)
                               .with_columns(pl.col('best_match').fill_null(False)))

    fail_rates: pl.DataFrame = _calculate_match_rate(match_log).collect()  # collect for iteration

    scores: list[pl.DataFrame] = []
    for accession, rate in zip(fail_rates['accession'].to_list(), fail_rates['fail_rate'].to_list()):
        if rate < (1 - min_overlap):
            df: pl.DataFrame = pl.DataFrame({'accession': [accession], 'score_pass': [True], 'match_rate': [1 - rate]})
            logger.debug(f"Score {accession} passes minimum matching threshold ({1 - rate:.2%}  variants match)")
            scores.append(df.with_column(pl.col('accession').cast(pl.Categorical)))
        else:
            df: pl.DataFrame = pl.DataFrame({'accession': [accession], 'score_pass': [False], 'match_rate': [1 - rate]})
            logger.error(f"Score {accession} fails minimum matching threshold ({1 - rate:.2%} variants match)")
            scores.append(df.with_column(pl.col('accession').cast(pl.Categorical)))

    score_summary: pl.LazyFrame = pl.concat(scores).lazy()
    filtered_scores: pl.DataFrame = (filtered_matches.join(score_summary, on='accession', how='left')
                                     .filter(pl.col('score_pass') == True))

    return filtered_scores, score_summary


def _calculate_match_rate(df: pl.DataFrame) -> pl.DataFrame:
    logger.debug("Calculating overlap between target genome and scoring file")
    return (df.groupby('accession')
            .agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
            .with_column((pl.col('no_match') / pl.col('count')).alias('fail_rate')))


def _filter_matches(df: pl.LazyFrame) -> pl.LazyFrame:
    logger.debug("Filtering variants with exclude flag")
    return df.filter((pl.col('best_match') == True) & (pl.col('exclude') == False))


def _join_filtered_matches(matches: pl.LazyFrame, scorefile: pl.LazyFrame, dataset: str) -> pl.LazyFrame:
    return (scorefile.join(matches, on=['row_nr', 'accession'], how='left')
            .with_column(pl.lit(dataset).alias('dataset'))
            .select(pl.exclude("^.*_right$")))


def _match_keys() -> list[str]:
    return ['chr_name', 'chr_position', 'effect_allele', 'other_allele',
            'accession', 'effect_type', 'effect_weight']
