import polars as pl
import logging

logger = logging.getLogger(__name__)


def get_all_matches(scorefile: pl.DataFrame, target: pl.DataFrame) -> pl.DataFrame:
    scorefile_cat, target_cat = _cast_categorical(scorefile, target)
    scorefile_oa = scorefile_cat.filter(pl.col("other_allele") != None)
    scorefile_no_oa = scorefile_cat.filter(pl.col("other_allele") == None)

    matches: list[pl.DataFrame] = []

    if scorefile_oa:
        logger.debug("Getting matches for scores with effect allele and other allele")
        matches.append(_match_variants(scorefile_cat, target_cat, effect_allele='REF', other_allele='ALT',
                                       match_type="refalt"))
        matches.append(_match_variants(scorefile_cat, target_cat, effect_allele='ALT', other_allele='REF',
                                       match_type="altref"))
        matches.append(_match_variants(scorefile_cat, target_cat, effect_allele='REF_FLIP',
                                       other_allele='ALT_FLIP',
                                       match_type="refalt_flip"))
        matches.append(_match_variants(scorefile_cat, target_cat, effect_allele='ALT_FLIP',
                                       other_allele='REF_FLIP',
                                       match_type="altref_flip"))

    if scorefile_no_oa:
        logger.debug("Getting matches for scores with effect allele only")
        matches.append(_match_variants(scorefile_no_oa, target_cat, effect_allele='REF', other_allele=None,
                                       match_type="no_oa_ref"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, effect_allele='ALT', other_allele=None,
                                       match_type="no_oa_alt"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, effect_allele='REF_FLIP',
                                       other_allele=None, match_type="no_oa_ref_flip"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, effect_allele='ALT_FLIP',
                                       other_allele=None, match_type="no_oa_alt_flip"))

    return pl.concat(matches)


def check_match_rate(scorefile: pl.DataFrame, matches: pl.DataFrame, min_overlap: float) -> None:
    scorefile: pl.DataFrame = scorefile.with_columns([
        pl.col('effect_type').cast(pl.Categorical),
        pl.col('accession').cast(pl.Categorical)])  # same dtypes for join
    match_log: pl.DataFrame = _join_matches(matches, scorefile)
    fail_rates: pl.DataFrame = (match_log.groupby('accession')
                                .agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
                                .with_column((pl.col('no_match') / pl.col('count')).alias('fail_rate'))
                                )

    for accession, rate in zip(fail_rates['accession'].to_list(), fail_rates['fail_rate'].to_list()):
        if rate < (1 - min_overlap):
            logger.debug(f"Score {accession} passes minimum matching threshold ({1-rate:.2%}  variants match)")
        else:
            logger.error(f"Score {accession} fails minimum matching threshold ({1-rate:.2%} variants match)")
            raise Exception


def _match_keys():
    return ['chr_name', 'chr_position', 'effect_allele', 'other_allele',
            'accession', 'effect_type', 'effect_weight']


def _join_matches(matches, scorefile):
    return scorefile.join(matches, on=_match_keys(), how='left')


def _match_variants(scorefile: pl.DataFrame,
                    target: pl.DataFrame,
                    effect_allele: str,
                    other_allele: str | None,
                    match_type: str) -> pl.DataFrame:
    logger.debug(f"Matching strategy: {match_type}")
    return (scorefile.join(target,
                           left_on=_scorefile_keys(other_allele),
                           right_on=_target_keys(effect_allele, other_allele),
                           how='inner')).pipe(_post_match, effect_allele, other_allele, match_type)


def _post_match(df: pl.DataFrame,
                effect_allele: str,
                other_allele: str,
                match_type: str) -> pl.DataFrame:
    """ Annotate matches with parameters """
    return df.with_columns([pl.col("*"),
                            pl.col("effect_allele").alias(effect_allele),
                            pl.col("other_allele").alias(other_allele),
                            pl.lit(match_type).alias("match_type")
                            ])[_matched_colnames()]


def _cast_categorical(scorefile, target) -> tuple[pl.DataFrame, pl.DataFrame]:
    """ Casting important columns to categorical makes polars fast """
    if scorefile:
        scorefile = scorefile.with_columns([
            pl.col("effect_allele").cast(pl.Categorical),
            pl.col("other_allele").cast(pl.Categorical),
            pl.col("effect_type").cast(pl.Categorical),
            pl.col("accession").cast(pl.Categorical)
        ])
    if target:
        target = target.with_columns([
            pl.col("REF").cast(pl.Categorical),
            pl.col("ALT").cast(pl.Categorical),
            pl.col("ALT_FLIP").cast(pl.Categorical),
            pl.col("REF_FLIP").cast(pl.Categorical)
        ])

    return scorefile, target


def _scorefile_keys(other_allele: str) -> list[str]:
    if other_allele:
        return ['chr_name', 'chr_position', 'effect_allele', 'other_allele']
    else:
        return ['chr_name', 'chr_position', 'effect_allele']


def _target_keys(effect_allele: str, other_allele: str) -> list[str]:
    if other_allele:
        return ['#CHROM', 'POS', effect_allele, other_allele]
    else:
        return ['#CHROM', 'POS', effect_allele]


def _matched_colnames() -> list[str]:
    return ['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type', 'accession',
            'ID', 'REF', 'ALT', 'REF_FLIP', 'ALT_FLIP', 'match_type']
