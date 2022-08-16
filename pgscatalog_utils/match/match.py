import polars as pl
import logging

from pgscatalog_utils.match.postprocess import postprocess_matches
from pgscatalog_utils.match.write import write_log

logger = logging.getLogger(__name__)


def get_all_matches(scorefile: pl.DataFrame, target: pl.DataFrame, remove_ambiguous: bool) -> pl.DataFrame:
    scorefile_cat, target_cat = _cast_categorical(scorefile, target)
    scorefile_oa = scorefile_cat.filter(pl.col("other_allele") != None)
    scorefile_no_oa = scorefile_cat.filter(pl.col("other_allele") == None)

    matches: list[pl.DataFrame] = []

    if scorefile_oa:
        logger.debug("Getting matches for scores with effect allele and other allele")
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="refalt"))
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="altref"))
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="refalt_flip"))
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="altref_flip"))

    if scorefile_no_oa:
        logger.debug("Getting matches for scores with effect allele only")
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_ref"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_alt"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_ref_flip"))
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_alt_flip"))

    return pl.concat(matches).pipe(postprocess_matches, remove_ambiguous)


def check_match_rate(scorefile: pl.DataFrame, matches: pl.DataFrame, min_overlap: float, dataset: str) -> None:
    scorefile: pl.DataFrame = scorefile.with_columns([
        pl.col('effect_type').cast(pl.Categorical),
        pl.col('accession').cast(pl.Categorical)])  # same dtypes for join
    match_log: pl.DataFrame = _join_matches(matches, scorefile, dataset)
    write_log(match_log, dataset)
    fail_rates: pl.DataFrame = (match_log.groupby('accession')
                                .agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
                                .with_column((pl.col('no_match') / pl.col('count')).alias('fail_rate'))
                                )

    for accession, rate in zip(fail_rates['accession'].to_list(), fail_rates['fail_rate'].to_list()):
        if rate < (1 - min_overlap):
            logger.debug(f"Score {accession} passes minimum matching threshold ({1 - rate:.2%}  variants match)")
        else:
            logger.error(f"Score {accession} fails minimum matching threshold ({1 - rate:.2%} variants match)")
            raise Exception


def _match_keys():
    return ['chr_name', 'chr_position', 'effect_allele', 'other_allele',
            'accession', 'effect_type', 'effect_weight']


def _join_matches(matches: pl.DataFrame, scorefile: pl.DataFrame, dataset: str):
    return scorefile.join(matches, on=_match_keys(), how='left').with_column(pl.lit(dataset).alias('dataset'))


def _match_variants(scorefile: pl.DataFrame, target: pl.DataFrame, match_type: str) -> pl.DataFrame:
    logger.debug(f"Matching strategy: {match_type}")
    match match_type:
        case 'refalt':
            score_keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
            target_keys = ["#CHROM", "POS", "REF", "ALT"]
        case 'altref':
            score_keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
        case 'refalt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF", "ALT"]
        case 'altref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
        case 'no_oa_ref':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "REF"]
        case 'no_oa_alt':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "ALT"]
        case 'no_oa_ref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF"]
        case 'no_oa_alt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT"]
        case _:
            logger.critical(f"Invalid match strategy: {match_type}")
            raise Exception

    return (scorefile.join(target, score_keys, target_keys, how='inner')
            .pipe(_post_match, target_keys, match_type))


def _post_match(df: pl.DataFrame,
                target_keys: list[str],
                match_type: str) -> pl.DataFrame:
    """ Annotate matches with parameters """
    if len(target_keys) == 3:
        logger.debug("Dropping missing other_allele during annotation")
        ref_key = target_keys[-1]
        alt_key = 'dummy'  # prevent trying to alias a column to None
    else:
        ref_key = target_keys[-2]
        alt_key = target_keys[-1]

    # aliases keep a copy of columns dropped during the join
    return df.with_columns([pl.col("*"),
                            pl.col("effect_allele").alias(ref_key),
                            pl.col("other_allele").alias(alt_key),
                            pl.lit(match_type).alias("match_type"),
                            ])[_matched_colnames()]


def _cast_categorical(scorefile, target) -> tuple[pl.DataFrame, pl.DataFrame]:
    """ Casting important columns to categorical makes polars fast """
    if scorefile:
        scorefile = scorefile.with_columns([
            pl.col("effect_allele").cast(pl.Categorical),
            pl.col("other_allele").cast(pl.Categorical),
            pl.col("effect_type").cast(pl.Categorical),
            pl.col("effect_allele_FLIP").cast(pl.Categorical),
            pl.col("other_allele_FLIP").cast(pl.Categorical),
            pl.col("accession").cast(pl.Categorical)
        ])
    if target:
        target = target.with_columns([
            pl.col("REF").cast(pl.Categorical),
            pl.col("ALT").cast(pl.Categorical)
        ])

    return scorefile, target


def _matched_colnames() -> list[str]:
    return ['chr_name', 'chr_position', 'effect_allele', 'effect_allele_FLIP', 'other_allele', 'other_allele_FLIP',
            'effect_weight', 'effect_type', 'accession', 'ID', 'REF', 'ALT', 'match_type', 'is_multiallelic']
