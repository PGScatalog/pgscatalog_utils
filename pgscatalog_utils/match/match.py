import logging

import polars as pl

from pgscatalog_utils.match.postprocess import postprocess_matches
from pgscatalog_utils.match.write import write_log

logger = logging.getLogger(__name__)


def get_all_matches(scorefile: pl.DataFrame, target: pl.DataFrame, remove_ambiguous: bool,
                    skip_flip: bool) -> pl.DataFrame:
    scorefile_cat, target_cat = _cast_categorical(scorefile, target)
    scorefile_oa = scorefile_cat.filter(pl.col("other_allele") != None)
    scorefile_no_oa = scorefile_cat.filter(pl.col("other_allele") == None)

    matches: list[pl.DataFrame] = []
    col_order = ['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
                 'accession', 'effect_allele_FLIP', 'other_allele_FLIP',
                 'ID', 'REF', 'ALT', 'is_multiallelic', 'matched_effect_allele', 'match_type']

    if scorefile_oa:
        logger.debug("Getting matches for scores with effect allele and other allele")
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="refalt").select(col_order))
        matches.append(_match_variants(scorefile_cat, target_cat, match_type="altref").select(col_order))
        if skip_flip is False:
            matches.append(_match_variants(scorefile_cat, target_cat, match_type="refalt_flip").select(col_order))
            matches.append(_match_variants(scorefile_cat, target_cat, match_type="altref_flip").select(col_order))

    if scorefile_no_oa:
        logger.debug("Getting matches for scores with effect allele only")
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_ref").select(col_order))
        matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_alt").select(col_order))
        if skip_flip is False:
            matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_ref_flip").select(col_order))
            matches.append(_match_variants(scorefile_no_oa, target_cat, match_type="no_oa_alt_flip").select(col_order))

    return pl.concat(matches).pipe(postprocess_matches, remove_ambiguous)


def check_match_rate(scorefile: pl.DataFrame, matches: pl.DataFrame, min_overlap: float, dataset: str) -> pl.DataFrame:
    scorefile: pl.DataFrame = scorefile.with_columns([
        pl.col('effect_type').cast(pl.Categorical),
        pl.col('accession').cast(pl.Categorical)])  # same dtypes for join
    match_log: pl.DataFrame = _join_matches(matches, scorefile, dataset)
    fail_rates: pl.DataFrame = (match_log.groupby('accession')
                                .agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
                                .with_column((pl.col('no_match') / pl.col('count')).alias('fail_rate'))
                                )
    pass_df: pl.DataFrame = pl.DataFrame()
    for accession, rate in zip(fail_rates['accession'].to_list(), fail_rates['fail_rate'].to_list()):
        if rate < (1 - min_overlap):
            df = pl.DataFrame({'accession': [accession], 'match_pass': [True], 'match_rate': [1 - rate]})
            pass_df = pl.concat([pass_df, df])
            logger.debug(f"Score {accession} passes minimum matching threshold ({1 - rate:.2%}  variants match)")
        else:
            df = pl.DataFrame({'accession': [accession], 'match_pass': [False], 'match_rate': [1 - rate]})
            pass_df = pl.concat([pass_df, df])
            logger.error(f"Score {accession} fails minimum matching threshold ({1 - rate:.2%} variants match)")

    # add match statistics to log and matches
    write_log((match_log.with_column(pl.col('accession').cast(str))
               .join(pass_df, on='accession', how='left')), dataset)

    return (matches.with_column(pl.col('accession').cast(str))
            .join(pass_df, on='accession', how='left'))



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
            effect_allele_column = "effect_allele"
        case 'altref':
            score_keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
            effect_allele_column = "effect_allele"
        case 'refalt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF", "ALT"]
            effect_allele_column = "effect_allele_FLIP"
        case 'altref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
            effect_allele_column = "effect_allele_FLIP"
        case 'no_oa_ref':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "REF"]
            effect_allele_column = "effect_allele"
        case 'no_oa_alt':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "ALT"]
            effect_allele_column = "effect_allele"
        case 'no_oa_ref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF"]
            effect_allele_column = "effect_allele_FLIP"
        case 'no_oa_alt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT"]
            effect_allele_column = "effect_allele_FLIP"
        case _:
            logger.critical(f"Invalid match strategy: {match_type}")
            raise Exception

    missing_cols = ['REF', 'ALT']
    if match_type.startswith('no_oa'):
        if match_type.startswith('no_oa_ref'):
            missing_cols = ['REF']
        else:
            missing_cols = ['ALT']
    join_cols = ['ID'] + missing_cols
    return (scorefile.join(target, score_keys, target_keys, how='inner')
            .with_columns([pl.col("*"),
                           pl.col(effect_allele_column).alias("matched_effect_allele"),
                           pl.lit(match_type).alias("match_type")])
            .join(target.select(join_cols), on="ID", how="inner"))  # get REF / ALT back after first join


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
            pl.col("ID").cast(pl.Categorical),
            pl.col("REF").cast(pl.Categorical),
            pl.col("ALT").cast(pl.Categorical)
        ])

    return scorefile, target
