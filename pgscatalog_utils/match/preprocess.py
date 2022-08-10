import polars as pl
import logging

logger = logging.getLogger(__name__)


def ugly_complement(df: pl.DataFrame) -> pl.DataFrame:
    """ Complementing alleles with a pile of regexes seems weird, but polars string functions are currently limited
    (i.e. no str.translate). This is fast, and I stole the regex idea from Scott.
    """
    logger.debug("Complementing target alleles")
    return df.with_columns([
        (pl.col("REF").str.replace_all("A", "V")
         .str.replace_all("T", "X")
         .str.replace_all("C", "Y")
         .str.replace_all("G", "Z")
         .str.replace_all("V", "T")
         .str.replace_all("X", "A")
         .str.replace_all("Y", "G")
         .str.replace_all("Z", "C"))
        .alias("REF_FLIP"),
        (pl.col("ALT").str.replace_all("A", "V")
         .str.replace_all("T", "X")
         .str.replace_all("C", "Y")
         .str.replace_all("G", "Z")
         .str.replace_all("V", "T")
         .str.replace_all("X", "A")
         .str.replace_all("Y", "G")
         .str.replace_all("Z", "C"))
        .alias("ALT_FLIP")
    ])


def handle_multiallelic(df: pl.DataFrame, remove_multiallelic: bool, pvar: bool) -> pl.DataFrame:
    # plink2 pvar multi-alleles are comma-separated
    df: pl.DataFrame = (df.with_column(
        pl.when(pl.col("ALT").str.contains(','))
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias('is_multiallelic')))

    if df['is_multiallelic'].sum() > 0:
        logger.debug("Multiallelic variants detected")
        if remove_multiallelic:
            if not pvar:
                logger.warning("--remove_multiallelic requested for bim format, which already contains biallelic "
                               "variant representations only")
            logger.debug('Dropping multiallelic variants')
            return df[~df['is_multiallelic']]
        else:
            logger.debug("Exploding dataframe to handle multiallelic variants")
            df.replace('ALT', df['ALT'].str.split(by=','))  # turn ALT to list of variants
            return df.explode('ALT')  # expand the DF to have all the variants in different rows
    else:
        logger.debug("No multiallelic variants detected")
        return df


def check_weights(df: pl.DataFrame) -> None:
    weight_count = df.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele']).count()['count']

    if any(weight_count > 1):
        logger.error("Multiple effect weights per variant per accession detected")
        raise Exception


def _annotate_multiallelic(df: pl.DataFrame) -> pl.DataFrame:
    df.with_column(pl.when(pl.col("ALT").str.contains(',')).then(pl.lit(True)).otherwise(pl.lit(False)).alias('is_multiallelic'))