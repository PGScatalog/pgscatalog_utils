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


def handle_multiallelic(df: pl.DataFrame, remove_multiallelic: bool) -> pl.DataFrame:
    is_ma: pl.Series = df['ALT'].str.contains(',')  # plink2 pvar multi-alleles are comma-separated
    if is_ma.sum() > 0:
        logger.debug("Multiallelic variants detected")
        if remove_multiallelic:
            logger.debug('Dropping multiallelic variants')
            return df[~is_ma]
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
