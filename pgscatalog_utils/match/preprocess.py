import logging

import polars as pl

logger = logging.getLogger(__name__)


def filter_target(df: pl.DataFrame) -> pl.DataFrame:
    """ Remove variants that won't be matched against the scorefile

    Chromosomes 1 - 22, X, and Y with an efficient join. Remmove variants with missing identifiers also
    """
    logger.debug("Filtering target to include chromosomes 1 - 22, X, Y")
    chroms = [str(x) for x in list(range(1, 23)) + ['X', 'Y']]
    return df.filter((pl.col('#CHROM').is_in(chroms)) & (pl.col('ID') != '.'))


def complement_valid_alleles(df: pl.DataFrame, flip_cols: list[str]) -> pl.DataFrame:
    """ Improved function to complement alleles. Will only complement sequences that are valid DNA.
    """
    for col in flip_cols:
        logger.debug(f"Complementing column {col}")
        new_col = col + '_FLIP'
        df = df.with_column(
            pl.when(pl.col(col).str.contains('^[ACGT]+$'))
            .then(pl.col(col).str.replace_all("A", "V")
                  .str.replace_all("T", "X")
                  .str.replace_all("C", "Y")
                  .str.replace_all("G", "Z")
                  .str.replace_all("V", "T")
                  .str.replace_all("X", "A")
                  .str.replace_all("Y", "G")
                  .str.replace_all("Z", "C"))
            .otherwise(pl.col(col))
            .alias(new_col)
        )
    return df


def annotate_multiallelic(df: pl.DataFrame) -> pl.DataFrame:
    """ Identify variants that are multiallelic with a column flag """
    # plink2 pvar multi-alleles are comma-separated
    df: pl.DataFrame = (df.with_column(
        pl.when(pl.col("ALT").str.contains(','))
        .then(pl.lit(True))
        .otherwise(pl.lit(False))
        .alias('is_multiallelic')))

    if (df.get_column('is_multiallelic')).any():
        logger.debug("Exploding dataframe to handle multiallelic variants")
        df.replace('ALT', df['ALT'].str.split(by=','))  # turn ALT to list of variants
        return df.explode('ALT')  # expand the DF to have all the variants in different rows
    else:
        logger.debug("No multiallelic variants detected")
        return df
