import polars as pl
import logging

logger = logging.getLogger(__name__)


def postprocess_matches(df: pl.DataFrame, remove_ambiguous: bool) -> pl.DataFrame:
    df = _label_biallelic_ambiguous(df)
    if remove_ambiguous:
        logger.debug("Removing ambiguous matches")
        return df.filter(pl.col("ambiguous") == False)
    else:
        logger.debug("Keeping best possible match from ambiguous matches")
        # pick the best possible match from the ambiguous matches
        # EA = REF and OA = ALT or EA = REF and OA = None
        ambiguous: pl.DataFrame = df.filter((pl.col("ambiguous") == True) & \
                                            (pl.col("match_type") == "refalt") |
                                            (pl.col("ambiguous") == True) & \
                                            (pl.col("match_type") == "no_oa_ref"))
        unambiguous: pl.DataFrame = df.filter(pl.col("ambiguous") == False)
        return pl.concat([ambiguous, unambiguous])


def _label_biallelic_ambiguous(df: pl.DataFrame) -> pl.DataFrame:
    # A / T or C / G may match multiple times
    df = df.with_columns([
        pl.col(["effect_allele", "other_allele", "REF", "ALT", "REF_FLIP", "ALT_FLIP"]).cast(str),
        pl.lit(True).alias("ambiguous")
    ])

    return (df.with_column(
        pl.when((pl.col("effect_allele") == pl.col("ALT_FLIP")) | (pl.col("effect_allele") == pl.col("REF_FLIP")))
        .then(pl.col("ambiguous"))
        .otherwise(False))).pipe(_get_distinct_weights)


def _get_distinct_weights(df: pl.DataFrame) -> pl.DataFrame:
    """ Get a single effect weight for each matched variant per accession """
    count: pl.DataFrame = df.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele']).count()
    singletons: pl.DataFrame = (count.filter(pl.col('count') == 1)[:, "accession":"effect_allele"]
                                .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele'], how='left'))

    # TODO: something more complex than .unique()?
    # TODO: prioritise unambiguous -> ref -> alt -> ref_flip -> alt_flip
    dups: pl.DataFrame = (count.filter(pl.col('count') > 1)[:, "accession":"effect_allele"]
                          .join(df, on=['accession', 'chr_name', 'chr_position', 'effect_allele'], how='left')
                          .distinct(subset=['accession', 'chr_name', 'chr_position', 'effect_allele']))
    distinct: pl.DataFrame = pl.concat([singletons, dups])

    assert all((distinct.groupby(['accession', 'chr_name', 'chr_position', 'effect_allele']).count()['count']) == 1), \
        "Duplicate effect weights for a variant"

    return distinct
