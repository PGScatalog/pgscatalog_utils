import gzip
import logging

import polars as pl

logger = logging.getLogger(__name__)


def write_log(df: pl.DataFrame, dataset: str) -> None:
    logger.debug("Compressing and writing log")
    with gzip.open(f"{dataset}_log.csv.gz", 'wb') as f:
        df.pipe(_prettify_log).write_csv(f)


def _prettify_log(df: pl.DataFrame) -> pl.DataFrame:
    keep_cols = ["chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight", "effect_type",
                 "accession", "row_nr", "ID", "REF", "ALT", "matched_effect_allele", "match_type", "is_multiallelic",
                 "ambiguous", "duplicate", "best_match", "dataset", "score_pass", "match_rate"]
    return df.select(keep_cols).select(pl.exclude("^.*_right"))
