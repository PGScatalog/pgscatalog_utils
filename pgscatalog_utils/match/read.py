
import logging

import polars as pl
import pgscatalog_utils.config as config
from pgscatalog_utils.match.preprocess import annotate_multiallelic, complement_valid_alleles, filter_target
from pgscatalog_utils.target import Target

logger = logging.getLogger(__name__)


def read_target(paths: list[str], low_memory: bool) -> pl.LazyFrame:
    targets: list[Target] = [Target.from_path(x, low_memory) for x in paths]

    logger.debug("Reading all target data complete")
    # handling multiallelic requires str methods, so don't forget to cast back or matching will break
    return (pl.concat([x.read() for x in targets])
            .pipe(filter_target)
            .pipe(annotate_multiallelic)
            .with_column(pl.col('ALT').cast(pl.Categorical))).lazy()


def read_scorefile(path: str) -> pl.LazyFrame:
    logger.debug("Reading scorefile")
    dtypes = {'chr_name': pl.Categorical,
              'chr_position': pl.UInt64,
              'effect_allele': pl.Utf8,  # str functions required to complement
              'other_allele': pl.Utf8,
              'effect_type': pl.Categorical,
              'accession': pl.Categorical}
    return (pl.read_csv(path, sep='\t', dtype=dtypes, n_threads=config.POLARS_MAX_THREADS)
            .lazy()
            .pipe(complement_valid_alleles, flip_cols=['effect_allele', 'other_allele'])).with_columns([
        pl.col("effect_allele").cast(pl.Categorical),
        pl.col("other_allele").cast(pl.Categorical),
        pl.col("effect_allele_FLIP").cast(pl.Categorical),
        pl.col("other_allele_FLIP").cast(pl.Categorical)])
