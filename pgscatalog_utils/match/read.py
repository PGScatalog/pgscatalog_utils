import glob
import logging

import polars as pl

from pgscatalog_utils.match.preprocess import handle_multiallelic, complement_valid_alleles, filter_target
from pgscatalog_utils.target import Target

logger = logging.getLogger(__name__)


def read_target(path: str, remove_multiallelic: bool) -> pl.DataFrame:
    """ Read one or more targets from a path (may contain a wildcard) """

    if '*' in path:
        logger.debug("Wildcard detected in target path: finding all matching files")
        paths: list[str] = glob.glob(path)
    else:
        logger.debug("")
        paths: list[str] = [path]

    targets: list[Target] = [Target.from_path(x) for x in paths]
    dfs: list[pl.DataFrame] = []
    for target in targets:
        assert target.file_format in ['bim', 'pvar']
        dfs.append(target.read())

    logger.debug("Reading all target data complete")
    # explicitly rechunk now, because reading is complete and the input data were read unchunked to save memory
    # only pipe functions once rechunking has happened to improve speed
    # handling multiallelic requires str methods, so don't forget to cast back or matching will break
    return (pl.concat(dfs, rechunk=True)
            .pipe(filter_target)
            .pipe(handle_multiallelic, remove_multiallelic=remove_multiallelic)
            .with_column(pl.col('ALT').cast(pl.Categorical)))


def read_scorefile(path: str) -> pl.DataFrame:
    logger.debug("Reading scorefile")
    dtypes = {'chr_name': pl.Categorical,
              'chr_position': pl.UInt64,
              'effect_allele': pl.Utf8,  # str functions required to complement
              'other_allele': pl.Utf8,
              'effect_type': pl.Categorical,
              'accession': pl.Categorical}
    return (pl.scan_csv(path, sep='\t', dtype=dtypes)
    .pipe(complement_valid_alleles, flip_cols=['effect_allele', 'other_allele'])
    .with_columns([
        pl.col("effect_allele").cast(pl.Categorical),
        pl.col("other_allele").cast(pl.Categorical),
        pl.col("effect_allele_FLIP").cast(pl.Categorical),
        pl.col("other_allele_FLIP").cast(pl.Categorical)
    ]))
