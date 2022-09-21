import glob
import logging

import polars as pl

from pgscatalog_utils.match.preprocess import handle_multiallelic, complement_valid_alleles
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
        dfs.append(target.read().pipe(handle_multiallelic, remove_multiallelic=remove_multiallelic,
                                      file_format=target.file_format))

    return pl.concat(dfs).filter(pl.col("ID") != '.')


def read_scorefile(path: str) -> pl.DataFrame:
    logger.debug("Reading scorefile")
    scorefile: pl.DataFrame = (pl.read_csv(path, sep='\t', dtype={'chr_name': str})
    .pipe(complement_valid_alleles, flip_cols=['effect_allele', 'other_allele'])
    .with_columns([
        pl.col('accession').cast(pl.Categorical),
        pl.col("effect_type").cast(pl.Categorical)]))

    return scorefile
