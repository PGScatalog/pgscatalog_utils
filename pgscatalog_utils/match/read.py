
import logging
import typing

import polars as pl
from pgscatalog_utils.match.tempdir import get_tmp_path

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
            .with_column(pl.col('ALT').cast(pl.Categorical)))


def read_scorefile(path: str, chrom: typing.Union[str, None]) -> pl.LazyFrame:
    logger.debug("Reading scorefile")
    dtypes = {'chr_name': pl.Categorical,
              'chr_position': pl.UInt64,
              'effect_allele': pl.Utf8,  # str functions required to complement
              'other_allele': pl.Utf8,
              'effect_type': pl.Categorical,
              'accession': pl.Categorical}

    # parse CSV and write to temporary feather file
    # enforce laziness! scanning is very fast and saves memory
    fout: str = get_tmp_path("scorefile", "scorefile.ipc.zst")
    (pl.read_csv(path, sep='\t', dtype=dtypes).write_ipc(fout, compression='zstd'))
    ldf: pl.LazyFrame = pl.scan_ipc(fout, memory_map=False)

    if chrom is not None:
        logger.debug(f"--chrom set, filtering scoring file to chromosome {chrom}")
        ldf = ldf.filter(pl.col('chr_name') == chrom)  # add filter to query plan
    else:
        logger.debug("--chrom parameter not set, using all variants in scoring file")

    return (ldf.pipe(complement_valid_alleles, flip_cols=['effect_allele', 'other_allele'])).with_columns([
        pl.col("effect_allele").cast(pl.Categorical),
        pl.col("other_allele").cast(pl.Categorical),
        pl.col("effect_allele_FLIP").cast(pl.Categorical),
        pl.col("other_allele_FLIP").cast(pl.Categorical)])


