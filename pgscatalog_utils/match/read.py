import glob
import logging
from typing import NamedTuple

import polars as pl

from pgscatalog_utils.match.preprocess import handle_multiallelic, complement_valid_alleles

logger = logging.getLogger(__name__)


def read_target(path: str, remove_multiallelic: bool, single_file: bool = False,
                chrom: str = "") -> pl.DataFrame:
    target: Target = _detect_target_format(path)
    d = {'column_1': str}  # column_1 is always CHROM. CHROM must always be a string

    if single_file:
        logger.debug(f"Scanning target genome for chromosome {chrom}")
        # scan target and filter to reduce memory usage on big files
        df: pl.DataFrame = (
            pl.scan_csv(path, sep='\t', has_header=False, comment_char='#', dtype=d)
            .filter(pl.col('column_1') == chrom)
            .collect())

        if df.is_empty():
            logger.warning(f"Chromosome missing from target genome: {chrom}")
            return df
    else:
        logger.debug(f"Reading target {path}")
        df: pl.DataFrame = pl.read_csv(path, sep='\t', has_header=False, comment_char='#', dtype=d)

    df.columns = target.header

    match target.file_format:
        case 'bim':
            return (df.select(_default_cols())
                    .filter(pl.col('ID') != '.')  # remove missing IDs
                    .pipe(handle_multiallelic, remove_multiallelic=remove_multiallelic, pvar=False))
        case 'pvar':
            return (df.select(_default_cols())
                    .filter(pl.col('ID') != '.')
                    .pipe(handle_multiallelic, remove_multiallelic=remove_multiallelic, pvar=True))
        case _:
            logger.error("Invalid file format detected")
            raise Exception


def read_scorefile(path: str) -> pl.DataFrame:
    logger.debug("Reading scorefile")
    scorefile: pl.DataFrame = (pl.read_csv(path, sep='\t', dtype={'chr_name': str})
                               .pipe(complement_valid_alleles, flip_cols=['effect_allele', 'other_allele'])
                               .with_columns([
        pl.col('accession').cast(pl.Categorical),
        pl.col("effect_type").cast(pl.Categorical)]))

    return scorefile


class Target(NamedTuple):
    """ Important summary information about a target genome. Cheap to compute (just reads the header). """
    file_format: str
    header: list[str]


def _detect_target_format(path: str) -> Target:
    file_format: str
    header: list[str]

    if "*" in path:
        logger.debug("Detecting target file format")
        path = glob.glob(path)[0]  # guess format from first file in directory

    with open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                logger.debug("pvar format detected")
                file_format = 'pvar'
                header = _pvar_header(path)
                break
            else:
                logger.debug("bim format detected")
                file_format = 'bim'
                header = _bim_header()
                break

    return Target(file_format, header)


def _default_cols() -> list[str]:
    return ['#CHROM', 'POS', 'ID', 'REF', 'ALT']  # only columns we want from a target genome


def _pvar_header(path: str) -> list[str]:
    """ Get the column names from the pvar file (not constrained like bim, especially when converted from VCF) """
    line: str = '#'
    with open(path, 'rt') as f:
        while line.startswith('#'):
            line: str = f.readline()
            if line.startswith('#CHROM'):
                return line.strip().split('\t')


def _bim_header() -> list[str]:
    return ['#CHROM', 'ID', 'CM', 'POS', 'REF', 'ALT']
