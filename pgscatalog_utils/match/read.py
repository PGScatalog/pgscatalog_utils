import polars as pl
import logging
from typing import NamedTuple

from pgscatalog_utils.match.preprocess import ugly_complement, handle_multiallelic, check_weights

logger = logging.getLogger(__name__)
log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(level=logging.DEBUG,
                            format=log_fmt,
                            datefmt='%Y-%m-%d %H:%M:%S')


def read_target(path: str, n_threads: int, remove_multiallelic: bool) -> pl.DataFrame:
    target: Target = _detect_target_format(path)
    d = {'column_1': str}  # column_1 is always CHROM. CHROM must always be a string
    df: pl.DataFrame = pl.read_csv(path, sep='\t', has_header=False, comment_char='#', dtype=d, n_threads=n_threads)
    df.columns = target.header

    match target.file_format:
        case 'bim':
            return (df[_default_cols()]
                    .pipe(ugly_complement))
        case 'pvar':
            return (df[_default_cols()]
                    .pipe(handle_multiallelic, remove_multiallelic=remove_multiallelic)
                    .pipe(ugly_complement))
        case _:
            logger.error("Invalid file format detected")
            raise Exception


def read_scorefile(path: str) -> pl.DataFrame:
    logger.debug("Reading scorefile")
    scorefile: pl.DataFrame = pl.read_csv(path, sep='\t', dtype={'chr_name': str})
    check_weights(scorefile)
    return scorefile


class Target(NamedTuple):
    """ Important summary information about a target genome. Cheap to compute (just reads the header). """
    file_format: str
    header: list[str]


def _detect_target_format(path: str) -> Target:
    file_format: str
    header: list[str]
    with open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                logging.debug("pvar format detected")
                file_format = 'pvar'
                header = _pvar_header(path)
                break
            else:
                logging.debug("bim format detected")
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


