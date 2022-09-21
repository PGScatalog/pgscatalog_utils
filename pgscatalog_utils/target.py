import zstandard
from dataclasses import dataclass
import io
import logging
import polars as pl

logger = logging.getLogger(__name__)


@dataclass
class Target:
    """ Class to detect and read a plink1/plink2 variant information file """
    file_format: str = None
    header: list[str] = None
    path: str = None
    compressed: bool = False

    @classmethod
    def from_path(cls, path):
        """ Create a Target object from a path. Cheaply detect file format and headers. """
        try:
            with open(path, 'r') as f:
                file_format, header = _get_header(f)
                compressed = False
        except UnicodeDecodeError:
            logger.error("Can't open target as a text file, so trying to read zstd compressed binary file")
            with open(path, 'rb') as f:
                dctx = zstandard.ZstdDecompressor()
                stream_reader = dctx.stream_reader(f)
                text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
                file_format, header = _get_header(text_stream)
                compressed = True

        return cls(file_format=file_format, path=path, header=header, compressed=compressed)

    def read(self) -> pl.DataFrame:
        """ Read variant information into a polars df (expensive operation). Automatically handle compressed data.  """
        # column_1 is always CHROM, which must always be a string or X/Y/MT/PAR will break inferred dtypes
        logger.debug("Reading target into memory")
        chrom_dtype = {'column_1': str}
        if self.compressed:
            with open(self.path, 'rb') as f:
                dctx = zstandard.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    df: pl.DataFrame = pl.read_csv(reader, sep='\t', has_header=False, comment_char='#', dtype=chrom_dtype)
                    df.columns = self.header
                    return df.select(_default_cols())
        else:
            df: pl.DataFrame = pl.read_csv(self.path, sep='\t', has_header=False, comment_char='#', dtype=chrom_dtype)
            df.columns = self.header
            return df.select(_default_cols())


def _get_header(fh) -> tuple[str, list[str]]:
    header = None
    file_format = None
    logger.debug(f"Scanning header to get file format and column names")
    for line in fh:
        if line.startswith('#'):
            logger.debug("pvar format detected")
            file_format = 'pvar'
            header = _pvar_header(fh)
            break
        else:
            logger.debug("bim format detected")
            file_format = 'bim'
            header = _bim_header()
            break

    return file_format, header


def _pvar_header(fh) -> list[str]:
    """ Get the column names from the pvar file (not constrained like bim, especially when converted from VCF) """
    line: str = '#'
    while line.startswith('#'):
        line: str = fh.readline()
        if line.startswith('#CHROM'):
            return line.strip().split('\t')


def _bim_header() -> list[str]:
    return ['#CHROM', 'ID', 'CM', 'POS', 'REF', 'ALT']


def _default_cols() -> list[str]:
    return ['#CHROM', 'POS', 'ID', 'REF', 'ALT']  # only columns we want from a target genome
