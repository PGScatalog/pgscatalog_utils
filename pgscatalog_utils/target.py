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

    # @profile
    def read(self):
        # this function is responsible for returning dfs allocated to contiguous memory, so manually rechunk
        if self.compressed:
            return self._read_compressed_chunks().rechunk()
        else:
            batch_size = 10000000
            n_rows_read = 0
            df_lst = []
            while True:
                df_lst.append(self._read_batch(batch_size=batch_size, n_skip=n_rows_read))
                n_rows_read = n_rows_read + batch_size

                if df_lst[-1].shape[0] < batch_size:
                    logger.debug("Finished reading final batch")
                    break

            return pl.concat(df_lst, rechunk=True)

    def _read_batch(self, batch_size, n_skip):
        logger.debug(f"{n_skip} target variants read, reading next batch")
        assert not self.compressed
        # TODO: lazy frame it
        logger.debug("Reading uncompressed data")
        return pl.read_csv(self.path, sep='\t', has_header=False, comment_char='#', n_threads=1,
                           dtype=_get_col_dtypes(self.file_format),
                           columns=_get_default_col_idx(self.file_format),
                           new_columns=_default_cols(),
                           rechunk=False,
                           n_rows=batch_size,
                           skip_rows_after_header=n_skip)

    def _read_compressed_chunks(self):
        logger.debug("Reading zstd compressed data")
        df_lst = []
        dtypes = _get_col_dtypes(self.file_format)
        columns = _get_default_col_idx(self.file_format)
        new_col_names = _default_cols()

        with open(self.path, 'rb') as fh:
            dctx = zstandard.ZstdDecompressor()
            chunk_buffer = b''

            # don't decode bytes stream to utf-8 with TextIOWrapper in python, polars + rust will be faster
            for chunk in dctx.read_to_iter(fh, read_size=int(1e+8)):  # read 100MB of compressed data per chunk
                if not chunk:
                    break

                end = chunk.rfind(b'\n') + 1  # only want to read complete rows
                if chunk_buffer:
                    row_chunk = b''.join([chunk_buffer, chunk[:end]])
                    chunk_buffer = b''
                else:
                    row_chunk = chunk[:end]

                df = pl.read_csv(row_chunk, sep='\t', has_header=False, comment_char='#', n_threads=1,
                                 dtype=dtypes,
                                 columns=columns,
                                 new_columns=new_col_names,
                                 rechunk=False)
                df_lst.append(df)
                chunk_buffer = b''.join([chunk_buffer, chunk[end:]])

        return pl.concat(df_lst, rechunk=False)


def _get_default_col_idx(file_format):
    # import default columns:
    #  ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
    match file_format:
        case 'bim':
            return [0, 1, 3, 4, 5]  # see _get_col_dtypes, dropping centimorgans
        case 'pvar':
            return [0, 1, 2, 3, 4]  # dropping QUAL FILTER INFO etc
        case _:
            logger.critical("Trying to get column idx for an invalid file format, TWENTY THREE NINETEEN")
            raise Exception


def _get_col_dtypes(file_format):
    """ Manually set up categorical dtypes """
    match file_format:
        case 'bim':
            # 1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
            # 2. Variant identifier
            # 3. Position in morgans or centimorgans (safe to use dummy value of '0')
            # 4. Base-pair coordinate (1-based; limited to 231-2)
            # 5. Allele 1 (corresponding to clear bits in .bed; usually minor)
            # 6. Allele 2 (corresponding to set bits in .bed; usually major)
            d = {'column_1': pl.Categorical, 'column_2': str, 'column_3': pl.Float64, 'column_4': pl.UInt64,
                 'column_5': pl.Categorical, 'column_6': pl.Categorical}
        case 'pvar':
            # 1. CHROM
            # 2. POS (base-pair coordinate)
            # 3. ID (variant ID; required)
            # 4. REF (reference allele)
            # 5. ALT (alternate alleles, comma-separated)
            # 6. QUAL (phred-scaled quality score for whether the locus is variable at all)
            # 7. FILTER ('PASS', '.', or semicolon-separated list of failing filter codes)
            # 8. INFO (semicolon-separated list of flags and key-value pairs, with types declared in header)
            d = {'column_1': pl.Categorical, 'column_2': pl.UInt64, 'column_3': pl.Utf8, 'column_4': pl.Categorical,
                 'column_5': pl.Utf8, 'column_6': pl.Float32, 'column_7': pl.Utf8, 'column_8': pl.Utf8}
            # can't cast ALT to cat yet, because of multiallelic variants!
        case _:
            logger.critical("Trying to set header dtypes for an invalid file format, time to explode")
            raise Exception
    return d


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
