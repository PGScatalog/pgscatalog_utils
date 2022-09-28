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
    path: str = None
    compressed: bool = False

    @classmethod
    def from_path(cls, path):
        """ Create a Target object from a path. Cheaply detect file format and headers. """
        try:
            with open(path, 'r') as f:
                file_format = _get_format(f)
                compressed = False
        except UnicodeDecodeError:
            logger.error("Can't open target as a text file, so trying to read zstd compressed binary file")
            with open(path, 'rb') as f:
                dctx = zstandard.ZstdDecompressor()
                stream_reader = dctx.stream_reader(f)
                text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
                file_format = _get_format(text_stream)
                compressed = True

        return cls(file_format=file_format, path=path, compressed=compressed)

    #@profile
    def read(self):
        if self.compressed:
            return self._read_compressed_chunks().lazy()
        else:
            return self._read_uncompressed_chunks().lazy()

    def _read_uncompressed_chunks(self):
        """ Read a CSV using a BufferedIOReader. This is a bit slower than pl.read_csv() (30s vs 5s).

        Lots of testing showed that lazy scanning and native polars reading used a lot of RAM, then freed a bunch.
        Plotting RAM usage against time looked like a spiky hedgehog.

        This function linearly consumes RAM in a more linear way by:
            1. Reading a batch of lines
            2. Dropping unused columns
            3. Setting categorical dtypes on read
            4. Don't rechunk until later
        """
        logger.debug("Started reading uncompressed chunks")

        df_lst = []
        dtypes = _get_col_dtypes(self.file_format)
        col_idxs = _get_default_col_idx(self.file_format)
        new_col_names = _default_cols()

        with open(self.path, "rb") as f:
            while True:
                buffer = b''.join(f.readlines(int(1e6)))

                if not buffer:
                    break

                df = (pl.read_csv(buffer, sep='\t', has_header=False, comment_char='#', n_threads=1,
                                  dtype=dtypes,
                                  columns=col_idxs,
                                  new_columns=new_col_names,
                                  rechunk=False))

                df_lst.append(df)

        logger.debug("Finished reading uncompressed chunks")
        logger.debug("Concatenating chunked data frames")
        return pl.concat(df_lst, rechunk=False)

    def _read_compressed_chunks(self):
        """ Like _read_uncompressed_chunks, but read chunks of bytes and handle incomplete rows

        zstd returns chunks of bytes, not lines, but encoding utf-8 will be faster in rust and polars
         """
        logger.debug("Started reading zstd compressed data")
        df_lst = []
        dtypes = _get_col_dtypes(self.file_format)
        columns = _get_default_col_idx(self.file_format)
        new_col_names = _default_cols()

        with open(self.path, 'rb') as fh:
            dctx = zstandard.ZstdDecompressor()
            chunk_buffer = b''

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

        logger.debug("Finished reading zstd compressed chunks")
        logger.debug("Concatenating chunked data frames")
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
    """ Manually set up dtypes. pl.Categorical saves a lot of RAM vs pl.Utf8 """
    match file_format:
        case 'bim':
            # 1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
            # 2. Variant identifier
            # 3. Position in morgans or centimorgans (safe to use dummy value of '0')
            # 4. Base-pair coordinate (1-based; limited to 231-2)
            # 5. Allele 1 (corresponding to clear bits in .bed; usually minor)
            # 6. Allele 2 (corresponding to set bits in .bed; usually major)
            d = {'column_1': pl.Categorical, 'column_2': pl.Categorical, 'column_3': pl.Float64, 'column_4': pl.UInt64,
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
            d = {'column_1': pl.Categorical, 'column_2': pl.UInt64, 'column_3': pl.Categorical, 'column_4': pl.Categorical,
                 'column_5': pl.Utf8, 'column_6': pl.Float32, 'column_7': pl.Utf8, 'column_8': pl.Utf8}
            # can't cast ALT to cat yet, because of multiallelic variants!
        case _:
            logger.critical("Trying to set header dtypes for an invalid file format, time to explode")
            raise Exception
    return d


def _get_format(fh) -> str:
    file_format = None
    logger.debug(f"Scanning header to get file format")
    for line in fh:
        if line.startswith('#'):
            logger.debug("pvar format detected")
            file_format = 'pvar'
            break
        else:
            logger.debug("bim format detected")
            file_format = 'bim'
            break

    return file_format


def _default_cols() -> list[str]:
    """ Standardise column names in a target genome """
    return ['#CHROM', 'POS', 'ID', 'REF', 'ALT']

