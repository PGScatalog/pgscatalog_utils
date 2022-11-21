import gc
import io
import logging
import os
import pathlib
from dataclasses import dataclass
from itertools import islice

import polars as pl
import zstandard

import pgscatalog_utils.config as config
from pgscatalog_utils.match.tempdir import get_tmp_path

logger = logging.getLogger(__name__)


@dataclass
class Target:
    """ Class to detect and read a plink1/plink2 variant information file """
    file_format: str = None
    path: str = None
    compressed: bool = False
    low_memory: bool = True  # targets can be big, and use a lot of RAM when reading

    @classmethod
    def from_path(cls, path, low_memory):
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

        return cls(file_format=file_format, path=path, compressed=compressed, low_memory=low_memory)

    # @profile  # decorator needed to annotate memory profiles, but will cause NameErrors outside of profiling
    def read(self):
        if self.low_memory:
            if self.compressed:
                logger.debug("Reading compressed chunks from target genome (slower, lower RAM usage)")
                return self._read_compressed_chunks()
            else:
                logger.debug("Reading uncompressed chunks from target genome (slower, lower RAM usage)")
                return self._read_uncompressed_chunks()
        else:
            if self.compressed:
                logger.debug("Reading compressed target genome (fast mode, high RAM usage)")
                return self._read_compressed()
            else:
                logger.debug("Reading uncompressed target genome (fast mode, high RAM usage)")
                return self._read_uncompressed()

    def _read_compressed(self) -> pl.LazyFrame:
        """ Read a zst compressed target as quickly as possible """
        with open(self.path, 'rb') as fh:
            dctx = zstandard.ZstdDecompressor()
            with dctx.stream_reader(fh) as reader:
                dtypes = _get_col_dtypes(self.file_format)
                col_idxs, new_col_names = _default_cols(self.file_format)

                fn: str = pathlib.Path(self.path).stem + ".ipc.zst"
                fout = get_tmp_path("input", fn)

                (pl.read_csv(reader, sep='\t', has_header=False, comment_char='#',
                             dtype=dtypes,
                             columns=col_idxs,
                             new_columns=new_col_names,
                             n_threads=config.N_THREADS)
                 .write_ipc(fout, compression='zstd'))
                return pl.scan_ipc(fout, memory_map=False)

    def _read_uncompressed(self) -> pl.LazyFrame:
        """ Read an uncompressed target as quickly as possible. Uses up to 16GB RAM on 1000 genomes pvar. """
        dtypes = _get_col_dtypes(self.file_format)
        col_idxs, new_col_names = _default_cols(self.file_format)

        fn: str = pathlib.Path(self.path).stem + ".ipc.zst"
        fout: str = get_tmp_path("input", fn)

        (pl.read_csv(self.path, sep='\t', has_header=False, comment_char='#',
                     dtype=dtypes,
                     columns=col_idxs,
                     new_columns=new_col_names,
                     n_threads=config.N_THREADS)
         .write_ipc(fout, compression='zstd'))
        return pl.scan_ipc(fout, memory_map=False)

    def _read_uncompressed_chunks(self) -> pl.LazyFrame:
        """ Read a CSV using a BufferedReader in batches to reduce memory usage.

        Reads 1 million variant chunks and immediately writes to feather format in a temporary directory.

        Read all temporary feather files and return a big pl.DataFrame. Reading feather is fast, and preserves dtypes.

        Uses ~ 2GB
        """
        dtypes = _get_col_dtypes(self.file_format)
        col_idxs, new_col_names = _default_cols(self.file_format)

        batch_n = 0
        batch_size = int(1e6)
        with open(self.path, 'rb') as f:
            while True:
                line_batch = b''.join(islice(f, batch_size))
                if not line_batch:
                    break

                fn: str = str(batch_n) + ".ipc.zst"
                fout: str = get_tmp_path("input", fn)

                (pl.read_csv(line_batch, sep='\t', has_header=False, comment_char='#',
                             dtype=dtypes,
                             columns=col_idxs,
                             new_columns=new_col_names,
                             n_threads=config.N_THREADS)
                 .write_ipc(fout, compression='zstd'))
                batch_n += 1

        gc.collect()  # just to be safe
        logger.debug(f"{batch_n} batches staged in temporary directory {config.TEMPDIR}")
        return pl.scan_ipc(os.path.join(config.TEMPDIR.name, "input", "*.ipc.zst"), memory_map=False)

    def _read_compressed_chunks(self) -> pl.LazyFrame:
        """ Like _read_uncompressed_chunks, but read chunks of bytes and handle incomplete rows

        zstd returns chunks of bytes, not lines, but encoding utf-8 will be faster in rust and polars
         """
        logger.debug("Started reading zstd compressed data")
        dtypes = _get_col_dtypes(self.file_format)
        columns, new_col_names = _default_cols(self.file_format)

        n_chunks = 0
        with open(self.path, 'rb') as fh:
            dctx = zstandard.ZstdDecompressor()
            chunk_buffer = b''

            for chunk in dctx.read_to_iter(fh, read_size=int(1e8), write_size=int(1e8)):
                if not chunk:
                    logger.debug("Finished reading zstd compressed chunks")
                    break

                end = chunk.rfind(b'\n') + 1  # only want to read complete rows, which end in \n
                if chunk_buffer:
                    row_chunk = b''.join([chunk_buffer, chunk[:end]])
                    chunk_buffer = b''
                else:
                    row_chunk = chunk[:end]

                fn: str = str(n_chunks) + ".ipc.zst"
                fout: str = get_tmp_path("input", fn)

                (pl.read_csv(row_chunk, sep='\t', has_header=False, comment_char='#',
                             dtype=dtypes,
                             columns=columns,
                             new_columns=new_col_names,
                             n_threads=config.N_THREADS)
                 .write_ipc(fout, compression='zstd'))

                chunk_buffer = b''.join([chunk_buffer, chunk[end:]])
                n_chunks += 1

            gc.collect()  # just to be safe
            logger.debug(f"{n_chunks} chunks")  # write_size will change n_chunks
            return pl.scan_ipc(os.path.join(config.TEMPDIR.name, "input", "*.ipc.zst"), memory_map=False)


def _get_col_dtypes(file_format):
    """ Manually set up dtypes to save memory. Repeated strings like REF / ALT / CHROM work best as pl.Categorical.

    ID shouldn't be pl.Categorical, or you'll create a massive string cache and waste RAM """
    match file_format:
        case 'bim':
            # 1. Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
            # 2. Variant identifier
            # 3. Position in morgans or centimorgans (safe to use dummy value of '0')
            # 4. Base-pair coordinate (1-based; limited to 231-2)
            # 5. Allele 1 (corresponding to clear bits in .bed; usually minor)
            # 6. Allele 2 (corresponding to set bits in .bed; usually major)
            d = {'column_1': pl.Categorical, 'column_2': pl.Utf8, 'column_3': pl.Float64, 'column_4': pl.UInt64,
                 'column_5': pl.Categorical, 'column_6': pl.Utf8}
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


def _default_cols(file_format) -> tuple[list[int], list[str]]:
    """ Return a list of column integers to keep, assuming plink default column sets """
    match file_format:
        case 'bim':
            idxs = [0, 1, 3, 4, 5]  # see _get_col_dtypes, dropping centimorgans
            names = ['#CHROM', 'ID', 'POS', 'REF', 'ALT']  # technically A1/A2, but it's ok
            return idxs, names
        case 'pvar':
            idxs = [0, 1, 2, 3, 4]  # dropping QUAL FILTER INFO etc
            names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
            return idxs, names
        case _:
            logger.critical("Trying to get column idx for an invalid file format, TWENTY THREE NINETEEN")
            raise Exception

