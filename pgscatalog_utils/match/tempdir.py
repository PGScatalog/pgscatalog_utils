import logging
import os

from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def get_tmp_path(subdir: str, fn: str) -> str:
    """ Create a subdirectory in the tempodir and return a full path to fn

     subdir: 'input', fn: 'test.txt' -> '/path/tp/tmpdir/input/test.txt'
     """
    path: str = os.path.join(config.TEMPDIR.name, subdir)
    if not os.path.exists(path):
        os.mkdir(path)

    return os.path.join(path, fn)


def cleanup():
    """ A temporary directory is used to store staged data in feather format.

    tempdir/
    ├── target
    ├── scorefile
    └── matched

    Data are staged to disk for a few different reasons:

    target/ and scorefile/:
    - Raw text data may be compressed with zstd or gzip, which can't be lazily read
    - Parsing a very large file causes big RAM spikes
        - To mitigate this, optionally read and parse in batches (compressed + uncompressed)
    - ipc are uncompressed to allow fast memory mapping
    - These files should always be cleaned up by python or the host OS (SIGTERM breaks atexit)

    matched/
    - Split - apply - combine means a common use case is to just write matches and exit
    - In this case, files are saved (moved to args.outdir) and used as input to combine_matches
    - In other cases, post-processing of matches makes complex query plans, which failed when collected
        - Re-scanning collected files on disk prevents this problem
    - ipc are compressed to save space. Further processing takes some time, so decompression is ok.

    This function is registered with atexit to run when the program ends.
    """
    logger.debug(f"Cleaning up tempdir path {config.TEMPDIR}")
    config.TEMPDIR.cleanup()
