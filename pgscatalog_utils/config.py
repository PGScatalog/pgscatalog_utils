import atexit
import logging
import os
import tempfile

import polars as pl

from pgscatalog_utils.match import tempdir

N_THREADS: int = 1  # dummy value, is reset by args.n_threads (default: 1)
OUTDIR: str = "."  # dummy value, reset by args.outdir
TEMPDIR: tempfile.TemporaryDirectory

logger = logging.getLogger(__name__)


def setup_tmpdir(outdir, combine=False):
    if combine:
        work_dir = "work_combine"
        dirs = [work_dir]
    else:
        work_dir = "work_match"
        dirs = [work_dir, "matches"]

    for d in dirs:
        if os.path.exists(d):
            logger.critical(f"{d} already exists, bailing out")
            logger.critical("Please choose a different --outdir or clean up")
            raise SystemExit(1)

    global TEMPDIR
    os.mkdir(os.path.join(outdir, work_dir))
    TEMPDIR = tempfile.TemporaryDirectory(dir=os.path.join(outdir, work_dir))


def setup_cleaning():
    logger.debug(F"Temporary directory set up: {TEMPDIR.name}")
    atexit.register(tempdir.cleanup)


def set_logging_level(verbose: bool):
    log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"

    if verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format=log_fmt,
                            datefmt='%Y-%m-%d %H:%M:%S')
        logging.debug("Verbose logging enabled")
    else:
        logging.basicConfig(level=logging.WARNING,
                            format=log_fmt,
                            datefmt='%Y-%m-%d %H:%M:%S')


def setup_polars_threads(n: int):
    global N_THREADS
    N_THREADS = n
    os.environ['POLARS_MAX_THREADS'] = str(N_THREADS)
    logger.debug(f"Using {N_THREADS} threads to read CSVs")
    logger.debug(f"polars threadpool size: {pl.threadpool_size()}")

    if pl.threadpool_size() != N_THREADS:
        logger.warning(f"polars threadpool doesn't match -n argument ({pl.threadpool_size()} vs {n})")
        logger.info("To silence this warning, set POLARS_MAX_THREADS to match -n before running combine_matches, e.g.:")
        logger.info("$ export POLARS_MAX_THREADS=x")
        logger.info("$ combine_matches ... -n x")
