import logging
import os

try:
    POLARS_MAX_THREADS: int = int(os.getenv('POLARS_MAX_THREADS'))
except TypeError:
    POLARS_MAX_THREADS = 1  # not defined, it's better to be slow than set to n_cores (polars default)


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
