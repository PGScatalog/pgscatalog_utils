import logging

POLARS_MAX_THREADS: int = 1  # dummy value, is reset by args.n_threads (default: 1)
OUTDIR: str = "."  # dummy value, reset by args.outdir


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
