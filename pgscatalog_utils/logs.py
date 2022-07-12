import logging


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