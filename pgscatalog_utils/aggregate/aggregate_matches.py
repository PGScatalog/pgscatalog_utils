import argparse
import logging
import os

import polars as pl

from pgscatalog_utils import config

from pgscatalog_utils.match.filter import filter_scores
from pgscatalog_utils.match.log import make_summary_log, check_log_count
from pgscatalog_utils.match.read import read_scorefile
from pgscatalog_utils.match.write import write_log, write_out

logger = logging.getLogger(__name__)


def aggregate_matches():
    args = _parse_args()
    config.set_logging_level(args.verbose)

    config.POLARS_MAX_THREADS = args.n_threads
    os.environ['POLARS_MAX_THREADS'] = str(config.POLARS_MAX_THREADS)
    # now the environment variable, parsed argument args.n_threads, and threadpool should agree
    logger.debug(f"Setting POLARS_MAX_THREADS environment variable: {os.getenv('POLARS_MAX_THREADS')}")
    logger.debug(f"Using {config.POLARS_MAX_THREADS} threads to read CSVs")
    logger.debug(f"polars threadpool size: {pl.threadpool_size()}")

    with pl.StringCache():
        scorefile = read_scorefile(path=args.scorefile, chrom=None)  # read all variants

        logs = pl.scan_ipc(args.logs)  # lazily read ipc to preserve dtypes
        valid_matches, filter_summary = filter_scores(scorefile=scorefile, matches=logs, dataset="test",
                                                      min_overlap=args.min_overlap)

        if valid_matches.fetch().is_empty():  # this can happen if args.min_overlap = 0
            logger.error("Error: no target variants match any variants in scoring files")
            raise Exception

        dataset = args.dataset.replace('_', '-')
        summary_log: pl.LazyFrame = make_summary_log(match_candidates=logs, filter_summary=filter_summary, dataset=dataset,
                                                     scorefile=scorefile)
        check_log_count(summary_log=summary_log, scorefile=scorefile)

        write_log(df=logs, prefix=dataset, chrom=None, outdir=args.outdir, file_format="csv")
        summary_log.collect().write_csv(f"{dataset}_summary.csv")
        write_out(valid_matches, args.split, args.outdir, dataset)


def _parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dataset', dest='dataset', required=True,
                        help='<Required> Label for target genomic dataset')
    parser.add_argument('-m', '--min_overlap', dest='min_overlap', required=True,
                        type=float, help='<Required> Minimum proportion of variants to match before error')
    parser.add_argument('-s', '--scorefile', dest='scorefile', required=True,
                        help='<Required> Path to scorefile')
    parser.add_argument('--split', dest='split', default=True, action='store_true',
                        help='<Optional> Split scorefile per chromosome?')
    parser.add_argument('-l', '--logs', dest='logs', required=True,
                        help='<Required> Glob of log files including quotation marks e.g. "*.ipc"')
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help='<Required> Output directory')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    parser.add_argument('-n', dest='n_threads', default=1, help='<Optional> n threads for matching', type=int)
    return parser.parse_args(args)


if __name__ == "__main__":
    aggregate_matches()
