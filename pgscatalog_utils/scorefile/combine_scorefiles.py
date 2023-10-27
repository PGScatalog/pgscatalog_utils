import argparse
import json
import logging
import logging.handlers
import multiprocessing
import pathlib
import sys
import textwrap

from pgscatalog_utils.scorefile.combine import combine
from pgscatalog_utils.scorefile.log import listener_process, listener_configurer, \
    log_format


def combine_scorefiles():
    args = _parse_args()

    # set up centralised logging for workers
    m = multiprocessing.Manager()
    queue = m.Queue(maxsize=-1)  # infinite queue size
    listener = multiprocessing.Process(target=listener_process,
                                       args=(queue, listener_configurer))
    listener.start()

    # set up logging for main process
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.DEBUG,
                        format=log_format())
    logger.debug("Starting combine scorefiles")

    if pathlib.Path(args.outfile).exists():
        logger.critical(f"Output file {args.outfile} already exists")
        raise FileNotFoundError

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")

    score_logs = []
    logger.debug(f"Setting up multiprocessing pool with {args.threads} workers")
    with multiprocessing.Pool(processes=args.threads) as pool:
        lock = multiprocessing.Manager().Lock()
        args_list = [(path, args, lock, queue) for path in paths]

        score_logs.append(pool.starmap(combine, args_list))

    with open(pathlib.Path(args.outfile).parent / "log_combined.json", "w") as f:
        json.dump(score_logs, f, indent=4)

    logger.debug("Finished :)")
    queue.put_nowait(None)
    listener.join()


def _description_text() -> str:
    return textwrap.dedent('''\
    Combine multiple scoring files in PGS Catalog format (see https://www.pgscatalog.org/downloads/ 
    for details) to a 'long' table of columns needed for variant matching and subsequent calculation. 

    Custom scorefiles in PGS Catalog format can be combined with PGS Catalog scoring files, and 
    optionally liftover genomic coordinates to GRCh37 or GRCh38. The script can accept a mix of
    unharmonised and harmonised PGS Catalog data. By default all variants are output (including 
    positions with duplicated data [often caused by rsID/liftover collions across builds]) and 
    variants with missing positions. 
    ''')


def _epilog_text() -> str:
    return textwrap.dedent('''\
    The long table is used to simplify intersecting variants in target genotyping datasets 
    and the scoring files with the match_variants program.
    ''')


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(),
                                     epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--scorefiles', dest='scorefiles', nargs='+',
                        help='<Required> Scorefile path (wildcard * is OK)',
                        required=True)
    parser.add_argument('--threads', dest='threads',
                        help='Number of threads to use',
                        default=1, type=int)
    parser.add_argument('--liftover', dest='liftover',
                        help='<Optional> Convert scoring file variants to target genome build?',
                        action='store_true')
    parser.add_argument('-t', '--target_build', dest='target_build',
                        choices=['GRCh37', 'GRCh38'],
                        help='<Required> Build of target genome',
                        required=True)
    parser.add_argument('-c', '--chain_dir', dest='chain_dir',
                        help='Path to directory containing chain files',
                        required="--liftover" in sys.argv)
    parser.add_argument('-m', '--min_lift', dest='min_lift',
                        help='<Optional> If liftover, minimum proportion of variants lifted over',
                        required="--liftover" in sys.argv, default=0.95, type=float)
    parser.add_argument('--drop_missing', dest='drop_missing', action='store_true',
                        help='<Optional> Drop variants with missing information (chr/pos) and '
                             'non-standard alleles (e.g. HLA=P/N) from the output file.')
    parser.add_argument('-o', '--outfile', dest='outfile', required=True,
                        default='combined.txt',
                        help='<Required> Output path to combined long scorefile '
                             '[ will compress output if filename ends with .gz ]')
    parser.add_argument('-l', '--logfile', dest='logfile', default='log_combined.json',
                        help='<Required> Name for the log file (score metadata) for combined scores.'
                             '[ will write to identical directory as combined scorefile]')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    combine_scorefiles()
