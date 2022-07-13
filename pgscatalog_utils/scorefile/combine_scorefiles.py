import argparse
import sys
import logging
import pandas as pd

from pgscatalog_utils.log_config import set_logging_level
from pgscatalog_utils.scorefile.read import load_scorefile
from pgscatalog_utils.scorefile.effect_type import set_effect_type
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.write import write_scorefile


def parse_args(args=None) -> argparse.Namespace:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Combine multiple scoring files')
    parser.add_argument('-s', '--scorefiles', dest='scorefiles', nargs='+',
                        help='<Required> Scorefile path (wildcard * is OK)', required=True)
    parser.add_argument('--liftover', dest='liftover',
                        help='<Optional> Convert scoring file variants to target genome build?', action='store_true')
    parser.add_argument('-t', '--target_build', dest='target_build', help='Build of target genome <GRCh37 / GRCh38>',
                        required='--liftover' in sys.argv)
    parser.add_argument('-c', '--chain_dir', dest='chain_dir', help='Path to directory containing chain files',
                        required="--liftover" in sys.argv)
    parser.add_argument('-m', '--min_lift', dest='min_lift',
                        help='If liftover, minimum proportion of variants lifted over',
                        required="--liftover" in sys.argv, default=0.95, type=float)
    parser.add_argument('-o', '--outfile', dest='outfile', required=True,
                        default='combined.txt',
                        help='<Required> Output path to combined long scorefile')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


def combine_scorefiles():
    args = parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")
    scorefiles: pd.DataFrame = pd.concat([_read_and_melt(x) for x in paths])

    if args.liftover:
        logger.debug("Annotating scorefiles with liftover parameters")
        scorefiles['target_build'] = args.target_build
        scorefiles = liftover(scorefiles, args.chain_dir, args.min_lift)

    write_scorefile(scorefiles, args.outfile)


def _read_and_melt(path):
    """ Load a scorefile, melt it, and set the effect types"""
    return (load_scorefile(path)
            .pipe(melt_effect_weights)
            .pipe(set_effect_type))


if __name__ == "__main__":
    combine_scorefiles()
