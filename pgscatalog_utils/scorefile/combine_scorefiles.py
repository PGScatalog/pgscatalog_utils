import argparse
import logging
import sys
import textwrap

import pandas as pd

from pgscatalog_utils.log_config import set_logging_level
from pgscatalog_utils.scorefile.read import load_scorefile
from pgscatalog_utils.scorefile.effect_type import set_effect_type
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.write import write_scorefile


def combine_scorefiles():
    args = _parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")
    scorefiles: pd.DataFrame = pd.concat([_read_and_melt(x, drop_missing=args.drop_missing) for x in paths])

    if args.liftover:
        logger.debug("Annotating scorefiles with liftover parameters")
        scorefiles = liftover(scorefiles, args.chain_dir, args.min_lift, args.target_build)

    write_scorefile(scorefiles, args.outfile)


def _read_and_melt(path, drop_missing: bool = False):
    """ Load a scorefile, melt it, and set the effect types"""
    return (load_scorefile(path, drop_missing=drop_missing)
            .pipe(melt_effect_weights)
            .pipe(set_effect_type))


if __name__ == "__main__":
    combine_scorefiles()


def _description_text() -> str:
    return textwrap.dedent('''\
    Combine multiple scoring files in PGS Catalog format (see 
    https://www.pgscatalog.org/downloads/ for details) to a 'long'
    table, and optionally liftover genomic coordinates to GRCh37 or
    GRCh38. Custom scorefiles in PGS Catalog format can be combined
    with PGS Catalog scoring files. The program can accept a mix of
    unharmonised and harmonised PGS Catalog data.     
    ''')


def _epilog_text() -> str:
    return textwrap.dedent('''\
    The long table is used to simplify intersecting variants in target
    genomes and the scoring files with the match_variants program.    
    ''')


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
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
    parser.add_argument('--drop_missing', dest='drop_missing', action='store_true',
                        help='Drop variants with missing information (chr/pos) and '
                             'non-standard alleles from the output file.')
    parser.add_argument('-o', '--outfile', dest='outfile', required=True,
                        default='combined.txt',
                        help='<Required> Output path to combined long scorefile')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)
