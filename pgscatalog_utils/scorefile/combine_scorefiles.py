import argparse
import logging
import sys
import textwrap

import pandas as pd

from pgscatalog_utils.log_config import set_logging_level
from pgscatalog_utils.scorefile.effect_type import set_effect_type
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.genome_build import build2GRC
from pgscatalog_utils.scorefile.harmonised import remap_harmonised
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.qc import quality_control
from pgscatalog_utils.scorefile.read import load_scorefile
from pgscatalog_utils.scorefile.write import write_scorefile


def combine_scorefiles():
    args = _parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")

    scorefiles = []
    for x in paths:
        # Read scorefile df and header
        h, score = load_scorefile(x)

        # Check if we should use the harmonized positions
        use_harmonised = False
        current_build = None
        if h.get('HmPOS_build') is not None:
            if h.get('HmPOS_build') == args.target_build:
                use_harmonised = True
                current_build = h.get('HmPOS_build')
            else:
                logger.error(
                    f"Cannot combine {x} (harmonized to {h.get('HmPOS_build')}) in target build {args.target_build}")
                raise Exception

        # Process/QC score and check variant columns
        score = (score.pipe(remap_harmonised, use_harmonised=use_harmonised)
                 .pipe(quality_control, drop_missing=args.drop_missing)
                 .pipe(melt_effect_weights)
                 .pipe(set_effect_type))

        # Annotate score with the genome_build (in GRCh notation)
        if current_build is None:
            current_build = build2GRC(h.get('genome_build'))
            if current_build is None:
                logger.error("Scorefile has no build information, "
                             "please add the build to the header with "
                             "('#genome_build=[insert variant build]")
                raise Exception

        score = score.assign(genome_build=current_build)

        if (current_build != args.target_build) and (args.liftover is False):
            logger.error(
                f"Cannot combine {x} (build={h.get('genome_build')}) with target build {args.target_build} without liftover")
            logger.error("Try running with --liftover and specifying the --chain_dir")
            raise Exception

        scorefiles.append(score)

    if len(scorefiles) > 0:
        scorefiles: pd.DataFrame = pd.concat(scorefiles)
    else:
        logger.error("No valid scorefiles could be combined")
        raise Exception

    if args.liftover:
        logger.debug("Annotating scorefiles with liftover parameters")
        scorefiles = liftover(scorefiles, args.chain_dir, args.min_lift, args.target_build)

    write_scorefile(scorefiles, args.outfile)


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
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--scorefiles', dest='scorefiles', nargs='+',
                        help='<Required> Scorefile path (wildcard * is OK)', required=True)
    parser.add_argument('--liftover', dest='liftover',
                        help='<Optional> Convert scoring file variants to target genome build?', action='store_true')
    parser.add_argument('-t', '--target_build', dest='target_build',
                        choices=['GRCh37', 'GRCh38'], help='<Required> Build of target genome',
                        required=True)
    parser.add_argument('-c', '--chain_dir', dest='chain_dir', help='Path to directory containing chain files',
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
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    combine_scorefiles()
