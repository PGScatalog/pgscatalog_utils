import argparse
import logging
import os
import sys
import textwrap
import json

from pgscatalog_utils.config import set_logging_level
from pgscatalog_utils.scorefile.effect_type import set_effect_type
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.genome_build import build2GRC
from pgscatalog_utils.scorefile.harmonised import remap_harmonised
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.qc import quality_control
from pgscatalog_utils.scorefile.read import load_scorefile, get_scorefile_basename
from pgscatalog_utils.scorefile.write import write_scorefile


headers2logs = [
    'pgs_id',
    'pgp_id',
    'pgs_name',
    'genome_build',
    'variants_number',
    'trait_reported',
    'trait_efo',
    'trait_mapped',
    'weight_type',
    'citation'
]
headers2logs_harmonisation = [
    'HmPOS_build',
    'HmPOS_date',
    'HmPOS_match_chr',
    'HmPOS_match_pos'
]

def combine_scorefiles():
    args = _parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")

    if os.path.exists(args.outfile):
        logger.critical(f"Output file {args.outfile} already exists")
        raise Exception

    # Score header logs - init
    score_logs = {}
    dir_output = os.path.dirname(args.outfile)
    if dir_output == '':
        dir_output = './'
    elif dir_output.endswith('/') is False:
        dir_output += '/'
    json_logs_file =  dir_output + args.logfile

    for x in paths:
        # Read scorefile df and header
        h, score = load_scorefile(x)
        score_shape_original = score.shape

        if score.empty:
            logger.critical(f"Empty scorefile {x} detected! Please check the input data")
            raise Exception

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

        if args.liftover:
            logger.debug("Annotating scorefile with liftover parameters")
            score = liftover(score, args.chain_dir, args.min_lift, args.target_build)

        if score.empty and (args.drop_missing is False):
            logger.critical("Empty output score detected, something went wrong while combining")
            raise Exception

        write_scorefile(score, args.outfile)

        # Build Score header logs
        score_id = get_scorefile_basename(x)
        score_header = score_logs[score_id] = {}
        # Scoring file header information
        for header in headers2logs:
            header_val = h.get(header)
            if (header in ['trait_efo', 'trait_mapped']) and (header_val is not None):
                header_val = header_val.split('|')
            score_header[header] = header_val
        # Other header information
        score_header['columns'] = list(score.columns)
        score_header['use_liftover'] = False
        if args.liftover:
             score_header['use_liftover'] = True
        # Harmonized header information
        score_header['use_harmonised'] = use_harmonised
        if use_harmonised:
            score_header['sources'] = sorted(score['hm_source'].unique().tolist())
            for hm_header in headers2logs_harmonisation:
                hm_header_val = h.get(hm_header)
                if hm_header_val:
                    if hm_header.startswith('HmPOS_match'):
                        hm_header_val = json.loads(hm_header_val)
                    score_header[hm_header] = hm_header_val
        if score_header['variants_number'] is None:
            score_header['variants_number'] = score_shape_original[0]

    # Write Score header logs file
    with open(json_logs_file, 'w') as fp:
        json.dump(score_logs, fp, indent=4)


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
    parser.add_argument('-l', '--logfile', dest='logfile', default='log_combined.json',
                        help='<Required> Name for the log file (score metadata) for combined scores.'
                             '[ will write to identical directory as combined scorefile]')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    combine_scorefiles()
