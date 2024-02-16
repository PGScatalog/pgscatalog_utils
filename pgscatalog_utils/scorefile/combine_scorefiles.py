import argparse
import json
import logging
import pathlib
import sys
import textwrap

from pgscatalog_utils.config import set_logging_level
from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.liftover import create_liftover
from pgscatalog_utils.scorefile.scoringfile import ScoringFile
from pgscatalog_utils.scorefile.write import write_combined


def combine_scorefiles():
    args = _parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    Config.batch_size = 100000
    Config.drop_missing = args.drop_missing
    Config.target_build = GenomeBuild.from_string(args.target_build)
    Config.liftover = args.liftover
    Config.min_lift = args.min_lift

    if args.chain_dir:
        Config.chain_dir = args.chain_dir
        Config.lo = create_liftover()

    if pathlib.Path(args.outfile).exists():
        raise FileExistsError(f"{args.outfile}")

    paths: list[str] = list(set(args.scorefiles))  # unique paths only
    logger.debug(f"Input scorefiles: {paths}")

    sfs = [ScoringFile.from_path(x) for x in paths]

    target_build = GenomeBuild.from_string(args.target_build)
    bad_builds = [x.accession for x in sfs if x.genome_build != target_build]

    if not args.liftover:
        for bad_file in bad_builds:
            logger.critical(f"{bad_file} doesn't match {target_build}, can't combine")
        if len(bad_builds) > 0:
            raise Exception

    # provide line counts when making the scoring files
    logs: dict[str, int] = write_combined(sfs, args.outfile)
    json_log = []
    for (k, v), sf in zip(logs.items(), sfs):
        json_log.append(sf.generate_log(v))

    log_out_path = pathlib.Path(args.outfile).parent / args.logfile
    with open(log_out_path, "w") as f:
        logger.info(f"Writing log to {f.name}")
        json.dump(json_log, f, indent=4)


def _description_text() -> str:
    return textwrap.dedent(
        """\
    Combine multiple scoring files in PGS Catalog format (see https://www.pgscatalog.org/downloads/ 
    for details) to a 'long' table of columns needed for variant matching and subsequent calculation. 

    Custom scorefiles in PGS Catalog format can be combined with PGS Catalog scoring files, and 
    optionally liftover genomic coordinates to GRCh37 or GRCh38. The script can accept a mix of
    unharmonised and harmonised PGS Catalog data. By default all variants are output (including 
    positions with duplicated data [often caused by rsID/liftover collions across builds]) and 
    variants with missing positions. 
    """
    )


def _epilog_text() -> str:
    return textwrap.dedent(
        """\
    The long table is used to simplify intersecting variants in target genotyping datasets 
    and the scoring files with the match_variants program.
    """
    )


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=_description_text(),
        epilog=_epilog_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--scorefiles",
        dest="scorefiles",
        nargs="+",
        help="<Required> Scorefile path (wildcard * is OK)",
        required=True,
    )
    parser.add_argument(
        "--liftover",
        dest="liftover",
        help="<Optional> Convert scoring file variants to target genome build?",
        action="store_true",
    )
    parser.add_argument(
        "-t",
        "--target_build",
        dest="target_build",
        choices=["GRCh37", "GRCh38"],
        help="<Required> Build of target genome",
        required=True,
    )
    parser.add_argument(
        "-c",
        "--chain_dir",
        dest="chain_dir",
        help="Path to directory containing chain files",
        required="--liftover" in sys.argv,
    )
    parser.add_argument(
        "-m",
        "--min_lift",
        dest="min_lift",
        help="<Optional> If liftover, minimum proportion of variants lifted over",
        default=0.95,
        type=float,
    )
    parser.add_argument(
        "--drop_missing",
        dest="drop_missing",
        action="store_true",
        help="<Optional> Drop variants with missing information (chr/pos) and "
        "non-standard alleles (e.g. HLA=P/N) from the output file.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        required=True,
        default="combined.txt",
        help="<Required> Output path to combined long scorefile "
        "[ will compress output if filename ends with .gz ]",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        default="log_combined.json",
        help="<Required> Name for the log file (score metadata) for combined scores."
        "[ will write to identical directory as combined scorefile]",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="<Optional> Extra logging information",
    )
    return parser.parse_args(args)


if __name__ == "__main__":
    combine_scorefiles()
