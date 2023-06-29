import argparse
import logging
import os
import pathlib
import textwrap
import typing

from pgscatalog_utils import config
from pgscatalog_utils.download.CatalogCategory import CatalogCategory
from pgscatalog_utils.download.Catalog import CatalogQuery, CatalogResult
from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.download.ScoringFileDownloader import ScoringFileDownloader


logger = logging.getLogger(__name__)


def download_scorefile() -> None:
    args = _parse_args()
    config.set_logging_level(args.verbose)
    _check_args(args)

    build: typing.Union[None, GenomeBuild]
    match args.build:
        case None:
            logger.info("Downloading scoring file(s) in the author-reported genome build")
            build = None
        case "GRCh37":
            build = GenomeBuild.GRCh37
        case "GRCh38":
            build = GenomeBuild.GRCh38
        case _:
            logger.critical(f"Invalid genome build specified: {args.build}")
            logger.critical("Only -b GRCh37 and -b GRCh38 are supported")
            raise Exception

    if build is not None:
        logger.info(f"Downloading harmonised scoring file(s) in build: {build}")

    if args.overwrite_existing_file:
        logger.debug("--overwrite, overwriting new version of the Scoring file, if available")
        logger.warning(
            "Existing Scoring files will be overwritten if new versions of the Scoring files are available for download.")

    if args.pgsc_calc:
        config.PGSC_CALC_VERSION = args.pgsc_calc
        logger.info(f"Setting user agent to {config.PGSC_CALC_VERSION} for PGS Catalog API queries")

    config.OUTDIR = pathlib.Path(args.outdir).resolve()
    logger.info(f"Download directory: {config.OUTDIR}")
    config.OUTDIR.mkdir(exist_ok=True)
    config.OVERWRITE = args.overwrite_existing_file

    results: list[list[CatalogResult]] = []
    if args.efo:
        inc_child = False
        if args.efo_include_children:
            logger.debug("--trait set, querying traits (including PGS for child terms)")
            inc_child = True
        else:
            logger.debug("--trait set, querying traits")
        for term in args.efo:
            results.append(CatalogQuery(CatalogCategory.TRAIT, term, include_children=inc_child,
                                    pgsc_calc_version=config.PGSC_CALC_VERSION).get())

    if args.pgp:
        logger.debug("--pgp set, querying publications")
        for term in args.pgp:
            results.append(CatalogQuery(CatalogCategory.PUBLICATION, term, pgsc_calc_version=config.PGSC_CALC_VERSION).get())

    if args.pgs:
        logger.debug("--id set, querying scores")
        results.append(
            CatalogQuery(CatalogCategory.SCORE, args.pgs,
                         pgsc_calc_version=config.PGSC_CALC_VERSION).get())  # pgs_lst: a list containing up to three flat lists

    flat_results = [element for sublist in results for element in sublist]

    ScoringFileDownloader(results=flat_results, genome_build=build, overwrite=config.OVERWRITE).download_files()

    # warn if missing PGS IDs in downloaded files
    requested_pgs: set[str]
    if args.pgs is None:
        requested_pgs = set()
    else:
        requested_pgs = set(args.pgs)
    missing_pgs: set[str] = requested_pgs.difference(flat_results[-1].pgs_ids)

    if missing_pgs:
        logger.warning(f"Requested PGS scoring file not downloaded: {missing_pgs}")
        logger.warning("Check if the accessions are valid")

    logger.info("Downloads complete")


def _check_args(args):
    if not args.efo:
        if not args.pgp:
            if not args.pgs:
                logger.critical(
                    "One of --trait, --pgp, or --id is required to download scorefiles"
                )
                raise Exception


def _description_text() -> str:
    return textwrap.dedent(
        """\
    Download a set of scoring files from the PGS Catalog using PGS
    Scoring IDs, traits, or publication IDs.
    
    The PGS Catalog API is queried to get a list of scoring file
    URLs. Scoring files are downloaded via FTP to a specified
    directory. PGS Catalog scoring files are staged with the name:

            {PGS_ID}.txt.gz

    If a valid build is specified harmonized files are downloaded as:
    
        {PGS_ID}_hmPOS_{genome_build}.txt.gz
    
    These harmonised scoring files contain genomic coordinates,
    remapped from author-submitted information such as rsids.
   """
    )


def _epilog_text() -> str:
    return textwrap.dedent(
        """\
    download_scorefiles will skip downloading a scoring file if it
    already exists in the download directory. This can be useful if
    the download process is interrupted and needs to be restarted
    later. You can track download progress with the verbose flag.    
   """
    )


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=_description_text(),
        epilog=_epilog_text(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i", "--pgs", nargs="+", dest="pgs", help="PGS Catalog ID(s) (e.g. PGS000001)"
    )
    parser.add_argument(
        "-t",
        "--efo",
        dest="efo",
        nargs="+",
        help="Traits described by an EFO term(s) (e.g. EFO_0004611)",
    )
    parser.add_argument(
        "-e",
        "--efo_direct",
        dest="efo_include_children",
        action="store_false",
        help="<Optional> Return only PGS tagged with exact EFO term "
             "(e.g. no PGS for child/descendant terms in the ontology)",
    )
    parser.add_argument(
        "-p",
        "--pgp",
        dest="pgp",
        help="PGP publication ID(s) (e.g. PGP000007)",
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--build",
        dest="build",
        choices=["GRCh37", "GRCh38"],
        help="Download Harmonized Scores with Positions in Genome build: GRCh37 or GRCh38",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        required=True,
        default="scores/",
        help="<Required> Output directory to store downloaded files",
    )
    parser.add_argument(
        "-w",
        "--overwrite",
        dest="overwrite_existing_file",
        action="store_true",
        help="<Optional> Overwrite existing Scoring File if a new version is available for download on the FTP",
    )
    parser.add_argument(
        "-c",
        "--pgsc_calc",
        dest="pgsc_calc",
        help="<Optional> Provide information about downloading scoring files via pgsc_calc",
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
    download_scorefile()
