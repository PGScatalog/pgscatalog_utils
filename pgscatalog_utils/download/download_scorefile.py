import logging
import argparse
import os
import shutil
from contextlib import closing
from functools import reduce
from urllib import request as request
import sys

from pgscatalog_utils.download.publication import query_publication
from pgscatalog_utils.download.score import get_url
from pgscatalog_utils.download.trait import query_trait
from pgscatalog_utils.log_config import set_logging_level

logger = logging.getLogger(__name__)


def parse_args(args=None) -> argparse.Namespace:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Download scoring files')
    parser.add_argument('-i', '--pgs', nargs='+', dest='pgs', help='PGS Catalog ID(s) (e.g. PGS000001)')
    parser.add_argument('-t', '--efo', dest='efo', nargs='+',
                        help='Traits described by an EFO term(s) (e.g. EFO_0004611)')
    parser.add_argument('-p', '--pgp', dest='pgp', help='PGP publication ID(s) (e.g. PGP000007)', nargs='+')
    parser.add_argument('-b', '--build', dest='build',
                        help='Download Harmonized Scores with Positions in Genome build: GRCh37 or GRCh38')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
                        default='scores/',
                        help='<Required> Output directory to store downloaded files')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='Extra logging information')
    return parser.parse_args(args)


def download_scorefile() -> None:
    args = parse_args()
    set_logging_level(args.verbose)
    _check_args(args)
    _mkdir(args.outdir)

    if args.build is None:
        logger.critical(f'Downloading scoring file(s) in the author-reported genome build')
    elif args.build in ['GRCh37', 'GRCh38']:
        logger.critical(f'Downloading harmonized scoring file(s) in build: {args.build}.')
    else:
        logger.critical(f'Invalid genome build specified: {args.build}. Only -b GRCh37 and -b GRCh38 are supported')
        raise Exception

    pgs_lst: list[list[str]] = []

    if args.efo:
        logger.debug("--trait set, querying traits")
        pgs_lst = pgs_lst + [query_trait(x) for x in args.efo]

    if args.pgp:
        logger.debug("--pgp set, querying publications")
        pgs_lst = pgs_lst + [query_publication(x) for x in args.pgp]

    if args.pgs:
        logger.debug("--id set, querying scores")
        pgs_lst.append(args.pgs)  # pgs_lst: a list containing up to three flat lists

    pgs_id: list[str] = list(set(reduce(lambda x, y: x + y, pgs_lst)))

    urls: dict[str, str] = get_url(pgs_id, args.build)

    for pgsid, url in urls.items():
        logger.debug(f"Downloading {pgsid} from {url}")
        if args.build is None:
            path: str = os.path.join(args.outdir, pgsid + '.txt.gz')
        else:
            path: str = os.path.join(args.outdir, pgsid + f'_hmPOS_{args.build}.txt.gz')
        _download_ftp(url, path)


def _mkdir(outdir: str) -> None:
    if not os.path.exists(outdir):
        logger.debug("Creating output directory")
        os.makedirs(outdir)


def _download_ftp(url: str, path: str) -> None:
    if os.path.exists(path):
        logger.warning(f"File already exists at {path}, skipping download")
        return
    else:
        with closing(request.urlopen(url)) as r:
            with open(path, 'wb') as f:
                shutil.copyfileobj(r, f)


def _check_args(args):
    if not args.efo:
        if not args.pgp:
            if not args.pgs:
                logger.critical("One of --trait, --pgp, or --id is required to download scorefiles")
                raise Exception


if __name__ == "__main__":
    download_scorefile()
