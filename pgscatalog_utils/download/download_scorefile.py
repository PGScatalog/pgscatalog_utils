import argparse
import logging
import os
import shutil
import textwrap
import time
from contextlib import closing
from functools import reduce
from urllib import request as request
from urllib.error import HTTPError, URLError

from pgscatalog_utils.download.publication import query_publication
from pgscatalog_utils.download.score import get_url
from pgscatalog_utils.download.trait import query_trait
from pgscatalog_utils.config import set_logging_level

logger = logging.getLogger(__name__)


def download_scorefile() -> None:
    args = _parse_args()
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

    pgsc_calc_info = None
    if args.pgsc_calc:
        pgsc_calc_info = args.pgsc_calc

    if args.efo:
        if args.efo_include_children:
            logger.debug("--trait set, querying traits (including PGS for child terms)")
        else:
            logger.debug("--trait set, querying traits")
        pgs_lst = pgs_lst + [query_trait(x, pgsc_calc_info, args.efo_include_children) for x in args.efo]


    if args.pgp:
        logger.debug("--pgp set, querying publications")
        pgs_lst = pgs_lst + [query_publication(x, pgsc_calc_info) for x in args.pgp]

    if args.pgs:
        logger.debug("--id set, querying scores")
        pgs_lst.append(args.pgs)  # pgs_lst: a list containing up to three flat lists

    pgs_id: list[str] = list(set(reduce(lambda x, y: x + y, pgs_lst)))

    urls: dict[str, str] = get_url(pgs_id, args.build, pgsc_calc_info)

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


def _download_ftp(url: str, path: str, retry:int = 0) -> None:
    if os.path.exists(path):
        logger.warning(f"File already exists at {path}, skipping download")
        return
    else:
        try:
            with closing(request.urlopen(url)) as r:
                with open(path, 'wb') as f:
                    shutil.copyfileobj(r, f)
        except (HTTPError, URLError) as error:
            max_retries = 5
            print(f'Download failed: {error.reason}')
            # Retry to download the file if the server is busy
            if '421' in error.reason and retry < max_retries:
                print(f'> Retry to download the file ... attempt {retry+1} out of {max_retries}.')
                retry += 1
                time.sleep(10)
                _download_ftp(url,path,retry)
            else:
                raise RuntimeError("Failed to download '{}'.\nError message: '{}'".format(url, error.reason))


def _check_args(args):
    if not args.efo:
        if not args.pgp:
            if not args.pgs:
                logger.critical("One of --trait, --pgp, or --id is required to download scorefiles")
                raise Exception


def _description_text() -> str:
    return textwrap.dedent('''\
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
   ''')


def _epilog_text() -> str:
    return textwrap.dedent('''\
    download_scorefiles will skip downloading a scoring file if it
    already exists in the download directory. This can be useful if
    the download process is interrupted and needs to be restarted
    later. You can track download progress with the verbose flag.    
   ''')


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--pgs', nargs='+', dest='pgs', help='PGS Catalog ID(s) (e.g. PGS000001)')
    parser.add_argument('-t', '--efo', dest='efo', nargs='+',
                        help='Traits described by an EFO term(s) (e.g. EFO_0004611)')
    parser.add_argument('-e', '--efo_direct', dest='efo_include_children', action='store_false',
                        help='<Optional> Return only PGS tagged with exact EFO term '
                             '(e.g. no PGS for child/descendant terms in the ontology)')
    parser.add_argument('-p', '--pgp', dest='pgp', help='PGP publication ID(s) (e.g. PGP000007)', nargs='+')
    parser.add_argument('-b', '--build', dest='build', choices=['GRCh37', 'GRCh38'],
                        help='Download Harmonized Scores with Positions in Genome build: GRCh37 or GRCh38')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
                        default='scores/',
                        help='<Required> Output directory to store downloaded files')
    parser.add_argument('-c', '--pgsc_calc', dest='pgsc_calc',
                        help='<Optional> Provide information about downloading scoring files via pgsc_calc')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    download_scorefile()
