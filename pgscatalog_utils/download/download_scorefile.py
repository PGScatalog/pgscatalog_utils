import argparse
import logging
import os
import shutil
import textwrap
import time
import hashlib
from contextlib import closing
from functools import reduce
from urllib import request as request
from urllib.error import HTTPError, URLError

from pgscatalog_utils.download.publication import query_publication
from pgscatalog_utils.download.score import get_url
from pgscatalog_utils.download.trait import query_trait
from pgscatalog_utils.config import set_logging_level

logger = logging.getLogger(__name__)

md5_sep = '  '
md5_ext = '.md5'

def download_scorefile() -> None:
    args = _parse_args()
    set_logging_level(args.verbose)
    _check_args(args)
    _mkdir(args.outdir)

    if args.build is None:
        logger.warning(f'Downloading scoring file(s) in the author-reported genome build')
    elif args.build in ['GRCh37', 'GRCh38']:
        logger.warning(f'Downloading harmonized scoring file(s) in build: {args.build}.')
    else:
        logger.critical(f'Invalid genome build specified: {args.build}. Only -b GRCh37 and -b GRCh38 are supported')
        raise Exception

    overwrite_existing_file = args.overwrite_existing_file
    if overwrite_existing_file:
        logger.debug("--overwrite, overwriting new version of the Scoring file, if available")
        logger.warning(f'Existing Scoring files will be overwritten if new versions of the Scoring files are available for download.')

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
        # Download Scoring file
        is_downloaded = _download_ftp(url, path, overwrite_existing_file)
        # Check the download
        if is_downloaded:
            # - Generate MD5 checksum for the downloaded Scoring File
            downloaded_file_md5 = _generate_md5_checksum(path)
            # - Download MD5 checksum from the FTP and compare it with the generated MD5 checksum for the downloaded Scoring File
            ftp_md5 = _get_md5_checksum_from_ftp(url)
            # - Compare MD5 checksums
            if downloaded_file_md5 != ftp_md5:
                raise RuntimeError("The download of the file wasn't done properly: the generated MD5 Checksums from the downloaded file differs with the MD5 Checksums on the PGS Catalog FTP")


def _mkdir(outdir: str) -> None:
    if not os.path.exists(outdir):
        logger.debug("Creating output directory")
        os.makedirs(outdir)


def _download_ftp(url: str, path: str, overwrite_existing_file: bool, retry: int = 0) -> None:
    """ Download the Scoring file via the PGS Catalog FTP """
    # Check if Scoring file already exist in the local directory
    if os.path.exists(path) and retry == 0:
        existing_file_md5 = _generate_md5_checksum(path)
        ftp_md5 = _get_md5_checksum_from_ftp(url)
        # Overwrite option for the Scoring file
        if overwrite_existing_file:
            if existing_file_md5 == ftp_md5:
                logger.warning(f"Identical Scoring file already exists at {path}, skipping download.")
                return False
            else:
                logger.warning(f"Existing Scoring file at {path}, will be overwritten by a new version of the file on the FTP.")
        # No overwrite option
        else:
            warning_msg = f"File already exists at {path}, skipping download."
            if existing_file_md5 != ftp_md5:
                warning_msg = warning_msg + " However a new version of the Scoring file is available on the FTP. You can overwrite the existing Scoring file using the parameter '--overwrite'."
            logger.warning(warning_msg)
            return False
    # Attempt to download the Scoring file
    try:
        with closing(request.urlopen(url)) as r:
            with open(path, 'wb') as f:
                shutil.copyfileobj(r, f)
        return True
    except (HTTPError, URLError) as error:
        max_retries = 5
        logger.warning(f'Download failed: {error.reason}')
        # Retry to download the file if the server is busy
        if '421' in error.reason and retry < max_retries:
            logger.warning(f'> Retry to download the file ... attempt {retry+1} out of {max_retries}.')
            retry += 1
            time.sleep(10)
            is_downloaded = _download_ftp(url,path,overwrite_existing_file,retry)
            if is_downloaded:
                return is_downloaded
        else:
            raise RuntimeError("Failed to download '{}'.\nError message: '{}'".format(url, error.reason))


def _generate_md5_checksum(filename,blocksize=4096):
    """ Returns MD5 checksum for the given file. """
    md5 = hashlib.md5()
    try:
        file = open(filename, 'rb')
        with file:
            for block in iter(lambda: file.read(blocksize), b""):
                md5.update(block)
    except IOError:
        logger.warning('File \'' + filename + '\' not found!')
        return None
    except:
        logger.warning("Error: the script couldn't generate a MD5 checksum for '" + filename + "'!")
        return None
    return md5.hexdigest()


def _get_md5_checksum_from_ftp(url:str, retry:int = 0) -> str:
    """ Download and extract MD5 checksum value from the FTP. """
    md5 = None
    try:
        url_md5 = url+md5_ext
        with closing(request.urlopen(url_md5)) as file:
            md5_content = file.read().decode()
            md5 = md5_content.split(md5_sep)[0]
    except (HTTPError, URLError) as error:
        max_retries = 5
        logger.warning(f'MD5 checksum download failed: {error.reason}')
        # Retry to download the MD% checksum file if the server is busy
        if '421' in error.reason and retry < max_retries:
            logger.warning(f'> Retry to download the file ... attempt {retry+1} out of {max_retries}.')
            retry += 1
            time.sleep(10)
            _get_md5_checksum_from_ftp(url,retry)
        else:
            raise RuntimeError(f"MD5 file download failed - Error message: '{error.reason}'")
    return md5


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
    parser.add_argument('-w', '--overwrite', dest='overwrite_existing_file', action='store_true',
                        help='<Optional> Overwrite existing Scoring File if a new version is available for download on the FTP')
    parser.add_argument('-c', '--pgsc_calc', dest='pgsc_calc',
                        help='<Optional> Provide information about downloading scoring files via pgsc_calc')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    download_scorefile()
