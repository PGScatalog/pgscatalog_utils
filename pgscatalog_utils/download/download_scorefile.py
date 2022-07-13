import logging
import argparse
import os
import shutil
from contextlib import closing
from urllib import request as request
from pgscatalog_utils.download.api import pgscatalog_result
from pgscatalog_utils.log_config import set_logging_level

logger = logging.getLogger(__name__)


def parse_args(args=None) -> argparse.Namespace:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='Download scoring files')
    parser.add_argument('-i', '--id', nargs='+', dest='pgs',
                        help='<Required> PGS Catalog ID', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
                        default='scores/',
                        help='<Required> Output directory to store downloaded files')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


def download_scorefile() -> None:
    args = parse_args()

    set_logging_level(args.verbose)

    _mkdir(args.outdir)

    urls: dict[str, str] = pgscatalog_result(args.pgs)

    for pgsid, url in urls.items():
        logger.debug(f"Downloading {pgsid} from {url}")
        path: str = os.path.join(args.outdir, pgsid + '.txt.gz')
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


if __name__ == "__main__":
    download_scorefile()
