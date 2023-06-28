import logging
import os
import pathlib
import time
import urllib.parse
from ftplib import FTP
from urllib.parse import urlsplit

import requests

from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def download_file(url: str, local_path: str, overwrite: bool, ftp_fallback: bool) -> None:
    if config.OUTDIR.joinpath(local_path).exists():
        if not overwrite:
            logger.warning(f"{config.OUTDIR.joinpath(local_path)} exists and overwrite is false, skipping download")
            return
        elif overwrite:
            logger.warning(f"Overwriting {config.OUTDIR.joinpath(local_path)}")

    logger.info(f"Downloading {local_path} from {url}")
    attempt: int = 0

    while attempt < config.MAX_RETRIES:
        response: requests.Response = requests.get(url)
        match response.status_code:
            case 200:
                with open(config.OUTDIR.joinpath(local_path), "wb") as f:
                    f.write(response.content)
                logger.info("HTTPS download complete")
                attempt = 0
                break
            case _:
                logger.warning(f"HTTP status {response.status_code} at download attempt {attempt}")
                attempt += 1
                time.sleep(5)

    if attempt > config.MAX_RETRIES:
        if ftp_fallback:
            logger.warning("Attempting FTP fallback")
            _ftp_fallback_download(url=url, local_path=local_path)
        else:
            raise Exception(f"Can't download {url} using HTTPS")


def _ftp_fallback_download(url: str, local_path: str) -> None:
    url = url.replace("https://", "ftp://")
    retries = 0

    while retries < config.MAX_RETRIES:
        try:
            spliturl: urllib.parse.SplitResult = urlsplit(url)
            ftp = FTP(spliturl.hostname)
            ftp.login()
            ftp.cwd(str(pathlib.Path(urlsplit(url).path).parent))
            with open(config.OUTDIR.joinpath(local_path), "wb") as file:
                ftp.retrbinary("RETR " + local_path, file.write)
                ftp.quit()
                logger.info("FTP download completed")
                return
        except Exception as e:
            if "421" in str(e):
                retries += 1
                logger.debug(f"FTP server is busy. Waiting and retrying. Retry {retries} of {config.MAX_RETRIES}")
                time.sleep(config.DOWNLOAD_WAIT_TIME)
            else:
                logger.critical(f"Download failed: {e}")
                raise Exception
