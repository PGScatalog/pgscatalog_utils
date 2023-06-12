import logging
import pathlib
import time
import urllib.parse
from ftplib import FTP
from urllib.parse import urlsplit

import requests

logger = logging.getLogger(__name__)


def download_file(url: str, local_path: str, overwrite: bool, ftp_fallback: bool) -> None:
    if pathlib.Path(local_path).exists():
        if not overwrite:
            logger.warning("File exists and overwrite is false, skipping download")
            return
        else:
            logger.warning("Overwriting existing scoring file")

    logger.info(f"Downloading {local_path} from {url}")
    max_retries: int = 5
    attempt: int = 0

    while attempt < max_retries:
        response: requests.Response = requests.get(url)
        match response.status_code:
            case 200:
                with open(local_path, "wb") as f:
                    f.write(response.content)
                logger.info("HTTPS download complete")
                attempt = 0
                break
            case _:
                logger.warning(f"HTTP status {response.status_code} at download attempt {attempt}")
                attempt += 1
                time.sleep(5)

    if max_retries < attempt:
        if ftp_fallback:
            logger.warning("Attempting FTP fallback")
            _ftp_fallback_download(url=url, local_path=local_path)
        else:
            raise Exception(f"Can't download {url} using HTTPS")


def _ftp_fallback_download(url, local_path) -> None:
    url = url.replace("https://", "ftp://")
    max_retries = 5
    wait_time = 30
    retries = 0

    while retries < max_retries:
        try:
            spliturl: urllib.parse.SplitResult = urlsplit(url)
            ftp = FTP(spliturl.hostname)
            ftp.login()
            ftp.cwd(str(pathlib.Path(urlsplit(url).path).parent))
            with open(local_path, "wb") as file:
                ftp.retrbinary("RETR " + local_path, file.write)
                ftp.quit()
                logger.info("FTP download completed")
                return
        except Exception as e:
            if "421" in str(e):
                retries += 1
                logger.debug(f"FTP server is busy. Waiting and retrying. Retry {retries} of {max_retries}")
                time.sleep(wait_time)
            else:
                logger.critical(f"Download failed: {e}")
                raise Exception
