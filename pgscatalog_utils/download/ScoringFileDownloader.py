import logging
import pathlib
import time
import typing
import urllib.parse
from dataclasses import dataclass
from enum import Enum, auto
from ftplib import FTP
from urllib.parse import urlsplit

import requests

from pgscatalog_utils.download.CatalogQuery import CatalogResult
from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.download.ScoringFile import ScoringFile

logger = logging.getLogger(__name__)


def _ftp_fallback_download(url, local_path):
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


@dataclass
class ScoringFileDownloader:
    """
    Simplifies downloading scoring files from the PGS Catalog using FTP or HTTPS.

    HTTPS is preferred but mysteriously fails sometimes. FTP gets busy.
    """
    class DownloadMethod(Enum):
        HTTPS = auto()
        FTP = auto()

    results: list[CatalogResult]
    genome_build: typing.Union[GenomeBuild, None]
    download_method: DownloadMethod = DownloadMethod.HTTPS
    ftp_fallback: bool = True
    overwrite: bool = True

    def _download_file(self, url, local_path):
        if pathlib.Path(local_path).exists():
            if not self.overwrite:
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
            if self.ftp_fallback:
                logger.warning("Attempting FTP fallback")
                _ftp_fallback_download(url=url, local_path=local_path)
            else:
                raise Exception(f"Can't download {url} using HTTPS")

    def download_files(self):
        url_dict = {}
        for result in self.results:
            url_dict = url_dict | result.get_download_urls()

        scoring_files: list[ScoringFile] = []
        for pgs_id, scoring_file_list in url_dict.items():
            scoring_files.append(list(filter(lambda x: x.build == self.genome_build, scoring_file_list))[0])

        for scoring_file in scoring_files:
            self._download_file(scoring_file.url, scoring_file.local_path)
