import logging
import typing
from dataclasses import dataclass

from pgscatalog_utils import config
from pgscatalog_utils.download.Catalog import CatalogResult
from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.download.ScoringFile import ScoringFile
from pgscatalog_utils.download.ScoringFileChecksum import ScoringFileChecksum
from pgscatalog_utils.download.download_file import download_file

logger = logging.getLogger(__name__)


@dataclass
class ScoringFileDownloader:
    """
    Simplifies downloading scoring files from the PGS Catalog using FTP or HTTPS.

    HTTPS is preferred but mysteriously fails sometimes. FTP gets busy.
    """
    results: list[CatalogResult]
    genome_build: typing.Union[GenomeBuild, None]
    ftp_fallback: bool = True
    overwrite: bool = True

    def download_files(self):
        url_dict = {}
        for result in self.results:
            url_dict = url_dict | result.get_download_urls()

        scoring_files: list[ScoringFile] = []
        for pgs_id, scoring_file_list in url_dict.items():
            scoring_files.append(list(filter(lambda x: x.build == self.genome_build, scoring_file_list))[0])

        for scoring_file in scoring_files:
            download_file(scoring_file.url, scoring_file.local_path, overwrite=self.overwrite,
                          ftp_fallback=self.ftp_fallback)
            checksum: ScoringFileChecksum = ScoringFileChecksum.from_scoring_file(scoring_file)

            if not checksum.matches:
                logger.warning(f"Scoring file {scoring_file.local_path} fails validation")
                logger.warning(f"Remote checksum: {checksum.remote_checksum}")
                logger.warning(f"Local checksum: {checksum.local_checksum}")
                attempt = 0
                while attempt < config.MAX_RETRIES:
                    download_file(scoring_file.url, scoring_file.local_path, ftp_fallback=self.ftp_fallback,
                                  overwrite=config.OVERWRITE)
                    checksum: ScoringFileChecksum = ScoringFileChecksum.from_scoring_file(scoring_file)

                    if checksum.matches:
                        break
                    else:
                        attempt += 1

            if checksum.matches:
                logger.info("Checksum matches")
                continue
