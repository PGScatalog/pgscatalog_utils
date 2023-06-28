import hashlib
import logging
import pathlib
import typing
from dataclasses import dataclass, field

from pgscatalog_utils import config
from pgscatalog_utils.download.ScoringFile import ScoringFile
from pgscatalog_utils.download.download_file import download_file

logger = logging.getLogger(__name__)


def _generate_md5_checksum(filename: str, blocksize=4096) -> typing.Union[str, None]:
    """ Returns MD5 checksum for the given file. """
    md5 = hashlib.md5()
    try:
        file = open(filename, 'rb')
        with file:
            for block in iter(lambda: file.read(blocksize), b""):
                md5.update(block)
    except IOError:
        logger.warning(f"File {filename} not found!")
        return None

    return md5.hexdigest()


@dataclass
class ScoringFileChecksum:
    local_path: str
    remote_url: str
    remote_checksum: str
    local_checksum: str
    matches: bool = field(init=False)

    def __post_init__(self):
        self.matches = self.remote_checksum == self.local_checksum

    @classmethod
    def from_scoring_file(cls, scoring_file: ScoringFile):
        assert config.OUTDIR.joinpath(scoring_file.local_path).exists(), "Scoring file must be downloaded first"

        logger.info("Downloading checksum file")
        md5_url: str = scoring_file.url + '.md5'
        md5_local_path: str = scoring_file.local_path + '.md5'
        md5_sep = "  "

        # calculate checksum of scoring file
        downloaded_checksum: str = _generate_md5_checksum(scoring_file.local_path)

        # grab checksum from pgs catalog and read it
        remote_checksum: str = ""
        download_file(url=md5_url, local_path=md5_local_path, overwrite=config.OVERWRITE, ftp_fallback=True)
        with open(config.OUTDIR.joinpath(md5_local_path), 'r') as f:
            remote_checksum = f.read().split(md5_sep)[0]

        return cls(local_path=md5_local_path, remote_url=md5_url, remote_checksum=remote_checksum,
                   local_checksum=downloaded_checksum)
