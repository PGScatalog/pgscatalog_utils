import logging
import pathlib
import typing
from dataclasses import dataclass, field

from pgscatalog_utils.download.GenomeBuild import GenomeBuild

logger = logging.getLogger(__name__)


@dataclass
class ScoringFile:
    """
    A scoring file stored remotely on the PGS Catalog
    """
    url: str
    harmonized: bool
    build: typing.Union[GenomeBuild, None]
    local_path: str = field(init=False)

    @classmethod
    def from_result(cls, result: dict):
        harmonized: bool
        url: str
        instances: list[ScoringFile] = []

        for build in [GenomeBuild.GRCh38, GenomeBuild.GRCh37, None]:
            match build:
                case GenomeBuild.GRCh37:
                    url = result.get("ftp_harmonized_scoring_files").get("GRCh37").get("positions")
                    harmonized = True
                case GenomeBuild.GRCh38:
                    url = result.get("ftp_harmonized_scoring_files").get("GRCh38").get("positions")
                    harmonized = True
                case None:
                    url = result.get("ftp_scoring_file")
                    harmonized = False
                case _:
                    raise Exception(f"Invalid genome build: {build}")
            instances.append(cls(url, harmonized, build))

        return instances

    def __post_init__(self):
        self.local_path = pathlib.Path(self.url).name
