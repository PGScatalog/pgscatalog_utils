import functools
import gzip
import pathlib
from dataclasses import dataclass

from pgscatalog_utils.scorefile.config import Config
from pgzip import pgzip

from pgscatalog_utils.download.GenomeBuild import GenomeBuild


@dataclass
class ScoringFileHeader:
    pgs_id: str
    pgp_id: str
    pgs_name: str
    genome_build: GenomeBuild
    variants_number: int
    trait_reported: str
    trait_efo: str
    trait_mapped: str
    weight_type: str
    citation: str
    HmPOS_build: GenomeBuild
    HmPOS_date: str
    format_version: str

    def __post_init__(self):
        if self.variants_number:
            self.variants_number = int(self.variants_number)

        self.genome_build = GenomeBuild.from_string(self.genome_build)
        if self.HmPOS_build:
            self.HmPOS_build = GenomeBuild.from_string(self.HmPOS_build)

    @classmethod
    def from_path(cls, path: pathlib.Path):
        raw_header: dict = raw_header_to_dict(read_header(path))
        # only keep keys needed by class but support partial headers with None values
        keep_keys = ScoringFileHeader.__annotations__.keys()
        header_dict = {k: raw_header.get(k) for k in keep_keys}
        # ... so we can unpack the dict into a dataclass

        if "HmPOS_build" not in header_dict:
            # working with pgs catalog formatted header but unharmonised data
            header_dict["HmPOS_build"] = None

        if not all([v is None for _, v in header_dict.items()]):
            return ScoringFileHeader(**header_dict)
        else:
            # no header available
            raise Exception("No header detected in scoring file")


def raw_header_to_dict(header):
    d = {}
    for item in header:
        key, value = item.split("=")
        d[key[1:]] = value  # drop # character from key
    return d


def read_header(path: pathlib.Path):
    """Parses the header of a PGS Catalog format scorefile into a dictionary"""
    open_function = auto_open(path)
    with open_function(path, "rt") as f:
        yield from _gen_header_lines(f)


def _gen_header_lines(f):
    for line in f:
        if line.startswith("#"):
            if "=" in line:
                yield line.strip()
        else:
            # stop reading lines
            break


def auto_open(filepath):
    with open(filepath, "rb") as test_f:
        if test_f.read(2) == b"\x1f\x8b":
            gzipped = True
        else:
            gzipped = False

    if gzipped and Config.threads > 1:
        return functools.partial(pgzip.open, thread=Config.threads)
    elif gzipped:
        return gzip.open
    elif not gzipped:
        return open
