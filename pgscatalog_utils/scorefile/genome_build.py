import gzip
import io
import logging
import re
from typing import TextIO
import pandas as pd

logger = logging.getLogger(__name__)


def annotate_build(df: pd.DataFrame) -> pd.DataFrame:
    """ Annotate the dataframe with genome build data """
    builds: pd.DataFrame = _get_builds(df['filename'].drop_duplicates())
    return df.merge(builds, how="left", on="filename")


def _read_header(f: TextIO) -> str:
    """ Extract genome build of scorefile from PGS Catalog header format """
    build_dict = {'GRCh37': 'hg19', 'GRCh38': 'hg38', 'hg19': 'hg19', 'hg38': 'hg38'}
    for line in f:
        if re.search("^#genome_build", line):
            # get #genome_build=GRCh37 from header
            header = line.replace('\n', '').replace('#', '').split('=')
            # and remap to liftover style
            try:
                build: str = build_dict[header[-1]]
                logger.debug("Valid genome build detected")
                return build
            except KeyError:
                raise Exception("Bad genome build detected in header")
        elif line[0] != '#':
            raise Exception("No genome build detected in header")


def _read_build(path: str) -> str:
    """ Open scorefiles and automatically handle compressed input """
    logger.debug(f'Reading header of {path}')
    try:
        with io.TextIOWrapper(gzip.open(path, 'r')) as f:
            return _read_header(f)
    except gzip.BadGzipFile:
        with open(path, 'r') as f:
            return _read_header(f)


def _get_builds(s: pd.Series) -> pd.DataFrame:
    """ Get genome builds for a series of scorefile paths
        | filename | -> | filename | genome_build |
        | x.txt.gz |    | x.txt.gz | hg19         |
    """
    return pd.concat([s, s.apply(_read_build).rename("genome_build")], axis=1)
