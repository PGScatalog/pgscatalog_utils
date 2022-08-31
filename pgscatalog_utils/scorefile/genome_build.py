import logging

import pandas as pd

logger = logging.getLogger(__name__)


def annotate_build(df: pd.DataFrame, target_build: str) -> pd.DataFrame:
    """ Annotate the dataframe with genome build data  """
    logger.debug(f"Annotating target build: {target_build}")
    build_dict: dict = {'GRCh37': 'hg19', 'GRCh38': 'hg38', 'hg19': 'hg19', 'hg38': 'hg38'}  # standardise build names
    df['chain_target_build'] = build_dict[target_build]
    df = df.assign(chain_genome_build=[build_dict[x] for x in df['genome_build']])
    return df


def build2GRC(build):
    """Map build names so they can be compared with GRCh37 and 38"""
    build_2_GRC_dict = {'GRCh37': 'GRCh37', 'GRCh38': 'GRCh38', 'hg19': 'GRCh37',
                        'hg38': 'GRCh38'}  # standardise build names
    if pd.isnull(build):
        return None
    else:
        return build_2_GRC_dict.get(build)
