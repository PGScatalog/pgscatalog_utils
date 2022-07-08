import os
import pandas as pd
import logging
from .harmonised import remap_harmonised
from .qc import quality_control

logger = logging.getLogger(__name__)


def load_scorefile(path: str, use_harmonised: bool = True) -> pd.DataFrame:
    logger.debug(f'Reading scorefile {path}')
    return (pd.read_table(path, dtype=_scorefile_dtypes(), comment='#', na_values=['None'])
            .pipe(remap_harmonised, use_harmonised=use_harmonised)
            .assign(filename_prefix=_get_basename(path),
                    filename=path)
            .pipe(quality_control))


def _scorefile_dtypes() -> dict[str]:
    """ Data types for columns that might be found in a scorefile """
    return {'rsID': str, 'chr_name': str, 'chr_position': pd.UInt64Dtype(), 'effect_allele': 'str',
            'effect_weight': float, 'locus_name': str, 'OR': float, 'hm_source': str, 'hm_rsID': str,
            'hm_chr': str, 'hm_pos': pd.UInt64Dtype(), 'hm_inferOtherAllele': str}


def _get_basename(path: str) -> str:
    """ Return the basename of a scoring file without extension """
    return os.path.basename(path).split('.')[0]

