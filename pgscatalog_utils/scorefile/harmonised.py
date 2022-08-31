import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def remap_harmonised(df: pd.DataFrame, use_harmonised) -> pd.DataFrame:
    """ Replace original columns with harmonised data, if available and appropriate """

    if any([re.match("hm_\\w+", x) for x in df.columns]) and use_harmonised:
        logger.debug("Harmonised columns detected and used")
        hm_colnames: dict[str: str] = {'hm_chr': 'chr_name', 'hm_pos': 'chr_position',
                                       'hm_inferOtherAllele': 'other_allele'}

        if 'other_allele' not in df or all(df['other_allele'].isnull()):
            logger.debug("other_allele column contains no information, replacing with hm_inferOtherAllele")
            return (df.drop(['chr_name', 'chr_position', 'other_allele'], axis=1, errors='ignore')
                    .rename(hm_colnames, axis=1))
        else:
            logger.debug("other_allele column contains information, dropping hm_inferOtherAllele")
            return (df.drop(['chr_name', 'chr_position', 'hm_inferOtherAllele'], axis=1, errors='ignore')
                    .rename(hm_colnames, axis=1))
    elif any([re.match("hm_\\w+", x) for x in df.columns]) and not use_harmonised:
        logger.debug(f"Harmonised columns detected but not used (use_harmonised={use_harmonised})")
        return df
    else:
        logger.debug("Harmonised columns not detected")
        return df
