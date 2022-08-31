import logging
import re

import pandas as pd

logger = logging.getLogger(__name__)


def melt_effect_weights(df: pd.DataFrame) -> pd.DataFrame:
    """ Ensure all dataframes are in long format, with one effect weight column and a score accession column """
    elongate = _detect_multiple_weight_columns(df)

    if elongate:
        logger.debug("Melting effect weights")
        return _melt(df)
    else:
        logger.debug("Skipping melt")
        df['accession'] = df['filename']
        return df


def _detect_multiple_weight_columns(df: pd.DataFrame) -> bool:
    """ Detect if multiple effect weight columns are present

    Single weight format:
    | chr_name | chr_pos | effect_allele | effect_weight

    Multiple weight format:
    | chr_name | chr_pos | effect_allele | effect_weight_score_1 | ... | effect_weight_score_n
    """
    columns: list[re.match | None] = [re.search("^effect_weight$", x) for x in df.columns.to_list()]
    columns_suffix: list[re.match | None] = [re.search("^effect_weight_[A-Za-z0-9]+$", x) for x
                                             in df.columns.to_list()]

    if any([col for col in columns]):
        logger.debug("Single effect weight column detected")
        return False
    elif any([col for col in columns_suffix]):
        logger.debug("Multiple weight weight columns detected")
        return True
    else:
        logger.error("ERROR: Missing valid effect weight columns")
        raise Exception("Bad effect weights")


def _melt(df: pd.DataFrame) -> pd.DataFrame:
    """ Melt a multiple effect weight format """
    ew_cols: list[str] = df.filter(regex="effect_weight_*").columns.to_list()
    return df.melt(value_vars=ew_cols, value_name="effect_weight", var_name="accession")
