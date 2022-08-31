import logging

import pandas as pd

logger = logging.getLogger(__name__)


def set_effect_type(df: pd.DataFrame) -> pd.DataFrame:
    if {'is_recessive', 'is_dominant'}.issubset(df.columns):
        _check_effect_types(df)
        return (df.assign(additive=lambda x: ~x["is_recessive"] & ~x["is_dominant"])
                .assign(effect_type=lambda x: x[["is_recessive", "is_dominant", "additive"]].idxmax(1)))
    else:
        return _set_default_effect_type(df)


def _check_effect_types(df: pd.DataFrame):
    """ Check that only one effect type is set per variant """
    bad_rows: pd.DataFrame = df[['is_dominant', 'is_recessive']].all(axis=1).any()

    error = ''' ERROR: Bad variants in scorefile
    is_recessive and is_dominant columns are both TRUE for a variant
    These columns are mutually exclusive (both can't be true)
    However, both can be FALSE for additive variant scores
    '''
    if bad_rows:
        logger.error(error)
        logger.error(bad_rows)
        raise Exception


def _set_default_effect_type(df: pd.DataFrame, effect_type: str = "additive") -> pd.DataFrame:
    logger.debug(f'No effect types set, using default ({effect_type})')
    return df.assign(effect_type=effect_type)
