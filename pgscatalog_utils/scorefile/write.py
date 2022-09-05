import logging

import pandas as pd

logger = logging.getLogger(__name__)


def write_scorefile(df: pd.DataFrame, path: str) -> None:
    cols: list[str] = ['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
                       'is_duplicated', 'accession', 'row_nr']

    if df.empty:
        logger.error("Empty scorefile output! Please check the input data")
        raise Exception
    else:
        out_df: pd.DataFrame = (df.drop('accession', axis=1)
                                .rename({'filename_prefix': 'accession'}, axis=1)
                                .pipe(_filter_failed_liftover))

        if 'other_allele' not in out_df:
            logger.warning("No other allele information detected, writing out as missing data")
            out_df['other_allele'] = None
        if path.endswith('.gz'):
            logger.debug("Writing out gzip-compressed combined scorefile")
            out_df[cols].to_csv(path, index=False, sep="\t", compression='gzip')
        else:
            logger.debug("Writing out combined scorefile")
            out_df[cols].to_csv(path, index=False, sep="\t")


def _filter_failed_liftover(df: pd.DataFrame) -> pd.DataFrame:
    if 'liftover' in df:
        logger.debug("Filtering variants that failed liftover")
        return df.query('liftover == True')
    else:
        return df
