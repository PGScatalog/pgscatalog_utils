import pandas as pd
import logging
import sqlite3

logger = logging.getLogger(__name__)


def write_scorefile(df: pd.DataFrame, path: str) -> None:
    cols: list[str] = ['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
                       'accession']

    if df.empty:
        logger.error("Empty scorefile output! Please check the input data")
        raise Exception
    else:
        logger.debug("Writing out combined scorefile")
        out_df: pd.DataFrame = (df.drop('accession', axis=1)
                                .rename({'filename_prefix': 'accession'}, axis=1)
                                .pipe(_filter_failed_liftover))
        _write_log(out_df)
        out_df[cols].to_csv(path, index=False, sep="\t")


def _filter_failed_liftover(df: pd.DataFrame) -> pd.DataFrame:
    if 'liftover' in df:
        logger.debug("Filtering variants that failed liftover")
        return df.query('liftover == True')
    else:
        return df


def _write_log(df: pd.DataFrame) -> None:
    logger.debug("Writing log to local database")
    conn: sqlite3.Connection = sqlite3.connect('scorefiles.db')

    if 'liftover' not in df:
        df = df.assign(liftover=None, lifted_chr=None, lifted_pos=None)

    cols: list[str] = ['chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
                       'accession', 'liftover', 'lifted_chr', 'lifted_pos']

    # change some column types for sqlite
    # nullable_ints: list[str] = ['liftover', 'lifted_chr', 'lifted_pos']
    # df[nullable_ints] = df[nullable_ints].astype(pd.Int64Dtype())
    df['other_allele'] = df['other_allele'].astype(str)
    df[cols].to_sql('scorefile', conn, if_exists='replace')
    conn.close()
