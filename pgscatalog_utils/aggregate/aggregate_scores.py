import argparse
import textwrap

import pandas as pd

from pgscatalog_utils.config import set_logging_level
import glob
import logging

logger = logging.getLogger(__name__)


def aggregate_scores():
    args = _parse_args()
    set_logging_level(args.verbose)
    df = aggregate(list(set(args.scores)))
    logger.debug("Compressing and writing combined scores")
    df.to_csv('aggregated_scores.txt.gz', sep='\t', compression='gzip')


def aggregate(scorefiles: list[str]):
    combined = pd.DataFrame()
    aggcols = set()

    for i, path in enumerate(scorefiles):
        logger.debug(f"Reading {path}")
        # pandas can automatically detect zst compression, neat!
        df = (pd.read_table(path)
              .assign(sampleset=path.split('_')[0])
              .set_index(['sampleset', '#IID']))

        df.index.names = ['sampleset', 'IID']

        # Subset to aggregatable columns
        df = df[_select_agg_cols(df.columns)]
        aggcols.update(set(df.columns))

        # Combine DFs
        if i == 0:
            logger.debug('Initialising combined DF')
            combined = df.copy()
        else:
            logger.debug('Adding to combined DF')
            combined = combined.add(df, fill_value=0)

    assert all([x in combined.columns for x in aggcols]), "All Aggregatable Columns are present in the final DF"

    return combined.pipe(_calculate_average)


def _calculate_average(combined: pd.DataFrame):
    logger.debug("Averaging data")
    avgs = combined.loc[:, combined.columns.str.endswith('_SUM')].divide(combined['DENOM'], axis=0)
    avgs.columns = avgs.columns.str.replace('_SUM', '_AVG')
    return pd.concat([combined, avgs], axis=1)


def _select_agg_cols(cols):
    keep_cols = ['DENOM']
    return [x for x in cols if (x.endswith('_SUM') and (x != 'NAMED_ALLELE_DOSAGE_SUM')) or (x in keep_cols)]


def _description_text() -> str:
    return textwrap.dedent('''
    Aggregate plink .sscore files into a combined TSV table.
    
    This aggregation sums scores that were calculated from plink
    .scorefiles. Scorefiles may be split to calculate scores over different
    chromosomes or effect types. The PGS Catalog calculator automatically splits
    scorefiles where appropriate, and uses this script to combine them.
    
    Input .sscore files can be optionally compressed with zstd or gzip. 
    
    The aggregated output scores are compressed with gzip.
   ''')


def _parse_args(args=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=_description_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--scores', dest='scores', required=True, nargs='+',
                        help='<Required> List of scorefile paths. Use a wildcard (*) to select multiple files.')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
                        default='scores/', help='<Required> Output directory to store downloaded files')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    aggregate_scores()

