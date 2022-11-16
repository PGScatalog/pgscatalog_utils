import argparse
import logging

import polars as pl

from pgscatalog_utils import config
from pgscatalog_utils.match.label import make_params_dict, label_matches
from pgscatalog_utils.match.match_variants import log_and_write, add_match_args
from pgscatalog_utils.match.read import read_scorefile

logger = logging.getLogger(__name__)


def combine_matches():
    args = _parse_args()
    config.set_logging_level(args.verbose)
    config.setup_polars_threads(args.n_threads)
    config.setup_tmpdir(args.outdir, combine=True)
    config.OUTDIR = args.outdir

    with pl.StringCache():
        scorefile = read_scorefile(path=args.scorefile, chrom=None)  # chrom=None to read all variants
        logger.debug("Reading matches")
        matches = pl.concat([pl.scan_ipc(x, memory_map=False, rechunk=False) for x in args.matches], rechunk=False)

        logger.debug("Labelling match candidates")
        params: dict[str, bool] = make_params_dict(args)
        matches = matches.pipe(label_matches, params)

        # make sure there's no duplicate variant_ids across matches in multiple pvars
        # processing batched chromosomes with overlapping variants might cause problems
        # e.g. chr1 1-100000, chr1 100001-500000
        _check_duplicate_vars(matches)

        dataset = args.dataset.replace('_', '-')  # _ used as delimiter in pgsc_calc
        log_and_write(matches=matches, scorefile=scorefile, dataset=dataset, args=args)


def _check_duplicate_vars(matches: pl.LazyFrame):
    max_occurrence: list[int] = (matches.filter(pl.col('match_status') == 'matched')
                                 .groupby(['accession', 'ID'])
                                 .agg(pl.count())
                                 .select('count')
                                 .max()
                                 .collect()
                                 .get_column('count')
                                 .to_list())
    assert max_occurrence == [1], "Duplicate IDs in final matches"


def _parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dataset', dest='dataset', required=True,
                        help='<Required> Label for target genomic dataset')
    parser.add_argument('-s', '--scorefile', dest='scorefile', required=True,
                        help='<Required> Path to scorefile')
    parser.add_argument('-m', '--matches', dest='matches', required=True, nargs='+',
                        help='<Required> List of match files')
    parser.add_argument('--min_overlap', dest='min_overlap', required=True,
                        type=float, help='<Required> Minimum proportion of variants to match before error')
    parser = add_match_args(parser) # params for labelling matches
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help='<Required> Output directory')
    parser.add_argument('--split', dest='split', default=False, action='store_true',
                        help='<Optional> Split scorefile per chromosome?')
    parser.add_argument('-n', dest='n_threads', default=1, help='<Optional> n threads for matching', type=int)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    combine_matches()
