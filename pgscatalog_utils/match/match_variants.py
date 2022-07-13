import argparse
import logging
import polars as pl
from pgscatalog_utils.match.postprocess import postprocess_matches

from pgscatalog_utils.log_config import set_logging_level
from pgscatalog_utils.match.match import get_all_matches, check_match_rate
from pgscatalog_utils.match.read import read_target, read_scorefile
from pgscatalog_utils.match.write import write_out


def match_variants():
    args = _parse_args()

    logger = logging.getLogger(__name__)
    set_logging_level(args.verbose)

    scorefile: pl.DataFrame = read_scorefile(path=args.scorefile)
    target: pl.DataFrame = read_target(path=args.target, n_threads=args.n_threads,
                                       remove_multiallelic=args.remove_multiallelic)

    with pl.StringCache():
        matches: pl.DataFrame = get_all_matches(scorefile, target).pipe(postprocess_matches, args.remove_ambiguous)
        check_match_rate(scorefile, matches, args.min_overlap)

    if matches.shape[0] == 0:  # this can happen if args.min_overlap = 0
        logger.error("Error: no target variants match any variants in scoring files")
        raise Exception

    write_out(matches, args.split, args.outdir)


def _parse_args(args=None):
    parser = argparse.ArgumentParser(description='Read and format scoring files')
    parser.add_argument('-s', '--scorefiles', dest='scorefile', required=True,
                        help='<Required> Combined scorefile path (output of read_scorefiles.py)')
    parser.add_argument('-t', '--target', dest='target', required=True,
                        help='<Required> A table of target genomic variants (.bim format)')
    parser.add_argument('--split', dest='split', default=False, action='store_true',
                        help='<Optional> Split scorefile per chromosome?')
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help='<Required> Output directory')
    parser.add_argument('-n', '--n_threads', dest='n_threads', default=1, type=int,
                        help='<Required> Number of threads used to match (default = 1)')
    parser.add_argument('-m', '--min_overlap', dest='min_overlap', required=True,
                        type=float, help='<Required> Minimum proportion of variants to match before error')
    parser.add_argument('--keep_ambiguous', dest='remove_ambiguous', action='store_false',
                        help='Flag to force the program to keep variants with ambiguous alleles, (e.g. A/T and G/C '
                             'SNPs), which are normally excluded (default: false). In this case the program proceeds '
                             'assuming that the genotype data is on the same strand as the GWAS whose summary '
                             'statistics were used to construct the score.'),
    parser.add_argument('--keep_multiallelic', dest='remove_multiallelic', action='store_false',
                        help='Flag to allow matching to multiallelic variants (default: false).')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    match_variants()


# join matches and scorefile with keys depending on liftover
# count match type column
# matches.groupby('accession').agg([pl.count(), (pl.col('match_type') == None).sum().alias('no_match')])
