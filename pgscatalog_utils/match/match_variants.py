import argparse
import logging


def parse_args(args=None):
    parser = argparse.ArgumentParser(description='Read and format scoring files')
    parser.add_argument('-d', '--dataset', dest='dataset', required=True,
                        help='<Required> Label for target genomic dataset (e.g. "-d thousand_genomes")')
    parser.add_argument('-s', '--scorefiles', dest='scorefile', required=True,
                        help='<Required> Combined scorefile path (output of read_scorefiles.py)')
    parser.add_argument('-t', '--target', dest='target', required=True,
                        help='<Required> A table of target genomic variants (.bim format)')
    parser.add_argument('--split', dest='split', default=False, action='store_true',
                        help='<Required> Split scorefile per chromosome?')
    parser.add_argument('-n', '--n_threads', dest='n_threads', default=1, type=int,
                        help='<Required> Number of threads used to match (default = 1)')
    parser.add_argument('--format', required=True, dest='plink_format', help='<Required> bim or pvar?')
    parser.add_argument('--db', dest='db', required=True, help='<Required> path to database')
    parser.add_argument('-m', '--min_overlap', dest='min_overlap', required=True,
                        type=float, help='<Required> Minimum proportion of variants to match before error')
    parser.add_argument('--keep_ambiguous', dest='remove_ambiguous', action='store_false',
                        help='Flag to force the program to keep variants with ambiguous alleles, (e.g. A/T and G/C '
                             'SNPs), which are normally excluded (default: false). In this case the program proceeds '
                             'assuming that the genotype data is on the same strand as the GWAS whose summary '
                             'statistics were used to construct the score.'),
    parser.add_argument('--keep_multiallelic', dest='remove_multiallelic', action='store_false',
                        help='Flag to allow matching to multiallelic variants (default: false).')
    return parser.parse_args(args)


def match_variants():
    args = parse_args()

    logger = logging.getLogger(__name__)
    log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"


if __name__ == "main":
    match_variants()