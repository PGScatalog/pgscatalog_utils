import argparse
import gzip
import io
import logging
import os
import shutil
from functools import reduce
import operator

from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def _parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Relabel the column values in one file based on a pair of columns in another",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--maps", help='mapping filenames', dest='map_files', nargs='+', required=True)
    parser.add_argument("--col_from", help='column to change FROM', dest='col_from', required=True)
    parser.add_argument("--col_to", help='column to change TO', dest='col_to', required=True)
    parser.add_argument("--target_file", help='target file', dest='target_file', required=True)
    parser.add_argument("--target_col", help='target column to revalue', dest='target_col', required=True)
    parser.add_argument("--out", help='output filename', dest='out_file', required=True)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


def _is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def _read_map(in_map, col_from, col_to):
    h = in_map.readline().strip().split()
    i_from = h.index(col_from)
    i_to = h.index(col_to)

    mapping = {}
    for line in in_map:
        line = line.strip().split()
        mapping[line[i_from]] = line[i_to]

    return mapping


def open_map(path, col_from, col_to):
    # Read the mapping file
    if _is_gz_file(path):
        logger.debug(f"Reading map file {path} with gzip.open")
        with gzip.open(path, 'rt') as in_map:
            return _read_map(in_map, col_from, col_to)
    else:
        logger.debug(f"Reading map file {path} with open")
        with open(path, 'r') as in_map:
            return _read_map(in_map, col_from, col_to)


def _open_output(path, header):
    logger.debug(f"Opening {path} and writing header")
    outf = gzip.open(path, 'wt')
    outf.write('\t'.join(header) + '\n')
    return outf


def _get_chrom(line, id_idx):
    return line[id_idx].split(':')[0]


def relabel_ids():
    args = _parse_args()
    config.set_logging_level(args.verbose)

    # more than one mapping file input -> assume split input, so produce split output
    if len(args.map_files) > 1:
        logger.debug(f"--maps n args {len(args.map_files)}, setting split_output = True")
        split_output = True
    else:
        logger.debug(f"--maps n args {len(args.map_files)}, setting split_output = False")
        split_output = False

    map_list = [open_map(x, args.col_from, args.col_to) for x in args.map_files]

    mapping = reduce(operator.ior, map_list, {})  # merge dicts quickly, ior is equivalent to | operator
    del map_list

    # Read, relabel and output file
    with open(args.target_file, 'r') as in_target:
        h = in_target.readline().strip().split()
        i_target_col = h.index(args.target_col)
        out_suffix = ".gz"

        if not split_output:
            current_chrom: str = 'ALL'
            outf = _open_output(f"{current_chrom}_{args.out_file}{out_suffix}", h)

        for i, line in enumerate(in_target):
            line = line.strip().split()
            # get the first column index that contains :
            # assume this contains a variant ID e.g. 1:1234:A:C
            # this column index can change across different types of files
            id_idx = [i for i, x in enumerate(line) if ':' in x][0]

            if split_output and i == 0:
                current_chrom = _get_chrom(line, id_idx)
                logger.debug(f"Creating split output, current chrom: {line[0]}")
                outf = _open_output(f"{current_chrom}_{args.out_file}{out_suffix}", h)

            if split_output and current_chrom != _get_chrom(line, id_idx):
                logger.debug(f"New chromosome {_get_chrom(line, id_idx)} detected in split mode, writing to new file")
                outf.close()

                current_chrom = _get_chrom(line, id_idx)
                outf = _open_output(f"{current_chrom}_{args.out_file}{out_suffix}", h)

            line[i_target_col] = mapping[line[i_target_col]]  # revalue column
            outf.write('\t'.join(line) + '\n')

        outf.close()
        logger.debug("Finished relabelling")


if __name__ == "__main__":
    relabel_ids()
