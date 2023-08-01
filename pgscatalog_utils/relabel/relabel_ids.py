import argparse
import gzip
import io
import logging
import operator
from functools import reduce

import zstandard

from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def _parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Relabel the column values in one file based on a pair of columns in another",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--dataset', dest='dataset', required=True,
                        help='<Required> Label for target genomic dataset')
    parser.add_argument("-m", "--maps", help='mapping filenames', dest='map_files', nargs='+', required=True)
    parser.add_argument("--col_from", help='column to change FROM', dest='col_from', required=True)
    parser.add_argument("--col_to", help='column to change TO', dest='col_to', required=True)
    parser.add_argument("--target_file", help='target file', dest='target_file', required=True)
    parser.add_argument("--target_col", help='target column to revalue', dest='target_col', required=True)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    parser.add_argument('--split', dest='split', action='store_true', required=False)
    parser.add_argument('--combined', dest='combined', action='store_true', required=False)
    parser.add_argument('-cc', '--comment_char', dest='comment_char', default='##')
    args = parser.parse_args()

    if not (args.split or args.combined):
        parser.error("At least one of --combined or --split is required")

    return args


def _is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def _is_zstd_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(4) == b'\x28\xb5\x2f\xfd'


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


def _detect_target(path):
    if _is_gz_file(path):
        return "gzip"
    elif _is_zstd_file(path):
        return "zstd"
    else:
        return "text"


def _open_output(path, header):
    logger.debug(f"Opening {path} and writing header")
    outf = gzip.open(path, 'wt')
    outf.write('\t'.join(header) + '\n')
    return outf


def _get_chrom(line, id_idx):
    return line[id_idx].split(':')[0]


def _get_outf_path(current_chrom, dataset):
    return f"{dataset}_{current_chrom}_relabelled.gz"


def _relabel_target(args, mapping, split_output):
    with open(args.target_file, 'rb') as f:
        match _detect_target(args.target_file):
            case 'zstd':
                dctx = zstandard.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    _relabel(in_target=io.TextIOWrapper(reader), mapping=mapping, split_output=split_output, args=args)
            case 'gzip':
                with gzip.open(f) as reader:
                    _relabel(in_target=io.TextIOWrapper(reader), mapping=mapping, split_output=split_output, args=args)
            case 'text':
                _relabel(in_target=io.TextIOWrapper(f), mapping=mapping, split_output=split_output, args=args)
            case _:
                raise Exception("Can't detect target format")


def _relabel(in_target, args, mapping, split_output):
    h = in_target.readline()
    if args.comment_char:
        while h.startswith(args.comment_char):
            h = in_target.readline()
    h = h.strip().split()

    i_target_col = h.index(args.target_col)

    if not split_output:
        current_chrom: str = 'ALL'
        outf_path = _get_outf_path(current_chrom=current_chrom, dataset=args.dataset)
        outf = _open_output(outf_path, h)

    for i, line in enumerate(in_target):
        line = line.strip().split()
        # get the first column index that contains :
        # assume this contains a variant ID e.g. 1:1234:A:C
        # this column index can change across different types of files
        id_idx = [i for i, x in enumerate(line) if ':' in x][0]

        if split_output and i == 0:
            current_chrom = _get_chrom(line, id_idx)
            logger.debug(f"Creating split output, current chrom: {_get_chrom(line, id_idx)}")
            outf_path = _get_outf_path(current_chrom=current_chrom, dataset=args.dataset)
            outf = _open_output(outf_path, h)

        if split_output and current_chrom != _get_chrom(line, id_idx):
            logger.debug(f"New chromosome {_get_chrom(line, id_idx)} detected in split mode, writing to new file")
            outf.close()

            current_chrom = _get_chrom(line, id_idx)
            outf_path = _get_outf_path(current_chrom=current_chrom, dataset=args.dataset)
            outf = _open_output(outf_path, h)

        line[i_target_col] = mapping[line[i_target_col]]  # revalue column
        outf.write('\t'.join(line) + '\n')

    outf.close()


def relabel_ids():
    args = _parse_args()
    config.set_logging_level(args.verbose)

    # sometimes we want to write the relabelled data out as split _and_ combined
    split_output = []
    if args.split:
        logger.debug("Writing split output enabled")
        split_output.append(True)

    if args.combined:
        logger.debug("Writing combined output enabled")
        split_output.append(False)

    map_list = [open_map(x, args.col_from, args.col_to) for x in args.map_files]
    mapping = reduce(operator.ior, map_list, {})  # merge dicts quickly, ior is equivalent to | operator
    del map_list

    if not len(mapping) > 1:
        logger.critical("Empty mapping file inputs, please check --maps files")
        raise Exception

    # Read, relabel and output file
    # note: if --split and --combined, reads the input file twice. not ideal but it's pretty quick
    [_relabel_target(args=args, mapping=mapping, split_output=x) for x in split_output]
    logger.debug("Finished relabelling")


if __name__ == "__main__":
    relabel_ids()
