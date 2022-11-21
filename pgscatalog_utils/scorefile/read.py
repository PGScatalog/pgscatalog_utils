import gzip
import io
import logging
import os

import pandas as pd

logger = logging.getLogger(__name__)


def load_scorefile(path: str) -> tuple[dict, pd.DataFrame]:
    logger.debug(f'Reading scorefile {path}')
    df = pd.read_table(path, dtype=_scorefile_dtypes(), comment='#', na_values=['None'], low_memory=False)
    return (_read_header(path),
            df.assign(filename_prefix=get_scorefile_basename(path), filename=path, row_nr=df.index))


def _read_header(path: str) -> dict:
    """Parses the header of a PGS Catalog format scorefle into a dictionary"""
    f = io.TextIOWrapper(gzip.open(path, 'r'))
    try:
        f.readline()
    except gzip.BadGzipFile:
        f = open(path, 'r')

    header = {}
    lastline = '#'
    while lastline.startswith('#'):
        lastline = f.readline()
        line = lastline.strip()
        if line.startswith('#'):
            if '=' in line:
                line = line[1:].split('=')
                field, val = [x.strip() for x in line]
                if field in remap_header:
                    header[remap_header[field]] = val
                else:
                    header[field] = val

    if ('genome_build' in header) and (header['genome_build'] == 'NR'):
        header['genome_build'] = None
    f.close()
    return header


def _scorefile_dtypes() -> dict[str]:
    """ Data types for columns that might be found in a scorefile """
    return {'rsID': str, 'chr_name': str, 'chr_position': pd.UInt64Dtype(), 'effect_allele': 'str',
            'effect_weight': float, 'locus_name': str, 'OR': float, 'hm_source': str, 'hm_rsID': str,
            'hm_chr': str, 'hm_pos': pd.UInt64Dtype(), 'hm_inferOtherAllele': str}


def get_scorefile_basename(path: str) -> str:
    """ Return the basename of a scoring file without extension """
    filename = os.path.basename(path)
    if filename.endswith('.txt.gz'):
        filename = filename.replace('.txt.gz', '')
    elif filename.endswith('.txt'):
        filename = filename.replace('.txt', '')
    return filename


remap_header = {
    'PGS ID': 'pgs_id',
    'PGS Name': 'pgs_name',
    'Reported Trait': 'trait_reported',
    'Original Genome Build': 'genome_build',
    'Number of Variants': 'variants_number',
    'PGP ID': 'pgp_id',
    'Citation': 'citation',
    'LICENSE': 'license',
    # Harmonization related
    'HmPOS Build': 'HmPOS_build',
    'HmPOS Date': 'HmPOS_date',
    'HmVCF Reference': 'HmVCF_ref',
    'HmVCF Date': 'HmVCF_date',
    'HmVCF N Matched Variants': 'HmVCF_n_matched',
    'HmVCF N Unmapped Variants': 'HmVCF_n_unmapped'
}  # Used to maintain reverse compatibility to old scoring files
