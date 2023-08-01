import argparse
import logging
import math
import os
import pathlib

import pandas as pd

from pathlib import Path
from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def _parse_args(args=None) -> argparse.Namespace:
    d: str = "Convert pgscatalog/pgsc_calc samplesheet file to JSON and check its contents."
    e: str = "Example usage: python check.py <FILE_IN> <FILE_OUT>"

    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=d, epilog=e)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def _truncate_chrom(chrom):
    try:
        return str(int(chrom))  # truncate numeric chromosomes 22.0 -> 22
    except ValueError:  # it's OK if chrom is a string e.g. MT / X / Y
        return chrom
    except TypeError:  # also OK if chrom is missing entirely
        return None


def _check_colnames(df: pd.DataFrame):
    mandatory: list[str] = ['sampleset', 'path_prefix', 'chrom', 'format']
    optional: list[str] = ['vcf_genotype_field']

    if not set(mandatory) == set(df.columns):
        if set(mandatory + optional) == set(df.columns):
            # this is fine
            return
        else:
            logger.critical("Samplesheet has invalid header row")
            logger.critical(f"Column names must only include: {mandatory}")
            [logger.critical(f"Invalid column name: {col}") for col in df if col not in mandatory]
            raise Exception


def _check_unique_paths(df: pd.DataFrame):
    """ Each row in a samplesheet should have a unique path """
    duplicated: pd.Series = df['path_prefix'].duplicated()
    for idx, duplicate in duplicated.items():
        if duplicate:
            bad_record = df.iloc[:idx]
            logger.critical(f"Duplicated path found in samplesheet:\n{bad_record}")


def _check_empty_paths(df: pd.DataFrame):
    """ Paths are mandatory """
    empty_paths: pd.Series = df['path_prefix'].isnull()
    for idx, empty in empty_paths.items():
        if empty:
            logger.critical(f"Empty path found in samplesheet:\n {df.iloc[[idx]]}")
            raise Exception


def _read_samplesheet(path: str) -> pd.DataFrame:
    csv: pd.DataFrame = pd.read_csv(path, sep=',', header=0)
    csv['chrom'] = csv['chrom'].apply(_truncate_chrom)
    return csv


def _check_paths(df: pd.DataFrame) -> None:
    _check_empty_paths(df)
    _check_unique_paths(df)


def _get_chrom_list(df: pd.DataFrame) -> dict[str, list[str | None]]:
    chrom_dict = {}
    for idx, row in df.iterrows():
        key = row['sampleset']
        value = row['chrom']
        try:
            if math.isnan(value):
                value = None
        except TypeError:
            pass
        chroms = chrom_dict.get(key, [])
        chroms.append(value)
        chrom_dict.update({key: chroms})

    return chrom_dict


def _check_chrom_duplicates(sampleset: str, chrom_list: dict) -> None:
    seen = set()
    duplicate_chromosomes: list[str] = [str(x) for x in chrom_list if x in seen or seen.add(x)]
    if duplicate_chromosomes:
        logger.critical(f"Duplicate chromosomes detected in sampleset {sampleset}")
        logger.critical(f"Duplicate chromosomes: {duplicate_chromosomes}")
        raise Exception


def _check_multiple_missing_chrom(sampleset: str, chrom_list: dict) -> None:
    for chrom in chrom_list:
        if chrom is None and len(chrom_list) != 1:
            logger.critical(f"Sampleset {sampleset} has rows with multiple missing chromosomes")
            logger.critical("If you have file with multiple chromosomes, delete the duplicate rows")
            logger.critical("If your data are split per chromosome, then chromosomes must be set for all rows")
            raise Exception


def _check_chrom(df: pd.DataFrame) -> None:
    # get a list of chroms per sampleset and check them for some basic errors
    chroms: dict = _get_chrom_list(df)

    for sampleset, chrom_list in chroms.items():
        _check_chrom_duplicates(sampleset, chrom_list)
        _check_multiple_missing_chrom(sampleset, chrom_list)


def _check_format(df: pd.DataFrame):
    """ Make sure the file format is a valid choice """
    for idx, row in df.iterrows():
        valid_formats: list[str] = ['vcf', 'pfile', 'bfile']
        if row['format'] not in valid_formats:
            logger.critical(f"Invalid format: {row['format']} must be one of {valid_formats}")
            logger.critical(f"\n{df.iloc[[idx]]}")
            raise Exception


def _setup_paths(df: pd.DataFrame) -> pd.DataFrame:
    """ Add suffix to path prefixes depending on file format / type """
    paths: list[pd.Series] = []
    for idx, row in df.iterrows():
        suffix: list[str]
        match row['format']:
            case 'vcf':
                logger.info("Setting VCF input")
                suffix = ['.vcf.gz']
            case 'bfile':
                logger.info("Setting plink1 binary fileset (bfile) input")
                suffix = ['.bed', '.bim', '.fam']
            case 'pfile':
                logger.info("Setting plink2 binary fileset (pfile) input")
                suffix = ['.pgen', '.pvar', '.psam']
            case _:
                raise Exception

        resolved_paths: list[str] = _resolve_paths([row['path_prefix'] + x for x in suffix], row['format'])
        paths.append(pd.Series(data=[resolved_paths], index=[idx]))

    df['path'] = pd.concat(paths)
    return df


def _resolve_compressed_variant_path(path: str) -> pathlib.Path:
    # .bim.zst | .bim -> OK
    # .pvar.zst | .pvar -> OK
    # anything else not OK
    zstd_ext: str = '.zst'
    compressed_path: pathlib.Path = pathlib.Path(path + zstd_ext).resolve()
    uncompressed_path: pathlib.Path = pathlib.Path(path).resolve()

    # prefer compressed data
    if compressed_path.exists():
        logger.info(f"Found compressed variant information file {compressed_path.name}")
        return compressed_path
    elif uncompressed_path.exists():
        logger.info(f"Couldn't find compressed variant information file, trying {uncompressed_path.name}")
        return uncompressed_path
    else:
        logger.critical(f"{compressed_path} doesn't exist")
        logger.critical(f"{uncompressed_path} doesn't exist")
        logger.critical("Couldn't find variant information files, please check samplesheet path_prefix and try again")
        raise Exception


def _resolve_paths(path_list: list[str], filetype: str) -> list[str]:
    resolved_list: list[str] = []
    for path in path_list:
        if not Path(path).is_absolute():
            logger.warning("Relative path detected in samplesheet. Set absolute paths to silence this warning.")
            logger.warning("Assuming program working directory is a nextflow work directory (e.g. work/4c/8585/...)")
            base_dir: Path = Path(os.getcwd()).parent.parent.parent
            logger.warning(f"Resolving paths relative to work directory parent {base_dir}")
            path = str(base_dir.joinpath(path))

        match filetype:
            case 'pfile' | 'bfile':
                if path.endswith('.bim') or path.endswith('.pvar'):
                    resolved = _resolve_compressed_variant_path(path)
                else:
                    # bed / pgen | fam / psam
                    resolved = pathlib.Path(path).resolve()
            case 'vcf':
                resolved = pathlib.Path(path).resolve()
            case _:
                logger.critical(f"Unsupported filetype {filetype}")
                raise Exception

        if resolved.exists():
            logger.info(f"{resolved} exists")
            resolved_list.append(str(resolved))
        else:
            logger.critical(f"{resolved} doesn't exist, please check samplesheet path_prefix and try again")
            raise FileNotFoundError

    return resolved_list


def _check_genotype_field(df: pd.DataFrame) -> pd.DataFrame:
    df['vcf_import_dosage'] = False  # (dosage off by default)
    if 'vcf_genotype_field' in df.columns:
        logger.debug("vcf_genotype_field detected")
        for index, row in df.iterrows():
            if row['vcf_genotype_field'] not in ['GT', 'DS']:
                missing: bool  # missing dosage is OK
                try:
                    missing = math.isnan(row['vcf_genotype_field'])
                except TypeError:
                    missing = False

                if not missing:
                    logger.critical(f"Invalid entry in vcf_genotype_field: {row['vcf_genotype_field']}")
                    logger.critical(f"\n {row}")
                    raise Exception

        df.loc[df['vcf_genotype_field'] == 'DS', 'vcf_import_dosage'] = True
    else:
        logger.info("no vcf_genotype_field detected")

    return df


def _check_reserved_names(df: pd.DataFrame):
    if any(df['sampleset'] == 'reference'):
        logger.critical("Samplesets must not be named 'reference', please rename in the sample sheet")
        raise Exception


def _check_one_sampleset(df: pd.DataFrame):
    samplesets = set(df['sampleset'].to_list())
    if len(samplesets) > 1:
        logger.critical(f"Multiple samplesets defined in the samplesheet {samplesets}")
        sampleset_error = """ Only one sampleset per samplesheet is supported
        Your genomic data should _only_ be split by chromosome
        pgsc_calc works best with cohorts
        Individual VCFs should be merged into a multi-sample VCF
        If you want to process multiple cohorts, please run pgsc_calc multiple times with different samplesheets. """
        [logger.critical(x.strip()) for x in sampleset_error.split('\n')]
        raise Exception("Multiple samplesets")


def check_samplesheet() -> None:
    """
    This function checks that the samplesheet follows the following structure:
    sampleset,vcf_path,bfile_path,chrom,chunk
    cineca_synthetic_subset,cineca_synthetic_subset.vcf.gz,,22,
    """
    args = _parse_args()
    config.set_logging_level(args.verbose)
    df = _read_samplesheet(args.FILE_IN)

    # check df for errors
    _check_one_sampleset(df)
    _check_reserved_names(df)
    _check_colnames(df)
    _check_paths(df)
    _check_chrom(df)
    _check_format(df)

    # add information to df
    df = _setup_paths(df)
    df = _check_genotype_field(df)  # dosages

    logger.info("Samplesheet checks complete")
    (df.drop(['path_prefix'], axis=1)
     .to_json(args.FILE_OUT, orient='records'))
    logger.info(f"JSON file successfully written to {args.FILE_OUT}")


if __name__ == "__main__":
    check_samplesheet()
