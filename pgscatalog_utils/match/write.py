import gzip
import logging
import os
import typing

import polars as pl

logger = logging.getLogger(__name__)


def write_log(df: pl.LazyFrame, prefix: str, chrom: typing.Union[str, None], file_format: str, outdir: str) -> None:
    # feather file preserves dtypes and is small
    # don't compress the feather file to allow memory mapping
    if chrom is None:
        log_name: str = os.path.join(os.path.abspath(outdir), f"{prefix}_log")
    else:
        log_name: str = os.path.join(os.path.abspath(outdir), f"{prefix}_chrom{chrom}_log")

    match file_format:
        case 'ipc':
            fout: str = ''.join([log_name, ".ipc.zst"])
            logger.debug(f"Writing {fout} in format: {file_format}")
            df.collect().write_ipc(fout, compression='zstd')  # gzip compression not supported
        case 'csv':
            fout: str = ''.join([log_name, ".csv.gz"])
            logger.debug(f"Writing {fout} in format: {file_format}")
            with gzip.open(fout, 'wb') as f:
                df.collect().write_csv(f)
        case _:
            logger.critical(f"Invalid format: {file_format}")
            raise Exception


def write_out(df: pl.LazyFrame, split: bool, outdir: str, dataset: str) -> None:
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    for effect_type in ['additive', 'dominant', 'recessive']:
        logger.debug(f"Splitting by effect type {effect_type}")
        for i, x in enumerate(_deduplicate_variants(effect_type, df)):
            effect_df: pl.LazyFrame = x.filter(pl.col('effect_type') == effect_type)
            if effect_df.fetch().shape[0] > 0:
                chroms: list[int] = effect_df.select("chr_name").unique().collect().get_column("chr_name").to_list()
                params = {'chroms': chroms, 'effect_type': effect_type, 'i': str(i)}
                _write_scorefile(params=params, scorefiles=effect_df, split=split, outdir=outdir, dataset=dataset)
            else:
                logger.debug(f"{effect_type} empty, skipping writing out")
                continue

    logger.debug("All scorefiles written, goodbye!")

def _write_scorefile(params: dict, scorefiles: pl.LazyFrame, split: bool, outdir: str, dataset: str) -> None:
    """ Write a list of scorefiles with the same effect type """
    effect_type = params.get('effect_type')
    i = params.get('i')
    chroms = params.get('chroms')

    dfs: list[pl.LazyFrame] = _format_scorefile(scorefiles, chroms)

    if not split:
        logger.debug("Writing combined scorefile")
        chroms: list[str] = ["ALL"]  # reset chroms list and merge into one df
        out_dfs: list[pl.LazyFrame] = [pl.concat(dfs)]
    else:
        out_dfs: list[pl.LazyFrame] = dfs
        logger.debug("Writing split scorefiles")

    for chrom, scorefile in zip(chroms, out_dfs):
        fout: str = os.path.join(outdir, f"{dataset}_{chrom}_{effect_type}_{i}.scorefile.gz")
        logger.debug(f"Writing matched scorefile to {fout}")
        with gzip.open(fout, 'wb') as f:
            scorefile.collect().write_csv(f, sep="\t")


def _format_scorefile(df: pl.LazyFrame, chroms: list[str]) -> list[pl.LazyFrame]:
    """ Format a dataframe to plink2 --score standard
    Minimum example:
    ID | effect_allele | effect_weight
    Multiple scores are OK too:
    ID | effect_allele | weight_1 | ... | weight_n
    """
    logger.debug("Formatting scorefile to plink2 standard")
    dfs = []
    for chrom in chroms:
        dfs.append(df.filter(pl.col("chr_name") == chrom).collect()
                   .pivot(index=["ID", "matched_effect_allele"], values="effect_weight", columns="accession")
                   .rename({"matched_effect_allele": "effect_allele"})
                   .fill_null(strategy="zero")
                   .lazy())
    return dfs


def _deduplicate_variants(effect_type: str, df: pl.LazyFrame) -> list[pl.LazyFrame]:
    """ Find variant matches that have duplicate identifiers
    When merging a lot of scoring files, sometimes a variant might be duplicated
    this can happen when the matched effect allele differs at the same position, e.g.:
        - chr1: chr2:20003:A:C A 0.3 NA
        - chr1: chr2:20003:A:C C NA 0.7
    where the last two columns represent different scores.  plink demands
    unique identifiers! so need to split, score, and sum later
    Parameters:
    df: A dataframe containing all matches, with columns ID, effect_allele, and
        effect_weight
    Returns:
        A list of dataframes, with unique ID - matched effect allele combinations
    """
    # 1. unique ID - EA is important because normal duplicates are already
    #   handled by pivoting, and it's pointless to split them unnecessarily
    # 2. use cumcount to number duplicate IDs
    # 3. join cumcount data on original DF, use this data for splitting
    ea_count: pl.LazyFrame = (df.select(["ID", "matched_effect_allele"])
    .unique()
    .with_columns([
        pl.col("ID").cumcount().over(["ID"]).alias("cumcount"),
        pl.col("ID").count().over(["ID"]).alias("count")
    ]))

    dup_label: pl.LazyFrame = df.join(ea_count, on=["ID", "matched_effect_allele"], how="left")

    # now split the matched variants, and make sure we don't lose any
    n_splits: int = ea_count.select("cumcount").max().collect()[0, 0] + 1  # cumcount = ngroup-1
    df_lst: list = []

    for i in range(0, n_splits):
        x: pl.LazyFrame = dup_label.filter(pl.col("cumcount") == i)
        df_lst.append(x)

    if len(df_lst) > 1:
        logger.debug(f"Duplicate variant identifiers split for effect type {effect_type}")
    else:
        logger.debug(f"No duplicate variant identifiers found for effect type {effect_type}")

    return df_lst
