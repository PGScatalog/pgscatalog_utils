import logging
import os
import typing
import pgzip

import polars as pl
from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def write_out(matches: pl.LazyFrame, split: bool, dataset: str):
    if split:
        # Loop through chromosomes with matched variants
        chroms: list[str] = matches.select("chr_name").unique().collect().get_column("chr_name").to_list()
        for chrom in chroms:
            # 1. filter by chromosome
            chrom_df: pl.LazyFrame = matches.filter(pl.col('chr_name') == chrom)
            # 2. split by effect type
            additive: pl.LazyFrame
            dominant: pl.LazyFrame
            recessive: pl.LazyFrame
            additive, dominant, recessive = _split_effect_type(chrom_df)

            # 3. deduplicate
            effect_types = ['additive', 'dominant', 'recessive']
            deduped = dict(zip(effect_types, [_deduplicate_variants(x) for x in [additive, dominant, recessive]]))

            # 4. pivot and write!
            _write_split(deduped, chrom, dataset)
    else:
        # 2. split by effect type
        additive: pl.LazyFrame
        dominant: pl.LazyFrame
        recessive: pl.LazyFrame
        additive, dominant, recessive = _split_effect_type(matches)

        # 3. deduplicate
        effect_types = ['additive', 'dominant', 'recessive']
        deduped = dict(zip(effect_types, [_deduplicate_variants(x) for x in [additive, dominant, recessive]]))

        # 4. pivot and write!
        _write_split(deduped, 'ALL', dataset)


def _write_split(deduplicated: dict[str: tuple[int, pl.LazyFrame]], chrom: str, dataset: str):
    for effect_type, df_lst in deduplicated.items():
        for i, et_df in df_lst:
            if i is False:
                # deduplication returned an empty dataframe, so skip (normally recessive or dominant)
                continue

            # pivoting is !! _expensive_ !! (it collects the lazyframe)
            pivoted: pl.LazyFrame = _pivot_score(et_df, chrom)
            fout = os.path.join(config.OUTDIR, f"{dataset}_{chrom}_{effect_type}_{i}.scorefile.gz")
            _write_scorefile(pivoted, fout)


def _write_scorefile(df, fout):
    logger.debug(f"Writing matched scorefile to {fout}")
    with pgzip.open(fout, 'wb', thread=config.POLARS_MAX_THREADS) as f:
        df.collect().write_csv(f)


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
            with pgzip.open(fout, 'wb', thread=config.POLARS_MAX_THREADS) as f:
                df.collect().write_csv(f)
        case _:
            logger.critical(f"Invalid format: {file_format}")
            raise Exception


def _pivot_score(df: pl.LazyFrame, chrom: str) -> pl.LazyFrame:
    """ Format a dataframe to plink2 --score standard
    Minimum example:
    ID | effect_allele | effect_weight
    Multiple scores are OK too:
    ID | effect_allele | weight_1 | ... | weight_n
    """
    logger.debug(f"Pivoting score for chromosome {chrom}")
    return (df.collect()
            .pivot(index=["ID", "matched_effect_allele", "effect_type"], values="effect_weight",
                   columns="accession")
            .rename({"matched_effect_allele": "effect_allele"})
            .fill_null(strategy="zero")
            .lazy())


def _deduplicate_variants(df: pl.LazyFrame) -> list[tuple[int, pl.LazyFrame]]:
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
    if df.select('ID').head().collect().is_empty():
        logger.info("Empty input: skipping deduplication")
        return [(False, df)]
    else:
        logger.debug("Deduplicating variants")

    # 1. unique ID - EA is important because normal duplicates are already
    #   handled by pivoting, and it's pointless to split them unnecessarily
    # 2. use cumcount to number duplicate IDs
    # 3. join cumcount data on original DF, use this data for splitting
    # note: effect_allele should be equivalent to matched_effect_allele
    ea_count: pl.LazyFrame = (df.select(['ID', 'matched_effect_allele'])
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
        x: pl.LazyFrame = (dup_label.filter(pl.col("cumcount") == i).drop(['cumcount', 'count']))
        df_lst.append((i, x))

    return df_lst


def _split_effect_type(df: pl.LazyFrame) -> tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]:
    additive = df.filter(pl.col('effect_type') == 'additive')
    dominant = df.filter(pl.col('effect_type') == 'dominant')
    recessive = df.filter(pl.col('effect_type') == 'recessive')
    return additive, dominant, recessive
