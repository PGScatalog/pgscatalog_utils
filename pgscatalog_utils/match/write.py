import logging
import os
import typing
import pgzip

import polars as pl
from pgscatalog_utils import config

logger = logging.getLogger(__name__)


def write_log(df: pl.LazyFrame, prefix: str, chrom: typing.Union[str, None], outdir: str) -> None:
    if chrom is None:
        log_name: str = os.path.join(os.path.abspath(outdir), f"{prefix}_log")
    else:
        log_name: str = os.path.join(os.path.abspath(outdir), f"{prefix}_chrom{chrom}_log")

    fout: str = ''.join([log_name, ".csv.gz"])
    if os.path.exists(fout):
        logger.warning(f"Overwriting log that already exists: {fout}")
        os.remove(fout)

    _write_text_pgzip(df=df, sep = ',', fout=fout)


def write_scorefiles(matches: pl.LazyFrame, split: bool, dataset: str):
    _check_column_types(matches)
    additive: pl.LazyFrame
    dominant: pl.LazyFrame
    recessive: pl.LazyFrame

    # collect and cache minimum required columns
    min_cols: list[str] = ['accession', 'effect_type', 'chr_name', 'ID', 'matched_effect_allele', 'effect_weight']
    matches: pl.LazyFrame = (matches.select(min_cols)
                             .collect()
                             .lazy())

    if split:
        chroms: list[str] = matches.select("chr_name").unique().collect().get_column("chr_name").to_list()
        for chrom in chroms:
            # 1. filter by chromosome
            chrom_df: pl.LazyFrame = matches.filter(pl.col('chr_name') == chrom)
            # 2. split by effect type
            additive, dominant, recessive = _split_effect_type(chrom_df)

            # 3. deduplicate
            effect_types = ['additive', 'dominant', 'recessive']
            deduped = dict(zip(effect_types, [_deduplicate_variants(x) for x in [additive, dominant, recessive]]))

            # 4. pivot and write!
            _write_split(deduped, chrom, dataset)
    else:
        # 1. split by effect type
        additive, dominant, recessive = _split_effect_type(matches)

        # 2. deduplicate
        effect_types = ['additive', 'dominant', 'recessive']
        deduped = dict(zip(effect_types, [_deduplicate_variants(x) for x in [additive, dominant, recessive]]))

        # 3. pivot and write!
        _write_split(deduped, 'ALL', dataset)


def _check_column_types(matches: pl.LazyFrame):
    logger.debug("Checking column types")
    # these columns are most important for writing out
    correct_schema = {'chr_name': pl.Categorical, 'chr_position': pl.UInt64, 'ID': pl.Utf8,
                      'matched_effect_allele': pl.Categorical, 'effect_weight': pl.Float64,
                      'effect_type': pl.Categorical, 'accession': pl.Categorical}
    col_types = {x: matches.schema.get(x) for x in list((matches.schema.keys() & correct_schema.keys()))}
    assert col_types == correct_schema


def _write_split(deduplicated: dict[str: tuple[int, pl.LazyFrame]], chrom: str, dataset: str):
    for effect_type, df_lst in deduplicated.items():
        for i, et_df in df_lst:
            if i is False:
                # deduplication returned an empty dataframe, so skip (normally recessive or dominant)
                continue

            # pivoting is !! _expensive_ !! (it collects the lazyframe)
            pivoted: pl.LazyFrame = _pivot_score(et_df, chrom)

            dout = os.path.abspath(config.OUTDIR)
            fout = os.path.join(dout, f"{dataset}_{chrom}_{effect_type}_{i}.scorefile.gz")
            _write_text_pgzip(pivoted, fout)


def _write_text_pgzip(df: pl.LazyFrame, fout: str, sep: str = '\t', append: bool = False):
    """ Write a df to a text file (e.g. CSV / TSV) using parallel gzip, optionally appending to an existing file

    Notes:
    - Compression performance isn't ideal when concatenating gzip streams (append = True)
    - Generally it's best to feed compression algorithms all data and write in one go
    - However, df will normally be very big
    - It's collected for the first time in this function, and joins _a lot_ of data (contains all match candidates)
    - The files created by this function must be human-readable text files, so feather / parquet isn't helpful
    - Hopefully appending gzip streams is a reasonable compromise to mitigate OOM errors
    """
    if append:
        logger.debug(f"Appending to {fout}")
        mode = 'ab'
    else:
        logger.debug(f"Writing to {fout}")
        mode = 'wb'

    with pgzip.open(fout, mode, thread=config.N_THREADS) as f:
        df.collect().write_csv(f, sep=sep)


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
            .drop("effect_type")
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
