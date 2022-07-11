import polars as pl
import logging

logger = logging.getLogger(__name__)


def _write_scorefile(effect_type: str, scorefiles: pl.DataFrame, split: bool) -> None:
    """ Write a list of scorefiles with the same effect type """
    fout: str = '{chr}_{et}_{split}.scorefile'

    # each list element contains a dataframe of variants
    # lists are split to ensure variants have unique ID - effect alleles
    for i, scorefile in enumerate(scorefiles):
        df_dict: dict[str, pl.DataFrame] = _format_scorefile(scorefile, split)  # may be split by chrom

        for k, v in df_dict.items():
            path: str = fout.format(chr=k, et=effect_type, split=i)
            v.write_csv(path, sep="\t")


def _format_scorefile(df: pl.DataFrame, split: bool) -> dict[str, pl.DataFrame]:
    """ Format a dataframe to plink2 --score standard
    Minimum example:
    ID | effect_allele | effect_weight
    Multiple scores are OK too:
    ID | effect_allele | weight_1 | ... | weight_n
    """
    if split:
        chroms: list[int] = df["chr_name"].unique().to_list()
        return {x: (df.filter(pl.col("chr_name") == x)
                    .pivot(index=["ID", "effect_allele"], values="effect_weight", columns="accession")
                    .fill_null(pl.lit(0)))
                for x in chroms}
    else:
        return {'false': (df.pivot(index=["ID", "effect_allele"], values="effect_weight", columns="accession")
                          .fill_null(pl.lit(0)))}


def _split_effect_type(df: pl.DataFrame) -> dict[str, pl.DataFrame]:
    logger.debug("Splitting matches by effect type")
    effect_types: list[str] = df["effect_type"].unique().to_list()
    return {x: df.filter(pl.col("effect_type") == x) for x in effect_types}


def _unduplicate_variants(df: pl.DataFrame) -> list[pl.DataFrame]:
    """ Find variant matches that have duplicate identifiers
    When merging a lot of scoring files, sometimes a variant might be duplicated
    this can happen when the effect allele differs at the same position, e.g.:
        - chr1: chr2:20003:A:C A 0.3 NA
        - chr1: chr2:20003:A:C C NA 0.7
    where the last two columns represent different scores.  plink demands
    unique identifiers! so need to split, score, and sum later
    Parameters:
    df: A dataframe containing all matches, with columns ID, effect_allele, and
        effect_weight
    Returns:
        A list of dataframes, with unique ID - effect allele combinations
    """
    # 1. unique ID - EA is important because normal duplicates are already
    #   handled by pivoting, and it's pointless to split them unnecessarily
    # 2. use cumcount to number duplicate IDs
    # 3. join cumcount data on original DF, use this data for splitting
    ea_count: pl.DataFrame = df.select(["ID", "effect_allele"]).unique().with_columns([
        pl.col("ID").cumcount().over(["ID"]).alias("cumcount"),
        pl.col("ID").count().over(["ID"]).alias("count")
    ])

    dup_label: pl.DataFrame = df.join(ea_count, on=["ID", "effect_allele"], how="left")

    # now split the matched variants, and make sure we don't lose any
    n_splits: int = ea_count.select("cumcount").max()[0, 0] + 1  # cumcount = ngroup-1
    df_lst: list = []
    n_var: int = 0

    for i in range(0, n_splits):
        x: pl.DataFrame = dup_label.filter(pl.col("cumcount") == i)
        n_var += x.shape[0]
        df_lst.append(x)

    if len(df_lst) > 1:
        logger.debug("Duplicate variant identifiers split")
    else:
        logger.debug("No duplicate variant identifiers found")

    assert n_var == df.shape[0]

    return df_lst
