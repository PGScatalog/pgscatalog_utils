import pandas as pd
import logging

logger = logging.getLogger(__name__)


def quality_control(df: pd.DataFrame) -> pd.DataFrame:
    """ Do quality control checks on a scorefile """
    _check_shape(df)
    _check_columns(df)
    logger.debug("Quality control: checking for bad variants")
    return (df.pipe(_drop_hla)
            .pipe(_drop_missing_variants)
            .pipe(_check_duplicate_identifiers)
            .pipe(_drop_multiple_oa))


def _drop_multiple_oa(df: pd.DataFrame) -> pd.DataFrame:
    """ Set alleles to None in hm_inferOtherAllele if they contain multiple alleles

    e.g. A / C / T -> None; A -> A; A / C -> None
    """
    if df['other_allele'].str.contains('/').any():
        logger.debug("Multiple inferred other alleles detected, dropping other alleles for ambiguous variants")

    df['other_allele'] = df['other_allele'].replace(regex='.+\\/.+', value=None)
    return df


def _drop_missing_variants(df: pd.DataFrame) -> pd.DataFrame:
    no_na: pd.DataFrame = df.dropna(subset=['chr_name', 'chr_position', 'effect_weight'])
    n_dropped = df.shape[0] - no_na.shape[0]

    if n_dropped > 0:
        logger.warning(f"{n_dropped} variants with missing values detected and dropped from scoring file")

    return no_na


def _drop_hla(df: pd.DataFrame) -> pd.DataFrame:
    """ Drop HLA effect alleles with present / absent encoding """

    no_hla: pd.DataFrame = df.query('effect_allele != "P" | effect_allele != "N"')

    if df.shape[0] > no_hla.shape[0]:
        logger.debug("HLA alleles detected and dropped")

    return no_hla


def _check_duplicate_identifiers(df: pd.DataFrame) -> pd.DataFrame:
    unique: pd.Series = df.groupby(['chr_name', 'chr_position', 'effect_allele', 'other_allele']).size() == 1

    if unique.all():
        return df
    else:
        raise Exception("Duplicate variants in scoring file")


def _check_shape(df: pd.DataFrame) -> None:
    assert len(df.columns) > 1, "ERROR: scorefile not formatted correctly (0 columns)"
    assert df.shape[0] > 0, "ERROR: No variants detected in input file (0 rows)"


def _check_columns(df: pd.DataFrame) -> None:
    assert {'chr_name', 'chr_position'}.issubset(df.columns), "If you're using rsids did you request harmonised data?"
    assert 'effect_allele' in df, "ERROR: Missing effect allele column"

