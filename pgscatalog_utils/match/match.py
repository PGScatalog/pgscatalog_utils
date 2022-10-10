import gc
import logging
import os
from tempfile import TemporaryDirectory

import polars as pl

logger = logging.getLogger(__name__)


# @profile  # decorator needed to annotate memory profiles, but will cause NameErrors outside of profiling
def get_all_matches(scorefile: pl.LazyFrame, target: pl.LazyFrame, low_memory: bool = True) -> pl.LazyFrame:
    scorefile_oa = scorefile.filter(pl.col("other_allele") != None)
    scorefile_no_oa = scorefile.filter(pl.col("other_allele") == None)

    matches: list[pl.LazyFrame()] = []
    col_order = ['row_nr', 'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type',
                 'accession', 'effect_allele_FLIP', 'other_allele_FLIP',
                 'ID', 'REF', 'ALT', 'is_multiallelic', 'matched_effect_allele', 'match_type']

    logger.debug("Getting matches for scores with effect allele and other allele")
    matches.append(_match_variants(scorefile=scorefile_oa, target=target, match_type="refalt").select(col_order))
    matches.append(_match_variants(scorefile_oa, target, match_type="altref").select(col_order))
    matches.append(_match_variants(scorefile_oa, target, match_type="refalt_flip").select(col_order))
    matches.append(_match_variants(scorefile_oa, target, match_type="altref_flip").select(col_order))

    logger.debug("Getting matches for scores with effect allele only")
    matches.append(_match_variants(scorefile_no_oa, target, match_type="no_oa_ref").select(col_order))
    matches.append(_match_variants(scorefile_no_oa, target, match_type="no_oa_alt").select(col_order))
    matches.append(_match_variants(scorefile_no_oa, target, match_type="no_oa_ref_flip").select(col_order))
    matches.append(_match_variants(scorefile_no_oa, target, match_type="no_oa_alt_flip").select(col_order))

    if low_memory:
        logger.debug("Batch collecting matches (low memory mode)")
        match_lf = _batch_collect(matches)
    else:
        logger.debug("Collecting all matches (parallel)")
        match_lf = pl.concat(pl.collect_all(matches))

    return match_lf.lazy()


def _batch_collect(matches: list[pl.LazyFrame]) -> pl.DataFrame:
    """ A slower alternative to pl.collect_all(), but this approach will use less peak memory

    This batches the .collect() and writes intermediate results to a temporary working directory

    IPC files are binary and remember column schema. Reading them can be extremely fast. """
    with TemporaryDirectory() as temp_dir:
        n_chunks = 0
        for i, match in enumerate(matches):
            out_path = os.path.join(temp_dir, str(i) + ".ipc")
            match.collect().write_ipc(out_path)
            n_chunks += 1
        logger.debug(f"Staged {n_chunks} match chunks to {temp_dir}")
        gc.collect()
        return pl.read_ipc(os.path.join(temp_dir, "*.ipc"))


def _match_variants(scorefile: pl.LazyFrame, target: pl.LazyFrame, match_type: str) -> pl.LazyFrame:
    logger.debug(f"Matching strategy: {match_type}")
    match match_type:
        case 'refalt':
            score_keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
            target_keys = ["#CHROM", "POS", "REF", "ALT"]
            effect_allele_column = "effect_allele"
        case 'altref':
            score_keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
            effect_allele_column = "effect_allele"
        case 'refalt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF", "ALT"]
            effect_allele_column = "effect_allele_FLIP"
        case 'altref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP", "other_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT", "REF"]
            effect_allele_column = "effect_allele_FLIP"
        case 'no_oa_ref':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "REF"]
            effect_allele_column = "effect_allele"
        case 'no_oa_alt':
            score_keys = ["chr_name", "chr_position", "effect_allele"]
            target_keys = ["#CHROM", "POS", "ALT"]
            effect_allele_column = "effect_allele"
        case 'no_oa_ref_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "REF"]
            effect_allele_column = "effect_allele_FLIP"
        case 'no_oa_alt_flip':
            score_keys = ["chr_name", "chr_position", "effect_allele_FLIP"]
            target_keys = ["#CHROM", "POS", "ALT"]
            effect_allele_column = "effect_allele_FLIP"
        case _:
            logger.critical(f"Invalid match strategy: {match_type}")
            raise Exception

    missing_cols = ['REF', 'ALT']
    if match_type.startswith('no_oa'):
        if match_type.startswith('no_oa_ref'):
            missing_cols = ['REF']
        else:
            missing_cols = ['ALT']
    join_cols = ['ID'] + missing_cols
    return (scorefile.join(target, score_keys, target_keys, how='inner')
            .with_columns([pl.col("*"),
                           pl.col(effect_allele_column).alias("matched_effect_allele"),
                           pl.lit(match_type).alias("match_type")])
            .join(target.select(join_cols), on="ID", how="inner"))  # get REF / ALT back after first join
