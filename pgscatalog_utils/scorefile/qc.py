import logging
import typing

from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.effecttype import EffectType
from pgscatalog_utils.scorefile.header import ScoringFileHeader
from pgscatalog_utils.scorefile.liftover import liftover

logger = logging.getLogger(__name__)


def quality_control(variants, header: ScoringFileHeader, harmonised: bool, wide: bool):
    # order is important for:
    # 1. liftover non-harmonised data (quite rare), failed lifts get None'd
    # 2. remap harmonised data, failed harmonisations get None'd
    # 3. check and optionally drop bad variants
    # where a bad variant has None in a mandatory ScoreVariant field
    # then continue with other QC

    if Config.liftover:
        variants = liftover(
            variants,
            harmonised=harmonised,
            current_build=header.genome_build,
            target_build=Config.target_build,
        )

    variants = remap_harmonised(variants, harmonised)
    variants = check_bad_variant(variants)

    if Config.drop_missing:
        variants = drop_hla(variants)

    variants = assign_effect_type(variants)
    variants = check_effect_weight(variants)
    variants = assign_other_allele(variants)

    if wide:
        # wide data must be sorted because check_duplicates requires sorted input
        variants = (x for x in sorted(variants, key=lambda x: x["accession"]))

    variants = check_duplicates(variants)

    return variants


def check_duplicates(variants):
    seen_ids: dict = {}
    current_accession: typing.Union[str, None] = None
    n_duplicates: int = 0
    n_variants: int = 0
    for variant in variants:
        accession: str = variant["accession"]

        if accession != current_accession:
            seen_ids = {}
            current_accession = accession

        # None other allele -> empty string
        id: str = ":".join(
            [
                str(variant[k] or "")
                for k in ["chr_name", "chr_position", "effect_allele", "other_allele"]
            ]
        )

        if id in seen_ids:
            variant["is_duplicated"] = True
            n_duplicates += 1
        else:
            variant["is_duplicated"] = False

        seen_ids[id] = True

        yield variant
        n_variants += 1

    if n_duplicates > 0:
        logger.warning(
            f"{n_duplicates} of {n_variants} variants are duplicated in: {current_accession}"
        )


def drop_hla(variants):
    n_dropped = 0
    for variant in variants:
        match variant:
            case {"effect_allele": "P"} | {"effect_allele": "N"}:
                n_dropped += 1
                continue
            case _:
                yield variant

    logger.warning(f"{n_dropped} HLA alleles detected and dropped")


def check_effect_weight(variants):
    for variant in variants:
        try:
            float(variant["effect_weight"])
            yield variant
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError


def assign_other_allele(variants):
    n_dropped = 0
    for variant in variants:
        if "/" in variant["other_allele"]:
            n_dropped += 1
            variant["other_allele"] = None

        yield variant

    if n_dropped > 0:
        logger.warning(f"Multiple other_alleles detected in {n_dropped} variants")
        logger.warning("Other allele for these variants is set to missing")


def assign_effect_type(variants):
    for variant in variants:
        match (variant.get("is_recessive"), variant.get("is_dominant")):
            case (None, None) | ("FALSE", "FALSE"):
                pass  # default value is additive
            case ("FALSE", "TRUE"):
                variant["effect_type"] = EffectType.DOMINANT
            case ("TRUE", "FALSE"):
                variant["effect_type"] = EffectType.RECESSIVE
            case _:
                logger.critical(f"Bad effect type setting: {variant}")
                raise Exception
        yield variant


def remap_harmonised(variants, harmonised: bool):
    if harmonised:
        for variant in variants:
            # using the harmonised field in the header to make sure we don't accidentally overwrite
            # positions with empty data (e.g. in an unharmonised file)
            # if harmonisation has failed we _always_ want to use that information
            variant["chr_name"] = variant["hm_chr"]
            variant["chr_position"] = variant["hm_pos"]
            if variant["other_allele"] is None:
                variant["other_allele"] = variant["hm_inferOtherAllele"]
            yield variant
    else:
        for variant in variants:
            # can't remap, so don't try
            yield variant


def check_bad_variant(variants):
    n_bad = 0
    for variant in variants:
        match variant:
            case {"chr_name": None} | {"chr_position": None} | {"effect_allele": None}:
                # (effect weight checked separately)
                n_bad += 1
                if not Config.drop_missing:
                    yield variant
            case _:
                yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} bad variants")
