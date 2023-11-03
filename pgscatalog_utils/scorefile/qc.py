import logging
import typing

from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.header import ScoringFileHeader
from pgscatalog_utils.scorefile.liftover import liftover

logger = logging.getLogger(__name__)


def quality_control(variants, header: ScoringFileHeader, harmonised: bool, wide: bool):
    variants = remap_harmonised(variants, harmonised)

    if Config.drop_missing:
        variants = drop_hla(variants)

    variants = assign_effect_type(variants)
    variants = check_effect_weight(variants)
    variants = assign_other_allele(variants)

    if wide:
        # wide data must be sorted because:
        # - check_duplicates requires sorted input
        # - output would be unsorted, which looks a little bit messy
        variants = (x for x in sorted(variants, key=lambda x: x["accession"]))

    variants = check_duplicates(variants)

    if Config.liftover:
        variants = liftover(
            variants,
            harmonised=harmonised,
            current_build=header.genome_build,
            target_build=Config.target_build,
        )

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
        if variant["effect_allele"] != "P" or variant["effect_allele"] != "N":
            yield variant
        else:
            n_dropped += 1

    logger.warning(f"{n_dropped} HLA alleles detected and dropped")


def check_effect_weight(variants):
    for variant in variants:
        try:
            variant["effect_weight"] = float(variant["effect_weight"])
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError
        yield variant


def assign_other_allele(variants):
    n_dropped = 0
    for variant in variants:
        if "other_allele" in variant:
            if "/" in variant["other_allele"]:
                # drop multiple other alleles
                n_dropped += 1
                variant["other_allele"] = None
        else:
            variant["other_allele"] = None

        yield variant

    if n_dropped > 0:
        logger.warning(f"Multiple other_alleles detected in {n_dropped} variants")
        logger.warning("Other allele for these variants is set to missing")


def assign_effect_type(variants):
    for variant in variants:
        if "is_recessive" not in variant and "is_dominant" not in variant:
            variant["effect_type"] = "additive"
        else:
            if variant["is_recessive"] == "TRUE":
                variant["effect_type"] = "recessive"
            elif variant["is_dominant"] == "TRUE":
                variant["effect_type"] = "dominant"
            elif variant["is_recessive"] == "TRUE" and variant["is_dominant"] == "TRUE":
                logger.critical(f"Bad effect type setting: {variant}")
                raise Exception

        yield variant


def remap_harmonised(variants, harmonised: bool):
    n_bad = 0
    if harmonised:
        for variant in variants:
            if variant["hm_chr"]:
                variant["chr_name"] = variant["hm_chr"]

            if variant["hm_pos"]:
                variant["chr_position"] = variant["hm_pos"]

            if "hm_inferOtherAllele" in variant and variant.get("other_allele") is None:
                variant["other_allele"] = variant["hm_inferOtherAllele"]

            if (
                "chr_name" in variant
                and "chr_position" in variant
                and "effect_weight" in variant
            ):
                yield variant
            elif Config.drop_missing:
                continue
                # (don't yield anything, filtering out missing variants)
            else:
                # assume a bad harmonisation with no genomic coordinates
                # these will get labelled as duplicates eventually (probably)
                variant["chr_name"] = None
                variant["chr_position"] = None
                yield variant
                n_bad += 1
    else:
        for variant in variants:
            yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} variants failed harmonisation")
