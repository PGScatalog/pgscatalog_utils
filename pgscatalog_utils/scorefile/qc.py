import logging
import typing


from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.effecttype import EffectType
from pgscatalog_utils.scorefile.scoringfileheader import ScoringFileHeader
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.scorevariant import ScoreVariant

logger = logging.getLogger(__name__)


def quality_control(
    variants: typing.Generator[ScoreVariant, None, None],
    header: ScoringFileHeader,
    harmonised: bool,
    wide: bool,
) -> typing.Generator[ScoreVariant, None, None]:
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
    variants = check_effect_allele(variants)
    variants = detect_complex(variants)

    if wide:
        # wide data must be sorted because check_duplicates requires sorted input
        variants = (x for x in sorted(variants, key=lambda x: x["accession"]))

    variants = check_duplicates(variants)

    return variants


def check_duplicates(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    seen_ids: dict = {}
    current_accession: typing.Union[str, None] = None
    n_duplicates: int = 0
    n_variants: int = 0
    for variant in variants:
        accession: str = variant.accession

        if accession != current_accession:
            seen_ids = {}
            current_accession = accession

        # None other allele -> empty string
        variant_id: str = ":".join(
            [
                str(getattr(variant, k) or "")
                for k in ["chr_name", "chr_position", "effect_allele", "other_allele"]
            ]
        )

        if variant_id in seen_ids:
            variant.is_duplicated = True
            n_duplicates += 1

        seen_ids[variant_id] = True

        yield variant
        n_variants += 1

    if n_duplicates > 0:
        logger.warning(
            f"{n_duplicates} of {n_variants} variants are duplicated in: {current_accession}"
        )


def drop_hla(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    n_dropped = 0
    for variant in variants:
        match variant:
            case {"effect_allele": "P"} | {"effect_allele": "N"}:
                n_dropped += 1
                continue
            case _:
                yield variant

    logger.warning(f"{n_dropped} HLA alleles detected and dropped")


def check_effect_weight(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    for variant in variants:
        try:
            float(variant.effect_weight)
            yield variant
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError


def assign_other_allele(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    n_dropped = 0
    for variant in variants:
        if "/" in variant.other_allele:
            n_dropped += 1
            variant.other_allele = None

        yield variant

    if n_dropped > 0:
        logger.warning(f"Multiple other_alleles detected in {n_dropped} variants")
        logger.warning("Other allele for these variants is set to missing")


def assign_effect_type(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    for variant in variants:
        match (variant.is_recessive, variant.is_dominant):
            case (None, None) | (False, False):
                pass  # default value is additive, pass to break match and yield
            case (False, True):
                variant.effect_type = EffectType.DOMINANT
            case (True, False):
                variant.effect_type = EffectType.RECESSIVE
            case _:
                logger.critical(f"Bad effect type setting: {variant}")
                raise Exception
        yield variant


def remap_harmonised(
    variants: typing.Generator[ScoreVariant, None, None], harmonised: bool
) -> typing.Generator[ScoreVariant, None, None]:
    if harmonised:
        for variant in variants:
            # using the harmonised field in the header to make sure we don't accidentally overwrite
            # positions with empty data (e.g. in an unharmonised file)
            # if harmonisation has failed we _always_ want to use that information
            variant.chr_name = variant.hm_chr
            variant.chr_position = variant.hm_pos
            if variant.other_allele is None:
                variant.other_allele = variant.hm_inferOtherAllele
            yield variant
    else:
        for variant in variants:
            # can't remap, so don't try
            yield variant


def check_bad_variant(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    n_bad = 0
    for variant in variants:
        match variant:
            case (
                ScoreVariant(chr_name=None)
                | ScoreVariant(chr_position=None)
                | ScoreVariant(effect_allele=None)
            ):
                # (effect weight checked separately)
                n_bad += 1
                if not Config.drop_missing:
                    yield variant
            case _:
                yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} bad variants")


def check_effect_allele(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    n_bad = 0
    for variant in variants:
        if not variant.effect_allele.is_snp:
            n_bad += 1

        yield variant

    if n_bad > 1:
        logger.warning(f"{n_bad} variants have invalid effect alleles (not ACTG)")


def detect_complex(
    variants: typing.Generator[ScoreVariant, None, None]
) -> typing.Generator[ScoreVariant, None, None]:
    """Some older scoring files in the PGS Catalog are complicated.
    They often require bespoke set up to support interaction terms, etc
    """
    is_complex = False

    for variant in variants:
        if not is_complex:
            if variant.is_complex:
                is_complex = True

        yield variant

    if is_complex:
        logger.warning("Complex scoring file detected")
        logger.warning(
            "Complex files are difficult to calculate properly and may require manual intervention"
        )
