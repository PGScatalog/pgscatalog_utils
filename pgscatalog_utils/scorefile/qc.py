import logging

logger = logging.getLogger(__name__)


def quality_control(variants, harmonised):
    variants = remap_harmonised(variants, harmonised)
    variants = drop_hla(variants)
    variants = assign_effect_type(variants)
    variants = check_effect_weight(variants)
    variants = assign_other_allele(variants)
    variants = check_duplicates(variants)
    return variants


def check_duplicates(variants):
    seen_ids = set()
    current_accession = None

    for variant in variants:
        accession = variant['accession']

        if accession != current_accession:
            seen_ids = set()
            current_accession = accession

        # None other allele -> empty string
        id = ":".join([str(variant[k] or "") for k in
                       ['chr_name', 'chr_position', 'effect_allele', 'other_allele']])

        if id in seen_ids:
            logger.warning(
                f"Duplicate variant found: {variant['accession']}: {id} {variant['row_nr']}")
            variant['is_duplicated'] = True
        else:
            variant['is_duplicated'] = False

        seen_ids.add(id)

        yield variant


def drop_hla(variants):
    for variant in variants:
        if variant['effect_allele'] != 'P' or variant['effect_allele'] != 'N':
            yield variant
        else:
            logger.warning("HLA alleles detected and dropped")


def check_effect_weight(variants):
    for variant in variants:
        try:
            variant['effect_weight'] = float(variant['effect_weight'])
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError
        yield variant


def assign_other_allele(variants):
    for variant in variants:
        if 'other_allele' in variant:
            if '/' in variant['other_allele']:
                # drop multiple other alleles
                variant['other_allele'] = None

        if 'other_allele' not in variant:
            variant['other_allele'] = None
        yield variant


def assign_effect_type(variants):
    for variant in variants:
        if 'is_recessive' not in variant and 'is_dominant' not in variant:
            variant['effect_type'] = 'additive'

        if 'is_recessive' in variant or 'is_dominant' in variant:
            logger.info("Recessive or dominant variant detected")
            if variant['is_recessive']:
                variant['effect_type'] = 'recessive'
            elif variant['is_dominant']:
                variant['effect_type'] = 'dominant'
            elif variant['is_recessive'] and variant['is_dominant']:
                logger.critical(f"Bad effect type setting: {variant}")
                raise Exception

        yield variant


def remap_harmonised(variants, harmonised: bool):
    for variant in variants:
        if harmonised:
            # if harmonised data are available, always overwrite
            if variant['hm_chr']:
                variant['chr_name'] = variant['hm_chr']

            if variant['hm_pos']:
                variant['chr_position'] = variant['hm_pos']

            if 'hm_inferOtherAllele' in variant and variant.get('other_allele') is None:
                logger.debug(
                    "Replacing missing other_allele with inferred other allele")
                variant['other_allele'] = variant['hm_inferOtherAllele']

            yield variant
        else:
            yield variant
