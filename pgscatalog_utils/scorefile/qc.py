import logging


logger = logging.getLogger(__name__)

def drop_hla(variants):
    logger.info("Checking for HLA alleles")
    for variant in variants:
        if variant['effect_allele'] != 'P' or variant['effect_allele'] != 'N':
            yield variant
        else:
            logger.warning("HLA alleles detected and dropped")


def check_effect_weight(variants):
    logger.info("Checking effect weights")
    for variant in variants:
        try:
            variant['effect_weight'] = float(variant['effect_weight'])
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError
        yield variant


def assign_other_allele(variants):
    for variant in variants:
        if 'other_allele' not in variant:
            variant['other_allele'] = None
        yield variant

def assign_effect_type(variants):
    logger.info("Assigning effect types")
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

            variant.pop('is_recessive')
            variant.pop('is_dominant')

        yield variant


def remap_harmonised(variants, harmonised: bool):
    if harmonised:
        logger.info("Using harmonised data fields")
    else:
        logger.info("Harmonised data fields not available")

    for variant in variants:
        if harmonised:
            variant['chr_name'] = variant['hm_chr']
            variant['chr_position'] = variant['hm_pos']

            if 'hm_inferOtherAllele' in variant and variant.get('other_allele') is None:
                logger.debug("Replacing missing other_allele with inferred other allele")
                variant['other_allele'] = variant['hm_inferOtherAllele']

            yield {k: v for k, v in variant.items() if not k.startswith("hm")}
        else:
            yield variant
