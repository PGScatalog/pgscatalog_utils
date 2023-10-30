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
    # don't actually use converted value
    for variant in variants:
        try:
            float(variant['effect_weight'])
        except ValueError:
            logger.critical(f"{variant} has bad effect weight")
            raise ValueError
        yield variant


def assign_other_allele(variants):
    for variant in variants:
        if 'other_allele' not in variant:
            variant['other_allele'] = None
        yield variant
