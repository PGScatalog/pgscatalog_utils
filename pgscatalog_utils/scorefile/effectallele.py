import logging

logger = logging.getLogger(__name__)


class EffectAllele:
    # (class attribute, so shared)
    _valid_bases = frozenset({"A", "C", "T", "G"})

    @classmethod
    def is_valid(cls, effect_allele: str) -> bool:
        return not frozenset(effect_allele) - cls._valid_bases
