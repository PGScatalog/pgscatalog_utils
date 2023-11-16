from enum import Enum


class EffectType(Enum):
    RECESSIVE = "recessive"
    DOMINANT = "dominant"
    ADDITIVE = "additive"

    def __str__(self):
        return str(self.value)
