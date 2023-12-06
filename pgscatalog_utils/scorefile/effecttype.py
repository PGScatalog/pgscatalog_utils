from enum import Enum


class EffectType(Enum):
    RECESSIVE = "recessive"
    DOMINANT = "dominant"
    ADDITIVE = "additive"

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        # pasting __repr__ output should be sufficient to construct the class
        return f"{type(self).__name__}.{self.name}"
