class EffectAllele:
    _valid_bases = frozenset({"A", "C", "T", "G"})
    __slots__ = ("allele", "is_valid")

    def __init__(self, allele: str):
        self.allele = allele
        self.is_valid = self.is_valid_allele()

    def __repr__(self):
        return f'{type(self).__name__}("{self.allele}")'

    def __str__(self):
        return self.allele

    def is_valid_allele(self) -> bool:
        return not frozenset(self.allele) - self._valid_bases
