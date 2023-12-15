class EffectAllele:
    """A class that represents an effect allele found in PGS Catalog scoring files

    The allele that's dosage is counted (e.g. {0, 1, 2}) and multiplied by the variant's
    weight (effect_weight) when calculating score. The effect allele is also known as
    the 'risk allele'.
    >>> simple_ea = EffectAllele("A")
    >>> simple_ea
    EffectAllele("A")
    >>> simple_ea.is_snp
    True
    >>> str(simple_ea)
    'A'
    >>> EffectAllele("AG")
    EffectAllele("AG")
    >>> hla_example = EffectAllele("+")
    >>> hla_example
    EffectAllele("+")
    >>> hla_example.is_snp
    False
    """

    _valid_snp_bases = frozenset({"A", "C", "T", "G"})
    __slots__ = ("allele", "is_snp")

    def __init__(self, allele):
        self.allele = str(allele)
        self.is_snp = self._is_snp()

    def __repr__(self):
        return f'{type(self).__name__}("{self.allele}")'

    def __str__(self):
        return self.allele

    def _is_snp(self) -> bool:
        """SNPs are the most common type of effect allele in PGS Catalog scoring
        files. More complex effect alleles, like HLAs or APOE genes, often require
        extra work to represent in genomes. Users should be warned about complex
        effect alleles.
        >>> EffectAllele("+")._is_snp()
        False
        >>> EffectAllele("A")._is_snp()
        True
        """
        return not frozenset(self.allele) - self._valid_snp_bases
