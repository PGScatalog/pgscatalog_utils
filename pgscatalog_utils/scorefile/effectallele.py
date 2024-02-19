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
    __slots__ = ("_allele", "_is_snp")

    def __init__(self, allele):
        self._allele = str(allele)
        self._is_snp = None  # computed when accessed

    def __repr__(self):
        return f'{type(self).__name__}("{self.allele}")'

    def __str__(self):
        return self.allele

    @property
    def allele(self):
        return self._allele

    @allele.setter
    def allele(self, value):
        self._allele = str(value)
        self._is_snp = None  # reset _is_snp when allele is changed

    @property
    def is_snp(self) -> bool:
        """SNPs are the most common type of effect allele in PGS Catalog scoring
        files. More complex effect alleles, like HLAs or APOE genes, often require
        extra work to represent in genomes. Users should be warned about complex
        effect alleles.
        >>> ea = EffectAllele("+")
        >>> ea.is_snp
        False
        >>> ea.allele = "A"
        >>> ea.is_snp
        True
        """
        if self._is_snp is None:
            self._is_snp = not frozenset(self.allele) - self._valid_snp_bases
        return self._is_snp
