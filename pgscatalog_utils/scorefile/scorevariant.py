from pgscatalog_utils.scorefile.effectallele import EffectAllele
from pgscatalog_utils.scorefile.effecttype import EffectType


class ScoreVariant:
    mandatory_fields: tuple[str] = (
        "effect_allele",
        "effect_weight",
        "accession",
        "row_nr",
    )
    optional_fields: tuple[str] = (
        "chr_name",
        "chr_position",
        "rsID",
        "other_allele",
        "hm_chr",
        "hm_pos",
        "hm_inferOtherAllele",
        "hm_source",
        "is_dominant",
        "is_recessive",
        "hm_rsID",
        "hm_match_chr",
        "hm_match_pos",
        "is_duplicated",
        "effect_type",
    )

    complex_fields: tuple[str] = ("is_haplotype", "is_diplotype", "is_interaction")

    # column names for output are used by __iter__ and when writing out
    output_fields: tuple[str] = (
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
        "effect_type",
        "is_duplicated",
        "accession",
        "row_nr",
    )

    # slots uses magic to improve speed and memory when making millions of objects
    __slots__ = mandatory_fields + optional_fields + ("is_complex",)

    # __init__ is intentionally verbose and avoids using loops or trickery to work:
    #   - attributes won't change often
    #   - class accepts keyword parameters only to init (not positional)
    #   - type hints are helpful in parameters
    #   - setting sensible defaults for optional fields is clear
    #   - being verbose helps prevent IDE warnings
    # extra kwargs are silently ignored
    # (yes, effect_weight is treated as a str, want to avoid rounding errors at this stage)
    def __init__(
        self,
        *,
        effect_allele: str,
        effect_weight: str,
        accession: str,
        row_nr: int,
        chr_name: str = None,
        chr_position: int = None,
        rsID: str = None,
        other_allele: str = None,
        hm_chr: str = None,
        hm_pos: int = None,
        hm_inferOtherAllele: str = None,
        hm_source: str = None,
        is_dominant: str = None,
        is_recessive: str = None,
        hm_rsID: str = None,
        hm_match_chr: str = None,
        hm_match_pos: str = None,
        is_duplicated: bool = False,
        effect_type: EffectType = EffectType.ADDITIVE,
        is_complex: bool = False,
        **kwargs,
    ):
        # start with mandatory attributes
        self.effect_allele: EffectAllele = EffectAllele(effect_allele)
        self.effect_weight: str = effect_weight
        self.accession = accession
        self.row_nr = row_nr

        # now set optional fields
        self.chr_name = chr_name
        self.chr_position = chr_position
        self.rsID = rsID
        self.other_allele = other_allele
        self.hm_chr = hm_chr
        self.hm_pos = hm_pos
        self.hm_inferOtherAllele = hm_inferOtherAllele
        self.hm_source = hm_source
        self.is_dominant = is_dominant
        self.is_recessive = is_recessive
        self.hm_rsID = hm_rsID
        self.hm_match_chr = hm_match_chr
        self.hm_match_pos = hm_match_pos
        self.is_duplicated = is_duplicated
        self.effect_type = effect_type

        # these fields are important to check if variants are complex
        if any([x in kwargs for x in self.complex_fields]):
            is_complex = True
        self.is_complex = is_complex

    def __repr__(self):
        class_name = type(self).__name__
        values = {}

        for key in ScoreVariant.__slots__:
            values[key] = getattr(self, key, None)

        return f"{class_name}({values})"

    def __iter__(self):
        for attr in self.output_fields:
            yield getattr(self, attr)
