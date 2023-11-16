"""
This module contains the class ScoreVariant, which is a custom dictionary used to consistently represent rows in a PGS Catalog scoring file
"""
import collections

from pgscatalog_utils.scorefile.effecttype import EffectType


class ScoreVariant(collections.UserDict):
    """A single variant from a scoring file structured to follow PGS Catalog standards,
    typically extracted from a row in a scoring file.

     See https://www.pgscatalog.org/downloads/#dl_scoring_files for field descriptions.

     This class is intentionally simple (a dict that checks for mandatory keys and fills
     optional keys) because a more complicated __init__ will be slow when lots of variants
     are read from a file. dicts use fast C magic, so try not to interfere too much.

     Some additional keys are included for quality control:
     - accession: a unique identifier to group variants in the same score)
     - row_nr: an incrementing integer, used to track the number of variants in an accession
     - is_duplicated: a label to mark variants with the same coordinates and alleles
     - effect_type: additive, recessive, or dominant

     >>> variant = ScoreVariant(**{"chr_name": "1", "chr_position": 1, "effect_allele": "A", "other_allele": "G", "effect_weight": 0.5, "accession": "PGS000822", "row_nr": 0})
     >>> variant
     {'chr_name': '1', 'chr_position': 1, 'effect_allele': 'A', 'other_allele': 'G', 'effect_weight': 0.5, 'accession': 'PGS000822', 'row_nr': 0, 'rsID': None, 'hm_chr': None, 'hm_pos': None, 'hm_inferOtherAllele': None, 'hm_source': None, 'is_dominant': None, 'is_recessive': None, 'hm_rsID': None, 'hm_match_chr': None, 'hm_match_pos': None, 'is_duplicated': None, 'effect_type': <EffectType.ADDITIVE: 'additive'>}

     Mandatory data fields match PGS Catalog harmonised data standards:

    >>> ScoreVariant(**{"chr_name": "1", "chr_position": 1})
    Traceback (most recent call last):
        ...
    ValueError: Mandatory field 'effect_allele' is missing.
    """

    mandatory_fields: tuple[str] = (
        "chr_name",
        "chr_position",
        "effect_allele",
        "effect_weight",
        "accession",
        "row_nr",
    )
    optional_fields: tuple[str] = (
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
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)  # creates the dict

        for field in self.mandatory_fields:
            if field not in self.data:
                raise ValueError(f"Mandatory field '{field}' is missing.")

        # set most optional fields to None...
        for field in self.optional_fields:
            self.data.setdefault(field, None)

        # ... except effect type, as the vast majority of variants are additive
        self.data.setdefault("effect_type", EffectType.ADDITIVE)
