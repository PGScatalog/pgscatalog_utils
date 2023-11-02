from enum import Enum


class GenomeBuild(Enum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"

    def __str__(self):
        return str(self.value)

    @classmethod
    def from_string(cls, build):
        match build:
            case "GRCh37" | "hg18":
                return cls(GenomeBuild.GRCh37)
            case "GRCh38" | "hg19":
                return cls(GenomeBuild.GRCh38)
            case "NR":
                return None
            case _:
                raise Exception
