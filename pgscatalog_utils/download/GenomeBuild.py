from enum import Enum


class GenomeBuild(Enum):
    GRCh37 = "GRCh37"
    GRCh38 = "GRCh38"
    # just included to handle older files, incompatible unless harmonised:
    NCBI36 = "NCBI36"  # ew

    def __str__(self):
        return str(self.value)

    @classmethod
    def from_string(cls, build):
        match build:
            case "GRCh37" | "hg19":
                return cls(GenomeBuild.GRCh37)
            case "GRCh38" | "hg38":
                return cls(GenomeBuild.GRCh38)
            case "NR":
                return None
            case "NCBI36" | "hg18":
                return cls(GenomeBuild.NCBI36)
            case _:
                raise Exception(f"Can't match {build=}")
