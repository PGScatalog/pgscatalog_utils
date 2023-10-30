from enum import Enum, auto


class GenomeBuild(Enum):
    GRCh37 = auto()
    GRCh38 = auto()

    @classmethod
    def from_string(cls, build):
        match build:
            case 'GRCh37' | 'hg18':
                return cls(GenomeBuild.GRCh37)
            case 'GRCh38' | 'hg19':
                return cls(GenomeBuild.GRCh38)
            case 'NR':
                return None
            case _:
                raise Exception