import copy

from pgscatalog_utils.scorefile.config import Config

from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.liftover import liftover, create_liftover


def test_liftover(hg38_coords, hg19_coords, chain_files):
    Config.chain_dir = chain_files
    Config.lo = create_liftover()
    Config.min_lift = 0.95

    hg38 = copy.deepcopy(hg38_coords)
    lifted = list(
        liftover(
            hg38,
            harmonised=False,
            current_build=GenomeBuild.GRCh38,
            target_build=GenomeBuild.GRCh37,
        )
    )

    assert [x["chr_position"] for x in lifted] == [
        x["chr_position"] for x in hg19_coords
    ]
    assert [x["chr_name"] for x in lifted] == [x["chr_name"] for x in hg19_coords]

    hg19 = copy.deepcopy(hg19_coords)
    lift_back = list(
        liftover(
            hg19,
            harmonised=False,
            current_build=GenomeBuild.GRCh37,
            target_build=GenomeBuild.GRCh38,
        )
    )
    assert [x["chr_position"] for x in lift_back] == [
        x["chr_position"] for x in hg38_coords
    ]
    assert [x["chr_name"] for x in lift_back] == [x["chr_name"] for x in hg38_coords]
