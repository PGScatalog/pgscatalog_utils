import logging
import os

import pyliftover

from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.config import Config

logger = logging.getLogger(__name__)


def liftover(
    variants, harmonised: bool, current_build: GenomeBuild, target_build: GenomeBuild
):
    if harmonised:
        skip_lo = True
    elif target_build == current_build:
        skip_lo = True
    else:
        skip_lo = False

    if skip_lo:
        for variant in variants:
            yield variant
    else:
        if current_build == GenomeBuild.GRCh37 and target_build == GenomeBuild.GRCh38:
            lo: pyliftover.LiftOver = Config.lo["hg19hg38"]
        elif current_build == GenomeBuild.GRCh38 and target_build == GenomeBuild.GRCh37:
            lo: pyliftover.LiftOver = Config.lo["hg19hg38"]
        else:
            raise Exception("Can't get pyliftover object")

        n_lifted = 0
        n = 0

        for variant in variants:
            chrom = "chr" + variant["chr_name"]
            pos = int(variant["chr_position"]) - 1  # VCF -> 1 based, UCSC -> 0 based
            lifted = lo.convert_coordinate(chrom, pos)
            if lifted:
                variant["chr_name"] = lifted[0][0][3:].split("_")[0]
                variant["chr_position"] = lifted[0][1] + 1  # reverse 0 indexing
                n_lifted += 1
            yield variant
            n += 1

        if (n_lifted / n) < Config.min_lift:
            logger.error(f"Liftover failed")
            raise Exception
        else:
            logger.info("Liftover successful")


def create_liftover() -> dict["str" : pyliftover.LiftOver]:
    """Create LiftOver objects that can remap genomic coordinates"""
    chain_dir: str = Config.chain_dir
    builds: list[str] = ["hg19hg38", "hg38hg19"]
    chains: list[str] = [
        os.path.join(chain_dir, x)
        for x in ["hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz"]
    ]
    lo: list[pyliftover.LiftOver] = [pyliftover.LiftOver(x) for x in chains]
    logger.debug("Chain files loaded for liftover")
    return dict(zip(builds, lo))
