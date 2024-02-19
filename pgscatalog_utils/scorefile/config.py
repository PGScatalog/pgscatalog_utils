from dataclasses import dataclass

import pyliftover

from pgscatalog_utils.download.GenomeBuild import GenomeBuild


@dataclass
class Config:
    drop_missing: bool
    liftover: bool
    lo: pyliftover.liftover
    chain_dir: str
    min_lift: float
    batch_size: int
    target_build: GenomeBuild
