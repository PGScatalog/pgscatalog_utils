from dataclasses import dataclass


@dataclass
class Config:
    threads: int
    drop_missing: bool
    liftover: bool
    chain_dir: str
    min_lift: float
    batch_size: int