from dataclasses import dataclass


@dataclass
class Config:
    drop_missing: bool
    liftover: bool
    chain_dir: str
    min_lift: float