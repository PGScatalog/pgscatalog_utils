import pandas as pd

from pgscatalog_utils.scorefile.liftover import liftover


def test_liftover(hg38_coords, hg19_coords, chain_files):
    lifted = liftover(hg38_coords, chain_files, min_lift=0.9, target_build='GRCh37')
    coords: pd.DataFrame = hg19_coords[['lifted_pos', 'lifted_chr']] == lifted[['lifted_pos', 'lifted_chr']]
    assert coords.all(axis=None)
