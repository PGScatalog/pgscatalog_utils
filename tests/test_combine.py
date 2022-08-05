import pandas as pd


def test_combine_scorefiles(combined_scorefile):
    df = pd.read_table(combined_scorefile)
    cols = {'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type', 'accession'}
    assert set(df.columns).issubset(cols)
    assert df.shape[0] == 51215  # combined number of variants


def test_liftover(lifted_scorefiles):
    df = pd.read_table(lifted_scorefiles)
    assert df.shape[0] > 50000  # approx size

