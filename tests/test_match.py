import os
from unittest.mock import patch
import pandas as pd
from pgscatalog_utils.match.match_variants import match_variants


def test_match(combined_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', combined_scorefile,
                       '-t', target_path,
                       '-m', 0.0,
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with patch('sys.argv', args):
        match_variants()

    df = pd.read_table(os.path.join(out_dir, "false_additive_0.scorefile"))

    # PGS000802 doesn't match any variants in the test dataset, so the column is missing
    assert set(df.columns).issubset({'ID', 'effect_allele', 'PGS001229'})
    assert df.shape[0] == 819


