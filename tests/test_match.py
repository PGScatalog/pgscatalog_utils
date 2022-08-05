import os
from unittest.mock import patch
import pandas as pd
import pytest

from pgscatalog_utils.match.match_variants import match_variants


def test_match_fail(combined_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', combined_scorefile,
                       '-t', target_path,
                       '-m', 0,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with pytest.raises(Exception):
        with patch('sys.argv', args):
            match_variants()


def test_match_pass(mini_scorefile, target_path, tmp_path):
    out_dir = str(tmp_path.resolve())

    args: list[str] = ['match_variants', '-s', mini_scorefile,
                       '-t', target_path,
                       '-m', 0,
                       '-d', 'test',
                       '--outdir', out_dir,
                       '--keep_ambiguous', '--keep_multiallelic']

    with patch('sys.argv', args):
        match_variants()

