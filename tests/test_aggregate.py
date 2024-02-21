import importlib.resources
import os
import pandas as pd
from unittest.mock import patch

from pgscatalog_utils.aggregate.aggregate_scores import aggregate_scores
from . import data


def test_aggregate(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("aggregated")
    score_path = importlib.resources.files(data) / "cineca_22_additive_0.sscore.zst"

    args = ["aggregate_scores", "-s", str(score_path), "-o", str(out_dir)]

    with patch("sys.argv", args):
        aggregate_scores()

    assert os.listdir(out_dir) == ["aggregated_scores.txt.gz"]
    df = pd.read_csv(out_dir / "aggregated_scores.txt.gz", delimiter="\t")
    assert list(df.columns) == ["sampleset", "IID", "PGS", "SUM", "DENOM", "AVG"]
    assert df.shape == (2504, 6)
