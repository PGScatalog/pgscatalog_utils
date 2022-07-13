import os
import pytest
from unittest.mock import patch
from pgscatalog_utils.download.api import pgscatalog_result
from pgscatalog_utils.download.download_scorefile import download_scorefile


@pytest.fixture(params=[["PGS000001"], ["PGS000001", "PGS000802"]])
def pgscatalog_api(request):
    return pgscatalog_result(request.param)


def test_pgscatalog_result(pgscatalog_api):
    # is the key a list of accessions?
    for k, v in pgscatalog_api.items():
        # make sure the key is a PGS ID
        assert k.startswith("PGS")
        # ensure sure ftp prefix is correctly used to avoid connection problems
        assert v.startswith("ftp://")
        # make sure the URL points to a reasonable looking compressed text file
        assert v.endswith(".txt.gz")


def test_download_scorefile(tmp_path):
    out_dir = str(tmp_path.resolve())
    args: list[str] = ['download_scorefiles', '-i', 'PGS000001', '-o', out_dir]

    with patch('sys.argv', args):
        download_scorefile()
        assert os.listdir(out_dir) == ['PGS000001.txt.gz']

