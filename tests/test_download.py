import os
import pytest
from unittest.mock import patch

from pgscatalog_utils.download.trait import query_trait
from pgscatalog_utils.download.publication import query_publication
from pgscatalog_utils.download.score import get_url
from pgscatalog_utils.download.download_scorefile import download_scorefile


@pytest.fixture(params=[["PGS000001"], ["PGS000001", "PGS000802"]])
def pgscatalog_api(request):
    return get_url(request.param, "GRCh37")


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
    args: list[str] = ['download_scorefiles', '-i', 'PGS000001', '-b', 'GRCh38', '-o', out_dir]

    with patch('sys.argv', args):
        download_scorefile()
        assert os.listdir(out_dir) == ['PGS000001.txt.gz']


def test_query_publication():
    # publications are relatively static
    assert not set(query_publication("PGP000001")).difference(['PGS000001', 'PGS000002', 'PGS000003'])


def test_query_trait():
    # new scores may be added to traits in the future
    assert {'PGS001901', 'PGS002115'}.issubset(set(query_trait("EFO_0004329")))
