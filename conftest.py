import pytest
from unittest.mock import patch

from pgscatalog_utils.match.read import read_target

from pgscatalog_utils.download.download_scorefile import download_scorefile
import os
import requests as req

from pgscatalog_utils.scorefile.combine_scorefiles import combine_scorefiles


@pytest.fixture(scope="session")
def pgs_accessions():
    return ['PGS001229', 'PGS000802']


@pytest.fixture(scope="session")
def scorefiles(tmp_path_factory, pgs_accessions):
    fn = tmp_path_factory.mktemp("scorefiles")
    args: list[str] = ['download_scorefiles', '-o', str(fn.resolve()), '-i'] + pgs_accessions

    with patch('sys.argv', args):
        download_scorefile()

    paths: list[str] = [os.path.join(fn.resolve(), x + '.txt.gz') for x in pgs_accessions]

    assert all([os.path.exists(x) for x in paths])

    return paths


@pytest.fixture(scope="session")
def target_path(tmp_path_factory):
    try:
        bim = req.get('https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/cineca_synthetic_subset.bim', timeout = 5)
    except (req.exceptions.ConnectionError, req.Timeout):
        bim = []

    if not bim:
        pytest.skip("Couldn't get test data from network")
    else:
        fn = tmp_path_factory.mktemp("target") / "data.bim"
        with open(fn, 'wb') as f:
            f.write(bim.content)

        return str(fn.resolve())


@pytest.fixture(scope="session")
def combined_scorefile(scorefiles, tmp_path_factory):
    out_path = tmp_path_factory.mktemp("scores") / "combined.txt"
    args: list[str] = ['combine_scorefiles', '-s'] + scorefiles + ['-o', str(out_path.resolve())]

    with patch('sys.argv', args):
        combine_scorefiles()

    return str(out_path.resolve())

