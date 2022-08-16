import pytest
from unittest.mock import patch
from pgscatalog_utils.download.download_scorefile import download_scorefile
import os
import requests as req
from pgscatalog_utils.scorefile.combine_scorefiles import combine_scorefiles
from pysqlar import SQLiteArchive
import pandas as pd
import glob


@pytest.fixture(scope="session")
def pgs_accessions():
    return ['PGS001229', 'PGS000922']


@pytest.fixture(scope="session")
def scorefiles(tmp_path_factory, pgs_accessions):
    fn = tmp_path_factory.mktemp("scorefiles")
    args: list[str] = ['download_scorefiles', '-b', 'GRCh37', '-o', str(fn.resolve()), '-i'] + pgs_accessions

    with patch('sys.argv', args):
        download_scorefile()

    return glob.glob(os.path.join(fn.resolve(), "*.txt.gz"))


@pytest.fixture(scope="session")
def target_path(tmp_path_factory):
    try:
        bim = req.get(
            'https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/cineca_synthetic_subset.bim',
            timeout=5)
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
def mini_score_path(tmp_path_factory):
    try:
        score = req.get('https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/PGS001229_22.txt',
                        timeout=5)
    except (req.exceptions.ConnectionError, req.Timeout):
        score = []

    if not score:
        pytest.skip("Couldn't get test data from network")
    else:
        fn = tmp_path_factory.mktemp("score") / "PGS001229_22.txt"
        with open(fn, 'wb') as f:
            f.write(score.content)

        return str(fn.resolve())


@pytest.fixture(scope="session")
def mini_scorefile(mini_score_path, tmp_path_factory):
    # The mini scorefile overlaps well with cineca synthetic subset
    out_path = tmp_path_factory.mktemp("scores") / "mini_score.txt"
    args: list[str] = ['combine_scorefiles', '-s'] + [mini_score_path] + ['-o', str(out_path.resolve())]

    with patch('sys.argv', args):
        combine_scorefiles()

    return str(out_path.resolve())


@pytest.fixture(scope="session")
def combined_scorefile(scorefiles, tmp_path_factory):
    # The combined scorefile overlaps poorly with cineca synthetic subset
    out_path = tmp_path_factory.mktemp("scores") / "combined.txt"
    args: list[str] = ['combine_scorefiles', '-s'] + scorefiles + ['-o', str(out_path.resolve())]

    with patch('sys.argv', args):
        combine_scorefiles()

    return str(out_path.resolve())


@pytest.fixture(scope="session")
def db(tmp_path_factory):
    database = _get_timeout(
        'https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/reference_data/pgsc_calc_ref.sqlar')
    db_path = tmp_path_factory.mktemp('database') / 'db.sqlar'
    if not database:
        pytest.skip("Couldn't get file from remote host")
    else:
        with open(db_path, 'wb') as f:
            f.write(database.content)
        return str(db_path.resolve())


@pytest.fixture(scope="session")
def chain_files(db, tmp_path_factory):
    chain_dir = tmp_path_factory.mktemp('chain_dir')

    with SQLiteArchive(db, 'ro') as ar:
        ar.extract("hg38ToHg19.over.chain.gz", path=chain_dir)
        ar.extract("hg19ToHg38.over.chain.gz", path=chain_dir)

    return str(chain_dir.resolve())


@pytest.fixture(scope="session")
def lifted_scorefiles(scorefiles, chain_files, tmp_path_factory):
    out_path = tmp_path_factory.mktemp("scores") / "lifted.txt"
    args: list[str] = ['combine_scorefiles', '-s'] + scorefiles + ['--liftover', '-c', chain_files, '-t', 'GRCh38',
                                                                   '-m', '0.8'] + ['-o', str(out_path.resolve())]

    with patch('sys.argv', args):
        combine_scorefiles()

    return str(out_path.resolve())


@pytest.fixture(scope="session")
def hg38_coords(tmp_path_factory):
    out_path = tmp_path_factory.mktemp("dummy") / "hg38.txt"
    d = {'rsid': ['rs11903757', 'rs6061231'], 'chr_name': ['2', '20'], 'chr_position': [191722478, 62381861]}
    df = pd.DataFrame(d)
    with open(out_path, 'w') as f:
        f.write('#genome_build=GRCh38\n')
    df.to_csv(out_path, mode='a', index=False)
    df['filename'] = str(out_path.resolve())
    df['accession'] = 'dummy'
    return df


@pytest.fixture(scope="session")
def hg19_coords(hg38_coords):
    # hg38_coords in GRCh37, from dbSNP
    d = {'lifted_chr': ['2', '20'], 'lifted_pos': [192587204, 60956917], 'liftover': [True, True]}
    return pd.DataFrame(d)


def _get_timeout(url):
    try:
        return req.get(url, timeout=5)
    except (req.exceptions.ConnectionError, req.Timeout):
        return []
