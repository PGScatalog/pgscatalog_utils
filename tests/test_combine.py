import pandas as pd
import requests as req
import pytest
from pysqlar import SQLiteArchive


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    database = get_timeout('https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/reference_data/pgsc_calc_ref.sqlar')
    db_path = tmp_path_factory.mktemp('database') / 'db.sqlar'
    if not database:
        pytest.skip("Couldn't get file from remote host")
    else:
        with open(db_path, 'wb') as f:
            f.write(database.content)
        return str(db_path.resolve())


@pytest.fixture(scope="module")
def chain_files(db, tmp_path_factory):
    chain_dir = tmp_path_factory.mktemp('chain_dir')

    # temporarily broken
    # with SQLiteArchive('/Users/bwingfield/Downloads/pgsc_calc_ref.sqlar', 'ro') as ar:
    #     ar.extract("hg38ToHg19.over.chain.gz", path=chain_dir)
    #     ar.extract("hg19ToHg38.over.chain.gz", path=chain_dir)
    #
    # return str(chain_dir.resolve())


def get_timeout(url):
    try:
        return req.get(url, timeout = 5)
    except (req.exceptions.ConnectionError, req.Timeout):
        return []


def test_combine_scorefiles(combined_scorefile):
    df = pd.read_table(combined_scorefile)
    cols = {'chr_name', 'chr_position', 'effect_allele', 'other_allele', 'effect_weight', 'effect_type', 'accession'}
    assert set(df.columns).issubset(cols)
    assert df.shape[0] == 51224  # combined number of variants


