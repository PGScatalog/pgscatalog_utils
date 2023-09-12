import gzip
import os
from unittest.mock import patch

from pgscatalog_utils.download.Catalog import CatalogQuery, CatalogResult
from pgscatalog_utils.download.CatalogCategory import CatalogCategory
from pgscatalog_utils.download.download_scorefile import download_scorefile


def test_download_scorefile_author(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-o', out_dir]

    with patch('sys.argv', args):
        # Test download
        download_scorefile()
        score_filename = f'{pgs_id}.txt.gz'
        assert score_filename in os.listdir(out_dir)


def test_download_scorefile_hmPOS(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-b', 'GRCh38', '-o', out_dir]

    with patch('sys.argv', args):
        # Test download
        download_scorefile()
        hm_score_filename = f'{pgs_id}_hmPOS_GRCh38.txt.gz'
        assert hm_score_filename in os.listdir(out_dir)


def test_download_existing_scorefile_no_overwrite(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    score_filename = f'{pgs_id}.txt.gz'
    local_score_file_path = f'{out_dir}/{score_filename}'
    # Create fake PGS scoring file
    with open(local_score_file_path, 'w') as file:
        file.write("Download test")
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-o', out_dir]

    with patch('sys.argv', args):
        # Attempt to download scoring file
        download_scorefile()

        # the existing file won't be overwritten even when the checksum fails to validate
        with open(local_score_file_path, 'r') as f:
            assert f.read() == "Download test"


def test_download_overwrite_existing_scorefile(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    score_filename = f'{pgs_id}.txt.gz'
    local_score_file_path = f'{out_dir}/{score_filename}'
    # Create fake PGS scoring file
    with open(local_score_file_path, 'w') as file:
        file.write("Download test")
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-o', out_dir, '-w']

    with patch('sys.argv', args):
        # Download scoring file (over file existing locally)
        download_scorefile()

        with gzip.open(local_score_file_path, 'rt') as f:
            # the test data will be overwritten with real content
            assert f.read() != "Download test"


def test_download_publication(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgp_id = 'PGP000001'
    args: list[str] = ['download_scorefiles', '-p', pgp_id, '-o', out_dir]
    with patch('sys.argv', args):
        download_scorefile()
        assert not {'PGS000001.txt.gz', 'PGS000002.txt.gz', 'PGS000003.txt.gz'}.difference(os.listdir(out_dir))


def test_download_trait(tmp_path):
    out_dir = str(tmp_path.resolve())
    trait_id = 'EFO_0004329'
    args: list[str] = ['download_scorefiles', '-t', trait_id, '-o', out_dir, '-e']
    with patch('sys.argv', args):
        download_scorefile()
        # seven scores in catalog for alcohol drinking, more may be added in the future
        assert not {'PGS001901.txt.gz', 'PGS002115.txt.gz', 'PGS002909.txt.gz',
                    'PGS002910.txt.gz', 'PGS002911.txt.gz', 'PGS002912.txt.gz',
                    'PGS002913.txt.gz'}.difference(os.listdir(out_dir))


def test_query_publication():
    # publications are relatively static
    query: list[CatalogResult] = CatalogQuery(CatalogCategory.PUBLICATION, accession="PGP000001").get()
    assert not query[0].pgs_ids.difference({'PGS000001', 'PGS000002', 'PGS000003'})


def test_query_trait():
    # new scores may be added to traits in the future
    query: list[CatalogResult] = CatalogQuery(CatalogCategory.TRAIT, accession="EFO_0004329").get()
    assert not {'PGS001901', 'PGS002115'}.difference(query[0].pgs_ids)
