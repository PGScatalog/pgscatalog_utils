import os
from unittest.mock import patch

import pytest

from pgscatalog_utils.download.download_scorefile import download_scorefile,generate_md5_checksum,get_md5_checksum_from_ftp
from pgscatalog_utils.download.publication import query_publication
from pgscatalog_utils.download.score import get_url
from pgscatalog_utils.download.trait import query_trait


ftp_url_root = 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores'

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


def test_download_scorefile_author(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-o', out_dir]

    with patch('sys.argv', args):
        # Test download
        download_scorefile()
        score_filename = f'{pgs_id}.txt.gz'
        assert os.listdir(out_dir) == [score_filename]
        # Test MD5
        ftp_md5 = get_md5_checksum_from_ftp(f'{ftp_url_root}/{pgs_id}/ScoringFiles/{score_filename}')
        downloaded_file_md5 = generate_md5_checksum(f'{out_dir}/{score_filename}')
        assert ftp_md5 == downloaded_file_md5


def test_download_scorefile_hmPOS(tmp_path):
    out_dir = str(tmp_path.resolve())
    pgs_id = 'PGS000001'
    args: list[str] = ['download_scorefiles', '-i', pgs_id, '-b', 'GRCh38', '-o', out_dir]

    with patch('sys.argv', args):
        # Test download
        download_scorefile()
        hm_score_filename = f'{pgs_id}_hmPOS_GRCh38.txt.gz'
        assert os.listdir(out_dir) == [hm_score_filename]
        # Test MD5
        ftp_md5 = get_md5_checksum_from_ftp(f'{ftp_url_root}/{pgs_id}/ScoringFiles/Harmonized/{hm_score_filename}')
        downloaded_file_md5 = generate_md5_checksum(f'{out_dir}/{hm_score_filename}')
        assert ftp_md5 == downloaded_file_md5


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
        # Test MD5
        ftp_md5 = get_md5_checksum_from_ftp(f'{ftp_url_root}/{pgs_id}/ScoringFiles/{score_filename}')
        downloaded_file_md5 = generate_md5_checksum(local_score_file_path)
        assert ftp_md5 != downloaded_file_md5


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
        # Test MD5
        ftp_md5 = get_md5_checksum_from_ftp(f'{ftp_url_root}/{pgs_id}/ScoringFiles/{score_filename}')
        downloaded_file_md5 = generate_md5_checksum(local_score_file_path)
        assert ftp_md5 == downloaded_file_md5


def test_query_publication():
    # publications are relatively static
    assert not set(query_publication("PGP000001")).difference(['PGS000001', 'PGS000002', 'PGS000003'])


def test_query_trait():
    # new scores may be added to traits in the future
    assert {'PGS001901', 'PGS002115'}.issubset(set(query_trait("EFO_0004329")))
