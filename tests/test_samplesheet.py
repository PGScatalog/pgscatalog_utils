import json
import os
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from pgscatalog_utils.samplesheet.check import check_samplesheet


@pytest.fixture
def existing_vcf_prefix(tmp_path):
    vcf_path = tmp_path / "test.vcf.gz"
    _touch(vcf_path)
    return str(vcf_path.parent.joinpath(Path(vcf_path.stem).stem))


@pytest.fixture
def samplesheet_df(existing_vcf_prefix):
    return pd.DataFrame(
        {"path_prefix": [existing_vcf_prefix], "format": ["vcf"], "sampleset": ["test"], "chrom": [None]})


@pytest.fixture
def good_samplesheet(samplesheet_df, tmp_path):
    path = tmp_path / "good_samplesheet.csv"
    samplesheet_df.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def bad_samplesheet(samplesheet_df, tmp_path):
    path = tmp_path / "bad_samplesheet.csv"
    bad_df = samplesheet_df.copy()
    bad_df['path_prefix'] = 'bad_path'  # path doesn't exist
    bad_df.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def multi_samplesets(samplesheet_df, tmp_path):
    path = tmp_path / "multi_samplesets.csv"
    multi_samplesets = pd.concat([samplesheet_df, samplesheet_df], ignore_index=True)
    multi_samplesets.loc[multi_samplesets.index == 1, 'sampleset'] = 'a_different_name'
    multi_samplesets.to_csv(path, index=False)
    return str(path)


@pytest.fixture
def vcf_dosage(samplesheet_df, tmp_path):
    path = tmp_path / "vcf_dosage.csv"
    dosage_samplesheet = samplesheet_df.copy()
    dosage_samplesheet["vcf_genotype_field"] = ["DS"]
    dosage_samplesheet.to_csv(path, index=False)
    return str(path)


def _touch(fname):
    if os.path.exists(fname):
        os.utime(fname, None)
    else:
        open(fname, 'a').close()


def test_good_samplesheet(good_samplesheet, tmp_path):
    out_path = str(tmp_path / "out.json")
    args = ['samplesheet_to_json', good_samplesheet, out_path]
    with patch('sys.argv', args):
        check_samplesheet()

    assert os.path.exists(out_path), "No file written"


def test_bad_samplesheet(bad_samplesheet, tmp_path):
    out_path = str(tmp_path / "out.json")
    args = ['samplesheet_to_json', bad_samplesheet, out_path]
    with patch('sys.argv', args):
        with pytest.raises(FileNotFoundError):
            check_samplesheet()


def test_multi_samplesets(multi_samplesets, tmp_path):
    out_path = str(tmp_path / "out.json")
    args = ['samplesheet_to_json', multi_samplesets, out_path]
    with patch('sys.argv', args):
        with pytest.raises(Exception, match="Multiple samplesets"):
            check_samplesheet()


def test_dosage_samplesheet(vcf_dosage, tmp_path):
    out_path = str(tmp_path / "out.json")
    args = ['samplesheet_to_json', vcf_dosage, out_path]
    with patch('sys.argv', args):
        check_samplesheet()

    assert os.path.exists(out_path), "Missing output file"

    with open(out_path, 'r') as f:
        converted = json.loads(f.read())
        assert converted[0]['vcf_import_dosage'], "Not importing dosage correctly"
