import csv
import importlib.resources
import json
from unittest.mock import patch

import pytest

from pgscatalog_utils.scorefile.combine_scorefiles import combine_scorefiles
from tests.data import combine


def test_pgscatalog_combine(pgscatalog_path, tmp_path, combine_output_header):
    out_path = tmp_path / "combined.txt"
    args: list[str] = (
        ["combine_scorefiles", "-t", "GRCh37", "-s"]
        + [str(pgscatalog_path)]
        + ["-o", str(out_path.resolve())]
    )

    with patch("sys.argv", args):
        combine_scorefiles()

    n = -1  # skip header line
    with open(out_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                cols = line.strip().split("\t")
                assert not set(cols).difference(set(combine_output_header))
            n += 1

    with open(out_path.parent / "log_combined.json") as f:
        header = json.load(f)[0]
        assert header["PGS001229_22"]["pgs_id"] == "PGS001229"
        assert header["PGS001229_22"]["pgs_name"] == "GBE_INI50"
        assert header["PGS001229_22"]["genome_build"] == "GRCh37"
        assert int(header["PGS001229_22"]["variants_number"]) == n
        assert not header["PGS001229_22"]["use_harmonised"]


def test_effect_type_combine(effect_type_path, tmp_path, combine_output_header):
    # these genomes are in build GRCh37, so combining with -t GRCh38 will raise an exception
    out_path = tmp_path / "combined.txt"
    args: list[str] = (
        ["combine_scorefiles", "-t", "GRCh37", "-s"]
        + [str(effect_type_path)]
        + ["-o", str(out_path.resolve())]
    )
    with patch("sys.argv", args):
        combine_scorefiles()

    with open(out_path) as f:
        n = 0
        for line in csv.DictReader(f, delimiter="\t"):
            cols = list(line.keys())

            if int(line["row_nr"]) == 0:
                assert line["effect_type"] == "dominant"

            if int(line["row_nr"]) == 1:
                assert line["effect_type"] == "recessive"

            n += 1

        assert not set(cols).difference(set(combine_output_header))

    with open(out_path.parent / "log_combined.json") as f:
        header = json.load(f)[0]
        assert (
            header["scorefile_dominant_and_recessive"]["pgs_name"]
            == "PGS001229_22_DominantRecessiveExample"
        )
        assert header["scorefile_dominant_and_recessive"]["genome_build"] == "GRCh37"
        assert header["scorefile_dominant_and_recessive"]["variants_number"] == n
        assert not header["scorefile_dominant_and_recessive"]["use_harmonised"]


def test_custom_combine(custom_score_path, tmp_path, combine_output_header):
    # these genomes are in build GRCh37, so combining with -t GRCh38 will raise an exception
    out_path = tmp_path / "combined.txt"
    args: list[str] = (
        ["combine_scorefiles", "-t", "GRCh37", "-s"]
        + [str(custom_score_path)]
        + ["-o", str(out_path.resolve())]
    )

    with patch("sys.argv", args):
        combine_scorefiles()

    # read combined file
    n = -1  # skip header line
    with open(out_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                cols = line.strip().split("\t")
                assert not set(cols).difference(set(combine_output_header))
            n += 1

    with open(out_path.parent / "log_combined.json") as f:
        header = json.load(f)[0]
        assert header["scorefile"]["pgs_name"] == "PGS001229_22"
        assert header["scorefile"]["genome_build"] == "GRCh37"
        assert header["scorefile"]["variants_number"] == n
        assert not header["scorefile"]["use_harmonised"]


@pytest.fixture
def pgscatalog_path(scope="session"):
    path = importlib.resources.files(combine) / "PGS001229_22.txt"
    return path


@pytest.fixture
def effect_type_path(scope="session"):
    path = importlib.resources.files(combine) / "scorefile_dominant_and_recessive.txt"
    return path


@pytest.fixture(scope="session")
def custom_score_path(tmp_path_factory):
    path = importlib.resources.files(combine) / "scorefile.txt"
    return path


@pytest.fixture(scope="session")
def combine_output_header():
    return [
        "chr_name",
        "chr_position",
        "effect_allele",
        "other_allele",
        "effect_weight",
        "effect_type",
        "is_duplicated",
        "accession",
        "row_nr",
    ]
