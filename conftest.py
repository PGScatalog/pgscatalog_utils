import glob
import importlib.resources
import os
import shutil
from unittest.mock import patch

import polars as pl
import pytest
import requests as req

from pgscatalog_utils.download.download_scorefile import download_scorefile
from pgscatalog_utils.match.preprocess import complement_valid_alleles
from pgscatalog_utils.scorefile.combine_scorefiles import combine_scorefiles

from tests.data import combine

pl.toggle_string_cache(True)


@pytest.fixture(scope="session")
def pgs_accessions():
    return ["PGS001229", "PGS000922"]


@pytest.fixture(scope="session")
def mini_score_path(tmp_path_factory):
    path = importlib.resources.files(combine) / "PGS001229_22.txt"
    return str(path)


@pytest.fixture(scope="session")
def mini_scorefile(mini_score_path, tmp_path_factory):
    # The mini scorefile overlaps well with cineca synthetic subset
    out_path = tmp_path_factory.mktemp("scores") / "mini_score.txt"
    args: list[str] = (
        ["combine_scorefiles", "-t", "GRCh37", "-s"]
        + [mini_score_path]
        + ["-o", str(out_path.resolve())]
    )

    with patch("sys.argv", args):
        combine_scorefiles()

    return str(out_path.resolve())


@pytest.fixture(scope="session")
def scorefiles(tmp_path_factory, pgs_accessions):
    fn = tmp_path_factory.mktemp("scorefiles")
    args: list[str] = [
        "download_scorefiles",
        "-b",
        "GRCh37",
        "-o",
        str(fn.resolve()),
        "-i",
    ] + pgs_accessions

    with patch("sys.argv", args):
        download_scorefile()

    return glob.glob(os.path.join(fn.resolve(), "*.txt.gz"))


@pytest.fixture(scope="session")
def target_path(tmp_path_factory):
    try:
        bim = req.get(
            "https://gitlab.ebi.ac.uk/nebfield/test-datasets/-/raw/master/pgsc_calc/cineca_synthetic_subset.bim",
            timeout=5,
        )
    except (req.exceptions.ConnectionError, req.Timeout):
        bim = []

    if not bim:
        pytest.skip("Couldn't get test data from network")
    else:
        fn = tmp_path_factory.mktemp("target") / "data.bim"
        with open(fn, "wb") as f:
            f.write(bim.content)

        return str(fn.resolve())


@pytest.fixture(scope="session")
def chain_files(tmp_path_factory):
    chain_dir = tmp_path_factory.mktemp("chain_dir")

    shutil.copy2("tests/data/hg19ToHg38.over.chain.gz", chain_dir)
    shutil.copy2("tests/data/hg38ToHg19.over.chain.gz", chain_dir)

    return str(chain_dir.resolve())


@pytest.fixture(scope="session")
def hg38_coords():
    rs11903757 = {"rsid": "rs11903757", "chr_name": "2", "chr_position": 191722478}
    rs6061231 = {"rsid": "rs6061231", "chr_name": "20", "chr_position": 62381861}
    return [rs11903757, rs6061231]


@pytest.fixture(scope="session")
def hg19_coords():
    # hg38_coords in GRCh37, from dbSNP
    rs11903757 = {"rsid": "rs11903757", "chr_name": "2", "chr_position": 192587204}
    rs6061231 = {"rsid": "rs6061231", "chr_name": "20", "chr_position": 60956917}
    return [rs11903757, rs6061231]


@pytest.fixture(scope="session")
def small_flipped_scorefile(small_scorefile):
    # simulate a scorefile on the wrong strand
    return (
        complement_valid_alleles(small_scorefile, ["effect_allele", "other_allele"])
        .drop(["effect_allele", "other_allele"])
        .rename(
            {"effect_allele_FLIP": "effect_allele", "other_allele_FLIP": "other_allele"}
        )
        .pipe(complement_valid_alleles, ["effect_allele", "other_allele"])
    )


@pytest.fixture(scope="session")
def small_target():
    return pl.DataFrame(
        {
            "#CHROM": [1, 2, 3],
            "POS": [1, 2, 3],
            "REF": ["A", "T", "T"],
            "ALT": ["C", "A", "G"],
            "ID": ["1:1:A:C", "2:2:T:A", "3:3:T:G"],
            "is_multiallelic": [False, False, False],
        }
    )


@pytest.fixture(scope="session")
def small_scorefile():
    df = pl.DataFrame(
        {
            "accession": ["test", "test", "test"],
            "row_nr": [1, 2, 3],
            "chr_name": [1, 2, 3],
            "chr_position": [1, 2, 3],
            "effect_allele": ["A", "A", "G"],
            "other_allele": ["C", "T", "T"],
            "effect_weight": [1, 2, 3],
            "effect_type": ["additive", "additive", "additive"],
        }
    )

    return complement_valid_alleles(df, ["effect_allele", "other_allele"])


@pytest.fixture(scope="session")
def small_scorefile_no_oa(small_scorefile):
    return small_scorefile.with_column(pl.lit(None).alias("other_allele"))


def _get_timeout(url):
    try:
        return req.get(url, timeout=5)
    except (req.exceptions.ConnectionError, req.Timeout):
        return []
