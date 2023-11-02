import csv
import logging
import os
import pathlib
import typing
from dataclasses import dataclass
from itertools import islice

from pgscatalog_utils.scorefile.config import Config

from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.header import ScoringFileHeader, auto_open
from pgscatalog_utils.scorefile.qc import quality_control

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class ScoringFile:
    path: pathlib.Path
    accession: str
    header: typing.Union[ScoringFileHeader, None]
    genome_build: typing.Union[GenomeBuild, None]
    harmonised: bool
    fields: list[str]
    variants: typing.Generator

    def __post_init__(self):
        if self.header.HmPOS_build:
            logger.info(
                f"{self.path} harmonised data detected: {self.header.HmPOS_build}"
            )
            self.genome_build = self.header.HmPOS_build

        mandatory_columns = {"chr_name", "effect_allele", "effect_weight"}
        if not mandatory_columns.issubset(self.fields) not in self.fields:
            err_msg = f"{self.path} missing fields"
            raise Exception(err_msg)

    @classmethod
    def from_path(cls, path: pathlib.Path):
        header = ScoringFileHeader.from_path(path)
        name = os.path.basename(path).split(".")[0]
        if header:
            if header.HmPOS_build:
                harmonised = True
                genome_build = header.HmPOS_build
            else:
                harmonised = False
                genome_build = header.genome_build
        else:
            harmonised = False
            genome_build = None

        start_line, cols = get_columns(path)
        is_wide = detect_wide(cols)

        logger.info(f"Lazily reading variants from {path}")
        variants = ScoringFile.read_variants(
            path=path, start_line=start_line, fields=cols, name=name, is_wide=is_wide
        )

        # note: these generator expressions aren't doing a bunch of iterations
        # it's just a data processing pipeline
        variants = quality_control(
            variants, header=header, harmonised=harmonised, wide=is_wide
        )

        return cls(
            path=path,
            header=header,
            genome_build=genome_build,
            harmonised=harmonised,
            fields=cols,
            variants=variants,
            accession=name,
        )

    @staticmethod
    def read_variants(path, fields, start_line, name: str, is_wide: bool):
        open_function = auto_open(path)
        row_nr = 0
        with open_function(path, mode="rt") as f:
            for _ in range(start_line + 1):
                # skip header
                next(f)

            while True:
                batch = list(islice(f, Config.batch_size))
                if not batch:
                    break

                csv_reader = csv.reader(batch, delimiter="\t")
                yield from read_rows(csv_reader, fields, name, row_nr, is_wide)


def read_rows(csv_reader, fields: list[str], name: str, row_nr: int, wide: bool):
    for row in csv_reader:
        variant = dict(zip(fields, row))

        if wide:
            ew_col_idxs: list[int] = [
                i for i, x in enumerate(["effect_weight_" in x for x in fields]) if x
            ]
            for i, weight_name in zip(ew_col_idxs, [fields[i] for i in ew_col_idxs]):
                keys = ["chr_name", "chr_position", "effect_allele", "other_allele"]
                yield {k: variant[k] for k in keys if k in variant} | {
                    "accession": weight_name,
                    "row_nr": row_nr,
                    "effect_weight": variant[weight_name],
                }
        else:
            keys = [
                "chr_name",
                "chr_position",
                "effect_allele",
                "other_allele",
                "effect_weight",
                "hm_chr",
                "hm_pos",
                "hm_inferOtherAllele",
                "is_dominant",
                "is_recessive",
                "accession",
                "row_nr",
            ]

            yield {k: variant[k] for k in keys if k in variant} | {
                "accession": name,
                "row_nr": row_nr,
            }

        row_nr += 1


def get_columns(path) -> tuple[int, list[str]]:
    open_function = auto_open(path)
    with open_function(path, mode="rt") as f:
        for i, line in enumerate(f):
            if line.startswith("#"):
                continue
            line_no, cols = i, line.strip().split("\t")
            if len(set(cols)) != len(cols):
                logger.critical(f"Duplicated column names: {cols}")
                raise ValueError

            return line_no, cols


def detect_wide(cols: list[str]) -> bool:
    """
    Check columns to see if multiple effect weights are present. Multiple effect weights must be present in the form:
    effect_weight_suffix1
    effect_weight_suffix2
    """
    if any(["effect_weight_" in x for x in cols]):
        logger.info("Wide scoring file detected with multiple effect weights")
        return True
    else:
        return False
