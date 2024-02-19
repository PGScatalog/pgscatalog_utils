import csv
import logging
import os
import pathlib
import typing
from dataclasses import dataclass
from itertools import islice

from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.scoringfileheader import ScoringFileHeader, auto_open
from pgscatalog_utils.scorefile.qc import quality_control
from pgscatalog_utils.scorefile.scorevariant import ScoreVariant

logger = logging.getLogger(__name__)


@dataclass
class ScoringFile:
    path: pathlib.Path
    accession: str
    header: typing.Union[ScoringFileHeader, None]
    genome_build: typing.Union[GenomeBuild, None]
    harmonised: bool
    fields: list[str]
    variants: typing.Generator[ScoreVariant, None, None]

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
        variants: typing.Generator[
            ScoreVariant, None, None
        ] = ScoringFile.read_variants(
            path=path, start_line=start_line, fields=cols, name=name, is_wide=is_wide
        )

        # the quality_control function normalises a list of variants to have a standard representation
        # attributes are overwritten using harmonised data, etc.
        variants: typing.Generator[ScoreVariant, None, None] = quality_control(
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

    def generate_log(self, counted: typing.Counter):
        log = {
            key: str(value) if value is not None else None
            for key, value in self.header.__dict__.items()
        }

        if log["variants_number"] is None:
            # custom scoring files might not have this information
            log["variants_number"] = counted["n_variants"]

        if (
            int(log["variants_number"]) != counted["n_variants"]
            and not Config.drop_missing
        ):
            logger.warning(
                f"Mismatch between header ({log['variants_number']}) and output row count ({counted['n_variants']}) for {self.accession}"
            )
            logger.warning(
                "This can happen with older scoring files in the PGS Catalog (e.g. PGS000028)"
            )

        # multiple terms may be separated with a pipe
        if log["trait_mapped"]:
            log["trait_mapped"] = log["trait_mapped"].split("|")

        if log["trait_efo"]:
            log["trait_efo"] = log["trait_efo"].split("|")

        log["columns"] = self.fields
        log["use_liftover"] = Config.liftover
        log["use_harmonised"] = self.harmonised
        log["sources"] = [k for k, v in counted.items() if k != "n_variants"]

        return {self.accession: log}

    @staticmethod
    def read_variants(
        path, fields, start_line, name: str, is_wide: bool
    ) -> typing.Generator[ScoreVariant, None, None]:
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
                yield from read_rows(csv_reader, fields, name, is_wide, row_nr)
                # this is important because row_nr resets for each batch
                row_nr += len(batch)


def read_rows(
    csv_reader, fields: list[str], name: str, wide: bool, row_nr: int
) -> typing.Generator[ScoreVariant, None, None]:
    for row in csv_reader:
        variant = dict(zip(fields, row))

        if wide:
            ew_col_idxs: list[int] = [
                i for i, x in enumerate(["effect_weight_" in x for x in fields]) if x
            ]
            for i, weight_name in zip(ew_col_idxs, [fields[i] for i in ew_col_idxs]):
                yield ScoreVariant(
                    **variant,
                    **{
                        "accession": weight_name,
                        "row_nr": row_nr,
                        "effect_weight": variant[weight_name],
                    },
                )
        else:
            yield ScoreVariant(**variant, **{"accession": name, "row_nr": row_nr})

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
