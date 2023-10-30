import csv
import gzip
import logging
import os
import pathlib
import typing
from dataclasses import dataclass

from pgscatalog_utils.download.GenomeBuild import GenomeBuild
from pgscatalog_utils.scorefile.header import ScoringFileHeader, auto_open
from pgscatalog_utils.scorefile.qc import drop_hla, check_effect_weight, \
    assign_other_allele, assign_effect_type, remap_harmonised

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class ScoringFile:
    path: pathlib.Path
    name: str
    header: typing.Union[ScoringFileHeader, None]
    genome_build: typing.Union[GenomeBuild, None]
    harmonised: bool
    fields: list[str]
    variants: typing.Generator

    def __post_init__(self):
        if self.header.HmPOS_build:
            logger.info(
                f"{self.path} harmonised data detected: {self.header.HmPOS_build}")
            self.genome_build = self.header.HmPOS_build

        mandatory_columns = {'chr_name', 'effect_allele', 'effect_weight'}
        if not mandatory_columns.issubset(self.fields) not in self.fields:
            err_msg = f"{self.path} missing fields"
            raise Exception(err_msg)

    @classmethod
    def from_path(cls, path: pathlib.Path):
        header = ScoringFileHeader.from_path(path)
        if header:
            name = header.pgs_id
            if header.HmPOS_build:
                harmonised = True
                genome_build = header.HmPOS_build
            else:
                harmonised = False
                genome_build = header.genome_build
        else:
            harmonised = False
            genome_build = None
            name = os.path.basename(path).split('.')[0]

        start_line, cols = get_columns(path)

        # generate variants (a list of dicts, one for each variants)
        logger.info(f"Lazily reading variants from {path}")
        variants = ScoringFile.read_variants(path=path, start_line=start_line,
                                             fields=cols, name=name)

        # note: these generator expressions aren't doing a bunch of iterations
        # it's just a data processing pipeline
        variants = remap_harmonised(variants, harmonised)
        variants = drop_hla(variants)
        variants = assign_effect_type(variants)
        variants = check_effect_weight(variants)
        variants = assign_other_allele(variants)

        return cls(path=path, header=header, genome_build=genome_build,
                   harmonised=harmonised,
                   fields=cols,
                   variants=variants,
                   name=name)

    @staticmethod
    def read_variants(path, fields, start_line, name: str):
        open_function = auto_open(path)
        with open_function(path, 'rt') as f:
            csv_reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(csv_reader):
                if i > start_line:
                    variant = dict(zip(fields, row)) | {'name': name}
                    keys = ["chr_name", "chr_position", "effect_allele", "other_allele",
                            "effect_weight", "hm_chr", "hm_pos", "hm_inferOtherAllele",
                            "name", "is_dominant", "is_recessive"]
                    yield {k: variant[k] for k in keys if k in variant}

    @staticmethod
    def write_combined(scoring_files, out_path):
        if out_path.endswith("gz"):
            open_function = gzip.open
        else:
            open_function = open

        with open_function(out_path, 'wt') as f:
            fieldnames = ["name", "chr_name", "chr_position", "effect_allele",
                          "other_allele", "effect_weight", "effect_type"]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            # write out in chunks for compression efficiency and speed
            chunk_size = 10000
            chunk = []
            for scoring_file in scoring_files:
                logger.info(f"Writing {scoring_file.name} variants")
                for variant in scoring_file.variants:
                    chunk.append(variant)
                    if len(chunk) == chunk_size:
                        writer.writerows(chunk)
                        chunk = []
                # handle last chunk
                if chunk:
                    writer.writerows(chunk)


def get_columns(path) -> tuple[int, list[str]]:
    open_function = auto_open(path)
    with open_function(path, 'rt') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            return i, line.strip().split('\t')


