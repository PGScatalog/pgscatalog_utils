import csv
import functools
import gzip
import logging
import os
import sqlite3
import typing
from collections import Counter
from itertools import islice

from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.scorevariant import ScoreVariant
from pgscatalog_utils.scorefile.scoringfile import ScoringFile

logger = logging.getLogger(__name__)


class DataWriter:
    def __init__(self, filename):
        self.filename = filename
        self.fieldnames = [
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
        logger.info(f"Output filename: {filename}")

    def write(self, batch):
        pass


class TextFileWriter(DataWriter):
    def __init__(self, compress, filename):
        super().__init__(filename)
        self.compress = compress

        if self.compress:
            logger.info("Writing with gzip")
            self.open_function = functools.partial(gzip.open, compresslevel=6)
        else:
            logger.info("Writing text file")
            self.open_function = open

    def write(self, batch):
        mode = "at" if os.path.exists(self.filename) else "wt"
        with self.open_function(self.filename, mode) as f:
            writer = csv.writer(
                f,
                delimiter="\t",
                lineterminator="\n",
            )
            if mode == "wt":
                writer.writerow(ScoreVariant.output_fields)

            writer.writerows(batch)


class SqliteWriter(DataWriter):
    def __init__(self, filename):
        super().__init__(filename)

    def write(self, batch):
        conn = sqlite3.connect(self.filename)
        cursor = conn.cursor()
        placeholders = ", ".join("?" for _ in self.fieldnames)

        values = [
            tuple(row[key] for key in self.fieldnames if key in row) for row in batch
        ]

        cursor.execute(
            f"CREATE TABLE IF NOT EXISTS variants ({', '.join(self.fieldnames)})"
        )
        cursor.executemany(f"INSERT INTO variants VALUES ({placeholders})", values)
        conn.commit()
        conn.close()


def write_combined(
    scoring_files: list[ScoringFile], out_path: str
) -> dict[str : typing.Counter]:
    # compresslevel can be really slow, default is 9
    if out_path.endswith("gz"):
        writer = TextFileWriter(compress=True, filename=out_path)
    elif out_path.endswith("txt"):
        writer = TextFileWriter(compress=False, filename=out_path)
    elif out_path.endswith(".sqlite"):
        writer = SqliteWriter(filename=out_path)
    else:
        raise Exception("Can't configure writer, please check out_path")

    counts = []
    log = {}
    for scoring_file in scoring_files:
        logger.info(f"Writing {scoring_file.accession} variants")
        while True:
            batch = list(islice(scoring_file.variants, Config.batch_size))
            if not batch:
                break
            writer.write(batch=batch)
            counts = calculate_log(batch, counts)

        log[scoring_file.accession] = sum(counts, Counter())
        counts = []

    return log


def calculate_log(batch: list[ScoreVariant], log: list[Counter]) -> list[Counter]:
    # these statistics can only be generated while iterating through variants
    n_variants = Counter("n_variants" for item in batch)
    hm_source = Counter(getattr(item, "hm_source") for item in batch)
    log.extend([n_variants + hm_source])
    return log
