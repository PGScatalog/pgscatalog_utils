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

try:
    import pyarrow as pa

    PYARROW_AVAILABLE = True
except ImportError:
    PYARROW_AVAILABLE = False

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


class PyarrowWriter(DataWriter):
    if PYARROW_AVAILABLE:
        schema = pa.schema(
            [
                pa.field("chr_name", pa.string()),
                pa.field("chr_position", pa.uint64()),
                pa.field("effect_allele", pa.string()),
                pa.field("other_allele", pa.string()),
                pa.field("effect_weight", pa.string()),
                pa.field("effect_type", pa.string()),
                pa.field("is_duplicated", pa.bool_()),
                pa.field("accession", pa.string()),
                pa.field("row_nr", pa.uint64()),
            ]
        )

    def __init__(self, filename):
        if not PYARROW_AVAILABLE:
            # TODO: provide a pip command
            raise ImportError(
                "pyarrow output not available, please install pyarrow as listed in the pyproject.toml extras section"
            )
        super().__init__(filename)

        self._sink = pa.OSFile(self.filename, "wb")
        self._writer: pa.RecordBatchFileWriter = pa.ipc.new_file(
            self._sink, self.schema
        )

    def write(self, batch: list[ScoreVariant]):
        batch_dict = {
            "chr_name": [x.chr_name for x in batch],
            "chr_position": [x.chr_position for x in batch],
            "effect_allele": [str(x.effect_allele) for x in batch],
            "other_allele": [x.other_allele for x in batch],
            "effect_weight": [x.effect_weight for x in batch],
            "effect_type": [str(x.effect_type) for x in batch],
            "is_duplicated": [x.is_duplicated for x in batch],
            "accession": [x.accession for x in batch],
            "row_nr": [x.row_nr for x in batch],
        }

        record_batch = pa.RecordBatch.from_pydict(batch_dict, schema=self.schema)
        self._writer.write(record_batch)

    def __del__(self):
        # it's very important to close the writer and file, or it gets corrupted
        # can't use a with statement, so close when the object gets deleted
        self._writer.close()
        if not self._sink.closed:
            self._sink.close()


def write_combined(
    scoring_files: list[ScoringFile], out_path: str
) -> dict[str : typing.Counter]:
    # compresslevel can be really slow, default is 9
    match fn := out_path.lower():
        case _ if fn.endswith("gz"):
            writer = TextFileWriter(compress=True, filename=out_path)
        case _ if fn.endswith("txt"):
            writer = TextFileWriter(compress=False, filename=out_path)
        case _ if fn.endswith("sqlite"):
            writer = SqliteWriter(filename=out_path)
        case _ if fn.endswith("ipc"):
            writer = PyarrowWriter(filename=out_path)
        case _:
            raise ValueError(f"Unsupported file extension: {out_path}")

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
