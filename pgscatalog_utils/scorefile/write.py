import csv
import functools
import gzip
import logging
from itertools import islice

import pgzip

from pgscatalog_utils.scorefile.config import Config
from pgscatalog_utils.scorefile.scoringfile import ScoringFile

logger = logging.getLogger(__name__)


def write_combined(scoring_files: list[ScoringFile], out_path: str):
    # compresslevel can be really slow, default is 9
    if out_path.endswith("gz") and Config.threads == 1:
        logger.info("Writing with gzip (slow)")
        open_function = functools.partial(gzip.open, compresslevel=6)
    elif Config.threads > 1:
        logger.info("Writing with pgzip (fast)")
        open_function = functools.partial(pgzip.open, compresslevel=6,
                                          thread=Config.threads, blocksize=2 * 10 ** 8)
    else:
        logger.info("Writing text file (fast)")
        open_function = open

    with open_function(out_path, mode='wt') as f:
        fieldnames = ["chr_name", "chr_position", "effect_allele",
                      "other_allele", "effect_weight", "effect_type", "accession", "row_nr"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()

        # write out in batches for compression efficiency and speed
        for scoring_file in scoring_files:
            logger.info(f"Writing {scoring_file.accession} variants")
            while True:
                batch = list(islice(scoring_file.variants, Config.batch_size))
                if not batch:
                    break
                writer.writerows(batch)
