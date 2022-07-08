import os

import pandas as pd

from pgscatalog_utils.scorefile.read import load_scorefile
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.effect_type import set_effect_type
import logging

logger = logging.getLogger(__name__)
log_fmt = "%(name)s: %(asctime)s %(levelname)-8s %(message)s"
logging.basicConfig(level=logging.DEBUG,
                    format=log_fmt,
                    datefmt='%Y-%m-%d %H:%M:%S')

paths = ["/Users/bwingfield/Documents/data/harmonised_scores/37/PGS000018_hmPOS_GRCh37.txt.gz",
         "/Users/bwingfield/Documents/data/harmonised_scores/37/PGS000049_hmPOS_GRCh37.txt.gz"]

scorefiles: pd.DataFrame = (pd.concat([load_scorefile(x) for x in paths])
                                .pipe(melt_effect_weights)
                                .pipe(set_effect_type))
