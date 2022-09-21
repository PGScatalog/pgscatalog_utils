import zstandard
from dataclasses import dataclass
import io
import logging
import polars as pl

logger = logging.getLogger(__name__)


@dataclass
class Score:
    """ A class that represents calculated scores (.sscore)"""
