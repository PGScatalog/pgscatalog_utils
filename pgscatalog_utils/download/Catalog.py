import logging
import time
import typing
from dataclasses import dataclass, field
from functools import reduce

import requests

from pgscatalog_utils import __version__ as pgscatalog_utils_version
from pgscatalog_utils import config
from pgscatalog_utils.download.CatalogCategory import CatalogCategory
from pgscatalog_utils.download.ScoringFile import ScoringFile

logger = logging.getLogger(__name__)


@dataclass
class CatalogResult:
    """
    PGS Catalog data returned from a CatalogQuery is stored as a CatalogResult

    This class mostly offers a convenient way to get a list of scoring file URLs

    CatalogResult will also automatically resolve PGS IDs associated with publications and traits
    """
    accession: typing.Union[str, list[str]]
    category: CatalogCategory
    response: typing.Union[dict, None]
    pgs_ids: set[str] = field(init=False)
    include_children: bool = False

    def _grab_pgs_ids(self) -> set[str]:
        logger.info(f"Valid response from PGS Catalog for term: {self.accession}")
        pgs: list[set[str]] = []

        match self.category:
            case CatalogCategory.TRAIT | CatalogCategory.PUBLICATION:
                keys: list[str] = ["associated_pgs_ids"]
                if self.include_children:
                    keys.append("child_associated_pgs_ids")
                for key in keys:
                    match self.response.get(key):
                        case list() as pgs_list:
                            pgs.append(set(pgs_list))
                        case dict() as pgs_dict:
                            pgs.append(reduce(lambda x, y: set(x).union(set(y)), pgs_dict.values()))

                return set().union(*pgs)
            case CatalogCategory.SCORE:
                for result in self.response.get("results"):
                    pgs.append({result.get("id")})
                return set().union(*pgs)
            case _:
                raise Exception(f"Invalid {self.category}")

    def __post_init__(self):
        if self.response is None or self.response == {}:
            raise Exception(f"Bad response from PGS Catalog for term: {self.accession}")
        else:
            self.pgs_ids = self._grab_pgs_ids()

    def get_download_urls(self) -> dict[str: ScoringFile]:
        urls = {}
        match self.category:
            case CatalogCategory.SCORE:
                # scores already have the scoring file URL information in the response
                for result in self.response.get("results"):
                    urls[result.get("id")] = ScoringFile.from_result(result)
            case CatalogCategory.TRAIT | CatalogCategory.PUBLICATION:
                # publications and traits have to query Catalog API again to grab score data
                results: list[CatalogResult] = CatalogQuery(CatalogCategory.SCORE,
                                                            accession=list(self.pgs_ids),
                                                            pgsc_calc_version=config.PGSC_CALC_VERSION).get()
                for result in results:
                    for pgs in result.response.get("results"):
                        urls[pgs["id"]] = ScoringFile.from_result(pgs)
        return urls


@dataclass
class CatalogQuery:
    """
    Efficiently batch query the PGS Catalog API using trait, score, or publication identifiers
    """
    category: CatalogCategory
    accession: typing.Union[str, list[str]]
    pgsc_calc_version: typing.Union[str, None]
    include_children: bool = False
    _rest_url_root: str = "https://www.pgscatalog.org/rest"
    _max_retries: int = 5
    _version: str = pgscatalog_utils_version
    _user_agent: dict[str: str] = field(init=False)

    def _resolve_query_url(self) -> typing.Union[str, list[str]]:
        child_flag: int = int(self.include_children)
        # TRAIT and PUBLICATION should have one accession (return many PGS IDs)
        # but SCORE queries should be preferably batched
        match (self.category, self.accession):
            case CatalogCategory.TRAIT, str():
                return f"{self._rest_url_root}/trait/{self.accession}?include_children={child_flag}"
            case CatalogCategory.SCORE, str():
                return f"{self._rest_url_root}/score/search?pgs_ids={self.accession}"
            case CatalogCategory.SCORE, list():
                urls: list[str] = []
                for chunk in self._chunk_query():
                    chunked_accession = ",".join(chunk)
                    urls.append(f"{self._rest_url_root}/score/search?pgs_ids={chunked_accession}")
                return urls
            case CatalogCategory.PUBLICATION, str():
                return f"{self._rest_url_root}/publication/{self.accession}"
            case _:
                raise Exception(f"Invalid CatalogCategory and accession type: {self.category}, type({self.accession})")

    def __post_init__(self):
        ua: str
        if self.pgsc_calc_version:
            ua = pgscatalog_utils_version
        else:
            ua = f"pgscatalog_utils/{self._version}"

        self._user_agent = {"User-Agent": ua}

    def _query_api(self, url: str):
        wait: int = 10
        retry: int = 0
        results_json = None

        while retry < self._max_retries:
            try:
                logger.info(f"Querying {url}")
                r: requests.models.Response = requests.get(url, headers=self._user_agent)
                r.raise_for_status()
                results_json = r.json()
                break
            except requests.exceptions.HTTPError as e:
                logging.warning(f"HTTP error {e}")
                if r.status_code in [421, 429]:
                    retry += 1
                    logger.warning(f"Retrying API call in {wait} seconds")
                    logger.warning(f"Attempt {retry} of {self._max_retries}")
                    time.sleep(wait)
            except requests.exceptions.RequestException as e:
                logger.critical(f"Request error. Can't retry, so bailing out: {e}")
                break

        return results_json

    def _chunk_query(self):
        size = 50  # /rest/score/{pgs_id} limit when searching multiple IDs
        return (self.accession[pos: pos + size] for pos in range(0, len(self.accession), size))

    def get(self) -> list[CatalogResult]:
        query_url: typing.Union[str, list[str]] = self._resolve_query_url()
        results: list[CatalogResult] = []
        match query_url:
            case str():
                results.append(CatalogResult(accession=self.accession,
                                             category=self.category,
                                             include_children=self.include_children,
                                             response=self._query_api(query_url)))
            case list():
                for url in query_url:
                    results.append(CatalogResult(accession=self.accession,
                                                 category=self.category,
                                                 response=self._query_api(url)))
            case _:
                raise Exception(f"Invalid query url type: {type(query_url)}")
        return results
