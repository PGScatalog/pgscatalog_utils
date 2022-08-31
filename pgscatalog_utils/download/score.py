import logging
import sys

import jq
import requests

logger = logging.getLogger(__name__)


def get_url(pgs: list[str], build: str) -> dict[str, str]:
    pgs_result: list[str] = []
    url_result: list[str] = []

    for chunk in _chunker(pgs):
        try:
            response = _parse_json_query(query_score(chunk), build)
            pgs_result = pgs_result + list(response.keys())
            url_result = url_result + list(response.values())
        except TypeError:
            logger.error(f"Bad response from PGS Catalog API. Is {pgs} a valid ID?")
            sys.exit(1)

    missing_pgs = set(pgs).difference(set(pgs_result))

    if missing_pgs:
        logger.warning(f"Some queries missing in PGS Catalog response: {missing_pgs}")

    return dict(zip(pgs_result, url_result))


def query_score(pgs_id: list[str]) -> dict:
    pgs: str = ','.join(pgs_id)
    api: str = f'https://www.pgscatalog.org/rest/score/search?pgs_ids={pgs}'
    r: requests.models.Response = requests.get(api)
    return r.json()


def _chunker(pgs: list[str]):
    size = 50  # /rest/score/{pgs_id} limit when searching multiple IDs
    return (pgs[pos: pos + size] for pos in range(0, len(pgs), size))


def _parse_json_query(json: dict, build: str | None) -> dict[str, str]:
    result = jq.compile(".results").input(json).first()
    if not result:
        logger.warning("No results in response from PGS Catalog API. Please check the PGS IDs.")
    else:
        return _extract_ftp_url(json, build)


def _extract_ftp_url(json: list[dict], build: str | None) -> dict[str, str]:
    id: list[str] = jq.compile('[.results][][].id').input(json).all()
    if build is None:
        result: list[str] = jq.compile(f'[.results][][].ftp_scoring_file').input(
            json).all()
    else:
        result: list[str] = jq.compile(f'[.results][][].ftp_harmonized_scoring_files.{build}.positions').input(
            json).all()
    return dict(zip(id, [x.replace('https', 'ftp') for x in result]))
