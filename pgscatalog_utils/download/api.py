import requests
import jq
import logging
import sys

logger = logging.getLogger(__name__)


def pgscatalog_result(pgs: list[str], build: str) -> dict[str, str]:
    result = _parse_json_query(_api_query(pgs), build)

    try:
        if len(pgs) > len(result):
            missing_pgs: set[str] = set(pgs).difference(set(result.keys()))
            logger.warning(f"Some queries missing in PGS Catalog response: {missing_pgs}")
    except TypeError:
        logger.error(f"Bad response from PGS Catalog API. Is {pgs} a valid ID?")
        sys.exit(1)

    return result


def _api_query(pgs_id: list[str]) -> dict:
    pgs: str = ','.join(pgs_id)
    api: str = f'https://www.pgscatalog.org/rest/score/search?pgs_ids={pgs}'
    r: requests.models.Response = requests.get(api)
    return r.json()


def _parse_json_query(json: dict, build: str) -> dict[str, str]:
    result = jq.compile(".results").input(json).first()
    if not result:
        logger.warning("No results in response from PS Catalog API. Please check the PGS IDs.")
    else:
        return _extract_ftp_url(json, build)


def _extract_ftp_url(json: list[dict], build: str) -> dict[str, str]:
    id: list[str] = jq.compile('[.results][][].id').input(json).all()
    result: list[str] = jq.compile(f'[.results][][].ftp_harmonized_scoring_files.{build}.positions').input(json).all()
    return dict(zip(id, [x.replace('https', 'ftp') for x in result]))


