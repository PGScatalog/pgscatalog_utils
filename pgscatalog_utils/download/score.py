import logging
import sys

import jq
import requests
import time
from pgscatalog_utils import __version__ as pgscatalog_utils_version

logger = logging.getLogger(__name__)


def get_url(pgs: list[str], build: str, user_agent:str = None) -> dict[str, str]:
    pgs_result: list[str] = []
    url_result: list[str] = []

    for chunk in _chunker(pgs):
        try:
            response = _parse_json_query(query_score(chunk,user_agent), build)
            pgs_result = pgs_result + list(response.keys())
            url_result = url_result + list(response.values())
        except (AttributeError, TypeError):
            logger.error(f"Bad response from PGS Catalog API. Is {pgs} a valid ID?")
            sys.exit(1)

    missing_pgs = set(pgs).difference(set(pgs_result))

    if missing_pgs:
        logger.warning(f"Some queries missing in PGS Catalog response: {missing_pgs}")

    return dict(zip(pgs_result, url_result))


def query_api(api: str, user_agent:str = None, retry:int = 0) -> dict:
    max_retries = 5
    wait = 60
    results_json = None
    rest_url_root = 'https://www.pgscatalog.org/rest'
    # Set pgscatalog_utils user agent if none provided
    if not user_agent:
        user_agent = 'pgscatalog_utils/'+pgscatalog_utils_version
    try:
        headers = {'User-Agent': user_agent}
        r: requests.models.Response = requests.get(rest_url_root+api, headers=headers)
        r.raise_for_status()
        results_json = r.json()
    except requests.exceptions.HTTPError as e:
        print(f'HTTP Error: {e}')
        if r.status_code in [421,429] and retry < 5:
            retry +=1
            print(f'> Retry to query the PGS Catalog REST API in {wait}s ... attempt {retry} out of {max_retries}.')
            time.sleep(wait)
            results_json = query_api(api,retry)
    except requests.exceptions.ConnectionError as e:
        print(f'Error Connecting: {e}')
    except requests.exceptions.Timeout as e:
        print(f'Timeout Error: {e}')
    except requests.exceptions.RequestException as e:
        print(f'Request Error: {e}')
    return results_json


def query_score(pgs_id: list[str], user_agent:str = None) -> dict:
    pgs: str = ','.join(pgs_id)
    api: str = f'/score/search?pgs_ids={pgs}'
    results_json = query_api(api, user_agent)
    return results_json


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
