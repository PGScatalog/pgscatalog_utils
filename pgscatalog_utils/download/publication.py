import logging
from functools import reduce

from pgscatalog_utils.download.score import query_api

logger = logging.getLogger(__name__)


def query_publication(pgp: str, user_agent:str = None) -> list[str]:
    logger.debug("Querying PGS Catalog with publication PGP ID")
    api: str = f'/publication/{pgp}'
    results_json = query_api(api, user_agent)

    if results_json == {} or results_json == None:
        logger.critical(f"Bad response from PGS Catalog for EFO term: {pgp}")
        raise Exception

    pgs: dict[str, list[str]] = results_json.get('associated_pgs_ids')
    logger.debug(f"Valid response from PGS Catalog for PGP ID: {pgp}")
    return list(reduce(lambda x, y: set(x).union(set(y)), pgs.values()))
