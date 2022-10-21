import logging
from functools import reduce

from pgscatalog_utils.download.score import query_api

logger = logging.getLogger(__name__)


def query_trait(trait: str, user_agent:str = None, include_children:bool = True) -> list[str]:
    logger.debug(f"Querying PGS Catalog with trait {trait}")
    api: str = f'/trait/{trait}?include_children=0'
    results_json = query_api(api, user_agent)

    if results_json == {} or results_json == None:
        logger.critical(f"Bad response from PGS Catalog for EFO term: {trait}")
        raise Exception

    keys: list[str] = ['associated_pgs_ids']
    if include_children:
        keys.append('child_associated_pgs_ids')

    pgs: list[str] = []
    for key in keys:
        pgs.append(results_json.get(key))

    logger.debug(f"Valid response from PGS Catalog for EFO term: {trait}")
    return list(reduce(lambda x, y: set(x).union(set(y)), pgs))
