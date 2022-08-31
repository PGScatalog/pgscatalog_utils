import logging
from functools import reduce

import requests

logger = logging.getLogger(__name__)


def query_trait(trait: str) -> list[str]:
    api: str = f'https://www.pgscatalog.org/rest/trait/{trait}?include_children=1'
    logger.debug(f"Querying PGS Catalog with trait {trait}")
    r: requests.models.Response = requests.get(api)

    if r.json() == {}:
        logger.critical(f"Bad response from PGS Catalog for EFO term: {trait}")
        raise Exception

    keys: list[str] = ['associated_pgs_ids', 'child_associated_pgs_ids']
    pgs: list[str] = []
    for key in keys:
        pgs.append(r.json().get(key))

    logger.debug(f"Valid response from PGS Catalog for EFO term: {trait}")
    return list(reduce(lambda x, y: set(x).union(set(y)), pgs))
