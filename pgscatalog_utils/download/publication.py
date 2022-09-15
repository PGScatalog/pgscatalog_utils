import logging
from functools import reduce

import requests

logger = logging.getLogger(__name__)


def query_publication(pgp: str) -> list[str]:
    api: str = f'https://www.pgscatalog.org/rest/publication/{pgp}'
    logger.debug("Querying PGS Catalog with publication PGP ID")
    r: requests.models.Response = requests.get(api)

    if r.json() == {}:
        logger.critical(f"Bad response from PGS Catalog for EFO term: {pgp}")
        raise Exception

    pgs: dict[str, list[str]] = r.json().get('associated_pgs_ids')
    logger.debug(f"Valid response from PGS Catalog for PGP ID: {pgp}")
    return list(reduce(lambda x, y: set(x).union(set(y)), pgs.values()))
