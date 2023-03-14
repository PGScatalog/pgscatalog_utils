import logging
import pandas as pd
import numpy as np
import os

logger = logging.getLogger(__name__)

def read_pcs(loc_pcs: list[str],dataset: str, loc_related_ids=None, nPCs=None):
    """
    Read the .pc file outputs of the fraposa_pgsc projection
    :param loc_pcs: list of locations for .pcs files
    :param dataset: name of the dataset being read (used for index)
    :param loc_related_ids: path to newline-delimited list of IDs for related samples that can be used to filter
    :return: pandas dataframe with PC information
    """
    proj = pd.DataFrame()

    for i, path in enumerate(loc_pcs):
        logger.debug("Reading PCA projection: {}".format(path))
        df = pd.read_csv(path, sep='\t')
        df['sampleset'] = dataset
        df.set_index(['sampleset', 'IID'], inplace=True)

        if i == 0:
            logger.debug('Initialising combined DF')
            proj = df.copy()
        else:
            logger.debug('Appending to combined DF')
            proj = pd.concat([proj, df])

    # Read/process IDs for unrelated samples (usually reference dataset)
    if loc_related_ids:
        logger.debug("Flagging related samples with: {}".format(loc_related_ids))
        proj['Unrelated'] = True
        with open(loc_related_ids, 'r') as infile:
            IDs_related = [x.strip() for x in infile.readlines()]
        proj.loc[proj.index.get_level_values(level=1).isin(IDs_related), 'Unrelated'] = False
    else:
        proj['Unrelated'] = np.nan

    # Drop PCs
    if nPCs:
        logger.debug('Filtering to relevant PCs')
        dropcols = []
        for x in proj.columns:
            if int(x[2:]) > nPCs:
                dropcols.append(x)
        proj = proj.drop(dropcols, axis=1)

    return proj

def read_projection(loc_sscores: list[str],dataset: str, loc_related_ids=None, nPCs=None):
    """
    Read PCA projection data from pgsc_calc pipeline (outputs of PLINK2_PROJECT)
    :param loc_sscores: list of paths to the result of PCA projection (.sscore format)
    :param dataset: name of the dataset (used to create multi-index)
    :param loc_related_ids: loc_related_ids: path to newline-delimited list of IDs for related samples that can be used to filter
    :param nPCs: maximum number of PCs to extract from the files
    :return: pandas dataframe with PC information
    """
    proj = pd.DataFrame()
    nvars = []

    for i, path in enumerate(loc_sscores):
        logger.debug("Reading PCA projection: {}".format(path))
        df = pd.read_csv(path, sep='\t')

        # Check nvars
        path_vars = path + '.vars'
        if os.path.isfile(path_vars):
            nvars.append(len(open(path_vars, 'r').read().strip().split('\n')))
        else:
            nvars.append(None)

        match (df.columns[0]):
            # handle case of #IID -> IID (happens when #FID is present)
            case '#IID':
                df.rename({'#IID': 'IID'}, axis=1, inplace=True)
            case '#FID':
                df.drop(['#FID'], axis=1,  inplace=True)
            case _:
                assert False, "Invalid columns"

        df['sampleset'] = dataset
        df.set_index(['sampleset', 'IID'], inplace=True)

        if i == 0:
            logger.debug('Initialising combined DF')
            proj = df.copy()
            aggcols = [x for x in df.columns if (x.startswith('PC') and x.endswith('_SUM'))]
        else:
            logger.debug('Adding to combined DF')
            proj = proj[aggcols].add(df[aggcols], fill_value=0)

    # Filter & rename PC columns
    if nPCs:
        logger.debug('Filtering to relevant PCs')
        dropcols = []
        for x in aggcols:
            if int(x.split('_')[0][2:]) > nPCs:
                dropcols.append(x)
        proj = proj.drop(dropcols, axis =1)

    proj.columns = [x.replace('_SUM', '') for x in proj.columns]

    # Read/process IDs for unrelated samples (usually reference dataset)
    if loc_related_ids:
        logger.debug("Flagging related samples with: {}".format(loc_related_ids))
        proj['Unrelated'] = True
        with open(loc_related_ids, 'r') as infile:
            IDs_related = [x.strip() for x in infile.readlines()]
        proj.loc[proj.index.get_level_values(level=1).isin(IDs_related), 'Unrelated'] = False
    else:
        proj['Unrelated'] = np.nan

    if None in nvars:
        return proj, None
    else:
        return proj, sum(nvars)


def read_ref_psam(loc_psam):
    psam = pd.read_csv(loc_psam, sep='\t', comment='##')


def read_pgs(loc_aggscore, onlySUM: bool):
    """
    Function to read the output of aggreagte_scores
    :param loc_aggscore: path to aggregated scores output
    :param onlySUM: whether to return only _SUM columns (e.g. not _AVG)
    :return:
    """
    logger.debug('Reading aggregated score data: {}'.format(loc_aggscore))
    df = pd.read_csv(loc_aggscore, sep='\t', index_col=['sampleset', 'IID'])
    if onlySUM:
        df = df[[x for x in df.columns if x.endswith('_SUM')]]
        rn = [x.rstrip('_SUM') for x in df.columns]
        df.columns = rn
    return df