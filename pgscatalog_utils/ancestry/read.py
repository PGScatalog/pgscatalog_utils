import pandas as pd
import numpy as np


def read_projection(loc_sscore, loc_related_ids=None):
    """
    Read PCA projection data from pgsc_calc pipeline
    :param loc_sscore: path to the result of PCA projection (.sscore format)
    :param loc_related_ids: path to newline-delimited list of IDs for related samples that can be used to filter
    :return: pandas dataframe with PC information
    """
    # Read combined data
    proj = pd.read_csv(loc_sscore, sep='\t')

    match (proj.columns[0]):
        # handle case of #IID -> IID (happens when #FID is present)
        case '#IID':
            pass
        case '#FID':
            proj.rename({'IID': '#IID'}, axis=g1, inplace = True)
            proj.drop(['#FID'], axis=1,  inplace = True)
        case _:
            assert False, "Invalid columns"

    proj.set_index('#IID')
    proj.columns = [x.replace('_SUM', '') for x in proj.columns]

    # Read/process IDs for unrelated samples (usually reference dataset)
    if loc_related_ids:
        proj['Unrelated'] = True
        with open(loc_related_ids, 'r') as infile:
            IDs_related = [x.strip() for x in infile.readlines()]
        proj.loc[proj.index.isin(IDs_related), 'Unrelated'] = False
    else:
        proj['Unrelated'] = np.nan

    return proj


def read_pgs(loc_sscore):
    pgs = pd.read_csv(loc_sscore, sep='\t', index_col='#IID')
    return pgs