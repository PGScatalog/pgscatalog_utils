import pandas as pd
from sklearn.covariance import MinCovDet, EmpiricalCovariance
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
from scipy.stats import chi2, percentileofscore

####
# Methods for assigning ancestry using PCA data & population labels
###

_assign_methods = ["Mahalanobis", "RF"]
_mahalanobis_methods = ["MinCovDet", "EmpiricalCovariance"]


def assign_ancestry(ref_df, ref_pop_col, target_df, ref_train_col=None, n_pcs=4, method='RF', covariance_method = 'MinCovDet', p_threshold=None):
    """
    :param ref_df:
    :param ref_pop_col:
    :param target_df:
    :param ref_train_col:
    :param n_pcs:
    :param method:
    :param covariance_method:
    :param p_threshold: suggest 0.9 for RF and 0.001 for Mahalanobis
    :return:
    """
    # Check that datasets have the correct columns
    assert method in _assign_methods, 'ancestry assignment method parameter must be Mahalanobis or RF'
    if method == 'Mahalanobis':
        assert covariance_method in _mahalanobis_methods, 'covariance estimation method must be MinCovDet or EmpiricalCovariance'

    cols_pcs = ['PC{}'.format(x + 1) for x in range(0, n_pcs)]
    assert all([col in ref_df.columns for col in cols_pcs]), \
        "Reference Dataset (ref_df) is missing some PC columns for population assignment (max:{})".format(n_pcs)
    assert all([col in target_df.columns for col in cols_pcs]), \
        "Target Dataset (target_df) is missing some PC columns for population assignment (max:{})".format(n_pcs)
    assert ref_pop_col in ref_df.columns, "Population label column ({}) is missing from reference dataframe".format(ref_pop_col)
    ref_populations = ref_df[ref_pop_col].unique()

    ## Extract columns for analysis
    ref_df = ref_df[cols_pcs + [ref_pop_col, ref_train_col]].copy()
    target_df = target_df[cols_pcs].copy()

    ## Create Training dfs
    if ref_train_col:
        assert ref_train_col in ref_df.columns, "Training index column({}) is missing from reference dataframe".format(ref_train_col)
        ref_train_df = ref_df.loc[ref_df[ref_train_col] == True,]
    else:
        ref_train_df = ref_df

    # Run Ancestry Assignment methods
    if method == 'Mahalanobis':
        # Calculate population distances
        pval_cols = []
        for pop in ref_populations:
            # Fit the covariance matrix for the current population
            colname_dist = 'Mahalanobis_dist_{}'.format(pop)
            colname_pval = 'Mahalanobis_P_{}'.format(pop)

            if covariance_method == 'MinCovDet':
                # Robust method
                covariance_fit = MinCovDet().fit(
                    ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs]
                )
            elif covariance_method == 'EmpiricalCovariance':
                covariance_fit = EmpiricalCovariance().fit(
                    ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs]
                )

            # Caclulate Mahalanobis distance of each sample to that population
            # Reference Samples
            ref_df[colname_dist] = covariance_fit.mahalanobis(ref_df[cols_pcs])
            ref_df[colname_pval] = chi2.sf(ref_df[colname_dist], n_pcs - 1)
            # Target Samples
            target_df[colname_dist] = covariance_fit.mahalanobis(target_df[cols_pcs])
            target_df[colname_pval] = chi2.sf(target_df[colname_dist], n_pcs - 1)

            pval_cols.append(colname_pval)

        # Assign population (maximum probability)
        ref_assign = ref_df[pval_cols].copy()
        ref_assign = ref_assign.assign(Predicted_Population=ref_assign.idxmax(axis=1))

        target_assign = target_df[pval_cols].copy()
        target_assign = target_assign.assign(Predicted_Population=target_assign.idxmax(axis=1))

        if p_threshold:
            ref_assign['Predicted_Population'] = [x.split('_')[-1] if (ref_assign[x][i] > p_threshold) else 'OTH' for i,x in enumerate(ref_assign['Predicted_Population'])]
            target_assign['Predicted_Population'] = [x.split('_')[-1] if (target_assign[x][i] > p_threshold) else 'OTH' for i, x in enumerate(target_assign['Predicted_Population'])]
        else:
            ref_assign['Predicted_Population'] = [x.split('_')[-1] for x in ref_assign['Predicted_Population']]
            target_assign['Predicted_Population'] = [x.split('_')[-1] for x in target_assign['Predicted_Population']]

    elif method == 'RF':
        # Assign SuperPop Using Random Forest (PCA loadings)
        clf_rf = RandomForestClassifier()
        clf_rf.fit(ref_train_df[cols_pcs],  # X (training PCs)
                   ref_train_df[ref_pop_col].astype(str))  # Y (pop label)
        # Assign population
        ref_assign = pd.DataFrame(clf_rf.predict_proba(ref_df[cols_pcs]), index=ref_df.index, columns=['RF_P_{}'.format(x) for x in clf_rf.classes_])
        ref_assign['Predicted_Population'] = clf_rf.predict(ref_df[cols_pcs])

        target_assign = pd.DataFrame(clf_rf.predict_proba(target_df[cols_pcs]), index=target_df.index, columns=['RF_P_{}'.format(x) for x in clf_rf.classes_])
        target_assign['Predicted_Population'] = clf_rf.predict(target_df[cols_pcs])

        if p_threshold:
            ref_assign['Predicted_Population'] = [x if (ref_assign['RF_P_{}'.format(x)][i] > p_threshold) else 'OTH'
                                                  for i, x in enumerate(clf_rf.predict(ref_df[cols_pcs]))]
            target_assign['Predicted_Population'] = [
                x if (target_assign['RF_P_{}'.format(x)][i] > p_threshold) else 'OTH' for i, x in
                enumerate(clf_rf.predict(target_df[cols_pcs]))]

    return ref_assign, target_assign


####
# Methods for adjusting/reporting polygenic score results that account for genetic ancestry
####

def pgs_adjust(ref_df, target_df, scorecols: list, ref_pop_col, target_pop_col, ref_train_col=None, n_pcs=5):
    # Check that datasets have the correct columns
    ## Check that score is in both dfs
    assert all(
        [x in ref_df.columns for x in scorecols]), "Reference Dataset (ref_df) is missing some PGS column(s)".format(
        scorecols)
    assert all(
        [x in target_df.columns for x in scorecols]), "Target Dataset (target_df) is missing some PGS column(s)".format(
        scorecols)

    ## Check that PCs is in both dfs
    cols_pcs = ['PC{}'.format(x + 1) for x in range(0, n_pcs)]
    assert all([col in ref_df.columns for col in cols_pcs]), \
        "Reference Dataset (ref_df) is missing some PC columns for PCA adjustment (max:{})".format(n_pcs)
    assert all([col in target_df.columns for col in cols_pcs]), \
        "Target Dataset (target_df) is missing some PC columns for PCA adjustment (max:{})".format(n_pcs)
    assert ref_pop_col in ref_df.columns, "Population label column ({}) is missing from reference dataframe".format(
        ref_pop_col)
    ref_populations = ref_df[ref_pop_col].unique()
    assert target_pop_col in target_df.columns, "Population label column ({}) is missing from target dataframe".format(
        target_pop_col)
    ## Extract columns for analysis
    ref_df = ref_df[cols_pcs + [ref_pop_col, ref_train_col]].copy()
    target_df = target_df[cols_pcs + [target_pop_col]].copy()

    ## Create Training dfs
    if ref_train_col:
        assert ref_train_col in ref_df.columns, "Training index column({}) is missing from reference dataframe".format(
            ref_train_col)
        ref_train_df = ref_df.loc[ref_df[ref_train_col] == True,]
    else:
        ref_train_df = ref_df

    # Empirical adjustment with reference population assignments
    for pop in ref_populations:
        if pop == 'OTH':
            # ToDo: implement handling of individuals who don't have population label (weighted average based on distance?)
            continue

        ref_pop = ref_train_df[ref_train_df[ref_pop_col] == pop]  # Reference dataset
        i_ref_pop = (ref_df[ref_col_pop] == pop)
        i_target_pop = (target_df[target_pop_col] == pop)

        for c_pgs in scorecols:
            # Score Distribution
            c_pgs_pop_dist = ref_pop[c_pgs]

            # Calculate Percentile
            percentile_col = '{}.adj_empirical_percentile'.format(c_pgs)
            ref_df.loc[i_ref_pop, percentile_col] = percentileofscore(c_pgs_pop_dist, ref_df.loc[i_ref_pop, c_pgs])
            target_df.loc[i_target_pop, percentile_col] = percentileofscore(c_pgs_pop_dist, target_df.loc[i_target_pop, c_pgs])

            # Calculate Z
            z_col = '{}.adj_empirical_Z'.format(c_pgs)
            c_pgs_mean = c_pgs_pop_dist.mean()
            c_pgs_std = c_pgs_pop_dist.std(ddof=0)

            ref_df.loc[i_ref_pop, z_col] = (ref_df.loc[i_ref_pop, c_pgs] - c_pgs_mean)/c_pgs_std
            target_df.loc[i_target_pop, z_col] = (target_df.loc[i_target_pop, c_pgs] - c_pgs_mean)/c_pgs_std

    # PCA-based adjustment
    for c_pgs in scorecols:
        ## Method 1 (Khera): normalize mean (doi:10.1161/CIRCULATIONAHA.118.035658)
        ### Fit to Reference Data
        pcs2pgs_fit = LinearRegression().fit(ref_train_df[cols_pcs], ref_train_df[c_pgs])
        ref_train_pgs_pred = pcs2pgs_fit.predict(ref_train_df[cols_pcs])
        ref_train_pgs_resid = ref_train_df[c_pgs] - ref_train_pgs_pred
        ref_train_pgs_resid_mean = ref_train_pgs_resid.mean()
        ref_train_pgs_resid_std = ref_train_pgs_resid.std(ddof=0)

        ref_pgs_resid = ref_df[c_pgs] - pcs2pgs_fit.predict(ref_df[cols_pcs])
        ref_df['{}.adj_1_Khera'.format(pgs_col)] = ref_pgs_resid / ref_train_pgs_resid_std
        ### Apply to Target Data
        target_pgs_pred = pcs2pgs_fit.predict(target_df[cols_pcs])
        target_pgs_resid = target_df[pgs_col] - target_pgs_pred
        target_df['{}.adj_1_Khera'.format(pgs_col)] = target_pgs_resid / ref_train_pgs_resid_std

        ## Method 2 (Khan): normalize variance (doi:10.1038/s41591-022-01869-1)
        ### Normalize based on residual deviation from mean of the distribution
        ### (e.g. reduce the effect of genetic ancestry on how far away you are from the mean [equalize population sds])
        pcs2var_fit = LinearRegression().fit(ref_train_df[cols_pcs], (ref_train_pgs_resid - ref_train_pgs_resid_mean)**2)

        ### Alternative (not adjusting for the mean... which should be 0 anyways because we're trying to fit it)
        #pcs2var_fit = LinearRegression().fit(ref_train_df[cols_pcs], ref_train_pgs_resid**2)

        ref_df['{}.adj_2_Khan'.format(pgs_col)] = ref_pgs_resid / np.sqrt(pcs2var_fit.predict(ref_df[cols_pcs]))
        ### Apply to Target
        target_df['{}.adj_2_Khan'.format(pgs_col)] = target_pgs_resid / np.sqrt(pcs2var_fit.predict(target_df[cols_pcs]))

    # Output adjustment columns in reference and target dataframes
    #return ref_df[x if ('.adj_' in x) for x in ref_df.columns], target_df[x if ('.adj_' in x) for x in target_df.columns]