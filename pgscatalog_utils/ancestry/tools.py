import pandas as pd
from sklearn.covariance import MinCovDet
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
from scipy.stats import chi2, zscore, percentileofscore

####
# Methods for assigning ancestry using PCA data & population labels
###

_assign_methods = ["Mahalanobis", "RF"]


def assign_ancestry(ref_df, ref_pop_col, target_df, ref_train_col=None, n_pcs=4, method='RF'):
    # Check that datasets have the correct columns
    assert method in _assign_methods, 'ancestry assignment method parameter must be Mahalanobis or RF'

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

            robust_cov = MinCovDet().fit(
                ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs]
            )

            # Caclulate Mahalanobis distance of each sample to that population
            # Reference Samples
            ref_df[colname_dist] = robust_cov.mahalanobis(ref_df[cols_pcs])
            ref_df[colname_pval] = chi2.sf(ref_df[colname_dist], n_pcs - 1)
            # Target Samples
            target_df[colname_dist] = robust_cov.mahalanobis(target_df[cols_pcs])
            target_df[colname_pval] = chi2.sf(target_df[colname_dist], n_pcs - 1)

            pval_cols.append(colname_pval)

        # Assign population (maximum probability)
        ref_assign = ref_df[pval_cols].copy()
        ref_assign = ref_assign.assign(Predicted_Population=ref_assign.idxmax(axis=1))
        ref_assign['Predicted_Population'] = [x.split('_')[-1] for x in ref_assign['Predicted_Population']]

        target_assign = target_df[pval_cols].copy()
        target_assign = target_assign.assign(Predicted_Population=target_assign.idxmax(axis=1))
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

    return ref_assign, target_assign


####
# Methods for adjusting/reporting polygenic score results that account for genetic ancestry
####
