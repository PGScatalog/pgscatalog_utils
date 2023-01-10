from sklearn.covariance import EmpiricalCovariance, MinCovDet
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression
from scipy.stats import chi2, zscore, percentileofscore

####
# Methods for assigning ancestry using PCA data & population labels
###

_assign_methods = ["Mahalanobis", "RF"]


def assign_ancestry(ref_df, ref_pop_col, target_df, ref_train_col=None, n_pcs=4, method='RF'):
    ref_df = ref_df.copy()
    target_df = target_df.copy()
    # Check that datasets have the correct columns
    assert method in _assign_methods, 'ancestry assignment method parameter must be Mahalanobis or RF'

    cols_pcs = ['PC{}'.format(x + 1) for x in range(0, n_pcs)]
    assert all([col in ref_df.columns for col in cols_pcs]), \
        "Reference Dataset (ref_df) is missing some PC columns for population assignment (max:{})".format(n_pcs)
    assert all([col in target_df.columns for col in cols_pcs]), \
        "Target Dataset (target_df) is missing some PC columns for population assignment (max:{})".format(n_pcs)

    assert ref_pop_col in ref_df.columns, "Population label column ({}) is missing from reference dataframe".format(ref_pop_col)
    ref_populations = ref_df[ref_pop_col].unique()

    if ref_train_col:
        assert ref_train_col in ref_df.columns, "Training index column({}) is missing from reference dataframe".format(ref_train_col)
        ref_train_df = ref_df.loc[ref_df[ref_train_col] == True,]
    else:
        ref_train_df = ref_df

    if method == 'Mahalanobis':
        # Calculate population distances
        d_cols = []
        for pop in ref_populations:
            # Fit the covariance matrix for the current population
            robust_cov = MinCovDet().fit(
                ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs]
            )

            # Caclulate Mahalanobis distance of each sample to that population
            # Reference Samples
            ref_df['MahalanobisD_{}'.format(pop)] = robust_cov.mahalanobis(ref_df[cols_pcs])
            # Target Samples
            target_df['MahalanobisD_{}'.format(pop)] = robust_cov.mahalanobis(target_df[cols_pcs])
            d_cols.append('MahalanobisD_{}'.format(pop))

        # Assign population (minimum distance)
        ref_assign = [x.split('_')[-1] for x in ref_df[d_cols].idxmin(axis=1)]
        target_assign = [x.split('_')[-1] for x in target_df[d_cols].idxmin(axis=1)]
    elif method == 'RF':
        # Assign SuperPop Using Random Forest (PCA loadings)
        clf_rf = RandomForestClassifier()
        clf_rf.fit(ref_train_df[cols_pcs],  # X (training PCs)
                   ref_train_df[ref_pop_col].astype(str))  # Y (pop label)
        # Assign population (ToDO: output the probability?)
        ref_assign = list(clf_rf.predict(ref_df[cols_pcs]))
        target_assign = list(clf_rf.predict(target_df[cols_pcs]))

    return ref_assign, target_assign




####
# Methods for adjusting/reporting polygenic score results that account for genetic ancestry
####
