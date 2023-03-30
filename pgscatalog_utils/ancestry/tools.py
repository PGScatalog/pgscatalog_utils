import logging
import pandas as pd
import numpy as np
from sklearn.covariance import MinCovDet, EmpiricalCovariance
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, GammaRegressor
from scipy.stats import chi2, percentileofscore, mannwhitneyu
import json
import gzip

# import seaborn as sns
# import matplotlib.pyplot as plt


logger = logging.getLogger(__name__)


####
# Methods for comparing sample ancestry to a reference dataset using PCA data & population labels
###

comparison_method_threshold = {"Mahalanobis": 1e-10,
                           "RandomForest": 0.5}  # default p-value thresholds to define suitable population matches
_mahalanobis_methods = ["MinCovDet", "EmpiricalCovariance"]


def choose_pval_threshold(args):
    set_threshold = comparison_method_threshold[args.method_compare]  # method default
    if args.pThreshold is not None:
        if (args.pThreshold > 0) and (args.pThreshold < 1):
            set_threshold = args.pThreshold
        else:
            logging.warning("p-value threshold out of range, assigning as method default: {}")
    return set_threshold


def compare_ancestry(ref_df: pd.DataFrame, ref_pop_col: str, target_df: pd.DataFrame, ref_train_col=None, n_pcs=4, method='RandomForest',
                     covariance_method='EmpiricalCovariance', p_threshold=None):
    """
    Function to compare target sample ancestry to a reference panel with PCA data
    :param ref_df: reference dataset
    :param ref_pop_col: training labels for population assignment in reference dataset
    :param target_df: target dataset
    :param ref_train_col: column name of TRUE/FALSE labels for inclusion in training ancestry assignments (e.g. unrelated)
    :param n_pcs: number of PCs to use in ancestry assignment
    :param method: One of Mahalanobis or RandomForest
    :param covariance_method: Used to calculate Mahalanobis distances One of EmpiricalCovariance or MinCovDet
    :param p_threshold: used to define LowConfidence population assignments
    :return: dataframes for reference (predictions on training set) and target (predicted labels) datasets
    """
    # Check that datasets have the correct columns
    assert method in comparison_method_threshold.keys(), 'comparison method parameter must be Mahalanobis or RF'
    if method == 'Mahalanobis':
        assert covariance_method in _mahalanobis_methods, 'ovariance estimation method must be MinCovDet or EmpiricalCovariance'

    cols_pcs = ['PC{}'.format(x + 1) for x in range(0, n_pcs)]
    assert all([col in ref_df.columns for col in cols_pcs]), \
        "Reference Dataset (ref_df) is missing some PC columns for ancestry comparison (max:{})".format(n_pcs)
    assert all([col in target_df.columns for col in cols_pcs]), \
        "Target Dataset (target_df) is missing some PC columns for ancestry comparison (max:{})".format(n_pcs)
    assert ref_pop_col in ref_df.columns, "Population label column ({}) is missing from reference dataframe".format(ref_pop_col)
    ref_populations = ref_df[ref_pop_col].unique()

    # Extract columns for analysis
    ref_df = ref_df[cols_pcs + [ref_pop_col, ref_train_col]].copy()
    target_df = target_df[cols_pcs].copy()

    # Create Training dfs
    if ref_train_col:
        assert ref_train_col in ref_df.columns, "Training index column({}) is missing from reference dataframe".format(ref_train_col)
        ref_train_df = ref_df.loc[ref_df[ref_train_col] == True,]
    else:
        ref_train_df = ref_df

    # Check if PCs only capture target/reference stratification
    compare_info = {}
    for col_pc in cols_pcs:
        mwu_pc = mannwhitneyu(ref_train_df[col_pc], target_df[col_pc])
        compare_info[col_pc] = {'U': mwu_pc.statistic, 'pvalue': mwu_pc.pvalue}
        if mwu_pc.pvalue < 1e-4:
            logger.warning("{} *may* be capturing target/reference stratification (Mann-Whitney p-value={}), "
                           "use visual inspection of PC plot to confirm".format(col_pc, mwu_pc.pvalue))

    # Run Ancestry Assignment methods
    if method == 'Mahalanobis':
        logger.debug("Calculating Mahalanobis distances")
        # Calculate population distances
        pval_cols = []
        for pop in ref_populations:
            logger.debug("Fitting Mahalanobis distances: {}".format(pop))
            # Fit the covariance matrix for the current population
            colname_dist = 'Mahalanobis_dist_{}'.format(pop)
            colname_pval = 'Mahalanobis_P_{}'.format(pop)

            match covariance_method:
                case 'MinCovDet':
                    covariance_model = MinCovDet()
                case 'EmpiricalCovariance':
                    covariance_model = EmpiricalCovariance()
                case _:
                    assert False, "Invalid covariance method"

            covariance_fit = covariance_model.fit(ref_train_df.loc[ref_train_df[ref_pop_col] == pop, cols_pcs])

            # Caclulate Mahalanobis distance of each sample to that population
            # Reference Samples
            ref_df[colname_dist] = covariance_fit.mahalanobis(ref_df[cols_pcs])
            ref_df[colname_pval] = chi2.sf(ref_df[colname_dist], n_pcs - 1)
            # Target Samples
            target_df[colname_dist] = covariance_fit.mahalanobis(target_df[cols_pcs])
            target_df[colname_pval] = chi2.sf(target_df[colname_dist], n_pcs - 1)

            pval_cols.append(colname_pval)

        # Assign population (maximum probability)
        logger.debug("Assigning Populations (max Mahalanobis probability)")
        ref_assign = ref_df[pval_cols].copy()
        ref_assign = ref_assign.assign(MostSimilarPop=ref_assign.idxmax(axis=1))

        target_assign = target_df[pval_cols].copy()
        target_assign = target_assign.assign(MostSimilarPop=target_assign.idxmax(axis=1))

        ref_assign['MostSimilarPop_LowConfidence'] = np.nan
        target_assign['MostSimilarPop_LowConfidence'] = np.nan

        if p_threshold:
            logger.debug("Comparing Population Similarity to p-value threshold (p < {})".format(p_threshold))
            ref_assign['MostSimilarPop_LowConfidence'] = [(ref_assign[x][i] < p_threshold)
                                                      for i,x in enumerate(ref_assign['MostSimilarPop'])]
            target_assign['MostSimilarPop_LowConfidence'] = [(target_assign[x][i] < p_threshold)
                                                         for i, x in enumerate(target_assign['MostSimilarPop'])]

        # Cleanup variable names
        ref_assign['MostSimilarPop'] = [x.split('_')[-1] for x in ref_assign['MostSimilarPop']]
        target_assign['MostSimilarPop'] = [x.split('_')[-1] for x in target_assign['MostSimilarPop']]

    elif method == 'RandomForest':
        # Assign SuperPop Using Random Forest (PCA loadings)
        logger.debug("Training RandomForest classifier")
        clf_rf = RandomForestClassifier(random_state=32)
        clf_rf.fit(ref_train_df[cols_pcs],  # X (training PCs)
                   ref_train_df[ref_pop_col].astype(str))  # Y (pop label)

        # Predict most similar population using RF classifier
        logger.debug("Find most similar Populations (max RF probability)")
        ref_assign = pd.DataFrame(clf_rf.predict_proba(ref_df[cols_pcs]), index=ref_df.index, columns=['RF_P_{}'.format(x) for x in clf_rf.classes_])
        ref_assign['MostSimilarPop'] = clf_rf.predict(ref_df[cols_pcs])

        target_assign = pd.DataFrame(clf_rf.predict_proba(target_df[cols_pcs]), index=target_df.index, columns=['RF_P_{}'.format(x) for x in clf_rf.classes_])
        target_assign['MostSimilarPop'] = clf_rf.predict(target_df[cols_pcs])

        # Define confidence using p-value thresholds
        ref_assign['MostSimilarPop_LowConfidence'] = np.nan
        target_assign['MostSimilarPop_LowConfidence'] = np.nan

        if p_threshold:
            logger.debug("Comparing Population Similarity to p-value threshold (p < {})".format(p_threshold))
            ref_assign['MostSimilarPop_LowConfidence'] = [(ref_assign['RF_P_{}'.format(x)][i] < p_threshold)
                                                      for i, x in enumerate(clf_rf.predict(ref_df[cols_pcs]))]
            target_assign['MostSimilarPop_LowConfidence'] = [(target_assign['RF_P_{}'.format(x)][i] < p_threshold)
                                                         for i, x in enumerate(clf_rf.predict(target_df[cols_pcs]))]

    return ref_assign, target_assign, compare_info


####
# Methods for adjusting/reporting polygenic score results that account for genetic ancestry
####
normalization_methods = ["empirical", "mean", "mean+var"]


def pgs_adjust(ref_df, target_df, scorecols: list, ref_pop_col, target_pop_col, use_method:list, ref_train_col=None, n_pcs=5):
    """
    Function to adjust PGS using population references and/or genetic ancestry (PCs)
    :param ref_df:
    :param target_df:
    :param scorecols:
    :param ref_pop_col:
    :param target_pop_col:
    :param use_method:
    :param ref_train_col:
    :param n_pcs:
    :return:
    """
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

    ## Create Training dfs
    if ref_train_col:
        assert ref_train_col in ref_df.columns, "Training index column({}) is missing from reference dataframe".format(
            ref_train_col)
        ref_train_df = ref_df.loc[ref_df[ref_train_col] == True,]
    else:
        ref_train_df = ref_df

    ## Create results structures
    results_ref = {}
    results_target = {}
    results_models = {}  # used to store regression information
    for c_pgs in scorecols:
        # Makes melting easier later
        sum_col = 'SUM|{}'.format(c_pgs)
        results_ref[sum_col] = ref_df[c_pgs]
        results_target[sum_col] = target_df[c_pgs]

    # Report PGS values with respect to distribution of PGS in the most similar reference population
    if 'empirical' in use_method:
        logger.debug("Adjusting PGS using most similar reference population distribution.")
        for c_pgs in scorecols:
            # Initialize Output
            percentile_col = 'adj_empirical_percentile|{}'.format(c_pgs)
            results_ref[percentile_col] = pd.Series(index=ref_df.index, dtype='float64')
            results_target[percentile_col] = pd.Series(index=target_df.index, dtype='float64')
            z_col = 'adj_empirical_Z|{}'.format(c_pgs)
            results_ref[z_col] = pd.Series(index=ref_df.index, dtype='float64')
            results_target[z_col] = pd.Series(index=target_df.index, dtype='float64')
            
            r_model = {}

            # Adjust for each population
            for pop in ref_populations:
                r_pop = {}
                i_ref_pop = (ref_df[ref_pop_col] == pop)
                i_target_pop = (target_df[target_pop_col] == pop)

                # Reference Score Distribution
                c_pgs_pop_dist = ref_train_df.loc[ref_train_df[ref_pop_col] == pop, c_pgs]

                # Calculate Percentile
                results_ref[percentile_col].loc[i_ref_pop] = percentileofscore(c_pgs_pop_dist, ref_df.loc[i_ref_pop, c_pgs])
                results_target[percentile_col].loc[i_target_pop] = percentileofscore(c_pgs_pop_dist, target_df.loc[i_target_pop, c_pgs])
                r_pop['percentiles'] = np.percentile(c_pgs_pop_dist, range(0,101,1))

                # Calculate Z
                r_pop['mean'] = c_pgs_pop_dist.mean()
                r_pop['std'] = c_pgs_pop_dist.std(ddof=0)

                results_ref[z_col].loc[i_ref_pop] = (ref_df.loc[i_ref_pop, c_pgs] - r_pop['mean'])/r_pop['std']
                results_target[z_col].loc[i_target_pop] = (target_df.loc[i_target_pop, c_pgs] - r_pop['mean'])/r_pop['std']

                r_model[pop] = r_pop

            results_models['adj_empirical|{}'.format(c_pgs)] = r_model
            # ToDo: explore handling of individuals who have low-confidence population labels
            #  -> Possible Soln: weighted average based on probabilities? Small Mahalanobis P-values will complicate this
    # PCA-based adjustment
    if any([x in use_method for x in ['mean', 'mean+var']]):
        logger.debug("Adjusting PGS using PCA projections")
        for c_pgs in scorecols:
            # Method 1 (Khera): normalize mean (doi:10.1161/CIRCULATIONAHA.118.035658)
            adj_col = 'adj_1_Khera|{}'.format(c_pgs)
            # Fit to Reference Data
            pcs2pgs_fit = LinearRegression().fit(ref_train_df[cols_pcs], ref_train_df[c_pgs])
            ref_train_pgs_pred = pcs2pgs_fit.predict(ref_train_df[cols_pcs])
            ref_train_pgs_resid = ref_train_df[c_pgs] - ref_train_pgs_pred
            ref_train_pgs_resid_mean = ref_train_pgs_resid.mean()
            ref_train_pgs_resid_std = ref_train_pgs_resid.std(ddof=0)

            ref_pgs_resid = ref_df[c_pgs] - pcs2pgs_fit.predict(ref_df[cols_pcs])
            results_ref[adj_col] = ref_pgs_resid / ref_train_pgs_resid_std
            # Apply to Target Data
            target_pgs_pred = pcs2pgs_fit.predict(target_df[cols_pcs])
            target_pgs_resid = target_df[c_pgs] - target_pgs_pred
            results_target[adj_col] = target_pgs_resid / ref_train_pgs_resid_std
            results_models[adj_col] = package_regression(pcs2pgs_fit)

            if 'mean+var' in use_method:
                # Method 2 (Khan): normalize variance (doi:10.1038/s41591-022-01869-1)
                # Normalize based on residual deviation from mean of the distribution [equalize population sds]
                # (e.g. reduce the correlation between genetic ancestry and how far away you are from the mean)
                adj_col = 'adj_2_Khan|{}'.format(c_pgs)
                pcs2var_fit = LinearRegression().fit(ref_train_df[cols_pcs], (ref_train_pgs_resid - ref_train_pgs_resid_mean)**2)

                # Alternative (not adjusting for the mean... which should be 0 already because we've tried to fit it)
                # ---> pcs2var_fit = LinearRegression().fit(ref_train_df[cols_pcs], ref_train_pgs_resid**2)

                # Handle NAs
                pred_var_ref = pcs2var_fit.predict(ref_df[cols_pcs])
                pred_var_ref[pred_var_ref < 0] = np.nan
                results_ref[adj_col] = ref_pgs_resid / np.sqrt(pred_var_ref)
                # Apply to Target
                pred_var_target = pcs2var_fit.predict(target_df[cols_pcs])
                pred_var_target[pred_var_target < 0] = np.nan
                results_target[adj_col] = target_pgs_resid / np.sqrt(pred_var_target)
                results_models[adj_col] = package_regression(pcs2var_fit)

                # Check for NAs
                # has_null = sum(results_ref[adj_col].isnull()) + sum(results_target[adj_col].isnull())
                # if has_null > 0:
                #     print('adj_2_Khan', c_pgs, sum(results_ref[adj_col].isnull()), sum(results_target[adj_col].isnull()))
                #     fig_outloc = 'resid/HasNA/{}.png'.format(c_pgs)
                #     fig = sns.histplot((ref_train_pgs_resid - ref_train_pgs_resid_mean) ** 2)
                #     plt.savefig(fig_outloc)
                #     plt.clf()

                # Attempt gamma distribution for predicted variance to constrain it to be positive
                # (b/c using linear regression we can get negative predictions for the sd)
                pcs2var_fit_gamma = GammaRegressor(max_iter=1000).fit(ref_train_df[cols_pcs], (ref_train_pgs_resid - ref_train_pgs_resid_mean) ** 2)
                adj_col = 'adj_2_Gamma|{}'.format(c_pgs)
                results_ref[adj_col] = ref_pgs_resid / np.sqrt(pcs2var_fit_gamma.predict(ref_df[cols_pcs]))
                results_target[adj_col] = target_pgs_resid / np.sqrt(pcs2var_fit_gamma.predict(target_df[cols_pcs]))
                results_models[adj_col] = package_regression(pcs2var_fit_gamma)

    # Only return results
    logger.debug("Outputting adjusted PGS & models")
    results_ref = pd.DataFrame(results_ref)
    results_target = pd.DataFrame(results_target)
    return results_ref, results_target, results_models


def package_regression(model):
    """Extract relevant details from sklearn regression model"""
    return {'params': model.get_params(),
            '_intercept': model.intercept_,
            '_coef': dict(zip(model.feature_names_in_, model.coef_))
            }


class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types (taken from:
    https://stackoverflow.com/questions/26646362/numpy-array-is-not-json-serializable)"""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def write_model(dict, outname):
    """Use numpy encoder to write json file for models"""
    logger.debug("Writing PGS adjustment models to: {}".format(outname))
    if outname.endswith('.gz'):
        outfile = gzip.open(outname, "wt")
    else:
        outfile = open(outname, "w")

    outfile.write(json.dumps(dict, indent=2, cls=NumpyEncoder))
    outfile.close()
