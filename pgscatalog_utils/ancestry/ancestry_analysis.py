import argparse
import textwrap
import logging
import os
import gzip

import pandas as pd

import pgscatalog_utils.config as config
from pgscatalog_utils.ancestry.read import read_pcs, read_pgs, extract_ref_psam_cols
from pgscatalog_utils.ancestry.tools import compare_ancestry, comparison_method_threshold, choose_pval_threshold, \
    pgs_adjust, normalization_methods, write_model

logger = logging.getLogger(__name__)


def ancestry_analysis():
    args = _parse_args()
    config.set_logging_level(args.verbose)
    config.OUTDIR = args.outdir

    # Load PCA data
    maxPCs = max([10, args.nPCs_popcomp, args.nPCs_normalization])  # save memory by not using all PCs
    loc_ref_pcs = args.ref_pcs
    reference_df = read_pcs(loc_pcs=loc_ref_pcs, dataset=args.d_ref,
                            loc_related_ids=args.ref_related, nPCs=maxPCs)
    loc_ref_psam = args.psam
    reference_df = extract_ref_psam_cols(loc_ref_psam, args.d_ref, reference_df, keepcols=[args.ref_label])

    loc_target_sscores = args.target_pcs
    target_df = read_pcs(loc_pcs=loc_target_sscores, dataset=args.d_target, nPCs=maxPCs)

    # Load PGS data & merge with PCA data
    pgs = read_pgs(args.scorefile, onlySUM=True)
    scorecols = list(pgs.columns)
    reference_df = pd.merge(reference_df, pgs, left_index=True, right_index=True)
    target_df = pd.merge(target_df, pgs, left_index=True, right_index=True)
    del pgs  # clear raw PGS from memory

    # Compare target sample ancestry/PCs to reference panel
    assignment_threshold_p = choose_pval_threshold(args)
    ancestry_ref, ancestry_target, compare_info = compare_ancestry(ref_df=reference_df,
                                                                   ref_pop_col=args.ref_label, ref_train_col='Unrelated',
                                                                   target_df=target_df,
                                                                   n_pcs=args.nPCs_popcomp,
                                                                   method=args.method_compare,
                                                                   p_threshold=assignment_threshold_p)

    reference_df = pd.concat([reference_df, ancestry_ref], axis=1)
    target_df = pd.concat([target_df, ancestry_target], axis=1)
    del ancestry_ref, ancestry_target

    # Adjust PGS values
    adjpgs_ref, adjpgs_target, adjpgs_models = pgs_adjust(reference_df, target_df, scorecols,
                                                          args.ref_label, 'MostSimilarPop',
                                                          use_method=args.method_normalization,
                                                          ref_train_col='Unrelated',
                                                          n_pcs=args.nPCs_normalization)
    adjpgs = pd.concat([adjpgs_ref, adjpgs_target], axis=0)
    del adjpgs_ref, adjpgs_target

    # Write outputs
    dout = os.path.abspath(config.OUTDIR)
    if os.path.isdir(dout) is False:
        os.mkdir(dout)
    reference_df['REFERENCE'] = True
    target_df['REFERENCE'] = False
    final_df = pd.concat([target_df, reference_df], axis=0)
    del reference_df, target_df

    # Write Models
    write_model({'pgs': adjpgs_models, 'compare_pcs': compare_info}, os.path.join(dout, f"{args.d_target}_info.json.gz"))

    # Melt & write PGS
    # Currently each PGS will have it's own row... but it might be more optimal for each normalization method
    #  to be on separate rows? My logic is that you might want to check correaltion between methods and it is easiest
    #  in this format.
    loc_pgs_out = os.path.join(dout, f"{args.d_target}_pgs.txt.gz")
    with gzip.open(loc_pgs_out, 'wt') as outf:
        logger.debug('Writing adjusted PGS values (long format) to: {}'.format(loc_pgs_out))
        for i, pgs_id in enumerate(scorecols):
            df_pgs = adjpgs.loc[:, adjpgs.columns.str.endswith(pgs_id)].melt(ignore_index=False)  # filter to 1 PGS
            df_pgs[['method', 'PGS']] = df_pgs.variable.str.split("|", expand=True)
            df_pgs = df_pgs.drop('variable', axis=1).reset_index().pivot(index=['sampleset', 'IID', 'PGS'],
                                                                         columns='method', values='value')
            if i == 0:
                logger.debug('{}/{}: Writing {}'.format(i+1, len(scorecols), pgs_id))
                colorder = list(df_pgs.columns)  # to ensure sort order
                df_pgs.to_csv(outf, sep='\t')
            else:
                logger.debug('{}/{}: Appending {}'.format(i+1, len(scorecols), pgs_id))
                df_pgs[colorder].to_csv(outf, sep='\t', header=False)

    # Write results of PCA & population similarity
    loc_popsim_out = os.path.join(dout, f"{args.d_target}_popsimilarity.txt.gz")
    logger.debug('Writing PCA and popsim results to: {}'.format(loc_popsim_out))
    final_df.drop(scorecols, axis=1).to_csv(loc_popsim_out, sep='\t')
    logger.info("Finished ancestry analysis")


def _description_text() -> str:
    return textwrap.dedent('Program to analyze ancestry outputs of the pgscatalog/pgsc_calc pipeline. Current inputs: '
                           '\n  - PCA projections from reference and target datasets (*.pcs)'
                           '\n  - calculated polygenic scores (e.g. aggregated_scores.txt.gz), '
                           '\n  - information about related samples in the reference dataset (e.g. '
                           'deg2_hg38.king.cutoff.out.id).')


def _parse_args(args=None):
    parser = argparse.ArgumentParser(description=_description_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--dataset', dest='d_target', required=True,
                        help='<Required> Label of the TARGET genomic dataset')
    parser.add_argument('-r', '--reference', dest='d_ref', required=True,
                        help='<Required> Label of the REFERENCE genomic dataset')
    parser.add_argument('--ref_pcs', dest='ref_pcs', required=True, nargs='+',
                        help='<Required> Principal components path (output from fraposa_pgsc)')
    parser.add_argument('--target_pcs', dest='target_pcs', required=True, nargs='+',
                        help='<Required> Principal components path (output from fraposa_pgsc)')
    parser.add_argument('--psam', dest='psam', required=True,
                        help='<Required> Reference sample information file path in plink2 psam format)')
    parser.add_argument('-x', '--reference_related', dest='ref_related',
                        help='File of related sample IDs (excluded from training ancestry assignments)')
    parser.add_argument('-p', '--pop_label', dest='ref_label', default='SuperPop',
                        help='Population labels in REFERENCE psam to use for assignment')
    parser.add_argument('-s', '--agg_scores', dest='scorefile', default='aggregated_scores.txt.gz',
                        help='Aggregated scores in PGS Catalog format ([sampleset, IID] indexed)')
    parser.add_argument('-a', '--ancestry_method', dest='method_compare',
                        choices=comparison_method_threshold.keys(), default='RandomForest',
                        help='Method used for population/ancestry assignment')
    parser.add_argument('--n_popcomp', dest='nPCs_popcomp', type=int, metavar="[1-20]", choices=range(1, 21),
                        default=5,
                        help='Number of PCs used for population comparison (default = 5)')
    parser.add_argument('-t', '--pval_threshold', dest='pThreshold', type=float,
                        help='p-value threshold used to identify low-confidence ancestry similarities')
    parser.add_argument('-n', '--normalization_method', nargs='+', dest='method_normalization',
                        choices=normalization_methods, default=["empirical", "mean", "mean+var"],
                        help='Method used for adjustment of PGS using genetic ancestry')
    parser.add_argument('--n_normalization', dest='nPCs_normalization', type=int, metavar="[1-20]",
                        choices=range(1, 21), default=4,
                        help='Number of PCs used for population NORMALIZATION (default = 4)')
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help='<Required> Output directory')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return parser.parse_args(args)


if __name__ == "__main__":
    ancestry_analysis()
