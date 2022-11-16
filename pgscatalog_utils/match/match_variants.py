import argparse
import logging
import os
import shutil
import sys
import textwrap

import polars as pl

import pgscatalog_utils.config as config
from pgscatalog_utils.match import tempdir
from pgscatalog_utils.match.filter import filter_scores
from pgscatalog_utils.match.label import label_matches, make_params_dict
from pgscatalog_utils.match.log import make_logs, make_summary_log, check_log_count
from pgscatalog_utils.match.match import get_all_matches
from pgscatalog_utils.match.read import read_target, read_scorefile
from pgscatalog_utils.match.write import write_log, write_scorefiles

logger = logging.getLogger(__name__)


def match_variants():
    args = _parse_args()
    config.set_logging_level(args.verbose)
    config.setup_polars_threads(args.n_threads)
    config.setup_tmpdir(args.outdir)
    config.setup_cleaning()
    config.OUTDIR = args.outdir

    with pl.StringCache():
        scorefile: pl.LazyFrame = read_scorefile(path=args.scorefile, chrom=args.chrom)
        target_paths = list(set(args.target))
        n_target_files = len(target_paths)
        matches: pl.LazyFrame

        if n_target_files == 0:
            logger.critical("No target genomes found, check the path")
            sys.exit(1)

        if n_target_files == 1 and not args.fast:
            low_memory: bool = True
            match_mode: str = 'single'
        elif n_target_files > 1 and not args.fast:
            low_memory: bool = True
            match_mode: str = 'multi'
        elif args.fast:
            low_memory: bool = False
            match_mode: str = 'fast'

        match match_mode:
            case "single":
                logger.debug(f"Match mode: {match_mode}")
                # _fast_match with low_memory = True reads one target in chunks
                matches: list[list[pl.LazyFrame]] = _fast_match(target_paths, scorefile, low_memory)
            case "multi":
                logger.debug(f"Match mode: {match_mode}")  # iterate over multiple targets, in chunks
                matches: list[list[pl.LazyFrame]] = _match_multiple_targets(target_paths, scorefile, low_memory)
            case "fast":
                logger.debug(f"Match mode: {match_mode}")
                # _fast_match with low_memory = False just read everything into memory for speed
                matches: list[list[pl.LazyFrame]] = _fast_match(target_paths, scorefile, low_memory)
            case _:
                logger.critical(f"Invalid match mode: {match_mode}")
                raise Exception

        dataset = args.dataset.replace('_', '-')  # underscores are delimiters in pgs catalog calculator
        match_dir, matches = _materialise_matches(matches, dataset, low_memory)

        if args.only_match:
            logger.debug(f"--only_match set, writing out match candidates {match_dir} and exiting")
            shutil.move(match_dir, args.outdir)
            logger.debug("Intermediate files can be processed with combine_matches")
            raise SystemExit(0)
        else:
            logger.debug("Labelling match candidates")
            params: dict[str, bool] = make_params_dict(args)
            matches = matches.pipe(label_matches, params)
            logger.debug("Filtering match candidates and making scoring files")
            log_and_write(matches=matches, scorefile=scorefile, dataset=dataset, args=args)


def log_and_write(matches: pl.LazyFrame, scorefile: pl.LazyFrame, dataset: str, args):
    """ Make match logs and write """
    valid_matches, filter_summary = filter_scores(scorefile=scorefile, matches=matches, dataset=dataset,
                                                  min_overlap=args.min_overlap)

    if filter_summary.filter(pl.col("score_pass") == True).collect().is_empty():
        # this can happen when args.min_overlap = 0
        logger.critical("Error: no target variants match any variants in scoring files")
        raise Exception("No valid matches found")

    write_scorefiles(valid_matches, args.split, dataset)

    big_log: pl.LazyFrame = make_logs(scorefile=scorefile, match_candidates=matches, dataset=dataset)
    summary_log: pl.LazyFrame = make_summary_log(match_candidates=matches, filter_summary=filter_summary,
                                                 dataset=dataset,
                                                 scorefile=scorefile)

    check_log_count(summary_log=summary_log, scorefile=scorefile)
    write_log(df=big_log, prefix=dataset, chrom=None, outdir=args.outdir)
    dout = os.path.abspath(config.OUTDIR)
    summary_log.collect().write_csv(os.path.join(dout, f"{dataset}_summary.csv"))


def _materialise_matches(matches: list[list[pl.LazyFrame]], dataset: str, low_memory: bool) -> tuple[str, pl.LazyFrame]:
    """ Collect query plan and store results in temporary files"""
    # outer list: [target_1, target_2]
    # inner list: [ match_1, match_2 ]
    for i, match in enumerate(matches):
        fout = tempdir.get_tmp_path("matches", f"{dataset}_match_{i}.ipc.zst")
        if low_memory:
            pl.concat([x.collect() for x in match]).write_ipc(fout, compression='zstd')
        else:
            pl.concat(pl.collect_all(match)).write_ipc(fout, compression='zstd')
    match_dir: str = tempdir.get_tmp_path("matches", "")
    ldf: pl.LazyFrame = pl.scan_ipc(match_dir + "*.ipc.zst", memory_map=False)
    return match_dir, ldf


def _check_target_chroms(target: pl.LazyFrame) -> None:
    chroms: list[str] = target.select(pl.col("#CHROM").unique()).collect().get_column("#CHROM").to_list()
    if len(chroms) > 1:
        logger.critical(f"Multiple chromosomes detected: {chroms}. Check input data.")
        raise Exception
    else:
        logger.debug("Split target genome contains one chromosome (good)")


def _fast_match(target_paths: list[str], scorefile: pl.LazyFrame, low_memory: bool) -> list[list[pl.LazyFrame]]:
    # fast match is fast because:
    #   1) all target files are read into memory without batching
    #   2) matching occurs without iterating through chromosomes
    target: pl.LazyFrame = read_target(paths=target_paths, low_memory=low_memory)
    return [get_all_matches(scorefile=scorefile, target=target)]


def _match_multiple_targets(target_paths: list[str], scorefile: pl.LazyFrame,
                            low_memory: bool) -> list[list[pl.LazyFrame]]:
    match_lst = []
    for i, loc_target_current in enumerate(target_paths):
        logger.debug(f'Matching scorefile(s) against target: {loc_target_current}')
        target: pl.LazyFrame = read_target(paths=[loc_target_current], low_memory=low_memory)
        if len(target_paths) > 1:
            _check_target_chroms(target)
        match_lst.append(get_all_matches(scorefile=scorefile, target=target))
    return match_lst


def _description_text() -> str:
    return textwrap.dedent('''\
    Match variants from a combined scoring file against a set of
    target genomes from the same fileset, and output scoring files
    compatible with the plink2 --score function.
    
    A combined scoring file is the output of the combine_scorefiles
    script. It has the following structure:
    
        | chr_name | chr_position | ... | accession |
        | -------- | ------------ | --- | --------- |
        | 1        | 1            | ... | PGS000802 |
    
    The combined scoring file is in long format, with one row per
    variant for each scoring file (accession). This structure is
    different to the PGS Catalog standard, because the long format
    makes matching faster and simpler.
    
    Target genomes can be in plink1 bim format or plink2 pvar
    format. Variant IDs should be unique so that they can be specified
    in the scoring file as: variant_id|effect_allele|[effect_weight column(s)...] 
    
    Only one set of target genomes should be matched at a time. Don't
    try to match target genomes from different plink filesets. Matching 
    against a set of chromosomes from the same fileset is OK (see --split). 
   ''')


def _epilog_text() -> str:
    return textwrap.dedent('''\
    match_variants will output at least one scoring file in a
    format compatible with the plink2 --score function. This
    output might be split across different files to ensure each
    variant ID, effect allele, and effect type appears only once
    in each file. Output files have the pattern:

        {dataset}_{chromosome}_{effect_type}_{n}.scorefile.

    If multiple chromosomes are combined into a single file (i.e. not
    --split), then {chromosome} is replaced with 'ALL'. Once the
    scorefiles are used to calculate a score with plink2, the .sscore
    files will need to be aggregated to calculate a single polygenic
    score for each dataset, sample, and accession (scoring file). The
    PGS Catalog Calculator does this automatically.
    ''')


def _parse_args(args=None):
    parser = argparse.ArgumentParser(description=_description_text(), epilog=_epilog_text(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--dataset', dest='dataset', required=True,
                        help='<Required> Label for target genomic dataset')
    parser.add_argument('-s', '--scorefiles', dest='scorefile', required=True,
                        help='<Required> Combined scorefile path (output of read_scorefiles.py)')
    parser.add_argument('-t', '--target', dest='target', required=True, nargs='+',
                        help='<Required> A list of paths of target genomic variants (.bim format)')
    parser.add_argument('-c', '--chrom', dest='chrom', required=False, type=str,
                        help='<Optional> Set which chromosome is in the target variant file to speed up matching ')
    parser.add_argument('-f', '--fast', dest='fast', action='store_true',
                        help='<Optional> Enable faster matching at the cost of increased RAM usage')
    parser.add_argument('--only_match', dest='only_match', action='store_true',
                        help="<Optional> Only match, then write intermediate files, don't make scoring files")
    parser.add_argument('--min_overlap', dest='min_overlap', required=False,
                        type=float, help='<Optional> Minimum proportion of variants to match before error')
    parser = add_match_args(parser) # params for labelling matches
    parser.add_argument('--outdir', dest='outdir', required=True,
                        help='<Required> Output directory')
    parser.add_argument('--split', dest='split', default=False, action='store_true',
                        help='<Optional> Split scorefile per chromosome?')
    parser.add_argument('-n', dest='n_threads', default=1, help='<Optional> n threads for matching', type=int)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='<Optional> Extra logging information')
    return _check_args(parser.parse_args(args))


def add_match_args(parser):
    parser.add_argument('--keep_ambiguous', dest='remove_ambiguous', action='store_false',
                        help='''<Optional> Flag to force the program to keep variants with
                        ambiguous alleles, (e.g. A/T and G/C SNPs), which are normally
                        excluded (default: false). In this case the program proceeds
                        assuming that the genotype data is on the same strand as the
                        GWAS whose summary statistics were used to construct the score.
                                ''')
    parser.add_argument('--keep_multiallelic', dest='remove_multiallelic', action='store_false',
                        help='<Optional> Flag to allow matching to multiallelic variants (default: false).')
    parser.add_argument('--ignore_strand_flips', dest='skip_flip', action='store_true',
                        help='''<Optional> Flag to not consider matched variants that may be reported 
                        on the opposite strand.  Default behaviour is to flip/complement unmatched variants and check if
                        they match.''')
    parser.add_argument('--keep_first_match', dest='keep_first_match', action='store_true',
                        help='''<Optional> If multiple match candidates for a variant exist that can't be prioritised,
                         keep the first match candidate (default: drop all candidates)''')
    return parser


def _check_args(args):
    if args.chrom is not None and not args.only_match:
        # filtering the scoring file will break overlap assumptions and calculations
        # e.g.:
        #   what if one chromosome matches well but another chromosome matches poorly?
        #   what if the poorly matching chromosome only has 5 variants to match?
        #
        # pgsc_calc uses global overlap % to decide if a score fails matching
        # --only_match skips overlap calculations (done in combine_matches instead)
        logger.critical("--chrom requires --only_match")
        sys.exit(1)
    if args.only_match and args.min_overlap is not None:
        # can't calculate min_overlap properly if just checking matches
        logger.critical("Invalid arguments: --only_match and --min_overlap (pick one!)")
        sys.exit(1)
    if not args.only_match and args.min_overlap is None:
        # need to calculate min_overlap before making scoring files
        logger.critical("Invalid arguments: set --min_overlap or --only_match")
        sys.exit(1)
    if args.split and args.only_match:
        # not writing scoring files, so split output doesn't make sense
        logger.critical("Invalid arguments: --only_match and --split (pick one!)")
        sys.exit(1)
    if any([x in sys.argv for x in ['--keep_first_match', '--ignore_strand_flips',
                                    '--keep_multiallelic', '--keep_ambiguous']]):
        logger.warning("Invalid arguments: --only_match and --keep_first_match, --ignore_strand_flips,"
                        "keep_multiallelic, or keep_ambiguous")
        logger.warning("Pass these arguments to combine_matches instead")

    return args


if __name__ == "__main__":
    match_variants()
