import logging

from pgscatalog_utils.scorefile.effect_type import set_effect_type
from pgscatalog_utils.scorefile.effect_weight import melt_effect_weights
from pgscatalog_utils.scorefile.genome_build import build2GRC
from pgscatalog_utils.scorefile.harmonised import remap_harmonised
from pgscatalog_utils.scorefile.liftover import liftover
from pgscatalog_utils.scorefile.log import make_log
from pgscatalog_utils.scorefile.qc import quality_control
from pgscatalog_utils.scorefile.read import load_scorefile
from pgscatalog_utils.scorefile.write import write_scorefile

logger = logging.getLogger(__name__)


def combine(scorefile_path, args, output_file_lock) -> dict:
    # Read scorefile df and header
    h, score = load_scorefile(scorefile_path)

    if score.empty:
        logger.critical(
            f"Empty scorefile {scorefile_path} detected! Please check the input data")
        raise Exception

    # Check if we should use the harmonized positions
    use_harmonised = False
    current_build = None
    if h.get('HmPOS_build') is not None:
        if h.get('HmPOS_build') == args.target_build:
            use_harmonised = True
            current_build = h.get('HmPOS_build')
        else:
            logger.error(
                f"Cannot combine {scorefile_path} (harmonized to {h.get('HmPOS_build')}) in target build {args.target_build}")
            raise Exception

    # Process/QC score and check variant columns
    score = (score.pipe(remap_harmonised, use_harmonised=use_harmonised)
             .pipe(quality_control, drop_missing=args.drop_missing)
             .pipe(melt_effect_weights)
             .pipe(set_effect_type))

    # Annotate score with the genome_build (in GRCh notation)
    if current_build is None:
        current_build = build2GRC(h.get('genome_build'))
        if current_build is None:
            logger.error("Scorefile has no build information, "
                         "please add the build to the header with "
                         "('#genome_build=[insert variant build]")
            raise Exception

    score = score.assign(genome_build=current_build)

    if (current_build != args.target_build) and (args.liftover is False):
        logger.error(
            f"Cannot combine {scorefile_path} (build={h.get('genome_build')}) with target build {args.target_build} without liftover")
        logger.error("Try running with --liftover and specifying the --chain_dir")
        raise Exception

    if args.liftover:
        logger.debug("Annotating scorefile with liftover parameters")
        score = liftover(score, args.chain_dir, args.min_lift, args.target_build)

    if score.empty and (args.drop_missing is False):
        logger.critical(
            "Empty output score detected, something went wrong while combining")
        raise Exception

    write_scorefile(df=score, path=args.outfile, lock=output_file_lock)

    return make_log(args, scorefile_path=scorefile_path, use_harmonised=use_harmonised,
                    score=score, h=h, score_shape_original=score.shape)
