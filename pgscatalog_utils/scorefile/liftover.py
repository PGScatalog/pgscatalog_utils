import logging
import os

import pandas as pd
import pyliftover

from pgscatalog_utils.scorefile.genome_build import annotate_build

logger = logging.getLogger(__name__)


def liftover(df: pd.DataFrame, chain_dir: str, min_lift: float, target_build: str) -> pd.DataFrame:
    """ Liftover genomic coordinates to a different genome build """
    df = annotate_build(df, target_build)  # get chain_target_build (e.g. in hg notation to match chain files)

    mapped, unmapped = pd.DataFrame(), pd.DataFrame()
    no_liftover: pd.DataFrame = df.query('chain_target_build == chain_genome_build')
    to_liftover: pd.DataFrame = df.query('chain_target_build != chain_genome_build')

    if no_liftover.empty:
        logger.debug("Liftover required for all scorefile variants")
    else:
        logger.debug("Skipping liftover for scorefiles with same build as target genome")
        no_liftover.loc[:, ['lifted_chr', 'lifted_pos']] = no_liftover[
            ['chr_name', 'chr_position']]  # assume col structure
        no_liftover.assign(liftover=None)

    if to_liftover.empty:
        logger.debug("Liftover skipped because no variants required it")
    else:
        lo: dict[str, pyliftover.LiftOver] = _create_liftover(chain_dir)  # loads chain files
        logger.debug("Lifting over scoring files")
        lifted: pd.DataFrame = to_liftover.apply(_convert_coordinates, axis=1, lo_dict=lo)
        to_liftover = pd.concat([to_liftover, lifted], axis=1)
        logger.debug("Liftover complete")

        mapped: pd.DataFrame = (to_liftover[~to_liftover[['lifted_chr', 'lifted_pos']].isnull().any(axis=1)]
                                .assign(liftover=True))
        unmapped: pd.DataFrame = (to_liftover[to_liftover[['lifted_chr', 'lifted_pos']].isnull().any(axis=1)] \
                                  .assign(liftover=False))
        _check_min_liftover(mapped, unmapped, min_lift)

    return pd.concat([mapped, unmapped, no_liftover])


def _check_min_liftover(mapped: pd.DataFrame, unmapped: pd.DataFrame, min_lift: float) -> None:
    """ Check that liftover process met minimum parameters"""
    df = pd.concat([mapped, unmapped])
    n_variants: pd.DataFrame = (pd.DataFrame(df.groupby('accession')['liftover'].count())
                                .reset_index()
                                .rename({'liftover': 'n_var'}, axis=1))
    lo_counts = (pd.DataFrame(df.groupby(['accession', 'liftover'])['liftover'].count()) \
                 .rename_axis(['accession', 'liftover_status'])
                 .reset_index())
    summary: pd.DataFrame = lo_counts.merge(n_variants, on='accession')
    summary['proportion'] = summary['liftover'] / summary['n_var']

    for row in summary.query('liftover_status == True')[['accession', 'proportion']].itertuples():
        if row.proportion < min_lift:
            logger.error(f'Liftover failed for scorefile {row.accession}')
            logger.error(f'{row.proportion} of variants lifted over, less than min_lift parameter ({min_lift})')
            raise Exception
        else:
            logger.debug(f'Minimum liftover threshold passed for scorefile {row.accession}')


def _convert_coordinates(df: pd.Series, lo_dict: dict[str, pyliftover.LiftOver]) -> pd.Series:
    """ Convert genomic coordinates to different build """
    converted: list[tuple[str, int, str, int]] | None

    if df[['chr_name', 'chr_position']].isnull().values.any():
        converted = None
    else:
        lo = lo_dict[df['chain_genome_build'] + df['chain_target_build']]  # extract lo object from dict
        chrom: str = 'chr' + str(df['chr_name'])
        pos: int = int(df['chr_position']) - 1  # liftOver is 0 indexed, VCF is 1 indexed
        # converted example: [('chr22', 15460378, '+', 3320966530)] or None
        converted = lo.convert_coordinate(chrom, pos)

    if converted:
        lifted_chrom: str = _parse_lifted_chrom(converted[0][0][3:])  # return first matching liftover
        lifted_pos: int = int(converted[0][1]) + 1  # reverse 0 indexing
        return pd.Series([lifted_chrom, lifted_pos], index=['lifted_chr', 'lifted_pos'])
    else:
        return pd.Series([None, None], index=['lifted_chr', 'lifted_pos'])


def _parse_lifted_chrom(i: str) -> str:
    """ Convert lifted chromosomes to tidy integers

    liftover needs chr suffix for chromosome input (1 -> chr1), and it also
    returns weird chromosomes sometimes (chr22 -> 22_KI270879v1_alt)
    """
    return i.split('_')[0]


def _create_liftover(chain_dir: str) -> dict['str': pyliftover.LiftOver]:
    """ Create LiftOver objects that can remap genomic coordinates """
    builds: list[str] = ["hg19hg38", "hg38hg19"]
    chains: list[str] = [os.path.join(chain_dir, x) for x in ["hg19ToHg38.over.chain.gz", "hg38ToHg19.over.chain.gz"]]
    lo: list[pyliftover.LiftOver] = [pyliftover.LiftOver(x) for x in chains]
    logger.debug("Chain files loaded for liftover")
    return dict(zip(builds, lo))
