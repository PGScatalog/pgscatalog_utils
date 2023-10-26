import json

from pgscatalog_utils.scorefile.read import get_scorefile_basename

headers2logs = [
    'pgs_id',
    'pgp_id',
    'pgs_name',
    'genome_build',
    'variants_number',
    'trait_reported',
    'trait_efo',
    'trait_mapped',
    'weight_type',
    'citation'
]
headers2logs_harmonisation = [
    'HmPOS_build',
    'HmPOS_date',
    'HmPOS_match_chr',
    'HmPOS_match_pos'
]


def make_log(args, scorefile_path, score, h, use_harmonised, score_shape_original):
    # Build Score header logs
    score_logs = {}
    score_id = get_scorefile_basename(scorefile_path)
    score_header = score_logs[score_id] = {}

    # Scoring file header information
    for header in headers2logs:
        header_val = h.get(header)
        if (header in ['trait_efo', 'trait_mapped']) and (header_val is not None):
            header_val = header_val.split('|')
        score_header[header] = header_val
    # Other header information
    score_header['columns'] = list(score.columns)
    score_header['use_liftover'] = False
    if args.liftover:
        score_header['use_liftover'] = True
    # Harmonized header information
    score_header['use_harmonised'] = use_harmonised
    if use_harmonised:
        score_header['sources'] = sorted(score['hm_source'].unique().tolist())
        for hm_header in headers2logs_harmonisation:
            hm_header_val = h.get(hm_header)
            if hm_header_val:
                if hm_header.startswith('HmPOS_match'):
                    hm_header_val = json.loads(hm_header_val)
                score_header[hm_header] = hm_header_val
    if score_header['variants_number'] is None:
        score_header['variants_number'] = score_shape_original[0]

    return score_logs
