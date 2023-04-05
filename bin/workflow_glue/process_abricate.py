#!/usr/bin/env python
"""Functions to help us parse files for flu."""
import json

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


to_replace = {
    'pdm09': {
        'HA': 'H1N1',
        'NA': 'H1N1',
        'other': 'H1N1'}}


def process_strain_result(input_list, segment):
    """Process strain results."""
    if len(input_list) == 0:
        return [None]
    result = list()
    for strain in input_list:
        if strain in to_replace:
            result.append(to_replace[strain][segment])
        else:
            result.append(strain)
    return list(set(result))


def parse_typing_file(typing_file):
    """Summarise abricate results."""
    abricate = pd.read_csv(typing_file, delimiter="\t")
    if abricate.empty:
        result = {
            'HA': None,
            'NA': None,
            'other': None,
            'type': None
        }
        return result

    # split up some columns to make them useful
    abricate[['segment', 'strain']] = abricate['GENE'].str.split('-', 1, expand=True)
    abricate[['type', 'x']] = abricate['SEQUENCE'].str.split('_', 1, expand=True)

    # get any results for HA and NA segments
    ha_strains = abricate[abricate['segment'] == 'HA']['strain'].tolist()
    na_strains = abricate[abricate['segment'] == 'NA']['strain'].tolist()

    # get any results for other segments/core genes
    other_strains = abricate[
                (abricate['segment'] != 'NA') &
                (abricate['segment'] != 'HA')]['strain'].tolist()

    ha_strains_result = process_strain_result(ha_strains, 'HA')
    na_strains_result = process_strain_result(na_strains, 'NA')
    other_strains_result = process_strain_result(other_strains, 'other')

    # output this result
    result = {
        'HA': ha_strains_result,
        'NA': na_strains_result,
        'other': other_strains_result,
        'type': list(abricate['type'].unique())
    }

    return result


def main(args):
    """Run the entry point."""
    logger = get_named_logger("process_abricate")
    result = parse_typing_file(args.typing)

    with open(args.output, 'w') as f:
        f.write(json.dumps(result, indent=4))

    logger.info(f"Typing result written to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("process_abricate")
    parser.add_argument(
        "--typing", default=None,
        help="abricate typing results.")
    parser.add_argument(
        "--output", default=None,
        help="processed abricate results.")
    return parser
