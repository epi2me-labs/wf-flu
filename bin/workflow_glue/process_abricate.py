#!/usr/bin/env python
"""Functions to help us parse files for flu."""
import json

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def parse_typing_file(typing_file):
    """Summarise abricate results."""
    abricate = pd.read_csv(typing_file, delimiter="\t", keep_default_na=False)
    if abricate.empty:
        result = {
            'HA': None,
            'NA': None,
            'type': None
        }
        return result

    # output this result
    result = {
        'HA': abricate[abricate['GENE'] == 'HA']['RESISTANCE'].unique().tolist(),
        'NA': abricate[abricate['GENE'] == 'NA']['RESISTANCE'].unique().tolist(),
        'type': abricate[abricate['GENE'] == 'M1']['RESISTANCE'].unique().tolist()
    }

    for k, v in result.items():
        if not v:
            result[k] = ['undetermined']
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
