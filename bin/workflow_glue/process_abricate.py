#!/usr/bin/env python
"""Turn raw Abricate output into the one typing call the workflow stores.

Abricate gives us a table of hits, but the workflow stores a compact JSON
typing result. This file runs the shared resolver logic and writes that JSON
for the rest of the pipeline.
"""
import json

from .typing_resolver import (  # noqa: ABS101
    as_processed_json,
    load_typing_table,
    resolve_typing_calls,
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def parse_typing_file(typing_file):
    """Summarise Abricate results as a single best typing call."""
    abricate = load_typing_table(typing_file)
    return as_processed_json(resolve_typing_calls(abricate))


def main(args):
    """Run the entry point."""
    logger = get_named_logger("process_abricate")
    result = parse_typing_file(args.typing)

    with open(args.output, 'w') as f:
        json.dump(result, f, indent=4)

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
