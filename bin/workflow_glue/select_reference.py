#!/usr/bin/env python
"""Pick the reference panel that the second pass should align against.

The workflow first makes a broad typing call, then this file uses that call to
write a smaller FASTA for the real polishing pass. Influenza A samples keep
the A internal segments plus the best HA and NA references, influenza B
samples keep the B segments, and uncertain samples fall back to the full panel
so we do not over-prune too early.
"""
import json

import pysam

from .typing_resolver import (  # noqa: ABS101
    load_typing_table,
    resolve_typing_calls,
    TYPE_A,
    TYPE_B,
)
from .util import get_named_logger, wf_parser  # noqa: ABS101


def choose_reference_headers(calls, headers):
    """Choose which records from the reference panel to keep."""
    if calls["type"] == TYPE_A:
        # Always keep the shared A internal segments, then narrow HA/NA to the
        # selected subtype when we have a confident call.
        selected = [
            header for header in headers
            if header.startswith("A_")
            and not header.startswith("A_HA_")
            and not header.startswith("A_NA_")
        ]

        selected_ha = (
            f"A_HA_{calls['HA']}" if calls["HA"] != "undetermined" else None
        )
        if selected_ha in headers:
            selected.append(selected_ha)
        else:
            # If the requested HA is absent from the panel, widen back out to
            # every A HA record rather than leaving HA unrepresented.
            selected.extend(header for header in headers if header.startswith("A_HA_"))

        selected_na = (
            f"A_NA_{calls['NA']}" if calls["NA"] != "undetermined" else None
        )
        if selected_na in headers:
            selected.append(selected_na)
        else:
            # Apply the same conservative fallback for NA when a precise match
            # is unavailable in the bundled panel.
            selected.extend(header for header in headers if header.startswith("A_NA_"))

        return selected

    if calls["type"] == TYPE_B:
        return [header for header in headers if header.startswith("B_")]

    return list(headers)


def write_reference(reference_panel, selected_headers, output):
    """Write a FASTA containing just the selected records."""
    with pysam.FastaFile(reference_panel) as fasta, open(output, "w") as handle:
        for header in selected_headers:
            sequence = fasta.fetch(header)
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def main(args):
    """Run the entry point."""
    logger = get_named_logger("select_ref")
    typing = load_typing_table(args.typing)
    with pysam.FastaFile(args.reference) as reference_panel:
        headers = list(reference_panel.references)
    calls = resolve_typing_calls(typing)
    fallback = calls["type"] == "undetermined"
    selected_headers = choose_reference_headers(calls, headers)

    write_reference(args.reference, selected_headers, args.output)

    summary = {
        "sample_alias": args.sample_alias,
        "selected_type": calls["type"],
        "selected_ha": calls["HA"],
        "selected_na": calls["NA"],
        "selected_lineage": calls["lineage"],
        "fallback_to_full_panel": fallback,
        "selected_headers": selected_headers,
    }
    with open(args.summary, "w") as handle:
        json.dump(summary, handle, indent=4)

    logger.info(
        "Selected reference for %s: %s", args.sample_alias, ",".join(selected_headers)
    )


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("select_reference")
    parser.add_argument(
        "--typing", required=True, help="First-pass abricate typing TSV."
    )
    parser.add_argument(
        "--reference", required=True, help="Reference panel FASTA used for selection."
    )
    parser.add_argument("--output", required=True, help="Selected reference FASTA.")
    parser.add_argument("--summary", required=True, help="Selection summary JSON.")
    parser.add_argument("--sample_alias", required=True, help="Sample alias.")
    return parser
