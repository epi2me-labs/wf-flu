#!/usr/bin/env python
"""Create insaflu db.

You can find insaflu here:
https://raw.githubusercontent.com/INSaFLU/INSaFLU/master/static/db/contigs2sequences/sequences_v10.fasta
You can then make the blast db:
makeblastdb -in sequences -dbtype nucl  -hash_index -title insaflu
"""

from pysam import FastaFile

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")
    out = list()
    sequences_object = FastaFile(args.insaflu_fasta)
    for reference in sequences_object.references:
        if "segment" in reference:
            name, segment, sequence = reference.split('~~~')
            seg, seg_count, seg_name = segment.split('_')
            sequence_stuff = sequence.split("_")
            flu_type = sequence_stuff[0]
            length = len(sequence_stuff)
            if flu_type == "A":
                strain = sequence_stuff[length-2]
            if flu_type == "B":
                strain = sequence_stuff[5]
            accession = sequence_stuff[length-1]
            header = f">insaflu~~~{seg_name}-{strain}~~~{accession}~~~{sequence}"
            out.append(f"{header}\n{sequences_object.fetch(reference)}")

    with open(args.output_fasta, "w") as f:
        f.write("\n".join(out))
        f.write("\n")\

    logger.info(f"Sequences written to {args.output_fasta}")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("prepare_insaflu")

    parser.add_argument("--insaflu_fasta", help="FASTA to prepare.")
    parser.add_argument("--output_fasta", help="Output FASTA file.")
    return parser
