#!/usr/bin/env python
"""Make a fasta suitable for blast db."""

import re

import pysam

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("create_ncbi_database")
    parser.add_argument("file", help="Genome file")
    return parser


def main(args):
    """Run the entry point."""
    genome_file = args.file

    fasta_open = pysam.Fastafile(genome_file)
    for reference in fasta_open.references:

        if 'reassortant' in reference:
            continue
        try:
            accession, genotype, long, country = reference.split("|")
        except ValueError:
            continue

        try:
            r1 = re.compile(r'(.*?)\s*\(')
            m1 = r1.match(long)
            virus_tmp = m1.group(1).replace("Influenza_", "")
            virus = virus_tmp.replace("_virus", "").replace("_", "")
        except AttributeError:
            continue

        accession = accession.replace("_", "")
        awful = re.findall(r'\((.*?)\)', long)

        if len(awful) == 0:
            continue
        try:
            if "(" in awful[0]:
                strain, geno = awful[0].split("(")
            else:
                strain = awful[0]
        except ValueError:
            continue

        if len(awful) > 1:
            gene = awful[1]

        if genotype == '':
            genotype = '-'

        if gene not in ['HA', 'M1', 'NA']:
            continue

        file = open("database.fasta", "w")
        file.write(f"""
            >epi2melabs~~~{gene}~~~{accession}~~~{virus}|
            {genotype} {strain.replace('_','')}\n""")

        file.write(f"{fasta_open.fetch(reference)}\n")
        file.close()
