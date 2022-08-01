#!/usr/bin/env python
"""Make a fasta sutibale for blast db."""

import argparse
import re

import pysam

parser = argparse.ArgumentParser()

parser.add_argument("file")

args = parser.parse_args()

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

    print(f"""
        >epi2melabs~~~{gene}~~~{accession}~~~{virus}|
        {genotype} {strain.replace('_','')}""")

    print(fasta_open.fetch(reference))
