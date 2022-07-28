#!/usr/bin/env python
"""Functions to help us parse files for flu."""


def parse_typing_file(typing_file):
    """
    Parse abricate file.

    Extract the sample typecode from an abricate typing file.
    From original ONT Applications code 2021.
    """
    try:

        with open(typing_file) as file:
            lines = file.readlines()

        # iterate over table, note this may be empty for failed runs
        genedict = {}
        for line in lines:
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            columns = line.strip('\r').strip('\n').split('\t')

            # gene code
            gene = columns[5]
            # type entry
            g_type = columns[-1]

            if gene not in genedict:
                genedict[gene] = g_type
            else:
                # flag genes appearing twice or more as 'mixed'
                genedict[gene] = 'mixed'

        # check for list of known genes
        core_genes = ['M1', 'HA']
        mixed = False
        unknown = False
        flu_type = []
        for gene in core_genes:
            if gene in genedict:
                g_type = genedict[gene]
                if g_type == 'mixed':
                    mixed = True
                    continue
                else:
                    flu_type.append(g_type)
            else:
                unknown = True
                continue
        extra_genes = ['NA']
        for gene in extra_genes:
            if gene in genedict:
                g_type = genedict[gene]
                if g_type == 'mixed':
                    mixed = True
                    continue
                else:
                    flu_type.append(g_type)

        # create typecode
        if mixed is False and unknown is False:
            typecode = ' '.join(flu_type)
        elif mixed is True:
            typecode = 'mixed'
        elif unknown is True:
            typecode = 'unknown'
    except ValueError:
        typecode = 'failure'

    return typecode
