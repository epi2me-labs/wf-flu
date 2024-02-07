#!/usr/bin/env python
"""Functions to help us parse files for flu."""
import csv
import json
import os

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def process_typing(typing_json):
    """Process abricate typing string."""
    typing = json.load(open(typing_json))
    # if it's a result that's mixed we can't do nextclade
    if typing['type'] is None:
        return None
    if len(typing['HA']) > 1:
        return None
    if len(typing['NA']) > 1:
        return None
    if len(typing['type']) > 1:
        return None
    return typing


def find_nextclade(typing, nextclade_datasets):
    """Find nextclade datasets suitable for sample."""
    datasets = list()
    with open(nextclade_datasets, "r") as n:
        csv_reader = csv.DictReader(n, delimiter=",")

        flu_type = typing['type'][0]

        if flu_type == "Type_A":
            strain = f'{typing["HA"][0]}{typing["NA"][0]}'
        elif flu_type == "Type_B":
            strain = f'{typing["HA"][0]}'
        elif flu_type == "undetermined":
            raise ValueError("Flu type is undetermined. \n \
                             If reads are from the RBK protocol, ensure \
                             --rbk has been set to prevent overfiltering.")
        else:
            raise ValueError(f'{flu_type} is not Type_A or Type_B')

        for record in csv_reader:
            if record["strain"] == strain:
                datasets.append(
                    {
                        "type": flu_type.replace("Type_", ""),
                        "strain": strain,
                        "dataset": record["dataset"],
                        "gene": record["gene"],
                    }
                )

    return datasets


def make_consensus(datasets, consensus, alias):
    """Make the consensus file for nextclade."""
    fasta = pysam.FastaFile(consensus)
    results = list()
    for dataset in datasets:
        result = dict(dataset=dataset["dataset"])

        if dataset['type'] == "A":
            if dataset["gene"] == "HA":
                strain = f'{dataset["type"]}_{dataset["gene"]}_{dataset["strain"][:2]}'
                sequence = fasta.fetch(strain)
            elif dataset["gene"] == "NA":
                strain = f'{dataset["type"]}_{dataset["gene"]}_{dataset["strain"][2:4]}'
                sequence = fasta.fetch(strain)
        elif dataset['type'] == "B":
            sequence = fasta.fetch(f'{dataset["type"]}_{dataset["gene"]}')
        sequence_id = f">{alias}"
        result['consensus'] = f'{sequence_id}\n{sequence}\n'
        results.append(result)

    return results


def main(args):
    """Run the entry point."""
    logger = get_named_logger("nextclade_helper")
    typing = process_typing(args.typing)
    if typing is not None:
        datasets = find_nextclade(typing, args.nextclade_datasets)
        consensus = make_consensus(datasets, args.consensus, args.sample_alias)

        if len(consensus) == 0:
            logger.warning("Flu strain not currently available in nextclade.")

        for i in consensus:
            dataset_dir = f'datasets/{i["dataset"]}'
            os.mkdir(dataset_dir)
            with open(f'{dataset_dir}/{args.sample_alias}.nextclade.fasta', "w") as f:
                f.write(i["consensus"])
    else:
        logger.warning("typing information mixed.")
    logger.info("nextclade helping done.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("nextclade_helper")
    parser.add_argument(
        "--typing",
        help="Typing json from abricate.")
    parser.add_argument(
        "--consensus",
        help="Consensus FASTA from sample.")
    parser.add_argument(
        "--nextclade_datasets",
        help="Dataset list from nextclade.")
    parser.add_argument(
        "--sample_alias",
        help="The sample alias.")
    return parser
