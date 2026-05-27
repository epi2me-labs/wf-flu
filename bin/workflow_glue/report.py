"""Build the final HTML report and results CSV for the workflow.

This file gathers the published per-sample outputs, turns the typing calls
into the summary table, adds the detailed typing and Nextclade sections, and
writes the report the user opens at the end of the run.
"""
import csv
import json
import os
import re
import sys


from dominate.tags import h5, p, span, table, tbody, td, th, thead, tr
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd


from .util import get_named_logger, wf_parser  # noqa: ABS101


def gather_sample_files(sample_details, args):
    """Create dictionary for sample with paths to data files."""
    sample_files = {}
    for sample in sample_details:
        sample_dir = os.path.join(args.data[0], sample["sample"])
        sample_files[sample["sample"]] = {
            "depth": os.path.join(sample_dir, "depth.txt"),
            "mapped_reads": os.path.join(sample_dir, "mapped_reads.txt"),
            "type_json": os.path.join(sample_dir, "processed_type.json"),
            "type_txt": os.path.join(sample_dir, "insaflu.typing.txt")
        }

    return sample_files


def typing(sample_details, sample_files):
    """Get typing results and add to a table and csv for export."""
    for sample in sample_details:
        file_path = sample_files[sample["sample"]]["type_json"]
        if not os.path.exists(file_path):
            continue
        with open(file_path) as handle:
            typing_calls = json.load(handle)

        for i in ['type', 'HA', 'NA']:
            if typing_calls[i] is not None:
                if len(typing_calls[i]) == 1:
                    sample[i] = typing_calls[i][0]
                else:
                    sample[i] = 'undetermined'
            else:
                sample[i] = 'undetermined'

    data = pd.DataFrame(sample_details).fillna("N/A")
    data.rename(columns={"alias": "sample"}, inplace=True)
    for_csv = data
    for_csv.to_csv("wf-flu-results.csv", sep=',', index=False)

    return data


def get_archetype(row):
    """Return archetype string from typing table."""
    if re.match(r"H\d+", row['HA']) and re.match(r"N\d+", row['NA']):
        return row['HA']+row['NA']
    elif row['HA'] == 'undetermined' and row['NA'] == 'undetermined':
        return 'undetermined'
    elif row['HA'] == 'undetermined' and row['NA'] != 'undetermined':
        return row['NA']
    elif row['HA'] != 'undetermined' and row['NA'] == 'undetermined':
        return row['HA']
    else:
        return 'undetermined'


def get_row_format(row_data):
    """Format typing table row."""
    cell = tr()
    badge_dict = {
        "Type_A": "badge bg-primary",
        "Type_B": "badge bg-info",
        "undetermined": "badge bg-warning"}
    badge = badge_dict.get(row_data["Type"], "badge bg-warning")
    for column, row in row_data.items():
        if column in ["Sample", "Barcode"]:
            cell.add(td(span(row)))
        else:
            cell.add(td(h5(span(row.replace("Type_", ""), cls=badge), cls="mb-0")))
    return cell


def selected_panel_label(sample):
    """Return the most specific panel label we can show in the report."""
    sample_type = sample.get("type", "undetermined").replace("Type_", "")
    archetype = get_archetype({
        "HA": sample.get("HA", "undetermined"),
        "NA": sample.get("NA", "undetermined"),
    })
    if sample_type == "undetermined":
        return archetype
    if archetype == "undetermined":
        return sample_type
    return f"{sample_type} / {archetype}"


def mapped_read_counts(sample_details, sample_files):
    """Build a dataframe of mapped reads for the final chosen panel."""
    rows = []
    for sample in sample_details:
        alias = sample["sample"]
        file_path = sample_files[alias]["mapped_reads"]
        if not os.path.exists(file_path):
            continue
        with open(file_path) as handle:
            mapped_reads = int(handle.read().strip() or 0)
        rows.append({
            "Sample": alias,
            "Mapped reads": mapped_reads,
            "Selected panel": selected_panel_label(sample),
        })

    return pd.DataFrame(rows)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Influenza Sequencing Report", "wf-flu",
        args.params, args.versions, args.workflow_version)
    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["sample"])
    sample_files = gather_sample_files(sample_details, args)
    with report.add_section("Typing", "Typing"):
        p(
            """
            This table gives the influenza type and strain for each sample. Samples are
            first aligned to IRMA to generate a consensus of alignments, and then
            typed with Abricate using the INSaFLU database. Please see the table in the
            section below ('Typing details') for full Abricate results. These results
            are especially useful if typing results are discordant.
            """
        )
        typing_df = typing(sample_details, sample_files)
        # add a new column 'Archetype'
        typing_df['Archetype'] = typing_df.apply(lambda row: get_archetype(row), axis=1)
        typing_df.columns = typing_df.columns.str.title()
        typing_df = typing_df[['Sample', 'Barcode', 'Type', 'Archetype']]

        with table(cls="table"):
            with thead():
                for columns in typing_df.columns:
                    th(f"{columns}")
            with tbody():
                for _, row in typing_df.iterrows():
                    get_row_format(row)

    with report.add_section("Typing details", "Typing details"):
        tabs = Tabs()
        with tabs.add_dropdown_menu('Typing details', change_header=True):
            for sample, files in sample_files.items():
                if not os.path.exists(files["type_txt"]):
                    continue
                df = pd.read_csv(files["type_txt"], sep="\t", header=0, index_col=0)
                df = df.drop([
                    'START',
                    'END',
                    'STRAND',
                    'DATABASE',
                    'ACCESSION',
                    'PRODUCT'], axis=1)
                df.rename(columns={'RESISTANCE': 'DETAILS'}, inplace=True)
                df.columns = [x.title() for x in df.columns]
                df.columns = df.columns.str.replace('_', ' ')
                with tabs.add_dropdown_tab(sample):
                    DataTable.from_pandas(df, use_index=False, export=True)

        p(
            """
            Select samples from the drop-down in this table to view detailed Abricate
            results.
            """
        )

    with report.add_section("Coverage", "Coverage"):
        p(
            """
            The bar chart below shows the number of reads mapped to the
            selected final reference panel for each sample, coloured by the
            resolved panel call for that sample. Use the per-sample
            `coverage/depth.txt` files for segment-level depth information.
            """
        )
        mapped_reads = mapped_read_counts(sample_details, sample_files)
        if not mapped_reads.empty:
            order = mapped_reads["Sample"].tolist()
            plot = ezc.barplot(
                data=mapped_reads,
                x="Sample",
                y="Mapped reads",
                hue="Selected panel",
                order=order)
            plot._fig.title.text = "Reads mapped to selected final reference panel"
            plot._fig.xaxis.major_label_orientation = 1.0
            EZChart(plot, 'epi2melabs')
        else:
            p("No mapped-read counts were available for the selected alignments.")

    with report.add_section('Nextclade results', 'Nextclade', True):
        nextclade_data = {}
        with open(args.nextclade_datasets, "r") as n:
            csv_reader = csv.DictReader(n, delimiter=",")
            for record in csv_reader:
                nxt_json = f"nextclade/{record['dataset']}.json"
                if nxt_json in args.nextclade_files:
                    output = {
                        "sample": [],
                        "strain": [],
                        "gene": [],
                        "coverage": [],
                        "clade": [],
                        "warnings": []
                        }
                    try:
                        with open(nxt_json) as handle:
                            nxt_results = json.load(handle)
                    except json.decoder.JSONDecodeError:
                        logger.error(
                            f"Unable to load JSON for {record['dataset']}"
                        )
                        sys.exit(1)
                    for nxt_result in nxt_results["results"]:
                        output['sample'].append(nxt_result['seqName'])
                        output['strain'].append(record['strain'])
                        output['gene'].append(record['gene'])
                        output['coverage'].append(
                            format(nxt_result['coverage'], ".2f")
                        )
                        output['clade'].append(nxt_result['clade'])
                        output['warnings'].append(
                            ",".join(nxt_result['warnings']))
                    nextclade_data[f"{record['strain']} - {record['gene']}"] = output
        if nextclade_data:
            tabs = Tabs()
            for tab_name, results in nextclade_data.items():
                with tabs.add_tab(tab_name):
                    df = pd.DataFrame.from_dict(results)
                    df.columns = df.columns.str.title()
                    DataTable.from_pandas(
                        df,
                        use_index=False,
                        export=True,
                        file_name='wf-flu-nextclade')
        else:
            p(
                """
                No sample met the criteria for Nextclade analysis.
                Please check the quality and strains of your samples.
                Nextclade analysis is available for the following strains:
                H1N1pdm, H3N2, Victoria, Yamagata
                """
            )

    if args.stats:
        with report.add_section("Read summary", "Read summary"):
            names = tuple(d['sample'] for d in sample_details)
            stats = tuple(args.stats)
            if len(stats) == 1:
                stats = stats[0]
                names = names[0]
            fastcat.SeqSummary(stats, sample_names=names)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*',
        help="Fastcat stats directories containing histogram files.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--data", nargs="+", required=True,
        help="Collected outputs per sample")
    parser.add_argument(
        "--nextclade_files", nargs='+', required=True,
        help="Outputs from nextclade")
    parser.add_argument(
        "--nextclade_datasets", required=True,
        help="Nextclade datasets.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--workflow_version", required=True,
        help="Workflow version")
    return parser
