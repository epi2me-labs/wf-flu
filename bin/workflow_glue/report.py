"""Create workflow report."""
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
        sample_dir = os.path.join(args.data[0], sample["alias"])
        sample_files[sample["alias"]] = {
            "depth": os.path.join(sample_dir, "depth.txt"),
            "type_json": os.path.join(sample_dir, "processed_type.json"),
            "type_txt": os.path.join(sample_dir, "insaflu.typing.txt")
        }

    return sample_files


def typing(sample_details, sample_files):
    """Get typing results and add to a table and csv for export."""
    for sample in sample_details:
        file_path = sample_files[sample["alias"]]["type_json"]
        if not os.path.exists(file_path):
            continue
        typing = json.load(open(file_path))

        for i in ['type', 'HA', 'NA']:
            if typing[i] is not None:
                if len(typing[i]) == 1:
                    sample[i] = typing[i][0]
                else:
                    sample[i] = 'undetermined'
            else:
                sample[i] = 'undetermined'

    data = pd.DataFrame(sample_details)
    for_csv = data.drop('type', axis=1)
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


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Influenza Sequencing Report", "wf-flu",
        args.params, args.versions)
    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'alias': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["alias"])
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
        typing_df.columns = typing_df.columns.str.title().str.replace("_", " ")
        typing_df.rename(columns={'Ha': 'HA', 'Na': 'NA'}, inplace=True)
        # add a new column 'Archetype'
        typing_df['Archetype'] = typing_df.apply(lambda row: get_archetype(row), axis=1)
        typing_df = typing_df[['Alias', 'Barcode', 'Type', 'Archetype']]

        with table(cls="table"):
            with thead():
                for columns in typing_df.columns:
                    th(f"{columns}")
            with tbody():
                for index, row in typing_df.iterrows():
                    with tr():
                        for i in range(4):
                            cell = td()
                            if row[i] == 'undetermined':
                                cell.add(h5(span(
                                    'Undetermined', cls="badge bg-warning"),
                                    cls="mb-0"))
                            elif row.index[i] == 'Type' and row[i] == 'Type_A':
                                cell.add(
                                    h5(span('A', cls="badge bg-primary"), cls="mb-0"))
                            elif row.index[i] == 'Type' and row[i] == 'Type_B':
                                cell.add(h5(span('B', cls="badge bg-info"), cls="mb-0"))
                            elif row.index[i] == 'Archetype' and row.Type == "Type_A":
                                cell.add(h5(span(
                                    row[i], cls="badge bg-primary"), cls="mb-0"))
                            elif row.index[i] == 'Archetype' and row.Type == "Type_B":
                                cell.add(
                                    h5(span(row[i], cls="badge bg-info"), cls="mb-0"))
                            elif row.index[i] == 'Barcode' and row[i] is None:
                                cell.add(span("None"))
                            else:
                                cell.add(span(row[i]))

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
        dfs = []
        for sample, files in sample_files.items():
            if not os.path.exists(files["depth"]):
                continue
            df = pd.read_csv(files["depth"], sep="\t", header=None, index_col=0)
            df.columns = ['position', sample]
            df.drop(columns='position', inplace=True)
            df.index.name = 'segment'
            median = df.groupby('segment').median()
            dfs.append(median)

        data = pd.concat(dfs, join='outer', axis=1)
        data = data.rename_axis("sample", axis=1)

        plot = ezc.heatmap(data, vmin=0, annot=False)

        EZChart(plot, 'epi2melabs')

        p(
            """
            The heatmap shows the median coverage per segment for each sample.
            Each box in the heatmap represents one segment in a sample and is
            colour-coded using the range of values in the slider (from zero to maximum
            median coverage across the whole batch). The slider can be manipulated
            to filter the heatmap by coverage levels, enabling a quick assessment of
            the coverage for each sample.
            """
        )

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
                        nxt_results = json.load(open(nxt_json))
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
                            f"{nxt_result['coverage']:.2f}"
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
            fastcat.SeqSummary(args.stats)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*',
        help="Fastcat per-read stats file(s).")
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
    return parser
