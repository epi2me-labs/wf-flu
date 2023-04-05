"""Create workflow report."""
import csv
import glob
import json

from dominate.tags import p
import ezcharts as ezc
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd


from .util import get_named_logger, wf_parser  # noqa: ABS101


def typing(sample_details, args):
    """Get typing results and add to a table and csv for export."""
    for sample in sample_details:

        typing = json.load(
            open(f"{args.processed_type}/{sample['alias']}.typing.json"))

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

    with report.add_section("Typing", "Typing"):
        p(
            """
            This table gives the influenza type and strain for each sample. Samples are
            first aligned to IRMA to generate a consensus of alignments and then
            typed with Abricate using the INSaFLU database.
            """
        )

        typing_df = typing(sample_details, args)
        typing_df.columns = typing_df.columns.str.title().str.replace("_", " ")
        DataTable.from_pandas(
            typing_df, use_index=False, export=True, file_name='wf-flu-types')

        p(
            """
            * If HA and NA differ this is because the segments were typed to a different
            flu strain. Full details in the Abricate table.
            """
        )

    with report.add_section("Coverage", "Coverage"):
        dfs = []
        for depth in glob.glob(f'{args.coverage}/*.txt'):
            df = pd.read_csv(depth, sep="\t", header=None, index_col=0)
            df.columns = ['position', depth.replace('.depth.txt', '')]
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
        tabs = Tabs()
        with open(args.nextclade_datasets, "r") as n:
            csv_reader = csv.DictReader(n, delimiter=",")
            for record in csv_reader:
                nxt_json = f"nextclade/{record['dataset']}.json"
                if nxt_json in args.nextclade_files:
                    with tabs.add_tab(f"{record['strain']} - {record['gene']}"):
                        output = dict(
                            sample=list(),
                            strain=list(),
                            gene=list(),
                            coverage=list(),
                            clade=list(),
                            warnings=list())
                        nxt_results = json.load(open(nxt_json))

                        for nxt_result in nxt_results["results"]:
                            output['sample'].append(nxt_result['seqName'])
                            output['strain'].append(record['strain'])
                            output['gene'].append(record['gene'])
                            output['coverage'].append(f"{nxt_result['coverage']:.2f}")
                            output['clade'].append(nxt_result['clade'])
                            output['warnings'].append(
                                ",".join(nxt_result['warnings']))
                        df = pd.DataFrame.from_dict(output)
                        df.columns = df.columns.str.title()
                        DataTable.from_pandas(
                            df,
                            use_index=False,
                            export=True,
                            file_name='wf-flu-nextclade')

    if args.fastqstats:
        with report.add_section("Read summary", "Read summary"):
            fastcat.SeqSummary(args.fastqstats)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument("--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--coverage", default='bed_files',
        help="depth of coverage files")
    parser.add_argument(
        "--typing", default='typing_files',
        help="abricate typing files")
    parser.add_argument(
        "--processed_type", default='processed_typing_files',
        help="processed_abricate typing files")
    parser.add_argument(
        "--fastqstats", default='fastqstats',
        help="fastqstats file from fastcat")
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
