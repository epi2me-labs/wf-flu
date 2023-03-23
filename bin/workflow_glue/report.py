#!/usr/bin/env python
"""Create workflow report."""

import json
import math

from aplanat import bars
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import ont_colors
from bokeh.models import (
    BasicTicker,  ColorBar, ColumnDataSource, LinearColorMapper)
from bokeh.plotting import figure
import pandas as pd

from .flu import parse_typing_file  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101


def qc(sample_details, fastqstats):
    """Make a QC megaplot."""
    df = pd.DataFrame(
        columns=[
            'sample',
            'barcode',
            'read_count',
            'mean_quality',
            'mean_length'])
    total_reads = 0

    fastqstats_df = pd.read_csv(fastqstats, sep='\t', header=0)

    for sample in sorted(sample_details):

        if sample == "unclassified":
            continue

        reads = fastqstats_df[fastqstats_df['sample_name'] == sample]

        df = df.append(
            {
                'sample': sample,
                'barcode': sample_details[sample]['barcode'],
                'read_count': len(reads.index),
                'mean_quality': reads["mean_quality"].mean(),
                'mean_length': reads["read_length"].mean()},
            ignore_index=True)

        total_reads += len(reads.index)

    read_count = bars.simple_bar(
        df['sample'].astype(str),
        df['read_count'],
        colors=[colors.BRAND_BLUE]*len(df.index),
        title=(
            'Number of reads per sample'),
        plot_width=1100
    )
    read_count.xaxis.major_label_orientation = math.pi/4

    mean_length = bars.simple_bar(
        df['sample'].astype(str),
        df['mean_length'],
        colors=[colors.BRAND_LIGHT_BLUE]*len(df.index),
        title=(
            'Mean read length for each sample'),
        plot_width=1100
    )
    mean_length.xaxis.major_label_orientation = math.pi/4

    return ([read_count, mean_length], total_reads)


def unclassified(files, total_reads):
    """Return a status bar of number of unclassified reads."""
    unclassified_reads = pd.read_csv(files['fastqstats'], sep='\t', header=0)

    count = len(unclassified_reads.index)

    percentage = round((count / (count+total_reads)) * 100, 1)

    progress_bar = f"""
    <div class="progress align-items-center" style="height:40px">
        <div class="progress-bar
                progress-bar-striped
                bg-brand-primary" role="progressbar"
                style="width: {percentage}%;height:40px;"
                aria-valuenow="{percentage}"
                aria-valuemin="0" aria-valuemax="100">
        </div>
        <span class="justify-content-center d-flex position-absolute w-100">
        <h6 style="margin-bottom: 0;">
            {percentage}% of reads are unclassified</h6>
        </span>
    </div>"""

    return progress_bar


def typing(sample_details):
    """Get typing results and add to a table and csv for export."""
    header = ['sample', 'barcode', 'type', 'subtype']
    typing_table = ['''<table class="table table-sm">
                        <thead>
                            <tr>
                                <th>Sample</th>
                                <th>Barcode</th>
                                <th>Typing</th>
                            </tr>
                        </thead>''']

    out = [','.join(header)]

    for sample in sorted(sample_details):

        if sample == "unclassified":
            continue

        typing_result = parse_typing_file(sample_details[sample]['typing'])
        typing = typing_result.split(" ")

        line = []
        line.append(sample)
        line.append(sample_details[sample]['barcode'])
        line.append(typing[0])

        if typing[0] == "Type_A":
            color = 'brand-primary'
        elif typing[0] == "Type_B":
            color = 'brand-info'
        elif typing[0] == "unknown":
            color = 'brand-light-grey'
        elif typing[0] == "mixed":
            color = 'brand-red'

        extra = ''
        safe = ''
        if len(typing) > 1:
            safe = typing.copy()
            line.append(''.join(safe[1:]))
            extra = f"""
                <span class="badge badge-{color}">{''.join(safe[1:])}</span>"""

        typing_table.append(f"""
            <tr>
                <td>{sample}</td>
                <td>{sample_details[sample]['barcode']}</td>
                <td>
                    <h6>
                        <span class="badge badge-{color}">{typing[0]}</span>
                        {extra}
                    </h6>
                </td>
            </tr>
        """)

        out.append(','.join(str(x) for x in line))

    typing_table.append("</table>")

    csv = '\n'.join(out)
    f = open("wf-flu-results.csv", "w")
    f.write(csv)
    f.close()
    return "\n".join(typing_table)


def coverage(sample_details):
    """Plot coverage of segments for each sample."""
    all_summaries = None
    for sample in sample_details:

        if sample == "unclassified":
            continue

        coverage = pd.read_csv(
            sample_details[sample]['coverage'], sep='\t', header=None)
        coverage = coverage.rename(
            columns={0: 'segment', 1: 'pos', 2: 'coverage'})

        summary = coverage.groupby('segment').agg(
            {'coverage': ['mean', 'median', 'min', 'max']})
        summary = summary.assign(sample=sample)
        if 'all_summaries' in locals():
            all_summaries = pd.concat([all_summaries, summary])
        else:
            all_summaries = summary

    p = figure(
        title=None,
        toolbar_location=None,
        x_axis_location="above",
        x_range=sorted(list(set(all_summaries.index.values))),
        y_range=sorted(list(set(all_summaries['sample']))),
        width=1100,
        height=100+(30*len(sample_details))
    )
    ont_gradient = [
        "#ffffff",
        "#e2edf2",
        "#c5dbe6",
        "#a8c9da",
        "#8ab7cd",
        "#6ba6c1",
        "#4695b5",
        "#0084a9"]

    source = ColumnDataSource(all_summaries)
    mapper = LinearColorMapper(
        palette=ont_gradient,
        low=all_summaries["coverage"]["median"].min(),
        high=all_summaries["coverage"]["median"].max())

    p.rect(
        source=source,
        width=1,
        height=1,
        x='segment',
        y='sample_',
        fill_color={'field': 'coverage_median', 'transform': mapper},
        line_color=None)

    color_bar = ColorBar(
        color_mapper=mapper,
        location=(0, 0),
        ticker=BasicTicker(desired_num_ticks=len(ont_gradient))
    )

    p.add_layout(color_bar, 'right')

    p.xaxis.major_label_orientation = math.pi/4
    p.xaxis.axis_label = 'Segment'
    p.yaxis.axis_label = 'Sample'
    p.toolbar.active_drag = None
    p.toolbar.active_scroll = None
    p.toolbar.active_tap = None

    return p


def main(args):
    """Run the entry point."""
    global colors
    colors = ont_colors

    with open(args.metadata) as metadata:
        sample_details = {
            d['alias']: {
                'type': d['type'],
                'barcode': d['barcode'],
                'coverage': f"{args.coverage}/{d['alias']}.depth.txt",
                'typing': f"{args.typing}/{d['alias']}.insaflu.typing.txt",
                'fastqstats': f"{args.fastqstats}/{d['alias']}.stats"
            } for d in json.load(metadata)
        }

    report = WFReport(
        "wf-flu Influenza Sequencing Report", "wf-flu",
        revision=args.revision, commit=args.commit)

    section = report.add_section()
    section._add_item("""<h3>Sample Typing</h3>

    <p>The table below details the typing for each sample using the insaflu
    blast db assigned using abricate.</p>""")

    section._add_item(typing(sample_details))

    section = report.add_section()
    section._add_item("""<h3>Quality Control</h3>

    <p>This section contains plots and tables that might be useful in
    determining the success of a run or samples on that run.</p>""")

    qc_results = qc(sample_details, args.fastqstats)
    for plot in qc_results[0]:

        section.plot(plot)

    total_reads = qc_results[1]

    if "unclassfied" in sample_details:
        progress_bar = unclassified(
            sample_details['unclassified'], total_reads)

        section._add_item(f"""
        <div class="card bg-light mt-4 mb-4">
      <div class="card-body">
        <div class="row align-items-center">
            <div class="col-md-6">
            <h5>Unclassified Reads</h5>
        <p>Unclassified reads are those reads which cannot be assigned
        a barcode with high confidence by guppy.</p>
        </div>
         <div class="col-md-6">{progress_bar}</div></div></div></div>""")

    section = report.add_section()
    section._add_item("""<h3>Segment Coverage</h3>

    <p>The heatmap below shows all of the samples on the run (y-axis) and
    the median coverage of each of the segments analysed.</p>""")

    section.plot(coverage(sample_details))

    section._add_item("""<br><br>""")

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
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
        "--coverage", default='bed_files',
        help="depth of coverage files")
    parser.add_argument(
        "--typing", default='typing_files',
        help="abricate typing files")
    parser.add_argument(
        "--fastqstats", default='fastqstats',
        help="fastqstats file from fastcat")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    return parser
