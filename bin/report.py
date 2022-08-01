#!/usr/bin/env python
"""Create workflow report."""

import argparse
import math

from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from bokeh.models import (
    BasicTicker,  ColorBar, ColumnDataSource, LinearColorMapper)
from bokeh.plotting import figure
from flu import parse_typing_file
import pandas as pd


def typing(sample_details):
    """Get typing results and add to a table."""
    typing_table = ['''<table class="table table-sm">
                        <thead>
                            <tr><th>Sample</th><th>Typing</th></tr>
                        </thead>''']
    for sample in sorted(sample_details):
        print(sample)
        typing_result = parse_typing_file(sample_details[sample]['typing'])
        typing = typing_result.split(" ")
        print(typing)
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

            extra = f"""
                <span class="badge badge-{color}">{''.join(safe[1:])}</span>"""

        typing_table.append(f"""
            <tr>
                <td>{sample}</td>
                <td>
                    <h5>
                        <span class="badge badge-{color}">{typing[0]}</span>
                        {extra}
                    </h5>
                </td>
            </tr>
        """)

    typing_table.append("</table>")

    return "\n".join(typing_table)


def coverage(sample_details):
    """Plot coverage of segments for each sample."""
    all_summaries = None
    for sample in sample_details:
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
        height=400
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

    p.xaxis.major_label_orientation = math.pi/2
    p.xaxis.axis_label = 'Segment'
    p.yaxis.axis_label = 'Sample'
    p.toolbar.active_drag = None
    p.toolbar.active_scroll = None
    p.toolbar.active_tap = None

    return p


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
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
        "--samples", nargs='+', default='unknown',
        help="git commit number")
    parser.add_argument(
        "--types", nargs='+', default='unknown',
        help="git commit number")

    args = parser.parse_args()

    sample_details = {
        sample: {
            'type': type,
            'coverage': f"{args.coverage}/{sample}.depth.txt",
            'typing': f"{args.typing}/{sample}.insaflu.typing.txt"
        } for sample, type in zip(
            args.samples, args.types
        )
    }
    print(sample_details)
    report = WFReport(
        "wf-flu Influenza Sequencing Report", "wf-flu",
        revision=args.revision, commit=args.commit)

    section = report.add_section()
    section._add_item("""<h3>Segment Coverage</h3>

    <p>The heatmap below shows all of the samples on the run (y-axis) and
    the median coverage of each of the segments analysed.</p>""")

    section.plot(coverage(sample_details))

    section._add_item("""<br><br>""")

    section = report.add_section()
    section._add_item("""<h3>Sample Typing</h3>

    <p>The table below details the typing for each sample using the insaflu
    blast db assigned using abricate.</p>""")

    section._add_item(typing(sample_details))

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
