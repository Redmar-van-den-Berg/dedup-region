#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go


def read_cov_depth(fname):
    depths = list()
    with open(fname) as fin:
        for line in fin:
            chrom, pos, depth = line.strip().split()
            depths.append(int(depth))
    return depths


def main(args):
    before = read_cov_depth(args.before)
    after = read_cov_depth(args.after)

    plot_coverage(before, after, args.output)


def plot_coverage(before, after, fname):
    fig = go.Figure(
            data=[go.Scatter(y=before, name="Before"),
                go.Scatter(y=after, name="After")
                ],
        )

    fig.update_layout(
            title=f"Exon coverage",
            xaxis_title="Position",
            yaxis_title="Coverage"
    )
    fig.write_html(fname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--before', required=True, help='Coverage before deduplication')
    parser.add_argument('--after', required=True, help='Coverage after deduplication')
    parser.add_argument('--output', required=True)

    arguments = parser.parse_args()
    main(arguments)
