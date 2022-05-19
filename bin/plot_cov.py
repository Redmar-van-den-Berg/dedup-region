#!/usr/bin/env python3

import argparse
import plotly.graph_objects as go

import gtf

def read_cov_depth(fname):
    depths = list()
    with open(fname) as fin:
        for line in fin:
            chrom, pos, depth = line.strip().split()
            depths.append(int(depth))
    return depths


def guess_gene_name(fname, transcript):
    with open(fname) as fin:
        for record in gtf.gtf_to_json(fin, '', transcript):
            return record['attribute']['gene_name']


def guess_transcript(fname):
    """ Guess the transcript name from the file path """
    print(fname)
    return fname.split('/')[1]


def main(args):
    before = read_cov_depth(args.before)
    after = read_cov_depth(args.after)

    # Guess the transcript name for the figure heading
    transcript_name = guess_transcript(args.before)

    # Get the gene name from the transcript name
    gene_name = guess_gene_name(args.gtf, transcript_name)

    display_name = f"{gene_name} ({transcript_name})"

    plot_coverage(before, after, args.output, display_name)


def plot_coverage(before, after, fname, display_name):
    fig = go.Figure(
            data=[go.Scatter(y=before, name="Before"),
                go.Scatter(y=after, name="After")
                ],
        )

    fig.update_layout(
            title=f"Exon coverage for {display_name}",
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
    parser.add_argument('--gtf', required=True, help='Used to determine the gene name')
    parser.add_argument('--output', required=True)

    arguments = parser.parse_args()
    main(arguments)
